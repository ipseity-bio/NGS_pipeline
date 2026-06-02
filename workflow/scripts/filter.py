import glob
import io
import os
import re

import numpy as np
import pandas as pd


def read_vcf_no_meta(path: str) -> pd.DataFrame:
    header = None
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                header = raw.lstrip("#").rstrip("\n").split("\t")
                continue
            if raw.strip():
                rows.append(raw.rstrip("\n").split("\t"))
    if header is None:
        raise ValueError(f"No #CHROM header found in {path}")
    return pd.DataFrame(rows, columns=header, dtype=str).rename(
        columns={"#CHROM": "CHROM"}
    )


def parse_sample_by_format(df: pd.DataFrame) -> pd.DataFrame:
    sample_col = df.columns[-1]
    fmt_keys = df["FORMAT"].fillna("").str.split(":")
    sample_vals = df[sample_col].fillna("").str.split(":")

    parsed = []
    for keys, vals in zip(fmt_keys, sample_vals):
        parsed.append({k: v for k, v in zip(keys, vals)})

    parsed_df = pd.DataFrame(parsed)
    for col in ["GT", "DP", "GQ", "AD", "AO", "RO"]:
        if col not in parsed_df.columns:
            parsed_df[col] = np.nan

    out = pd.concat([df.reset_index(drop=True), parsed_df], axis=1)
    ad_alt = out["AD"].astype(str).str.split(",").str[1]
    out["AD_ALT"] = pd.to_numeric(ad_alt, errors="coerce")
    out["DP_num"] = pd.to_numeric(out["DP"], errors="coerce")
    out["VAF_fb"] = out["AD_ALT"] / out["DP_num"]
    return out


def add_pos_candidates(df: pd.DataFrame, location_col: str = "Location") -> pd.DataFrame:
    start_pos = (
        df[location_col].astype(str).str.split(r"[:|-]", regex=True).str[1]
    )
    start_pos = pd.to_numeric(start_pos, errors="coerce")
    df["POS_START"] = start_pos
    df["POS_ANCHOR"] = start_pos - 1
    return df


def load_freebayes_full(freebayes_path: str) -> pd.DataFrame:
    fb = parse_sample_by_format(read_vcf_no_meta(freebayes_path))
    fb["POS"] = pd.to_numeric(fb["POS"], errors="coerce")
    fb["ALT_LIST"] = fb["ALT"].astype(str).str.split(",")
    fb["AD_LIST"] = fb["AD"].astype(str).str.split(",").apply(
        lambda xs: [pd.to_numeric(x, errors="coerce") for x in xs]
        if isinstance(xs, list)
        else []
    )
    fb = fb.rename(columns={"GT": "GT_fb", "DP": "DP_fb", "GQ": "GQ_fb"})
    return fb[
        ["CHROM", "POS", "REF", "ALT", "ALT_LIST", "AD_LIST", "GT_fb", "DP_fb", "GQ_fb"]
    ]


def read_vep_tab(path: str) -> pd.DataFrame:
    header_cols = None
    data_rows = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.lstrip().startswith("##"):
                continue
            line = raw.lstrip().rstrip("\n")
            if header_cols is None and line.startswith("#"):
                header_cols = re.split(r"\t+", line[1:].strip())
                continue
            if line:
                data_rows.append(line.split("\t"))

    if header_cols is None:
        raise ValueError(f"No VEP header found in {path}")

    width = len(header_cols)
    normalized_rows = []
    for row in data_rows:
        if len(row) < width:
            row = row + [""] * (width - len(row))
        elif len(row) > width:
            row = row[:width]
        normalized_rows.append(row)

    return pd.DataFrame(normalized_rows, columns=header_cols, dtype=str)


def first_non_null(series: pd.Series):
    series = series.dropna()
    series = series[series.astype(str).str.strip() != ""]
    if not series.empty:
        return series.iloc[0]
    return np.nan


def fill_with_anchor_fallback(
    left_df: pd.DataFrame,
    lookup_df: pd.DataFrame,
    left_on_start,
    left_on_anchor,
    right_on,
    value_columns,
) -> pd.DataFrame:
    out = left_df.merge(lookup_df, left_on=left_on_start, right_on=right_on, how="left")
    missing = out[value_columns[0]].isna()
    if missing.any():
        fallback = left_df.loc[missing].copy()
        fallback["_row_id"] = fallback.index
        fallback_merge = fallback.merge(
            lookup_df, left_on=left_on_anchor, right_on=right_on, how="left"
        )
        fallback_merge = fallback_merge.drop_duplicates(subset="_row_id", keep="first")
        fallback_merge = fallback_merge.set_index("_row_id")
        target_index = out.index[missing]
        aligned = fallback_merge.reindex(target_index)
        out.loc[target_index, value_columns] = aligned[value_columns].to_numpy()
    return out


def attach_sample_metrics(
    df: pd.DataFrame, partner_path: str, freebayes_path: str
) -> pd.DataFrame:
    with open(partner_path, "r", encoding="utf-8", errors="replace") as handle:
        lines_1 = [line for line in handle if not line.startswith("##")]

    haplocall = pd.read_csv(
        io.StringIO("".join(lines_1)), dtype=str, sep="\t"
    ).rename(columns={"#CHROM": "CHROM"})

    sample_col = haplocall.columns[-1]
    df = add_pos_candidates(df, location_col="Location")
    haplocall["POS"] = pd.to_numeric(haplocall["POS"], errors="coerce")

    out = fill_with_anchor_fallback(
        left_df=df,
        lookup_df=haplocall,
        left_on_start="POS_START",
        left_on_anchor="POS_ANCHOR",
        right_on="POS",
        value_columns=list(haplocall.columns),
    )

    out["Genotype (GT)"] = out[sample_col].str.split(":").str[0]
    out["Depth of Coverage (DP)"] = out[sample_col].str.split(":").str[2]
    out["Genotype Quality (GQ)"] = out[sample_col].str.split(":").str[3]
    out["AD_ALT"] = out[sample_col].str.split(":").str[1].str.split(",").str[1]
    out["AD_ALT"] = pd.to_numeric(out["AD_ALT"], errors="coerce")
    out["Depth of Coverage (DP)"] = pd.to_numeric(
        out["Depth of Coverage (DP)"], errors="coerce"
    )
    out["VAF (Sample)"] = out["AD_ALT"] / out["Depth of Coverage (DP)"]

    if os.path.exists(freebayes_path):
        fb_full = load_freebayes_full(freebayes_path)
        out["CHROM"] = out["Location"].astype(str).str.split(r"[:|-]").str[0]
        out = add_pos_candidates(out, location_col="Location")

        out = fill_with_anchor_fallback(
            left_df=out,
            lookup_df=fb_full,
            left_on_start=["CHROM", "POS_START"],
            left_on_anchor=["CHROM", "POS_ANCHOR"],
            right_on=["CHROM", "POS"],
            value_columns=list(fb_full.columns),
        )

        def pick_alt_idx(allele, alt_list):
            if not isinstance(alt_list, list) or not alt_list:
                return np.nan
            allele = str(allele)
            for i, alt in enumerate(alt_list):
                if allele == alt:
                    return i + 1
            for i, alt in enumerate(alt_list):
                if allele in alt or alt in allele:
                    return i + 1
            return int(np.argmin([abs(len(alt) - len(allele)) for alt in alt_list])) + 1

        out["ALT_IDX"] = out.apply(
            lambda row: pick_alt_idx(row["Allele"], row["ALT_LIST"]), axis=1
        )
        out["AD_ALT_fb"] = out.apply(
            lambda row: row["AD_LIST"][int(row["ALT_IDX"])]
            if isinstance(row["AD_LIST"], list)
            and pd.notna(row["ALT_IDX"])
            and int(row["ALT_IDX"]) < len(row["AD_LIST"])
            else np.nan,
            axis=1,
        )
        out["VAF_fb"] = pd.to_numeric(
            out["AD_ALT_fb"], errors="coerce"
        ) / pd.to_numeric(out["DP_fb"], errors="coerce")

        out["Genotype (GT)"] = out["Genotype (GT)"].where(
            out["Genotype (GT)"].notna() & (out["Genotype (GT)"] != ""), out["GT_fb"]
        )
        out["Depth of Coverage (DP)"] = out["Depth of Coverage (DP)"].where(
            out["Depth of Coverage (DP)"].notna(),
            pd.to_numeric(out["DP_fb"], errors="coerce"),
        )
        out["Genotype Quality (GQ)"] = out["Genotype Quality (GQ)"].where(
            out["Genotype Quality (GQ)"].notna(),
            pd.to_numeric(out["GQ_fb"], errors="coerce"),
        )
        out["VAF (Sample)"] = out["VAF (Sample)"].where(
            out["VAF (Sample)"].notna(), out["VAF_fb"]
        )

    return out


def build_output_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df[
        [
            "Location",
            "REF_ALLELE",
            "Allele",
            "BIOTYPE",
            "Existing_variation",
            "SYMBOL",
            "GeneSymbol",
            "CLIN_SIG",
            "ClinicalSignificance",
            "SPDI",
            "HGVSg",
            "HGVSc",
            "HGVSp",
            "HGVS_OFFSET",
            "SIFT",
            "PolyPhen",
            "VAF (Sample)",
            "MAX_AF",
            "Genotype (GT)",
            "ZYG",
            "Depth of Coverage (DP)",
            "Genotype Quality (GQ)",
            "PhenotypeIDS",
            "PhenotypeList",
        ]
    ].drop_duplicates().reset_index(drop=True)

    out = out.rename(
        columns={
            "ClinicalSignificance": "ClinicalSignificance (ClinVar)",
            "MAX_AF": "MAX_AF (Population)",
        }
    )
    return out


def process_vcf(
    file_pattern: str,
    variant_summary_path: str,
    output_dir: str,
    haplotype_dir: str,
    freebayes_dir: str,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(glob.glob(file_pattern))
    clinvar = pd.read_csv(variant_summary_path, sep="\t", low_memory=False)

    for file in files:
        vcf = read_vep_tab(file)
        if "#CHROM" in vcf.columns:
            vcf = vcf.rename(columns={"#CHROM": "CHROM"})

        variant_key = [
            col
            for col in [
                "Location",
                "REF_ALLELE",
                "Allele",
                "BIOTYPE",
                "SYMBOL",
                "SPDI",
                "HGVSg",
                "HGVSc",
                "HGVSp",
                "HGVS_OFFSET",
                "SIFT",
                "PolyPhen",
                "MAX_AF",
                "ZYG",
                "CLIN_SIG",
            ]
            if col in vcf.columns
        ]

        vcf_lookup = vcf.copy()
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].fillna("").astype(str)
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].str.split(",")
        vcf_lookup = vcf_lookup.explode("Existing_variation")
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].fillna("").astype(str).str.strip()
        vcf_lookup["rs_id"] = pd.to_numeric(
            vcf_lookup["Existing_variation"].str.extract(r"^rs(\d+)$")[0],
            errors="coerce",
            downcast="integer",
        )

        clinvar_left = pd.merge(
            vcf_lookup,
            clinvar,
            left_on="rs_id",
            right_on="RS# (dbSNP)",
            how="left",
        )
        clinvar_left = clinvar_left[
            clinvar_left["Assembly"].isna()
            | clinvar_left["Assembly"].astype(str).str.startswith("GRCh37")
        ]

        all_candidates = clinvar_left.groupby(variant_key, dropna=False, as_index=False).agg(
            {
                "Existing_variation": first_non_null,
                "GeneSymbol": first_non_null,
                "ClinicalSignificance": first_non_null,
                "PhenotypeIDS": first_non_null,
                "PhenotypeList": first_non_null,
                "rs_id": first_non_null,
            }
        )
        all_candidates["GeneSymbol"] = all_candidates["GeneSymbol"].fillna(
            all_candidates.get("SYMBOL")
        )

        x = vcf_lookup[vcf_lookup["Existing_variation"].str.startswith("rs")].copy()
        x["rs_id"] = pd.to_numeric(
            x["Existing_variation"].str.replace("rs", "", regex=False),
            errors="coerce",
            downcast="integer",
        )

        result_df = pd.merge(
            x, clinvar, left_on="rs_id", right_on="RS# (dbSNP)", how="inner"
        )
        result_df = result_df[
            result_df["Assembly"].astype(str).str.startswith("GRCh37")
        ].drop_duplicates()

        result_df["CLIN_SIG"] = result_df["CLIN_SIG"].str.split(",")
        result_df = result_df.explode("CLIN_SIG").drop_duplicates().reset_index(drop=True)
        final = result_df.drop_duplicates(
            subset=[
                "SPDI",
                "CLIN_SIG",
                "REF_ALLELE",
                "Allele",
                "Existing_variation",
                "PhenotypeList",
            ],
            keep="first",
        ).reset_index(drop=True)

        fname = os.path.basename(file)
        sample = fname.replace("_merged_vep.vcf", "").replace("_vep.vcf", "")

        partner_path = os.path.join(haplotype_dir, f"{sample}_variants.vcf")
        freebayes_path = os.path.join(freebayes_dir, f"{sample}_freebayes.vcf")

        all_candidates = attach_sample_metrics(all_candidates, partner_path, freebayes_path)
        all_out = build_output_table(all_candidates)
        all_out.to_csv(os.path.join(output_dir, f"{sample}_all_candidates.csv"), index=False)

        dp_gq_zygo = attach_sample_metrics(final, partner_path, freebayes_path)
        out = build_output_table(dp_gq_zygo)

        out.to_csv(os.path.join(output_dir, f"{sample}_raw_output.csv"), index=False)

        out_f1 = out[
            out["CLIN_SIG"].isin(
                [
                    "pathogenic",
                    "pathogenic/established_risk_allele",
                    "likely_pathogenic",
                    "pathogenic/likely_pathogenic",
                ]
            )
        ]
        out_f1 = out_f1.drop_duplicates(
            subset=out_f1.columns.difference(
                ["CLIN_SIG", "ClinicalSignificance (ClinVar)", "PhenotypeIDS", "PhenotypeList"]
            )
        )
        out_f1.to_csv(os.path.join(output_dir, f"{sample}_output_p.csv"), index=False)
