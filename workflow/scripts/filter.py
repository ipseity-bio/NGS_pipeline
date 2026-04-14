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


def process_vcf(
    file_pattern: str,
    variant_summary_path: str,
    output_dir: str,
    haplotype_dir: str,
    freebayes_dir: str,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(glob.glob(file_pattern))
    clinvar = pd.read_csv(variant_summary_path, sep="\t")

    for file in files:
        vcf = read_vep_tab(file)
        if "#CHROM" in vcf.columns:
            vcf = vcf.rename(columns={"#CHROM": "CHROM"})

        vcf["Existing_variation"] = vcf["Existing_variation"].str.split(",")
        vcf = vcf.explode("Existing_variation").dropna(subset=["Existing_variation"])
        vcf["Existing_variation"] = vcf["Existing_variation"].astype(str)

        x = vcf[vcf["Existing_variation"].str.startswith("rs")].copy()
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

        with open(partner_path, "r", encoding="utf-8", errors="replace") as handle:
            lines_1 = [line for line in handle if not line.startswith("##")]

        haplocall = pd.read_csv(
            io.StringIO("".join(lines_1)), dtype=str, sep="\t"
        ).rename(columns={"#CHROM": "CHROM"})

        sample_col = haplocall.columns[-1]
        final = add_pos_candidates(final, location_col="Location")
        haplocall["POS"] = pd.to_numeric(haplocall["POS"], errors="coerce")

        m1 = final.merge(haplocall, left_on="POS_START", right_on="POS", how="left")
        need2 = m1[sample_col].isna()
        m2 = final.loc[need2].merge(haplocall, left_on="POS_ANCHOR", right_on="POS", how="left")
        m1.loc[need2, haplocall.columns] = m2[haplocall.columns].values
        dp_gq_zygo = m1

        dp_gq_zygo["Genotype (GT)"] = dp_gq_zygo[sample_col].str.split(":").str[0]
        dp_gq_zygo["Depth of Coverage (DP)"] = dp_gq_zygo[sample_col].str.split(":").str[2]
        dp_gq_zygo["Genotype Quality (GQ)"] = dp_gq_zygo[sample_col].str.split(":").str[3]
        dp_gq_zygo["AD_ALT"] = dp_gq_zygo[sample_col].str.split(":").str[1].str.split(",").str[1]
        dp_gq_zygo["AD_ALT"] = pd.to_numeric(dp_gq_zygo["AD_ALT"], errors="coerce")
        dp_gq_zygo["Depth of Coverage (DP)"] = pd.to_numeric(
            dp_gq_zygo["Depth of Coverage (DP)"], errors="coerce"
        )
        dp_gq_zygo["VAF (Sample)"] = (
            dp_gq_zygo["AD_ALT"] / dp_gq_zygo["Depth of Coverage (DP)"]
        )

        if os.path.exists(freebayes_path):
            fb_full = load_freebayes_full(freebayes_path)
            dp_gq_zygo["CHROM"] = (
                dp_gq_zygo["Location"].astype(str).str.split(r"[:|-]").str[0]
            )
            dp_gq_zygo = add_pos_candidates(dp_gq_zygo, location_col="Location")

            m1 = dp_gq_zygo.merge(
                fb_full,
                left_on=["CHROM", "POS_START"],
                right_on=["CHROM", "POS"],
                how="left",
                suffixes=("", "_fb"),
            )
            need2 = m1["GT_fb"].isna()
            m2 = m1.loc[need2].drop(columns=fb_full.columns.difference(["CHROM"]))
            m2 = m2.merge(
                fb_full,
                left_on=["CHROM", "POS_ANCHOR"],
                right_on=["CHROM", "POS"],
                how="left",
                suffixes=("", "_fb"),
            )
            m1.loc[need2, fb_full.columns] = m2[fb_full.columns].values
            dp_gq_zygo = m1

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

            dp_gq_zygo["ALT_IDX"] = dp_gq_zygo.apply(
                lambda row: pick_alt_idx(row["Allele"], row["ALT_LIST"]), axis=1
            )
            dp_gq_zygo["AD_ALT_fb"] = dp_gq_zygo.apply(
                lambda row: row["AD_LIST"][int(row["ALT_IDX"])]
                if isinstance(row["AD_LIST"], list)
                and pd.notna(row["ALT_IDX"])
                and int(row["ALT_IDX"]) < len(row["AD_LIST"])
                else np.nan,
                axis=1,
            )
            dp_gq_zygo["VAF_fb"] = pd.to_numeric(
                dp_gq_zygo["AD_ALT_fb"], errors="coerce"
            ) / pd.to_numeric(dp_gq_zygo["DP_fb"], errors="coerce")

            dp_gq_zygo["Genotype (GT)"] = dp_gq_zygo["Genotype (GT)"].where(
                dp_gq_zygo["Genotype (GT)"].notna()
                & (dp_gq_zygo["Genotype (GT)"] != ""),
                dp_gq_zygo["GT_fb"],
            )
            dp_gq_zygo["Depth of Coverage (DP)"] = dp_gq_zygo["Depth of Coverage (DP)"].where(
                dp_gq_zygo["Depth of Coverage (DP)"].notna(),
                pd.to_numeric(dp_gq_zygo["DP_fb"], errors="coerce"),
            )
            dp_gq_zygo["Genotype Quality (GQ)"] = dp_gq_zygo["Genotype Quality (GQ)"].where(
                dp_gq_zygo["Genotype Quality (GQ)"].notna(),
                pd.to_numeric(dp_gq_zygo["GQ_fb"], errors="coerce"),
            )
            dp_gq_zygo["VAF (Sample)"] = dp_gq_zygo["VAF (Sample)"].where(
                dp_gq_zygo["VAF (Sample)"].notna(), dp_gq_zygo["VAF_fb"]
            )

        out = dp_gq_zygo[
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
