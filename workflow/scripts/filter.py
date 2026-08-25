import glob
import gzip
import io
import os
import re

import numpy as np
import pandas as pd


def read_vcf_no_meta(path: str) -> pd.DataFrame:
    header = None
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
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
    start_pos = df[location_col].astype(str).str.split(r"[:|-]", regex=True).str[1]
    start_pos = pd.to_numeric(start_pos, errors="coerce")
    df["POS_START"] = start_pos
    df["POS_ANCHOR"] = start_pos - 1
    return df


def normalize_chromosome(value) -> str:
    text = str(value).strip() if not pd.isna(value) else ""
    if text.lower().startswith("chr"):
        text = text[3:]
    return text


def add_vcf_match_keys(df: pd.DataFrame, location_col: str = "Location") -> pd.DataFrame:
    out = add_pos_candidates(df.copy(), location_col=location_col)
    out["CHROM_KEY"] = out[location_col].astype(str).str.split(":").str[0].apply(normalize_chromosome)
    out["REF_KEY"] = out.get("REF_ALLELE", pd.Series("", index=out.index)).fillna("").astype(str).str.strip()
    out["ALT_KEY"] = out.get("Allele", pd.Series("", index=out.index)).fillna("").astype(str).str.strip()
    return out


def prepare_clinvar_lookup(clinvar: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "RS# (dbSNP)",
        "GeneSymbol",
        "ClinicalSignificance",
        "PhenotypeIDS",
        "PhenotypeList",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    available = [col for col in columns if col in clinvar.columns]
    out = clinvar[available].copy()
    if "Assembly" in out.columns:
        out = out[out["Assembly"].fillna("").astype(str).str.startswith("GRCh37")].copy()
    out["rs_id"] = pd.to_numeric(out.get("RS# (dbSNP)", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["CLINVAR_CHROM_KEY"] = out.get("Chromosome", pd.Series("", index=out.index)).apply(normalize_chromosome)
    out["CLINVAR_POS"] = pd.to_numeric(out.get("PositionVCF", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["CLINVAR_REF_KEY"] = out.get("ReferenceAlleleVCF", pd.Series("", index=out.index)).fillna("").astype(str).str.strip()
    out["CLINVAR_ALT_KEY"] = out.get("AlternateAlleleVCF", pd.Series("", index=out.index)).fillna("").astype(str).str.strip()
    return out


def merge_clinvar_annotations(vcf_lookup: pd.DataFrame, clinvar_lookup: pd.DataFrame) -> pd.DataFrame:
    keyed_vcf = add_vcf_match_keys(vcf_lookup, location_col="Location")
    rsid_lookup = clinvar_lookup[
        clinvar_lookup["rs_id"].notna() & (clinvar_lookup["rs_id"] > 0)
    ].copy()

    rsid_matches = keyed_vcf.merge(
        rsid_lookup,
        on="rs_id",
        how="left",
    )

    coordinate_lookup = clinvar_lookup[
        clinvar_lookup["CLINVAR_CHROM_KEY"].ne("")
        & clinvar_lookup["CLINVAR_POS"].notna()
        & clinvar_lookup["CLINVAR_REF_KEY"].ne("")
        & clinvar_lookup["CLINVAR_ALT_KEY"].ne("")
    ].drop(columns=["rs_id", "RS# (dbSNP)"], errors="ignore")

    start_matches = keyed_vcf.merge(
        coordinate_lookup,
        left_on=["CHROM_KEY", "POS_START", "REF_KEY", "ALT_KEY"],
        right_on=["CLINVAR_CHROM_KEY", "CLINVAR_POS", "CLINVAR_REF_KEY", "CLINVAR_ALT_KEY"],
        how="left",
    )
    anchor_matches = keyed_vcf.merge(
        coordinate_lookup,
        left_on=["CHROM_KEY", "POS_ANCHOR", "REF_KEY", "ALT_KEY"],
        right_on=["CLINVAR_CHROM_KEY", "CLINVAR_POS", "CLINVAR_REF_KEY", "CLINVAR_ALT_KEY"],
        how="left",
    )

    return pd.concat([rsid_matches, start_matches, anchor_matches], ignore_index=True, sort=False)


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


def split_terms(value) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if text in {"", "-", "."}:
        return []
    return [
        part.strip()
        for part in re.split(r"[,|;]", text)
        if part.strip() and part.strip() not in {"-", "."}
    ]


def normalize_term(value: str) -> str:
    return re.sub(r"\s+", "_", str(value).strip().lower())


def aggregate_unique(series: pd.Series):
    values = []
    seen = set()
    for value in series.dropna():
        text = str(value).strip()
        if not text or text in {"-", "."}:
            continue
        for part in split_terms(text):
            if part not in seen:
                seen.add(part)
                values.append(part)
    if values:
        return "|".join(values)
    return np.nan


def normalize_caller_value(value) -> str:
    text = str(value).strip() if not pd.isna(value) else ""
    lowered = text.lower()
    if lowered in {"hc", "hc_only", "concordant"}:
        return "HC"
    if lowered in {"fb", "fb_only"}:
        return "FB"
    terms = {term.lower() for term in split_terms(text)}
    if {"hc", "hc_only", "concordant"} & terms:
        return "HC"
    if {"fb", "fb_only"} & terms:
        return "FB"
    return ""


def aggregate_caller_support(series: pd.Series):
    values = [normalize_caller_value(value) for value in series.dropna()]
    if "HC" in values:
        return "HC"
    if "FB" in values:
        return "FB"
    return np.nan


def canonical_rank(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.upper().eq("YES").astype(int)


PATHOGENIC_CLIN_SIG = {
    "pathogenic",
    "pathogenic/established_risk_allele",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
    "pathogenic/likely pathogenic",
    "pathogenic_low_penetrance",
    "pathogenic/low_penetrance",
}

HIGH_IMPACT_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
}

MODERATE_IMPACT_CONSEQUENCES = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "coding_sequence_variant",
    "regulatory_region_ablation",
}


def is_pathogenic_clinvar(df: pd.DataFrame) -> pd.Series:
    def has_pathogenic_term(row) -> bool:
        terms = split_terms(row.get("CLIN_SIG", ""))
        terms.extend(split_terms(row.get("ClinicalSignificance", "")))
        return any(normalize_term(term) in PATHOGENIC_CLIN_SIG for term in terms)

    return df.apply(has_pathogenic_term, axis=1)


def is_rare_high_impact_candidate(
    df: pd.DataFrame,
    require_protein_coding: bool = True,
    reportable_max_af: float = 0.001,
) -> pd.Series:
    impact = df.get("IMPACT", pd.Series("", index=df.index)).fillna("").astype(str).str.upper()
    consequence = df.get("Consequence", pd.Series("", index=df.index)).fillna("").astype(str)
    has_high_impact = impact.apply(lambda value: "HIGH" in re.split(r"[,&|;]", value))
    has_high_consequence = consequence.apply(
        lambda value: any(term in HIGH_IMPACT_CONSEQUENCES for term in re.split(r"[,&|;]", value))
    )
    max_af = pd.to_numeric(df.get("MAX_AF", pd.Series(np.nan, index=df.index)), errors="coerce")
    biotype = df.get("BIOTYPE", pd.Series("", index=df.index)).fillna("").astype(str)
    rare = max_af.isna() | (max_af <= reportable_max_af)
    if require_protein_coding:
        biotype_ok = biotype.eq("protein_coding")
    else:
        biotype_ok = pd.Series(True, index=df.index)
    return rare & biotype_ok & (has_high_impact | has_high_consequence)


def is_rare_moderate_candidate(
    df: pd.DataFrame,
    reportable_max_af: float = 0.001,
    require_protein_coding: bool = True,
) -> pd.Series:
    impact = df.get("IMPACT", pd.Series("", index=df.index)).fillna("").astype(str).str.upper()
    consequence = df.get("Consequence", pd.Series("", index=df.index)).fillna("").astype(str)
    max_af = pd.to_numeric(df.get("MAX_AF", pd.Series(np.nan, index=df.index)), errors="coerce")
    biotype = df.get("BIOTYPE", pd.Series("", index=df.index)).fillna("").astype(str)
    rare = max_af.isna() | (max_af <= reportable_max_af)
    has_moderate_impact = impact.apply(lambda value: "MODERATE" in re.split(r"[,&|;]", value))
    has_moderate_consequence = consequence.apply(
        lambda value: any(term in MODERATE_IMPACT_CONSEQUENCES for term in re.split(r"[,&|;]", value))
    )
    if require_protein_coding:
        biotype_ok = biotype.eq("protein_coding")
    else:
        biotype_ok = pd.Series(True, index=df.index)
    return rare & biotype_ok & (has_moderate_impact | has_moderate_consequence)


def add_report_category(
    df: pd.DataFrame,
    require_protein_coding_for_rare_high_impact: bool = True,
    include_rare_moderate_candidates: bool = False,
    require_protein_coding_for_rare_moderate: bool = True,
    reportable_max_af: float = 0.001,
) -> pd.DataFrame:
    out = df.copy()
    pathogenic = is_pathogenic_clinvar(out)
    high_impact = is_rare_high_impact_candidate(
        out,
        require_protein_coding=require_protein_coding_for_rare_high_impact,
        reportable_max_af=reportable_max_af,
    )
    moderate = (
        is_rare_moderate_candidate(
            out,
            reportable_max_af=reportable_max_af,
            require_protein_coding=require_protein_coding_for_rare_moderate,
        )
        if include_rare_moderate_candidates
        else pd.Series(False, index=out.index)
    )
    out["ReportCategory"] = np.select(
        [pathogenic, high_impact, moderate],
        ["clinvar_pathogenic", "rare_high_impact_candidate", "rare_moderate_candidate"],
        default="candidate",
    )
    return out


def vcf_variant_keys(path: str) -> set[tuple[str, int, str, str]]:
    if not os.path.exists(path):
        return set()
    variants = read_vcf_no_meta(path)
    keys = set()
    for _, row in variants.iterrows():
        chrom = str(row.get("CHROM", "")).strip()
        pos = pd.to_numeric(row.get("POS"), errors="coerce")
        ref = str(row.get("REF", "")).strip()
        alts = str(row.get("ALT", "")).split(",")
        if not chrom or pd.isna(pos) or not ref:
            continue
        for alt in alts:
            alt = alt.strip()
            if alt:
                keys.add((chrom, int(pos), ref, alt))
    return keys


def add_caller_support(
    df: pd.DataFrame, hc_filtered_path: str, fb_filtered_path: str
) -> pd.DataFrame:
    out = df.copy()

    def support_from_vcf(row, hc_keys, fb_keys) -> str:
        chrom = str(row.get("Location", "")).split(":")[0]
        ref = str(row.get("REF_ALLELE", "")).strip()
        alt = str(row.get("Allele", "")).strip()
        positions = [
            pd.to_numeric(row.get("POS_START"), errors="coerce"),
            pd.to_numeric(row.get("POS_ANCHOR"), errors="coerce"),
        ]
        candidate_keys = {
            (chrom, int(pos), ref, alt)
            for pos in positions
            if pd.notna(pos) and ref and alt
        }
        in_hc = any(key in hc_keys for key in candidate_keys)
        in_fb = any(key in fb_keys for key in candidate_keys)
        if in_hc:
            return "HC"
        if in_fb:
            return "FB"
        raise ValueError(
            "CallerSupport could not be resolved for "
            f"{row.get('Location')} {row.get('REF_ALLELE')}>{row.get('Allele')}. "
            "Regenerate merged VEP files with the current merge_vep.py so each row is tagged as HC or FB."
        )

    if "CallerSupport" in out.columns:
        normalized = out["CallerSupport"].apply(normalize_caller_value)
        needs_fallback = normalized.eq("")
        out["CallerSupport"] = normalized
        if needs_fallback.any():
            hc_keys = vcf_variant_keys(hc_filtered_path)
            fb_keys = vcf_variant_keys(fb_filtered_path)
            out = add_pos_candidates(out, location_col="Location")
            out.loc[needs_fallback, "CallerSupport"] = out.loc[needs_fallback].apply(
                lambda row: support_from_vcf(row, hc_keys, fb_keys), axis=1
            )
    else:
        hc_keys = vcf_variant_keys(hc_filtered_path)
        fb_keys = vcf_variant_keys(fb_filtered_path)
        out = add_pos_candidates(out, location_col="Location")
        out["CallerSupport"] = out.apply(
            lambda row: support_from_vcf(row, hc_keys, fb_keys), axis=1
        )
    return out


def fill_missing_columns(df: pd.DataFrame, columns) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = np.nan
    return out


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

    if "FORMAT" in out.columns:
        fmt_keys = out["FORMAT"].fillna("").astype(str).str.split(":")
        sample_vals = out[sample_col].fillna("").astype(str).str.split(":")
        parsed = []
        for keys, vals in zip(fmt_keys, sample_vals):
            parsed.append({k: v for k, v in zip(keys, vals)})
        parsed_df = pd.DataFrame(parsed)
        for col in ["GT", "DP", "GQ", "AD"]:
            if col not in parsed_df.columns:
                parsed_df[col] = np.nan
        out["Genotype (GT)"] = parsed_df["GT"]
        out["Depth of Coverage (DP)"] = pd.to_numeric(parsed_df["DP"], errors="coerce")
        out["Genotype Quality (GQ)"] = pd.to_numeric(parsed_df["GQ"], errors="coerce")
        ad_alt = parsed_df["AD"].astype(str).str.split(",").str[1]
        out["AD_ALT"] = pd.to_numeric(ad_alt, errors="coerce")
        out["VAF (Sample)"] = out["AD_ALT"] / out["Depth of Coverage (DP)"]
    else:
        out["Genotype (GT)"] = out[sample_col].str.split(":").str[0]
        out["Depth of Coverage (DP)"] = pd.to_numeric(out[sample_col].str.split(":").str[2], errors="coerce")
        out["Genotype Quality (GQ)"] = pd.to_numeric(out[sample_col].str.split(":").str[3], errors="coerce")
        ad_alt = out[sample_col].str.split(":").str[1].str.split(",").str[1]
        out["AD_ALT"] = pd.to_numeric(ad_alt, errors="coerce")
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
    columns = [
        "Location",
        "REF_ALLELE",
        "Allele",
        "Feature",
        "Feature_type",
        "Consequence",
        "IMPACT",
        "CANONICAL",
        "BIOTYPE",
        "Existing_variation",
        "SYMBOL",
        "GeneSymbol",
        "CLIN_SIG",
        "ClinicalSignificance",
        "ReportCategory",
        "CallerSupport",
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
    out = fill_missing_columns(df, columns)
    out["_canonical_rank"] = canonical_rank(out["CANONICAL"])
    out = out.sort_values("_canonical_rank", ascending=False)
    out = out[columns].drop_duplicates().reset_index(drop=True)

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
    require_protein_coding_for_rare_high_impact: bool = True,
    include_rare_moderate_candidates: bool = False,
    require_protein_coding_for_rare_moderate: bool = True,
    reportable_max_af: float = 0.001,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(glob.glob(file_pattern))
    clinvar = pd.read_csv(variant_summary_path, sep="\t", low_memory=False)
    clinvar_lookup = prepare_clinvar_lookup(clinvar)

    for file in files:
        vcf = read_vep_tab(file)
        if "#CHROM" in vcf.columns:
            vcf = vcf.rename(columns={"#CHROM": "CHROM"})

        requires_biotype = require_protein_coding_for_rare_high_impact or (
            include_rare_moderate_candidates and require_protein_coding_for_rare_moderate
        )
        if "BIOTYPE" not in vcf.columns and requires_biotype:
            raise ValueError(
                f"BIOTYPE column missing in {file}; cannot enforce protein-coding rare candidate gate."
            )

        variant_key = [
            col
            for col in [
                "Location",
                "REF_ALLELE",
                "Allele",
                "Feature",
                "Feature_type",
                "Consequence",
                "IMPACT",
                "CANONICAL",
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
            ]
            if col in vcf.columns
        ]

        vcf_lookup = vcf.copy()
        if "Existing_variation" not in vcf_lookup.columns:
            vcf_lookup["Existing_variation"] = ""
        if "CallerSupport" not in vcf_lookup.columns:
            vcf_lookup["CallerSupport"] = ""
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].fillna("").astype(str)
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].str.split(",")
        vcf_lookup = vcf_lookup.explode("Existing_variation")
        vcf_lookup["Existing_variation"] = vcf_lookup["Existing_variation"].fillna("").astype(str).str.strip()
        vcf_lookup["rs_id"] = pd.to_numeric(
            vcf_lookup["Existing_variation"].str.extract(r"^rs(\d+)$")[0],
            errors="coerce",
            downcast="integer",
        )

        clinvar_left = merge_clinvar_annotations(vcf_lookup, clinvar_lookup)

        all_candidates = clinvar_left.groupby(variant_key, dropna=False, as_index=False).agg(
            {
                "Existing_variation": aggregate_unique,
                "CLIN_SIG": aggregate_unique,
                "CallerSupport": aggregate_caller_support,
                "Feature": first_non_null,
                "Feature_type": first_non_null,
                "Consequence": first_non_null,
                "IMPACT": first_non_null,
                "CANONICAL": first_non_null,
                "GeneSymbol": first_non_null,
                "ClinicalSignificance": aggregate_unique,
                "PhenotypeIDS": aggregate_unique,
                "PhenotypeList": aggregate_unique,
                "rs_id": first_non_null,
            }
        )
        all_candidates["GeneSymbol"] = all_candidates["GeneSymbol"].fillna(
            all_candidates.get("SYMBOL")
        )

        fname = os.path.basename(file)
        sample = fname.replace("_merged_vep.vcf", "").replace("_vep.vcf", "")

        partner_path = os.path.join(haplotype_dir, f"{sample}_variants.vcf")
        freebayes_path = os.path.join(freebayes_dir, f"{sample}_freebayes.vcf")
        hc_filtered_path = os.path.join(haplotype_dir, f"{sample}_filtered.vcf.gz")
        fb_filtered_path = os.path.join(freebayes_dir, f"{sample}_freebayes_filtered.vcf.gz")

        all_candidates = add_report_category(
            add_caller_support(
                attach_sample_metrics(all_candidates, partner_path, freebayes_path),
                hc_filtered_path,
                fb_filtered_path,
            ),
            require_protein_coding_for_rare_high_impact=require_protein_coding_for_rare_high_impact,
            include_rare_moderate_candidates=include_rare_moderate_candidates,
            require_protein_coding_for_rare_moderate=require_protein_coding_for_rare_moderate,
            reportable_max_af=reportable_max_af,
        )
        all_out = build_output_table(all_candidates)
        all_out.to_csv(os.path.join(output_dir, f"{sample}_all_candidates.csv"), index=False)

        report_categories = [
            "clinvar_pathogenic",
            "rare_high_impact_candidate",
        ]
        if include_rare_moderate_candidates:
            report_categories.append("rare_moderate_candidate")

        max_af_population = pd.to_numeric(
            all_out.get("MAX_AF (Population)", pd.Series(np.nan, index=all_out.index)),
            errors="coerce",
        )
        rare_or_absent = max_af_population.isna() | (max_af_population <= reportable_max_af)
        reportable = all_out["ReportCategory"].isin(report_categories) & rare_or_absent
        out = all_out[reportable].copy()
        out["_canonical_rank"] = canonical_rank(out["CANONICAL"])
        out = out.sort_values("_canonical_rank", ascending=False)
        out = out.drop_duplicates(
            subset=["Location", "REF_ALLELE", "Allele", "HGVSc", "HGVSp", "ReportCategory"],
            keep="first",
        ).drop(columns=["_canonical_rank"], errors="ignore").reset_index(drop=True)

        out.to_csv(os.path.join(output_dir, f"{sample}_raw_output.csv"), index=False)
        out.to_csv(os.path.join(output_dir, f"{sample}_output_p.csv"), index=False)
