import argparse
from pathlib import Path

import pandas as pd


def gt_to_zyg(gt):
    if pd.isna(gt):
        return "-"
    gt = str(gt).strip().replace("|", "/")
    if gt in ["./.", ".|."]:
        return "-"
    alleles = gt.split("/")
    if len(alleles) != 2 or "." in alleles:
        return "-"
    a, b = alleles
    if a == "0" and b == "0":
        return "HOM_REF"
    if a == b and a != "0":
        return "HOM_ALT"
    if a != b:
        return "HET"
    return "-"


def format_phenotype_list(value):
    if pd.isna(value):
        return "-"
    generic_terms = {"", "-", ".", "not provided", "not specified"}
    phenotypes = []
    seen = set()
    for part in str(value).split("|"):
        phenotype = part.strip()
        if phenotype.lower() in generic_terms:
            continue
        if phenotype not in seen:
            seen.add(phenotype)
            phenotypes.append(phenotype)
    return "|".join(phenotypes) if phenotypes else "-"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument(
        "--reportable-max-af",
        type=float,
        default=0.001,
        help="Maximum population allele frequency for report-facing variant retention.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    coverage_df = pd.read_csv(args.coverage_summary)

    for input_file in sorted(input_dir.glob("*raw_output.csv")):
        base_name = input_file.stem
        sample = base_name.split("_raw")[0]
        output_file = output_dir / f"{base_name}_filtered_tp.csv"
        report_output_file = output_dir / f"{base_name}_filtered_tp_report.csv"

        data = pd.read_csv(input_file)
        zyg_missing = data["ZYG"].astype(str).str.strip().isin(["-", ".", "nan", "NaN", ""])
        data.loc[zyg_missing, "ZYG"] = data.loc[zyg_missing, "Genotype (GT)"].apply(gt_to_zyg)
        data["MAX_AF (Population)"] = pd.to_numeric(
            data["MAX_AF (Population)"].replace("-", None), errors="coerce"
        )

        data["Avg Read Depth"] = data["Depth of Coverage (DP)"]
        data["Gene"] = data["GeneSymbol"]
        data["Transcript"] = data["HGVSc"].str.extract(r"^(NM_\d+\.\d+)")
        data["Variant"] = data["HGVSc"].str.extract(r"(c\.\S+)")[0] + data["HGVSp"].str.extract(r"(p\.\S+)")[0].radd(" ; ").fillna("")
        data["Inheritance"] = data["ZYG"]
        data["PhenotypeList"] = data["PhenotypeList"].fillna("").astype(str)
        data["Phenotype"] = data["PhenotypeList"].apply(format_phenotype_list)
        if "ReportCategory" not in data.columns:
            data["ReportCategory"] = "candidate"
        data["ReportCategory"] = data["ReportCategory"].fillna("candidate")
        data["CLIN_SIG"] = data["CLIN_SIG"].fillna("-")
        data["ClinicalSignificance (ClinVar)"] = data["ClinicalSignificance (ClinVar)"].fillna("-")
        data["Classification"] = "(" + data["CLIN_SIG"] + " ; " + data["ClinicalSignificance (ClinVar)"] + ")"
        high_impact_label = data["ReportCategory"].eq("rare_high_impact_candidate")
        moderate_label = data["ReportCategory"].eq("rare_moderate_candidate")
        data.loc[high_impact_label, "Classification"] = "Rare high-impact candidate"
        data.loc[moderate_label, "Classification"] = "Rare moderate candidate"
        data["Allele State"] = data["ZYG"]
        data["Allelic Read Depths"] = (
            "Ref(" + data["REF_ALLELE"] + "), Alt(" + data["Allele"] + ") VAF:" + (data["VAF (Sample)"] * 100).astype(str) + "%"
        )
        data["Genomic Position"] = data["HGVSg"]
        data["Variant Frequency"] = data["MAX_AF (Population)"].apply(
            lambda value: "Not identified in a large population studies"
            if pd.isna(value) or value == ""
            else f"{value * 100}% Max frequency observed in Annotated 1000 genomes/ESP/gnomAD."
        )

        rare_or_absent = (
            data["MAX_AF (Population)"].isna()
            | (data["MAX_AF (Population)"] <= args.reportable_max_af)
        )
        clinvar_reportable = (
            data["CLIN_SIG"].isin(["pathogenic", "pathogenic/likely_pathogenic"])
            | data["ClinicalSignificance (ClinVar)"].isin(["Pathogenic", "Pathogenic/Likely pathogenic"])
            | data["ReportCategory"].eq("clinvar_pathogenic")
        )
        rare_candidate_reportable = data["ReportCategory"].isin([
            "rare_high_impact_candidate",
            "rare_moderate_candidate",
        ])
        filtered_data = data[(clinvar_reportable | rare_candidate_reportable) & rare_or_absent].copy()
        if "CANONICAL" not in filtered_data.columns:
            filtered_data["CANONICAL"] = ""
        filtered_data["_canonical_rank"] = (
            filtered_data["CANONICAL"].fillna("").astype(str).str.upper().eq("YES")
        ).astype(int)
        filtered_data = filtered_data.sort_values("_canonical_rank", ascending=False)

        report_out = filtered_data[
            [
                "_canonical_rank",
                "Avg Read Depth",
                "Gene",
                "Transcript",
                "Variant",
                "Inheritance",
                "Phenotype",
                "Classification",
                "Consequence",
                "IMPACT",
                "CallerSupport",
                "Location",
                "Allele State",
                "Allelic Read Depths",
                "Genomic Position",
                "Variant Frequency",
            ]
        ].copy()

        sample_str = str(sample).strip()
        matching_rows = coverage_df[coverage_df["Sample"].astype(str).str.strip() == sample_str]
        if matching_rows.empty:
            matching_rows = coverage_df[coverage_df["Sample"].astype(str).str.contains(sample_str, na=False, regex=False)]

        if not matching_rows.empty and "Mean depth of coverage" in matching_rows.columns:
            value_md = matching_rows["Mean depth of coverage"].values[0]
            report_out["Avg Read Depth"] = f"{value_md}X"
        else:
            report_out["Avg Read Depth"] = "-"

        if not matching_rows.empty and "Percentage of bases covered at 30X" in matching_rows.columns:
            value_p30 = matching_rows["Percentage of bases covered at 30X"].values[0]
            report_out["Panel Coverage"] = f"{value_p30}%"
        else:
            report_out["Panel Coverage"] = "-"

        report_out = report_out[["Panel Coverage"] + [col for col in report_out.columns if col != "Panel Coverage"]]
        tp_dedup_subset = ["REF_ALLELE", "Allele", "HGVSg"]
        report_dedup_subset = ["Allele State", "Allelic Read Depths", "Genomic Position"]

        filtered_data.drop_duplicates(subset=tp_dedup_subset).drop(
            columns=["_canonical_rank"], errors="ignore"
        ).to_csv(output_file, index=False)
        report_out.drop_duplicates(subset=report_dedup_subset).drop(
            columns=["_canonical_rank"], errors="ignore"
        ).to_csv(report_output_file, index=False)


if __name__ == "__main__":
    main()
