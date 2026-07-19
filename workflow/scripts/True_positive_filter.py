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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--coverage-summary", required=True)
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
        data["Phenotype"] = data["PhenotypeList"].str.split("|").apply(
            lambda values: next((p for p in values if p not in ["not provided", "not specified"]), None)
            if isinstance(values, list)
            else None
        )
        if "ReportCategory" not in data.columns:
            data["ReportCategory"] = "candidate"
        data["ReportCategory"] = data["ReportCategory"].fillna("candidate")
        data["CLIN_SIG"] = data["CLIN_SIG"].fillna("-")
        data["ClinicalSignificance (ClinVar)"] = data["ClinicalSignificance (ClinVar)"].fillna("-")
        data["Classification"] = "(" + data["CLIN_SIG"] + " ; " + data["ClinicalSignificance (ClinVar)"] + ")"
        high_impact_label = data["ReportCategory"].eq("rare_high_impact_candidate")
        data.loc[high_impact_label, "Classification"] = "Rare high-impact candidate"
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
            | (data["MAX_AF (Population)"] <= 0.001)
        )
        clinvar_reportable = (
            data["CLIN_SIG"].isin(["pathogenic", "pathogenic/likely_pathogenic"])
            | data["ClinicalSignificance (ClinVar)"].isin(["Pathogenic", "Pathogenic/Likely pathogenic"])
            | data["ReportCategory"].eq("clinvar_pathogenic")
        )
        high_impact_reportable = data["ReportCategory"].eq("rare_high_impact_candidate")
        filtered_data = data[(clinvar_reportable | high_impact_reportable) & rare_or_absent].copy()
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
                "Location",
                "Allele State",
                "Allelic Read Depths",
                "Genomic Position",
                "Variant Frequency",
            ]
        ].copy()

        value_md = coverage_df.loc[
            coverage_df["Sample"].str.contains(sample, na=False), "Mean depth of coverage"
        ].values[0]
        report_out["Avg Read Depth"] = f"{value_md}X"

        value_p30 = coverage_df.loc[
            coverage_df["Sample"].str.contains(sample, na=False), "Percentage of bases covered at 30X"
        ].values[0]
        report_out["Panel Coverage"] = f"{value_p30}%"
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
