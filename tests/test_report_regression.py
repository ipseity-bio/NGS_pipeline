from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import filter  # noqa: E402


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_report_category_regression_for_af_and_biotype_gates() -> None:
    df = pd.DataFrame(
        {
            "CLIN_SIG": ["pathogenic,drug_response", "", "", "", ""],
            "ClinicalSignificance": ["", "", "", "", "Pathogenic/Likely pathogenic"],
            "MAX_AF": [0.01, 0.0001, None, 0.0001, 0.0001],
            "IMPACT": ["MODERATE", "HIGH", "MODERATE", "MODERATE", "LOW"],
            "Consequence": [
                "missense_variant",
                "frameshift_variant",
                "missense_variant",
                "missense_variant",
                "synonymous_variant",
            ],
            "BIOTYPE": ["protein_coding", "protein_coding", "protein_coding", "lncRNA", "protein_coding"],
        }
    )

    strict = filter.add_report_category(
        df,
        include_rare_moderate_candidates=True,
        require_protein_coding_for_rare_moderate=True,
        reportable_max_af=0.001,
    )
    assert list(strict["ReportCategory"]) == [
        "clinvar_pathogenic",  # Category assignment preserves ClinVar evidence; report retention gates AF later.
        "rare_high_impact_candidate",
        "rare_moderate_candidate",
        "candidate",
        "clinvar_pathogenic",
    ]

    broad = filter.add_report_category(
        df,
        include_rare_moderate_candidates=True,
        require_protein_coding_for_rare_moderate=False,
        reportable_max_af=0.001,
    )
    assert broad.loc[3, "ReportCategory"] == "rare_moderate_candidate"

    no_moderate = filter.add_report_category(
        df,
        include_rare_moderate_candidates=False,
        require_protein_coding_for_rare_moderate=False,
        reportable_max_af=0.001,
    )
    assert no_moderate.loc[2, "ReportCategory"] == "candidate"


def test_aggregation_and_caller_support_regression() -> None:
    assert filter.aggregate_unique(pd.Series(["Pathogenic,drug_response", "Pathogenic|Likely benign", "-"])) == (
        "Pathogenic|drug_response|Likely benign"
    )
    assert filter.aggregate_caller_support(pd.Series(["FB_only", "HC"])) == "HC"
    assert filter.aggregate_caller_support(pd.Series(["FB_only", "FB"])) == "FB"
    assert filter.normalize_caller_value("concordant") == "HC"


def test_process_vcf_report_outputs_are_candidate_preserving_and_reportable_af_gated(tmp_path: Path) -> None:
    annotation_dir = tmp_path / "annotation"
    variant_dir = tmp_path / "variants"
    output_dir = tmp_path / "reports"
    annotation_dir.mkdir()
    variant_dir.mkdir()

    vep_header = [
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
        "CLIN_SIG",
        "SPDI",
        "HGVSg",
        "HGVSc",
        "HGVSp",
        "HGVS_OFFSET",
        "SIFT",
        "PolyPhen",
        "MAX_AF",
        "ZYG",
        "CallerSupport",
    ]
    vep_rows = [
        ["chr1:100", "A", "T", "NM_000001.1", "Transcript", "missense_variant", "MODERATE", "YES", "protein_coding", "rs1001", "COMMONCLIN", "", "", "chr1:g.100A>T", "NM_000001.1:c.1A>T", "NP_000001.1:p.Lys1Asn", "-", "-", "-", "0.01", "HET", "HC"],
        ["chr1:200", "C", "G", "NM_000002.1", "Transcript", "missense_variant", "MODERATE", "YES", "protein_coding", "rs1002", "RARECLIN", "", "", "chr1:g.200C>G", "NM_000002.1:c.2C>G", "NP_000002.1:p.Ala2Gly", "-", "-", "-", "0.0001", "HET", "HC"],
        ["chr1:300", "G", "A", "NM_000003.1", "Transcript", "frameshift_variant", "HIGH", "YES", "protein_coding", "", "NOVELHIGH", "", "", "chr1:g.300G>A", "NM_000003.1:c.3del", "NP_000003.1:p.Gly3fsTer9", "-", "-", "-", "", "HET", "HC"],
        ["chr1:400", "T", "C", "NM_000004.1", "Transcript", "missense_variant", "MODERATE", "YES", "protein_coding", "", "MODPC", "", "", "chr1:g.400T>C", "NM_000004.1:c.4T>C", "NP_000004.1:p.Val4Ala", "-", "-", "-", "0.0001", "HET", "HC"],
        ["chr1:500", "G", "C", "NR_000005.1", "Transcript", "regulatory_region_ablation", "MODERATE", "YES", "lncRNA", "", "MODRNA", "", "", "chr1:g.500G>C", "NR_000005.1:n.5G>C", "-", "-", "-", "-", "0.0001", "HET", "HC"],
        ["chr1:600", "A", "C", "NM_000006.1", "Transcript", "synonymous_variant", "LOW", "YES", "protein_coding", "", "NORSIDCLIN", "", "", "chr1:g.600A>C", "NM_000006.1:c.6A>C", "NP_000006.1:p.Gly2=", "-", "-", "-", "0.0001", "HET", "HC"],
    ]
    write_tsv(annotation_dir / "S1_merged_vep.vcf", ["#" + vep_header[0], *vep_header[1:]], vep_rows)

    write_tsv(
        variant_dir / "S1_variants.vcf",
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"],
        [
            ["chr1", "100", ".", "A", "T", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:15,15:30:99"],
            ["chr1", "200", ".", "C", "G", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:20,20:40:99"],
            ["chr1", "300", ".", "G", "A", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:25,25:50:99"],
            ["chr1", "400", ".", "T", "C", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:30,30:60:99"],
            ["chr1", "500", ".", "G", "C", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:35,35:70:99"],
            ["chr1", "600", ".", "A", "C", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:40,40:80:99"],
        ],
    )

    clinvar = pd.DataFrame(
        {
            "RS# (dbSNP)": [1001, 1002, -1],
            "Assembly": ["GRCh37", "GRCh37", "GRCh37"],
            "GeneSymbol": ["COMMONCLIN", "RARECLIN", "NORSIDCLIN"],
            "ClinicalSignificance": ["Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic"],
            "PhenotypeIDS": ["MedGen:C1", "MedGen:C2", "MedGen:C3"],
            "PhenotypeList": ["condition common", "condition rare", "condition no rsid"],
            "Chromosome": ["1", "1", "1"],
            "PositionVCF": [100, 200, 600],
            "ReferenceAlleleVCF": ["A", "C", "A"],
            "AlternateAlleleVCF": ["T", "G", "C"],
        }
    )
    clinvar_path = tmp_path / "variant_summary.txt"
    clinvar.to_csv(clinvar_path, sep="\t", index=False)

    filter.process_vcf(
        file_pattern=str(annotation_dir / "*_merged_vep.vcf"),
        variant_summary_path=str(clinvar_path),
        output_dir=str(output_dir),
        haplotype_dir=str(variant_dir),
        freebayes_dir=str(variant_dir),
        require_protein_coding_for_rare_high_impact=True,
        include_rare_moderate_candidates=True,
        require_protein_coding_for_rare_moderate=True,
        reportable_max_af=0.001,
    )

    all_candidates = pd.read_csv(output_dir / "S1_all_candidates.csv")
    output_p = pd.read_csv(output_dir / "S1_output_p.csv")

    assert set(all_candidates["SYMBOL"]) >= {"COMMONCLIN", "RARECLIN", "NOVELHIGH", "MODPC", "MODRNA", "NORSIDCLIN"}
    assert set(output_p["SYMBOL"]) == {"RARECLIN", "NOVELHIGH", "MODPC", "NORSIDCLIN"}
    assert "COMMONCLIN" not in set(output_p["SYMBOL"])  # ClinVar P/LP is now reportable_max_af gated.
    assert "MODRNA" not in set(output_p["SYMBOL"])  # MODERATE protein-coding gate is enabled.
    assert output_p.set_index("SYMBOL").loc["NOVELHIGH", "ReportCategory"] == "rare_high_impact_candidate"
    assert output_p.set_index("SYMBOL").loc["MODPC", "ReportCategory"] == "rare_moderate_candidate"
    assert output_p.set_index("SYMBOL").loc["NORSIDCLIN", "ReportCategory"] == "clinvar_pathogenic"
    assert output_p.set_index("SYMBOL").loc["NORSIDCLIN", "ClinicalSignificance (ClinVar)"] == "Pathogenic"


def test_true_positive_filter_reportable_af_gate(tmp_path: Path) -> None:
    rows = []
    base = {
        "Location": "chr1:100",
        "REF_ALLELE": "A",
        "Allele": "T",
        "Feature": "NM_000001.1",
        "Feature_type": "Transcript",
        "Consequence": "missense_variant",
        "IMPACT": "MODERATE",
        "CANONICAL": "YES",
        "BIOTYPE": "protein_coding",
        "Existing_variation": "rs1",
        "SYMBOL": "GENE",
        "GeneSymbol": "GENE",
        "CLIN_SIG": "-",
        "ClinicalSignificance (ClinVar)": "-",
        "ReportCategory": "candidate",
        "CallerSupport": "HC",
        "SPDI": "",
        "HGVSg": "chr1:g.100A>T",
        "HGVSc": "NM_000001.1:c.1A>T",
        "HGVSp": "NP_000001.1:p.Lys1Asn",
        "HGVS_OFFSET": "-",
        "SIFT": "-",
        "PolyPhen": "-",
        "VAF (Sample)": 0.5,
        "MAX_AF (Population)": 0.0001,
        "Genotype (GT)": "0/1",
        "ZYG": "HET",
        "Depth of Coverage (DP)": 50,
        "Genotype Quality (GQ)": 99,
        "PhenotypeIDS": "MedGen:C1",
        "PhenotypeList": "not provided|Charcot-Marie-Tooth disease|type I|not specified",
    }
    for gene, category, clin_sig, clinvar, max_af in [
        ("COMMONCLIN", "clinvar_pathogenic", "pathogenic", "Pathogenic", 0.01),
        ("RARECLIN", "clinvar_pathogenic", "pathogenic", "Pathogenic", 0.0001),
        ("NOVELHIGH", "rare_high_impact_candidate", "-", "-", None),
        ("COMMONHIGH", "rare_high_impact_candidate", "-", "-", 0.1),
        ("MODPC", "rare_moderate_candidate", "-", "-", 0.0001),
    ]:
        row = base.copy()
        row.update(
            {
                "SYMBOL": gene,
                "GeneSymbol": gene,
                "ReportCategory": category,
                "CLIN_SIG": clin_sig,
                "ClinicalSignificance (ClinVar)": clinvar,
                "MAX_AF (Population)": "" if max_af is None else max_af,
                "HGVSc": f"NM_000001.1:c.{len(rows)+1}A>T",
                "HGVSp": f"NP_000001.1:p.Lys{len(rows)+1}Asn",
                "HGVSg": f"chr1:g.{100 + len(rows)}A>T",
                "Location": f"chr1:{100 + len(rows)}",
            }
        )
        rows.append(row)

    input_dir = tmp_path / "reports"
    output_dir = tmp_path / "out"
    input_dir.mkdir()
    write_csv(input_dir / "S1_raw_output.csv", rows)
    coverage = pd.DataFrame(
        {
            "Sample": ["S1"],
            "Mean depth of coverage": [100],
            "Percentage of bases covered at 30X": [99.5],
        }
    )
    coverage_path = tmp_path / "coverage_summary.csv"
    coverage.to_csv(coverage_path, index=False)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "True_positive_filter.py"),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--coverage-summary",
            str(coverage_path),
            "--reportable-max-af",
            "0.001",
        ],
        check=True,
    )

    filtered = pd.read_csv(output_dir / "S1_raw_output_filtered_tp.csv")
    assert set(filtered["GeneSymbol"]) == {"RARECLIN", "NOVELHIGH", "MODPC"}
    assert "COMMONCLIN" not in set(filtered["GeneSymbol"])
    assert "COMMONHIGH" not in set(filtered["GeneSymbol"])
    report = pd.read_csv(output_dir / "S1_raw_output_filtered_tp_report.csv")
    assert set(report["Phenotype"]) == {"Charcot-Marie-Tooth disease|type I"}



def test_clinvar_blank_or_negative_rsid_only_matches_by_coordinate() -> None:
    vcf = pd.DataFrame(
        [
            {"Location": "chr1:100", "REF_ALLELE": "A", "Allele": "T", "Existing_variation": "", "rs_id": pd.NA},
            {"Location": "chr1:200", "REF_ALLELE": "C", "Allele": "G", "Existing_variation": "rs123", "rs_id": 123},
        ]
    )
    clinvar = pd.DataFrame(
        {
            "RS# (dbSNP)": ["-1", "", "123"],
            "Assembly": ["GRCh37", "GRCh37", "GRCh37"],
            "GeneSymbol": ["NO_RSID_COORD", "BLANK_RSID_OTHER", "RSID_MATCH"],
            "ClinicalSignificance": ["Pathogenic", "Pathogenic", "Pathogenic"],
            "PhenotypeIDS": ["P1", "P2", "P3"],
            "PhenotypeList": ["D1", "D2", "D3"],
            "Chromosome": ["1", "1", "1"],
            "PositionVCF": [100, 999, 200],
            "ReferenceAlleleVCF": ["A", "G", "C"],
            "AlternateAlleleVCF": ["T", "A", "G"],
        }
    )

    merged = filter.merge_clinvar_annotations(vcf, filter.prepare_clinvar_lookup(clinvar))
    assert set(merged["GeneSymbol"].dropna()) == {"NO_RSID_COORD", "RSID_MATCH"}


def test_process_vcf_fails_when_required_biotype_column_is_missing(tmp_path: Path) -> None:
    annotation_dir = tmp_path / "annotation"
    variant_dir = tmp_path / "variants"
    output_dir = tmp_path / "reports"
    annotation_dir.mkdir()
    variant_dir.mkdir()

    write_tsv(
        annotation_dir / "S1_merged_vep.vcf",
        ["#Location", "REF_ALLELE", "Allele", "Feature", "Consequence", "IMPACT", "Existing_variation", "SYMBOL", "MAX_AF"],
        [["chr1:100", "A", "T", "NM_000001.1", "frameshift_variant", "HIGH", "", "GENE", "0.0001"]],
    )
    write_tsv(
        variant_dir / "S1_variants.vcf",
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"],
        [["chr1", "100", ".", "A", "T", ".", "PASS", ".", "GT:AD:DP:GQ", "0/1:10,20:30:99"]],
    )
    clinvar_path = tmp_path / "variant_summary.txt"
    pd.DataFrame({"RS# (dbSNP)": [], "Assembly": []}).to_csv(clinvar_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="BIOTYPE column missing"):
        filter.process_vcf(
            file_pattern=str(annotation_dir / "*_merged_vep.vcf"),
            variant_summary_path=str(clinvar_path),
            output_dir=str(output_dir),
            haplotype_dir=str(variant_dir),
            freebayes_dir=str(variant_dir),
            require_protein_coding_for_rare_high_impact=True,
        )


def test_attach_sample_metrics_parses_format_keys_independent_of_order(tmp_path: Path) -> None:
    partner = tmp_path / "S1_variants.vcf"
    write_tsv(
        partner,
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"],
        [["chr1", "100", ".", "A", "T", ".", "PASS", ".", "DP:GT:GQ:AD", "55:0/1:88:20,35"]],
    )
    df = pd.DataFrame(
        [{"Location": "chr1:100", "REF_ALLELE": "A", "Allele": "T", "Feature": "NM_000001.1"}]
    )

    out = filter.attach_sample_metrics(df, str(partner), str(tmp_path / "missing_freebayes.vcf"))
    assert out.loc[0, "Genotype (GT)"] == "0/1"
    assert out.loc[0, "Depth of Coverage (DP)"] == 55
    assert out.loc[0, "Genotype Quality (GQ)"] == 88
    assert round(out.loc[0, "VAF (Sample)"], 6) == round(35 / 55, 6)


def test_true_positive_filter_prefers_canonical_dedup_and_handles_missing_coverage(tmp_path: Path) -> None:
    input_dir = tmp_path / "reports"
    output_dir = tmp_path / "out"
    input_dir.mkdir()
    base = {
        "Location": "chr1:100",
        "REF_ALLELE": "A",
        "Allele": "T",
        "Feature": "NM_000001.1",
        "Feature_type": "Transcript",
        "Consequence": "frameshift_variant",
        "IMPACT": "HIGH",
        "BIOTYPE": "protein_coding",
        "Existing_variation": "",
        "SYMBOL": "GENE",
        "GeneSymbol": "GENE",
        "CLIN_SIG": "-",
        "ClinicalSignificance (ClinVar)": "-",
        "ReportCategory": "rare_high_impact_candidate",
        "CallerSupport": "HC",
        "SPDI": "",
        "HGVSg": "chr1:g.100A>T",
        "HGVSc": "NM_000001.1:c.1A>T",
        "HGVSp": "NP_000001.1:p.Lys1Asn",
        "HGVS_OFFSET": "-",
        "SIFT": "-",
        "PolyPhen": "-",
        "VAF (Sample)": 0.5,
        "MAX_AF (Population)": 0.0001,
        "Genotype (GT)": "0/1",
        "ZYG": "HET",
        "Depth of Coverage (DP)": 50,
        "Genotype Quality (GQ)": 99,
        "PhenotypeIDS": "MedGen:C1",
        "PhenotypeList": "not specified|condition",
    }
    rows = []
    for canonical in ["NO", "YES"]:
        row = base.copy()
        row["CANONICAL"] = canonical
        row["Feature"] = "NM_CANONICAL.1" if canonical == "YES" else "NM_NONCANONICAL.1"
        rows.append(row)
    write_csv(input_dir / "SAMPLE_raw_output.csv", rows)
    coverage_path = tmp_path / "coverage_summary.csv"
    pd.DataFrame(
        {"Sample": ["OTHER"], "Mean depth of coverage": [100], "Percentage of bases covered at 30X": [99.5]}
    ).to_csv(coverage_path, index=False)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "True_positive_filter.py"),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--coverage-summary",
            str(coverage_path),
            "--reportable-max-af",
            "0.001",
        ],
        check=True,
    )

    filtered = pd.read_csv(output_dir / "SAMPLE_raw_output_filtered_tp.csv")
    report = pd.read_csv(output_dir / "SAMPLE_raw_output_filtered_tp_report.csv")
    assert len(filtered) == 1
    assert filtered.loc[0, "CANONICAL"] == "YES"
    assert report.loc[0, "Avg Read Depth"] == "-"
    assert report.loc[0, "Panel Coverage"] == "-"
