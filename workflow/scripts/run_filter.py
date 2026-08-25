import argparse
from pathlib import Path

import filter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--variant-dir", required=True)
    parser.add_argument("--variant-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--require-protein-coding-for-rare-high-impact",
        default="true",
        help="Require protein_coding BIOTYPE for rare high-impact candidate reporting. Accepts true/false, case-insensitive.",
    )
    parser.add_argument(
        "--include-rare-moderate-candidates",
        default="false",
        help="Include rare MODERATE-impact candidates in report-facing outputs. Accepts true/false, case-insensitive.",
    )
    parser.add_argument(
        "--require-protein-coding-for-rare-moderate",
        default="true",
        help="Require protein_coding BIOTYPE for rare MODERATE-impact candidate reporting. Accepts true/false, case-insensitive.",
    )
    parser.add_argument(
        "--reportable-max-af",
        type=float,
        default=0.001,
        help="Maximum population allele frequency for report-facing variant retention.",
    )
    args = parser.parse_args()

    annotation_dir = Path(args.annotation_dir)
    filter.process_vcf(
        file_pattern=str(annotation_dir / "*_merged_vep.vcf"),
        variant_summary_path=args.variant_summary,
        output_dir=args.output_dir,
        haplotype_dir=args.variant_dir,
        freebayes_dir=args.variant_dir,
        require_protein_coding_for_rare_high_impact=(
            str(args.require_protein_coding_for_rare_high_impact).lower() == "true"
        ),
        include_rare_moderate_candidates=(
            str(args.include_rare_moderate_candidates).lower() == "true"
        ),
        require_protein_coding_for_rare_moderate=(
            str(args.require_protein_coding_for_rare_moderate).lower() == "true"
        ),
        reportable_max_af=args.reportable_max_af,
    )


if __name__ == "__main__":
    main()
