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
    )


if __name__ == "__main__":
    main()
