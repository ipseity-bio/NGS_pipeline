import argparse
from pathlib import Path

import filter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--variant-dir", required=True)
    parser.add_argument("--variant-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    annotation_dir = Path(args.annotation_dir)
    filter.process_vcf(
        file_pattern=str(annotation_dir / "*_merged_vep.vcf"),
        variant_summary_path=args.variant_summary,
        output_dir=args.output_dir,
        haplotype_dir=args.variant_dir,
        freebayes_dir=args.variant_dir,
    )


if __name__ == "__main__":
    main()
