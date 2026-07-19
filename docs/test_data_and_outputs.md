# Test Data and Outputs

## Test Dataset and Reproducibility Data

Representative Coriell reference FASTQ files used in this study are available at Zenodo:

`https://zenodo.org/records/17802399`

These data can be used as example validation inputs or as an installation-verification dataset after configuring the required reference resources.

For a minimal validation-style test, users may run the workflow on:

- `Coriell_NA23721_S2_R1_001.fastq.gz`
- `Coriell_NA23721_S2_R2_001.fastq.gz`

and verify recovery of the expected `MPZ c.418T>A (p.Ser140Thr)` finding in the final report.

## Outputs

All outputs are written below `output_dir`, including:

- `fastp/`
- `alignment/`
- `variants/`
- `annotation/`
- `compare/`
- `logs/`
- `reports/`

Representative outputs include:

- `alignment/<sample>_recal.bam`
- `variants/<sample>_filtered.vcf.gz`
- `annotation/<sample>_vep.vcf`
- `annotation/<sample>_fb_vep.vcf`
- `annotation/<sample>_merged_vep.vcf`
- `reports/<sample>_all_candidates.csv`
- `reports/<sample>_raw_output.csv`
- `reports/<sample>_output_p.csv`
- `reports/*_tp.csv`
- `reports/*_tp_report.csv`
- `reports/coverage/`
- `reports/.postprocess.done`
