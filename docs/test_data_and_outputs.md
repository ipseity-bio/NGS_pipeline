# Test Data and Outputs

## Test Dataset and Reproducibility Data

Representative Coriell reference FASTQ files used in this study are available at Zenodo:

`https://zenodo.org/records/17802399`

These data can be used as example validation inputs or as an installation-verification dataset after configuring the required reference resources.

For a minimal validation-style test, users may run the workflow on:

- `Coriell_NA23721_S2_R1_001.fastq.gz`
- `Coriell_NA23721_S2_R2_001.fastq.gz`

and verify recovery of the expected `MPZ c.418T>A (p.Ser140Thr)` finding in the final report.


## Report Output Tiers

- `reports/<sample>_all_candidates.csv`: comprehensive retained annotated variants after upstream filtering and ClinVar left-join annotation. This file preserves variants without dbSNP rsID or ClinVar matches, including rare `MODERATE`, `LOW`, and `MODIFIER` records.
- `reports/<sample>_raw_output.csv`: report-facing candidate set before final formatting. By default this includes ClinVar pathogenic/likely pathogenic records and rare high-impact candidates when `MAX_AF` is missing or within the configured `reporting.reportable_max_af`; rare moderate candidates are included only when configured, with protein-coding restriction controlled separately.
- `reports/<sample>_output_p.csv`: same report-facing candidate tier used by downstream TP/report formatting.
- `reports/<sample>_raw_output_filtered_tp.csv`: structured final review table generated from `raw_output.csv`.
- `reports/<sample>_raw_output_filtered_tp_report.csv`: formatted review report generated from `raw_output.csv`.

Report-facing files include `ReportCategory` and `CallerSupport`. `ReportCategory` distinguishes `clinvar_pathogenic`, `rare_high_impact_candidate`, optional `rare_moderate_candidate`, and background `candidate` records. `CallerSupport` records whether a variant is observed in the filtered HaplotypeCaller callset only (`HC`) or FreeBayes-only (`FB`); if both callers support the same reported event, `HC` is retained as the primary caller label.

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
