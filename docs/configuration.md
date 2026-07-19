# Configuration

Edit `config/config.yaml` before running the workflow.

## Main Settings

- `input_dir`: directory containing FASTQ files
- `output_dir`: root directory for all pipeline outputs
- `references.ref`: reference FASTA
- `references.bwa_index`: BWA index prefix
- `references.known_sites`: known-sites VCF for BQSR
- `references.bed_file`: target BED file
- `references.vep_cache`: local VEP cache directory
- `references.pseudogene_bed`: pseudogene BED file
- `filter_thresholds.DP`: depth threshold
- `filter_thresholds.GQ`: genotype-quality threshold
- `gatk_hard_filters.snps.*`: configurable SNP hard-filter thresholds
- `gatk_hard_filters.indels.*`: configurable indel hard-filter thresholds
- `reporting.require_protein_coding_for_rare_high_impact`: require `BIOTYPE=protein_coding` for rare high-impact candidate reporting (`true` by default)
- `resources.gatk_mem_mb`: Java memory for GATK rules
- `resources.tmpdir`: temporary directory
- `postprocess.variant_summary`: ClinVar variant summary file
- `postprocess.qc_min_30x`: QC threshold for 30X coverage
- `postprocess.positive_control_min_30x`: positive-control threshold
- `postprocess.ntc_max_30x`: NTC threshold
- `containers.*`: rule-level container images

## FASTQ Naming Convention

FASTQs must follow:

```text
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
```

## Rare High-Impact Candidate Reporting

Rare high-impact candidate reporting is controlled by VEP impact/consequence fields, population frequency, and optional biotype restriction.

By default:

```yaml
reporting:
  require_protein_coding_for_rare_high_impact: true
```

Set this to `false` to surface rare high-impact candidates regardless of VEP biotype:

```yaml
reporting:
  require_protein_coding_for_rare_high_impact: false
```

## Modularity and Customization

The workflow is organized as modular Snakemake rules with configuration-driven paths, thresholds, resource settings, and container definitions. Users can adapt reference resources, filtering thresholds, execution settings, and downstream reporting behavior to their own validated assay requirements.

To change common settings:

- Edit `config/config.yaml` under `filter_thresholds` to modify DP and GQ cutoffs.
- Edit `config/config.yaml` under `gatk_hard_filters` to change GATK-style hard-filter thresholds for SNPs and indels, including `QD`, `FS`, `SOR`, `MQ`, `MQRankSum`, and `ReadPosRankSum` where applicable.
- Edit `config/config.yaml` under `reporting.require_protein_coding_for_rare_high_impact` to control whether rare high-impact candidate reporting is restricted to `protein_coding` biotypes.
- Edit `config/config.yaml` under `resources` to change Java memory and temporary-directory settings.
- Edit `config/config.yaml` under `references` to change reference FASTA, known-sites VCF, BED file, pseudogene BED, and VEP cache paths.
- Edit `config/config.yaml` under `containers` to change or pin alternative container images.

To add or adjust tool-specific parameters:

- Modify the relevant rule in `workflow/rules/` for alignment, calling, filtering, or annotation parameters.
- Modify the relevant script in `workflow/scripts/` for post-processing and report-generation behavior.

Examples:

- To change the bcftools filtering thresholds, edit `workflow/rules/calling.smk`.
- To change the configured GATK-style SNP/indel hard filters, update `gatk_hard_filters` in `config/config.yaml`.
- To change VEP flags or annotation fields, edit `workflow/rules/annotation.smk`.
- To broaden rare high-impact candidate reporting beyond protein-coding biotypes, set `reporting.require_protein_coding_for_rare_high_impact: false` in `config/config.yaml`.
- To modify final report logic beyond the exposed configuration options, edit `workflow/scripts/filter.py`, `workflow/scripts/run_filter.py`, or `workflow/scripts/True_positive_filter.py`.
