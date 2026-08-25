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
- `duplicate_handling.remove_all_duplicates`: controls whether `MarkDuplicatesSpark` only marks duplicates (`false`, default) or removes duplicates (`true`)
- `read_groups.enabled`: enable per-sample read-group metadata from `read_groups.metadata_file`
- `read_groups.metadata_file`: YAML file containing per-sample read-group fields
- `gatk_hard_filters.snps.*`: configurable SNP hard-filter thresholds
- `gatk_hard_filters.indels.*`: configurable indel hard-filter thresholds
- `reporting.reportable_max_af`: population-frequency threshold for report-facing candidate retention, including ClinVar pathogenic/likely pathogenic, rare high-impact, and optional rare MODERATE categories (`0.001` by default; missing MAX_AF is retained)
- `reporting.require_protein_coding_for_rare_high_impact`: require `BIOTYPE=protein_coding` for rare high-impact candidate reporting (`true` by default)
- `reporting.include_rare_moderate_candidates`: optionally promote rare `MODERATE` VEP-impact candidates into report-facing outputs (`false` by default)
- `reporting.require_protein_coding_for_rare_moderate`: require `BIOTYPE=protein_coding` for rare `MODERATE` candidate promotion (`true` by default)
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


## Duplicate Handling and Read Groups

Duplicate handling is configurable:

```yaml
duplicate_handling:
  remove_all_duplicates: false
```

With the default `remove_all_duplicates: false`, duplicate reads are retained and marked rather than removed. This is generally preferable for germline capture workflows unless duplicate removal has been validated for the assay context.

Per-sample read-group metadata can be enabled:

```yaml
read_groups:
  enabled: true
  metadata_file: "config/read_groups.yaml"
```

Edit `config/read_groups.yaml` so sample names match the FASTQ-derived sample names. Supported fields are `ID`, `SM`, `LB`, `PL`, `PU`, `CN`, and `DT`. If read-group metadata is disabled, the workflow uses fallback values equivalent to `ID=<sample>`, `SM=<sample>`, `LB=lib1`, `PL=ILLUMINA`, and `PU=unit1`.

## Rare Candidate Reporting

Report-facing retention is controlled by population frequency for all prioritized categories. ClinVar pathogenic/likely pathogenic records, rare high-impact candidates, and optional rare MODERATE candidates are retained in report-facing outputs only when `MAX_AF` is missing or `<= reporting.reportable_max_af`.

By default:

```yaml
reporting:
  reportable_max_af: 0.001
  require_protein_coding_for_rare_high_impact: true
  include_rare_moderate_candidates: false
  require_protein_coding_for_rare_moderate: true
```

Report-facing candidates are retained only when `MAX_AF` is missing or `<= reportable_max_af`. Rare high-impact candidates additionally require `BIOTYPE=protein_coding` when the high-impact protein-coding gate is enabled, and VEP reports `IMPACT=HIGH` or one of the configured high-impact consequence terms. These candidates do not require a dbSNP rsID or ClinVar annotation. Rare `MODERATE` candidates can be promoted separately when `include_rare_moderate_candidates` is enabled; by default, this promotion also requires `BIOTYPE=protein_coding` unless `require_protein_coding_for_rare_moderate` is set to `false`.

Rare protein-coding `MODERATE` candidates, such as missense and in-frame indel records, are preserved in `<sample>_all_candidates.csv` by default but are not promoted into report-facing outputs unless this option is enabled:

```yaml
reporting:
  include_rare_moderate_candidates: true
  require_protein_coding_for_rare_moderate: true
```

Set `require_protein_coding_for_rare_high_impact` to `false` only if rare high-impact candidate surfacing should not be restricted by VEP biotype:

```yaml
reporting:
  require_protein_coding_for_rare_high_impact: false
```

The frequency threshold is configurable and should be adjusted for the validated assay context, ancestry context, and local review policy. It is not an ACMG/AMP pathogenicity classification rule.


## Report Term Sets

The report-category logic uses explicit term sets in `workflow/scripts/filter.py`. Users can extend these sets to match their local review policy or future annotation vocabulary changes.

Edit these variables in `workflow/scripts/filter.py`:

```python
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
```

After changing these sets, rerun the report-layer regression tests:

```bash
python3 -m pytest -q tests/test_report_regression.py
```

Any local expansion of these term sets should be documented and validated for the assay context before use in report-facing review.

## Modularity and Customization

The workflow is organized as modular Snakemake rules with configuration-driven paths, thresholds, resource settings, and container definitions. Users can adapt reference resources, filtering thresholds, execution settings, and downstream reporting behavior to their own validated assay requirements.

To change common settings:

- Edit `config/config.yaml` under `filter_thresholds` to modify DP and GQ cutoffs.
- Edit `config/config.yaml` under `gatk_hard_filters` to change GATK-style hard-filter thresholds for SNPs and indels, including `QD`, `FS`, `SOR`, `MQ`, `MQRankSum`, and `ReadPosRankSum` where applicable.
- Edit `config/config.yaml` under `reporting.reportable_max_af` to change the report-facing population-frequency threshold.
- Edit `config/config.yaml` under `reporting.require_protein_coding_for_rare_high_impact` to control whether rare high-impact candidate reporting is restricted to `protein_coding` biotypes.
- Edit `config/config.yaml` under `reporting.include_rare_moderate_candidates` to optionally promote rare `MODERATE` candidates into report-facing outputs.
- Edit `config/config.yaml` under `reporting.require_protein_coding_for_rare_moderate` to control whether rare `MODERATE` candidate promotion is restricted to `protein_coding` biotypes.
- Edit `config/config.yaml` under `resources` to change Java memory and temporary-directory settings.
- Edit `config/config.yaml` under `references` to change reference FASTA, known-sites VCF, BED file, pseudogene BED, and VEP cache paths.
- Edit `config/config.yaml` under `containers` to change or pin alternative container images.

To add or adjust tool-specific parameters:

- Modify the relevant rule in `workflow/rules/` for alignment, calling, filtering, or annotation parameters.
- Modify the relevant script in `workflow/scripts/` for post-processing and report-generation behavior.

Examples:

- To change the bcftools filtering thresholds, edit `workflow/rules/calling.smk`.
- To change the configured GATK-style SNP/indel hard filters without editing rule code, update `gatk_hard_filters` in `config/config.yaml`.
- To change VEP flags or annotation fields, edit `workflow/rules/annotation.smk`.
- To broaden rare high-impact candidate reporting beyond protein-coding biotypes, set `reporting.require_protein_coding_for_rare_high_impact: false` in `config/config.yaml`.
- To promote rare `MODERATE` candidates into report-facing outputs, set `reporting.include_rare_moderate_candidates: true` in `config/config.yaml`. By default this requires `BIOTYPE=protein_coding`; set `reporting.require_protein_coding_for_rare_moderate: false` to broaden MODERATE candidate promotion beyond protein-coding biotypes.
- To modify final report logic beyond the exposed configuration options, edit `workflow/scripts/filter.py`, `workflow/scripts/run_filter.py`, or `workflow/scripts/True_positive_filter.py`.
- To add ClinVar clinical-significance terms or VEP consequence terms used for report categorization, update `PATHOGENIC_CLIN_SIG`, `HIGH_IMPACT_CONSEQUENCES`, or `MODERATE_IMPACT_CONSEQUENCES` in `workflow/scripts/filter.py`, then rerun the regression tests.
