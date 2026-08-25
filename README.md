# NGS Variant Calling and Annotation Pipeline

Snakemake workflow for targeted germline small-variant calling, annotation, and report-oriented filtering in a GRCh37-based analysis setting.

## Workflow Summary

The pipeline performs:

- `fastp` read trimming
- `bwa` alignment
- `GATK` duplicate marking, BQSR, and HaplotypeCaller
- `FreeBayes` secondary calling
- `bcftools` normalization and filtering
- `Ensembl VEP` annotation
- coverage analysis, pseudogene-region flagging, and final report generation

## Scope

The workflow is intended for:

- targeted-panel germline analysis
- short-read paired-end Illumina data
- single-sample SNV and short-indel calling
- GRCh37-based reference resources

The workflow does not support structural variants, copy-number variants, repeat expansions, mitochondrial variant analysis, joint genotyping, trio-aware calling, or population-scale calling. Phenotype-related information is report-oriented; the workflow does not perform automated patient-phenotype matching or phenotype-driven ranking.

## Repository Layout

```text
README.md
run.sh
config/
  config.yaml
containers/
  postprocess/
docs/
tests/
workflow/
  Snakefile
  rules/
  scripts/
```

## Requirements

For standard containerized execution, the host system requires:

- `Snakemake` compatible with `--use-singularity`
- `Singularity` or `Apptainer` for rule-level container execution
- `bash`

The workflow tools are executed through pinned container images defined in `config/config.yaml`. Host installation of `bwa`, `samtools`, `gatk`, `bcftools`, `freebayes`, `vep`, `bedtools`, `pandas`, or `numpy` is not required for standard containerized workflow execution.

Duplicate handling and read groups are configurable. By default, `duplicate_handling.remove_all_duplicates: false` marks duplicates without removing them; set it to `true` only if duplicate removal has been validated for the assay context. Set `read_groups.enabled: true` and edit `config/read_groups.yaml` to provide per-sample `ID`, `SM`, `LB`, `PL`, `PU`, `CN`, or `DT` fields.

## Quick Start

1. Place paired-end FASTQs in `input_dir` using this naming convention:

```text
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
```

2. Configure paths and thresholds in `config/config.yaml`.

3. Run the workflow:

```bash
CORES=16 bash run.sh
```

If input data or reference resources are outside the repository directory, bind those paths into the container runtime:

```bash
SINGULARITY_ARGS="--bind /home --bind /mnt" CORES=16 bash run.sh
```

## Main Outputs

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
- `annotation/<sample>_merged_vep.vcf`
- `reports/<sample>_all_candidates.csv`
- `reports/<sample>_raw_output.csv`
- `reports/<sample>_output_p.csv`
- `reports/<sample>_raw_output_filtered_tp.csv`
- `reports/<sample>_raw_output_filtered_tp_report.csv`
- `reports/coverage/`
- `reports/.postprocess.done`

## Report Tiers

The reporting layer uses a candidate-preserving path. `all_candidates.csv` retains all annotated variants after upstream filtering, including variants without dbSNP rsID or ClinVar matches. ClinVar `variant_summary.txt` enrichment uses rsID matching where available and GRCh37 exact CHROM/POS/REF/ALT matching through `Chromosome`, `PositionVCF`, `ReferenceAlleleVCF`, and `AlternateAlleleVCF`, allowing ClinVar records without dbSNP rsID to be annotated when the VCF representation matches. Report-facing outputs include ClinVar pathogenic/likely pathogenic records and rare high-impact candidates by default when `MAX_AF` is missing or `<= reporting.reportable_max_af`. Rare `MODERATE` candidates can be promoted to report-facing outputs with `reporting.include_rare_moderate_candidates: true`; by default, this promotion is restricted to `protein_coding` biotypes through `reporting.require_protein_coding_for_rare_moderate: true`. Report tables include `CallerSupport` labels (`HC` or `FB`; variants seen in both callers are represented as `HC`).


Report-category term sets are defined in `workflow/scripts/filter.py`. Users can add locally validated ClinVar clinical-significance terms or VEP consequence terms by updating `PATHOGENIC_CLIN_SIG`, `HIGH_IMPACT_CONSEQUENCES`, or `MODERATE_IMPACT_CONSEQUENCES`, then rerunning the regression tests.

Report-layer regression tests and expected outputs are documented in [Reproducibility](docs/reproducibility.md).

## Documentation

Detailed setup and usage information is organized under `docs/`:

- [Reference resources](docs/reference_resources.md): required GRCh37 resources, ClinVar, VEP cache, pseudogene BED, and reference indexing.
- [Configuration](docs/configuration.md): `config/config.yaml` keys, thresholds, containers, and customization points.
- [Running with containers](docs/running_with_containers.md): `run.sh`, direct Snakemake execution, bind mounts, and post-processing-only execution.
- [Test data and outputs](docs/test_data_and_outputs.md): Zenodo test data, minimal Coriell test, output folders, and representative files.
- [Reproducibility](docs/reproducibility.md): pinned container images, post-processing container, licensing, and GRCh37 build note.
- [Troubleshooting](docs/troubleshooting.md): FASTQ detection, bind paths, VEP cache, and reference path issues.

## Modularity

The workflow is organized as modular Snakemake rules with configuration-driven paths, thresholds, resource settings, and container definitions. Users can adapt reference resources, filtering thresholds, execution settings, and downstream reporting behavior to their own validated assay requirements.

## Test Dataset

Representative Coriell reference FASTQ files used in this study are available at Zenodo:

`https://zenodo.org/records/17802399`

See [Test data and outputs](docs/test_data_and_outputs.md) for included files and a minimal validation-style run.

## Reference Build Note

The workflow is evaluated in a GRCh37-based configuration because the assay context, CAP proficiency materials, Coriell reference benchmarking set, and downstream interpretation resources used in this study were aligned to GRCh37. Future GRCh38 support should be treated as a separately validated configuration.
