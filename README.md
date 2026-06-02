# NGS Variant Calling and Annotation Pipeline

Snakemake workflow for targeted germline small-variant calling, annotation, and report-oriented prioritization in a GRCh37-based analysis setting.

## Repository Layout

```text
README.md
run.sh
config/
  config.yaml
workflow/
  Snakefile
  rules/
  scripts/
```

## Workflow Summary

The pipeline performs:

- `fastp` read trimming
- `bwa` alignment
- `GATK` duplicate marking, BQSR, and HaplotypeCaller
- `FreeBayes` secondary calling
- `bcftools` normalization and filtering
- `Ensembl VEP` annotation
- coverage, pseudogene-region flagging, and final report generation

## Requirements

### Software

- `Snakemake` compatible with `--use-singularity`
- `Singularity` or `Apptainer` for rule-level container execution
- `bash`
- `python3`
- `pandas`
- `numpy`
- `samtools`
- `bcftools`
- `bedtools`
- `zip`

The main analysis tools are executed through pinned container images defined in `config/config.yaml`. Host installation of `bwa`, `samtools`, `gatk`, `bcftools`, `freebayes`, or `vep` is not required for the core Snakemake workflow.

### Host-side reporting dependencies

The report-generation stage is now tracked by Snakemake, but the current release invokes the repository post-processing scripts directly from the host environment. For that stage, the host must provide:

- `python3`
- `pandas`
- `numpy`
- `bash`
- `samtools`
- `bcftools`
- `bedtools`
- `zip`

## Scope

The current workflow is intended for:

- targeted-panel germline analysis
- short-read paired-end Illumina data
- single-sample SNV and short-indel calling
- GRCh37-based reference resources

The current workflow does not support:

- structural variants
- copy-number variants
- repeat expansions
- mitochondrial variant analysis
- joint genotyping
- trio-aware or population-scale calling

Phenotype-related information in the current implementation is report-oriented. The workflow does not perform automated patient-phenotype matching or phenotype-driven ranking.

## Data Requirements

The workflow will not run without GRCh37-matched reference and annotation resources.

Required inputs and resources:

- paired-end FASTQ files in `input_dir`
- reference FASTA in `references.ref`
- BWA index generated from the same FASTA in `references.bwa_index`
- target BED file in `references.bed_file`
- known-sites VCF for BQSR in `references.known_sites`
- local VEP cache directory in `references.vep_cache`
- pseudogene BED file in `references.pseudogene_bed`
- ClinVar `variant_summary.txt` in `postprocess.variant_summary`

ClinVar variant summary reference:

- download from `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`
- use `variant_summary.txt.gz` from ClinVar
- decompress it to `variant_summary.txt`
- set the path in `config/config.yaml` under `postprocess.variant_summary`

Example:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
```

Reference FASTA prerequisites:

- FASTA file
- FASTA index (`.fai`)
- sequence dictionary (`.dict`)
- BWA index sidecar files

If these are not already available, create them with:

```bash
samtools faidx GRCh37.p13.genome.fa
gatk CreateSequenceDictionary -R GRCh37.p13.genome.fa
bwa index GRCh37.p13.genome.fa
```

Expected VEP cache layout:

```text
VEP_data/
  homo_sapiens/
    115_GRCh37/
```

## Configuration

Edit `config/config.yaml`.

Main settings:

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
- `resources.gatk_mem_mb`: Java memory for GATK rules
- `resources.tmpdir`: temporary directory
- `postprocess.variant_summary`: ClinVar variant summary file
- `postprocess.qc_min_30x`: QC threshold for 30X coverage
- `postprocess.positive_control_min_30x`: positive-control threshold
- `postprocess.ntc_max_30x`: NTC threshold
- `containers.*`: rule-level container images

FASTQs must follow:

```text
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
```

## Running

Run the full workflow:

```bash
CORES=16 bash run.sh
```

If your input, reference, or cache paths live outside the project directory, bind them into the container runtime.

Typical example:

```bash
SINGULARITY_ARGS="--bind /home --bind /mnt" CORES=16 bash run.sh
```

If all required paths share a common parent, binding that parent is usually sufficient.

Common cases that require bind mounts:

- FASTQs in `/home/...` or `/mnt/...`
- reference FASTA and BWA index outside the repository
- VEP cache outside the repository
- output directory outside the repository

Direct Snakemake execution:

```bash
snakemake \
  --use-singularity \
  --singularity-args "--bind /home --bind /mnt" \
  --cores 16 \
  --keep-going \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml
```

## Modularity and Customization

The workflow is organized as modular Snakemake rules with configuration-driven paths, thresholds, resource settings, and container definitions. Users can adapt reference resources, filtering thresholds, execution settings, and downstream reporting behavior to their own validated assay requirements.

To change common settings:

- Edit `config/config.yaml` under `filter_thresholds` to modify DP and GQ cutoffs.
- Edit `config/config.yaml` under `gatk_hard_filters` to change GATK-style hard-filter thresholds for SNPs and indels, including `QD`, `FS`, `SOR`, `MQ`, `MQRankSum`, and `ReadPosRankSum` where applicable.
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
- To modify final report logic, edit `workflow/scripts/filter.py`, `workflow/scripts/run_filter.py`, or `workflow/scripts/True_positive_filter.py`.

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

## Reproducibility

Core workflow dependencies are pinned through rule-level container images declared in `config/config.yaml`. The workflow uses Docker-compatible image URIs, including BioContainers images for `fastp`, `bwa`, `samtools`, `GATK4`, `bcftools`, `FreeBayes`, and `Ensembl VEP`. Execution settings that vary by system, such as GATK memory and temporary-directory location, are configured through `config/config.yaml` rather than hard-coded paths.

## Troubleshooting

### Container bind issues

If tools inside the container cannot see the reference, FASTQ, or cache paths, add bind mounts through `SINGULARITY_ARGS`.

Example:

```bash
SINGULARITY_ARGS="--bind /home --bind /mnt" CORES=16 bash run.sh
```

### FASTQs not detected

Check both:

- the files are located under `input_dir`
- the filenames follow `<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`

### VEP cache not found

Confirm that:

- `references.vep_cache` points to the cache directory
- the cache path is visible inside the container via `SINGULARITY_ARGS`

### Known-sites or reference path errors

Confirm that:

- `references.ref`, `references.bwa_index`, and `references.known_sites` are correct in `config/config.yaml`
- those paths are bind-mounted if they are outside the project directory

## License and Third-Party Components

The repository code is released under the MIT License. Third-party tools, container images, and reference resources used by the workflow retain their own licenses and terms of use, and users are responsible for complying with those requirements.
