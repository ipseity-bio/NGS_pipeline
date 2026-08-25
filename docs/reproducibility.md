# Reproducibility

Workflow dependencies are pinned through rule-level container images declared in `config/config.yaml`. The workflow uses Docker-compatible image URIs, including BioContainers images for `fastp`, `bwa`, `samtools`, `GATK4`, `bcftools`, `FreeBayes`, and `Ensembl VEP`, plus a dedicated post-processing image for downstream reporting. Execution settings that vary by system, such as GATK memory and temporary-directory location, are configured through `config/config.yaml` rather than hard-coded paths.

## Post-Processing Container

The reporting layer is part of the Snakemake DAG and runs in the `ipseity-ngs-postprocess` container. This avoids requiring users to manage Python packages, htslib-based utilities, or bedtools manually on the host system.

The downstream reporting stage uses a dedicated post-processing container:

```yaml
postprocess: "docker://ghcr.io/viswajitmulpuru/ipseity-ngs-postprocess:1.0.1"
```

This image contains Python, pandas, numpy, samtools, bcftools, bedtools, bash/coreutils, and zip. A separate container is used because downstream report generation combines Python table processing with htslib-based and interval-processing command-line tools in the same execution environment. The single-tool containers used for earlier workflow steps do not provide this complete combined runtime.

## Container Images

Example image URIs currently used:

- `docker://quay.io/biocontainers/fastp:0.24.0--heae3180_1`
- `docker://quay.io/biocontainers/bwa:0.7.19--h577a1d6_1`
- `docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0`
- `docker://quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0`
- `docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1`
- `docker://quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0`
- `docker://quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1`
- `docker://ghcr.io/viswajitmulpuru/ipseity-ngs-postprocess:1.0.1`

## Licensing

The repository source code is distributed under the MIT License. Third-party tools and containers retain their own licenses. Users are responsible for reviewing and complying with the licenses of external tools, reference resources, databases, and container images used in their deployment.


## Regression tests

The repository includes lightweight report-layer regression tests that do not require FASTQ files, reference genomes, VEP cache files, or Snakemake execution. These tests use synthetic VEP, VCF, and ClinVar-like fixtures to verify candidate-preserving reporting behavior, `reportable_max_af` filtering, ClinVar pathogenic/likely pathogenic AF gating, no-rsID ClinVar coordinate enrichment, rare high-impact and rare MODERATE candidate promotion, protein-coding gates, and caller-support aggregation. The tests also cover ClinVar `RS# (dbSNP) = -1` and blank-rsID behavior, prevention of false rsID matching, fail-closed behavior when `BIOTYPE` is missing under protein-coding gates, VCF sample-metric parsing by `FORMAT` keys rather than fixed field order, preservation of multi-part phenotype terms in TP reports, canonical-preferred TP deduplication, and safe handling of missing coverage-summary rows.

Run from the repository root:

```bash
python3 -m pytest -q tests/test_report_regression.py
```

Expected result:

```text
8 passed
```
