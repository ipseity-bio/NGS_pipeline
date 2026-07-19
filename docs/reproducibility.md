# Reproducibility

Workflow dependencies are pinned through rule-level container images declared in `config/config.yaml`. The workflow uses Docker-compatible image URIs, including BioContainers images for `fastp`, `bwa`, `samtools`, `GATK4`, `bcftools`, `FreeBayes`, and `Ensembl VEP`, plus a dedicated post-processing image for downstream reporting. Execution settings that vary by system, such as GATK memory and temporary-directory location, are configured through `config/config.yaml` rather than hard-coded paths.

The downstream reporting stage uses a dedicated post-processing container:

```yaml
postprocess: "docker://ghcr.io/viswajitmulpuru/ipseity-ngs-postprocess:1.0.1"
```

This image contains Python, pandas, numpy, samtools, bcftools, bedtools, bash/coreutils, and zip.

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
