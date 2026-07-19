# Reference Resources

The workflow will not run without GRCh37-matched reference and annotation resources.

## Required Inputs and Resources

- paired-end FASTQ files in `input_dir`
- reference FASTA in `references.ref`
- BWA index generated from the same FASTA in `references.bwa_index`
- target BED file in `references.bed_file`
- known-sites VCF for BQSR in `references.known_sites`
- local VEP cache directory in `references.vep_cache`
- pseudogene BED file in `references.pseudogene_bed`
- ClinVar `variant_summary.txt` in `postprocess.variant_summary`

## ClinVar Variant Summary Reference

ClinVar variant summary is used during filtering and reporting.

Download source:

`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`

Use `variant_summary.txt.gz`, decompress it to `variant_summary.txt`, and set the path in `config/config.yaml` under `postprocess.variant_summary`.

Example:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
```

## Reference FASTA Prerequisites

Required files:

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

## VEP Cache Layout

Expected VEP cache layout:

```text
VEP_data/
  homo_sapiens/
    115_GRCh37/
```

Set this path in `config/config.yaml` under `references.vep_cache`.

## Reference Build Note

The workflow is evaluated in a GRCh37-based configuration because the assay context, CAP proficiency materials, Coriell reference benchmarking set, and downstream interpretation resources used in this study were aligned to GRCh37. Future GRCh38 support should be treated as a separately validated configuration.
