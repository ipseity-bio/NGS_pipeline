# NGS Variant Calling & Annotation Pipeline

*A complete Snakemake-based workflow for GRCh37*

<p align="center">
  <img src="https://img.shields.io/badge/Workflow-Snakemake-blue?style=for-the-badge">
  <img src="https://img.shields.io/badge/Annotation-VEP%20(Docker)-orange?style=for-the-badge">
  <img src="https://img.shields.io/badge/Genome-GRCh37-green?style=for-the-badge">
</p>

---

## 📘 Overview

This repository provides a full NGS analysis workflow for:

* FASTQ QC trimming
* Alignment using **bwa**
* Duplicate marking & BQSR (GATK)
* Variant calling (GATK + FreeBayes)
* Filtering & normalization
* Annotation using **Ensembl VEP via Docker**
* Coverage analysis
* Pseudogene filtering
* Automated true-positive variant reporting

The primary workflow is implemented in **Snakemake**, with downstream processing in **Python** and **Bash**.

---

## ⚙️ Requirements

### **Software**

| Tool                               | Purpose                                    |
| ---------------------------------- | ------------------------------------------ |
| Snakemake                          | Workflow management                        |
| fastp                              | Read trimming                              |
| bwa                                | Alignment                                  |
| samtools                           | Sorting, indexing, depth                   |
| GATK4                              | MarkDuplicatesSpark, BQSR, HaplotypeCaller |
| freebayes                          | Alternative variant caller                 |
| bcftools                           | Variant filtering/normalization            |
| bedtools                           | Pseudogene region intersections            |
| Docker                             | Required for VEP annotation                |
| Python 3                           | Required for post-processing scripts       |
| Python packages: `pandas`, `numpy` | CSV/VCF processing                         |

---

## 🧬 Data Requirements (Must Provide)

The pipeline **will not run without these files**.  
All coordinates must match **GRCh37**.

### 1) Reference Genome (GRCh37)
You need:
- `GRCh37.fa` (FASTA)
- `GRCh37.fa.fai` (FASTA index)
- `GRCh37.dict` (sequence dictionary)

Create indexes if missing:

```bash
samtools faidx GRCh37.fa
gatk CreateSequenceDictionary -R GRCh37.fa
bwa index GRCh37.fa
````

`config.yaml` keys:

```yaml
ref: "./GRCh37.p13.genome.fa"
bwa_index: "./GRCh37.p13.genome.fa"
```

---

### 2) Known Sites VCF (for BQSR)

Any GRCh37-known-sites VCF, e.g. Mills/1000G indels.

`config.yaml`:

```yaml
known_sites: "./Mills_and_1000G_gold_standard.indels.b37.vcf"
```

---

### 4) VEP Cache

Download and store the Ensembl VEP cache **locally**.

Example cache directory:

```
VEP_data/
└── homo_sapiens/
    └── 115_GRCh37/ 
```

`config.yaml`:

```yaml
docker_data: "./VEP_data"
```

> The pipeline mounts this directory into the VEP docker container.

---

### 5) ClinVar / Variant Summary Reference

Used during filtering + reporting.

Required:

* `variant_summary.txt` (download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ and update the path in run_filter.py)

---

## 🔧 Configuration

Modify `config.yaml` to match your reference paths and thresholds:

```yaml
bed_file: "regions.bed"

filter_thresholds:
  DP: 30
  GQ: 99

known_sites: "./Mills_and_1000G_gold_standard.indels.b37.vcf"
ref: "./GRCh37.p13.genome.fa"
docker_data: "./VEP_data"
bwa_index: "./GRCh37.p13.genome.fa"
```

---

### FASTQ Naming Convention

FASTQs must be in the form:

```
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
```

Example:

```
Patient01_R1_001.fastq.gz → sample = Patient01
```

---

## ▶️ Running the Pipeline

### **🔹 Option A — One-Click Execution**

This runs the entire Snakemake pipeline + post-processing:

```bash
bash run.sh
```

---

### **🔹 Option B — Manual Execution**

#### **1. Run Snakemake**

```bash
snakemake \
  --cores 32 \
  --keep-going \
  --snakefile snakemake.smk \
  --configfile config.yaml
```

---

#### **2. Run Post-Processing**

```bash
bash run_post_process.sh
```

This will run:

1. `run_filter.py`
2. `depth.py`
3. `pseudogene_check.sh`
4. `coverage_table.py`
5. `coverage_summary.py`
6. `True_positive_filter.py` 

---

## 📤 Output Summary

### **After Snakemake**

| Output                     | Description                  |
| -------------------------- | ---------------------------- |
| `<sample>_recal.bam`       | Final BAM after BQSR         |
| `<sample>_variants.vcf`    | GATK HaplotypeCaller output  |
| `<sample>_freebayes.vcf`   | FreeBayes output             |
| `<sample>_filtered.vcf.gz` | Variant-level filtered VCF   |
| `<sample>_vep.vcf`         | GATK variants annotated      |
| `<sample>_fb_vep.vcf`      | FreeBayes variants annotated |
| `<sample>_merged_vep.vcf`  | Combined annotation          |

---

### **After Post-Processing**

| Output                             | Purpose                                       |
| ---------------------------------- | --------------------------------------------- |
| `<sample>_raw_output.csv`          | VEP + ClinVar merged dataset                  |
| `<sample>_output_p.csv`            | Pathogenic/Likely Pathogenic filtered list    |
| `<sample>_all_coverage.txt`        | Per-base depth across regions                 |
| `<folder>_coverage_summary.csv`    | Combined coverage summary                     |
| `<sample>_pseudogene_variants.txt` | Intersections with pseudogene regions         |
| `<sample>_filtered_tp.csv`         | True-positive filtered table (from TP script) |
| `<sample>_filtered_tp_report.csv`  | Final clinical-style TP report                |
| `true_positive.zip`                | Packaged TP files                             |
| `true_positive_report.zip`         | Packaged report files                         |

---

## 🧩 What Each Script Does

### ✔ `run_filter.py` + `filter.py`

* Joins VEP annotations with ClinVar data
* Produces raw & pathogenic-only CSVs

---

### ✔ `depth.py`

Computes per-base depth:

```bash
samtools depth -b regions.bed <sample>_recal.bam
```

---

### ✔ `coverage_table.py`

Aggregates all sample-level depth files → single summary CSV.

---

### ✔ `coverage_summary.py`

Adds QC fields:

* Mean depth
* % bases ≥30×
* QC status

---

### ✔ `pseudogene_check.sh`

Uses bedtools to detect pseudogene-overlapping variants.

---

### ✔ `True_positive_filter.py`

* Enforces population frequency thresholds
* Creates formatted TP report
* Adds panel coverage + mean depth
* Outputs

  * `<sample>_filtered_tp.csv`
  * `<sample>_filtered_tp_report.csv`

---

## ⚠️ Troubleshooting

### ❗ VEP Errors

Ensure your Docker VEP cache is downloaded and the path exists:

```
docker_data: "./VEP_data"
```

---

### ❗ FASTQs Not Detected

Check naming:

```
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
```

---

## 📜 Citation

If you use this pipeline in research or reporting, please cite:

* GATK
* FreeBayes
* VEP
* Snakemake
* fastp
* bwa
* samtools
* bcftools
