import os
import glob

# ---- Discover samples from FASTQ ----
fastq_files = glob.glob("*.fastq.gz")
sample_names = sorted(set(file.split("_R")[0] for file in fastq_files))

# ---- Resolve host paths from config ----
REF_PATH = config["ref"]
REF = config["ref"]
KNOWN_SITES = config["known_sites"]
BED_PATH = config["bed_file"]
VEP_CACHE = config["docker_data"]
BWA_INDEX = config["bwa_index"]

# ---- Mount points for Docker ----
REF_DIR = os.path.dirname(REF_PATH)
KNOWN_DIR = os.path.dirname(KNOWN_SITES)
BED_DIR = os.path.dirname(BED_PATH)

REF_BASENAME = os.path.basename(REF_PATH)
KNOWN_BASENAME = os.path.basename(KNOWN_SITES)
BED_BASENAME = os.path.basename(BED_PATH)

# ---- Output naming conventions ----
FASTP_OUT_1 = "{sample}_R1_fp.fastq"
FASTP_OUT_2 = "{sample}_R2_fp.fastq"

BAM_GROUP = "{sample}_grouped.bam"
BAM_INDEX = "{sample}_grouped.bam.bai"
BAM_MARKED = "{sample}_marked.bam"

BAM_RECAL = "{sample}_recal.bam"
RECAL_FILE = "{sample}_recal.table"
VCF_OUT = "{sample}_variants.vcf"
VCF_FREEBAYES = "{sample}_freebayes.vcf"
VCF_NO_STRAND_BIAS = "{sample}_no_strand_bias.vcf"
VCF_FILTERED = "{sample}_filtered.vcf.gz"
VEP_OUT = "{sample}_vep.vcf" 

VCF_FREEBAYES_FILTERED = "{sample}_freebayes_filtered.vcf.gz"
VEP_OUT_FB = "{sample}_fb_vep.vcf"

UNIQUE_VCF = "dir_{sample}/0000.vcf"
MERGED_VEP = "{sample}_merged_vep.vcf"



ruleorder: ensembl_vep_hc > ensembl_vep_fb > merge_vep


rule all:
    input:
        expand(MERGED_VEP, sample=sample_names)

rule fastp:
    input:
        r1=lambda wc: f"{wc.sample}_R1_001.fastq.gz",
        r2=lambda wc: f"{wc.sample}_R2_001.fastq.gz"
    output:
        fo1=FASTP_OUT_1,
        fo2=FASTP_OUT_2
    threads: 8
    log:
        "fastp_{sample}.log"
    shell:
        """
        fastp \
            -w {threads} \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.fo1} \
            -O {output.fo2} \
            -h {wildcards.sample}_fastp.html \
            > {log} 2>&1
        """


rule bwa_mem_sort:
    input:
        r1=FASTP_OUT_1,
        r2=FASTP_OUT_2
    output:
        BAM_GROUP
    threads: 24
    log:
        "bwa_mem_sort_{sample}.log"
    shell:
        """
        bwa mem -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1" \
            {BWA_INDEX} {input.r1} {input.r2} 2> {log} \
        | samtools sort -@ {threads} -o {output}
        """


rule samtools_index:
    input:
        BAM_GROUP
    output:
        BAM_INDEX
    threads: 8
    log:
        "samtools_index_{sample}.log"
    shell:
        """
        samtools index {input} -@ {threads} > {log} 2>&1
        """


rule gatk_MarkDuplicatesSpark:
    input:
        BAM_GROUP
    output:
        BAM_MARKED
    threads: 8
    log:
        "gatk_mark_duplicates_{sample}.log"
    shell:
        """
        /root/anaconda3/bin/java -Xmx32g \
            -jar /root/anaconda3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar \
            MarkDuplicatesSpark \
            -I {input} \
            -O {output} \
            --remove-all-duplicates true \
            --tmp-dir /tmp \
            > {log} 2>&1
        samtools index {output} -@ {threads}
        """


rule gatk_BaseRecalibratorSpark:
    input:
        bam=BAM_MARKED,
        ref=REF
    output:
        RECAL_FILE
    threads: 16
    log:
        "gatk_bqsr_{sample}.log"
    shell:
        """
        /root/anaconda3/bin/java -Xmx32g \
            -jar /root/anaconda3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar \
            BaseRecalibratorSpark \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {KNOWN_SITES} \
            -O {output} \
            --tmp-dir /tmp \
            --spark-master local[{threads}] \
            > {log} 2>&1
        """


rule gatk_ApplyBQSRSpark:
    input:
        bam=BAM_MARKED,
        ref=REF,
        table=RECAL_FILE
    output:
        BAM_RECAL
    threads: 16
    log:
        "gatk_applybqsr_{sample}.log"
    shell:
        """
        /root/anaconda3/bin/java -Xmx32g \
            -jar /root/anaconda3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar \
            ApplyBQSRSpark \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output} \
            --tmp-dir /tmp \
            --spark-master local[{threads}] \
            > {log} 2>&1
        """

rule gatk_HaplotypeCaller:
    input:
        bam=BAM_RECAL
    output:
        vcf=VCF_OUT
    threads: 16
    log:
        "gatk_HaplotypeCaller_{sample}.log"
    shell:
        """
        gatk HaplotypeCaller \
            -R {REF} \
            -I {input.bam} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            -L {BED_PATH} \
            > {log} 2>&1
        """


rule freebayes:
    input:
        bam=BAM_RECAL
    output:
        vcf=VCF_FREEBAYES
    threads: 8
    log:
        "freebayes_{sample}.log"
    shell:
        """
        freebayes \
            -f {REF} \
            --genotype-qualities \
            {input.bam} \
            > {output.vcf} 2> {log}
        """

rule filter_strand_bias:
    input:
        VCF_OUT
    output:
        VCF_NO_STRAND_BIAS
    log:
        "filter_strand_bias_{sample}.log"
    shell:
        r"""
        # Left-align + trim (no atomize), then apply an SOR window
        bcftools norm -f {config[ref]} --check-ref w {input} \
        | bcftools norm -f {config[ref]} -d none \
        | bcftools filter -i 'INFO/SOR > 0.5 & INFO/SOR < 2.0' -O v -o {output} -
        """

rule bcftools_filter:
    input:
        vcf=VCF_NO_STRAND_BIAS
    output:
        vcf_filtered=VCF_FILTERED
    params:
        DP=config["filter_thresholds"]["DP"],
        GQ=config["filter_thresholds"]["GQ"]
    threads: 8
    log:
        "bcftools_filter_{sample}.log"
    shell:
        """
        bcftools filter -i '(FORMAT/DP > {params.DP} & FORMAT/GQ >= {params.GQ}) & \
        (GT="0/1" || GT="1/0" || GT="0/2" || GT="2/0" || GT="1/2" || GT="2/1" || GT="1/1" || GT="2/2")' \
        {input.vcf} -Oz -o {output.vcf_filtered}
        bcftools index -f {output.vcf_filtered}
        """

rule bcftools_filter_fb:
    input:
        vcf=VCF_FREEBAYES
    output:
        vcf_filtered=VCF_FREEBAYES_FILTERED
    params:
        DP=config["filter_thresholds"]["DP"],  
        GQ=config["filter_thresholds"]["GQ"]  
    threads: 8
    log:
        "bcftools_filter_fb_{sample}.log"
    shell:
        r"""
        bcftools norm -f {config[ref]} --check-ref w {input.vcf} \
        | bcftools norm -f {config[ref]} -d none \
        | bcftools filter -i '(FMT/DP >= {params.DP}) & (FMT/GQ >= {params.GQ}) & (GT!="0/0")' \
        -Oz -o {output.vcf_filtered} -
        bcftools index -f {output.vcf_filtered}
        """

rule compare_vcf:
    input:
        fb=VCF_FREEBAYES_FILTERED,
        hc=VCF_FILTERED
    output:
        UNIQUE_VCF
    log:
        "compare_vcf_{sample}.log"
    threads: 8
    shell:
        """
        bcftools isec -n-1 -c all {input.fb} {input.hc} -p dir_{wildcards.sample}
        """

rule ensembl_vep_hc:
    input:
        vep=VCF_FILTERED,
        done=UNIQUE_VCF
    output:
        vep_out=VEP_OUT
    log:
        "ensembl_vep_{sample}.log"
    threads: 8
    shell:
        """
        docker run --rm -u $(id -u):$(id -g) \
            -v $(pwd):/work:Z \
            -v {config[docker_data]}:/data:Z \
            -w /work \
            ensemblorg/ensembl-vep \
            vep --fork 50 --offline --sift b --polyphen b --ccds --symbol --canonical \
            --protein --biotype --af --pubmed --uniprot --variant_class --gene_phenotype \
            --mirna --check_existing --hgvs --fasta /data/{REF_BASENAME} \
            --spdi --hgvsg --individual all --allele_number --dir_cache /data \
            --buffer_size 5000 --force_overwrite --assembly GRCh37 \
            --input_file {input.vep} --output_file {output.vep_out} \
            --tab --fields "Uploaded_variation,Location,Allele,Gene,Existing_variation,CLIN_SIG,SIFT,PolyPhen,CCDS,SYMBOL,ENSP,BIOTYPE,AF,MAX_AF,PUBMED,VARIANT_CLASS,GENE_PHENO,SOMATIC,PHENO,HGVSg,HGVSc,HGVSp,HGVS_OFFSET,SPDI,ZYG,IND,ALLELE_NUM,REF_ALLELE" \
            --cache_version 115 --refseq --show_ref_allele --max_af
        """

rule ensembl_vep_fb:
    input:
        vcf=UNIQUE_VCF
    output:
        vep_out=VEP_OUT_FB
    log:
        "ensembl_vep_rsid_{sample}.log"
    threads: 8
    shell:
        """
        docker run --rm -u $(id -u):$(id -g) \
            -v $(pwd):/work:Z \
            -v {config[docker_data]}:/data:Z \
            -w /work \
            ensemblorg/ensembl-vep \
            vep -i {input.vcf} \
            --assembly GRCh37 --fasta /data/{REF_BASENAME} \
            --tab --fields "Uploaded_variation,Location,Allele,Gene,Existing_variation,CLIN_SIG,SIFT,PolyPhen,CCDS,SYMBOL,ENSP,BIOTYPE,AF,MAX_AF,PUBMED,VARIANT_CLASS,GENE_PHENO,SOMATIC,PHENO,HGVSg,HGVSc,HGVSp,HGVS_OFFSET,SPDI,ZYG,IND,ALLELE_NUM,REF_ALLELE" \
            --force_overwrite --dir_cache /data --cache_version 115 \
            --output_file {output.vep_out} --refseq --cache --offline --hgvs --spdi --hgvsg --show_ref_allele --max_af --af --symbol
        """

rule merge_vep:
    input:
        hc=VEP_OUT,
        fb=VEP_OUT_FB
    output:
        merged=MERGED_VEP
    log:
        "merge_vep_{sample}.log"
    threads: 8
    run:
        import pandas as pd
        from pathlib import Path
        import io, re

        def read_vep(path: str) -> pd.DataFrame:
            header_cols = None
            data_rows = []
            with open(path, 'r', encoding='utf-8-sig', errors='replace') as fh:
                for raw in fh:
                    if raw.lstrip().startswith('##'):
                        continue
                    line = raw.rstrip('\n')
                    if header_cols is None and line.lstrip().startswith('#'):
                        hdr = line.lstrip()[1:].strip()
                        header_cols = re.split(r'\t+', hdr)
                        continue
                    if line.strip() == '':
                        continue
                    data_rows.append(line.split('\t'))

            if header_cols is None:
                raise ValueError(f"No single-# VEP header in {path}")

            n = len(header_cols)
            fixed = []
            for r in data_rows:
                if len(r) < n:
                    r = r + [''] * (n - len(r))
                elif len(r) > n:
                    r = r[:n]
                fixed.append(r)

            df = pd.DataFrame(fixed, columns=[(c or '').lstrip('\ufeff').strip() for c in header_cols], dtype=str)
            return df

        df1 = read_vep(input.hc) 
        df2 = read_vep(input.fb)  


        for c in df1.columns:
            if c not in df2.columns:
                df2[c] = ''   
        df2 = df2[df1.columns]

        merged_df = pd.concat([df1, df2], ignore_index=True)


        meta = []
        with open(input.hc, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('##'):
                    meta.append(line)
                elif line.lstrip().startswith('#'):
                    break
        header_line = '#' + '\t'.join(df1.columns) + '\n'

        print('HC col0 name:', repr(df1.columns[0])) 
        print('FB col0 name:', repr(df2.columns[0]))   
        print('First HC value:', df1.iloc[0, 0])       
        print('First FB value:', df2.iloc[0, 0])    

        with open(output.merged, 'w', encoding='utf-8') as out:
            out.writelines(meta)
            out.write(header_line)
            merged_df.to_csv(out, sep='\t', index=False, header=False)

        Path(log[0]).write_text(
            f"Merged {len(df1)} HC + {len(df2)} FB rows = {len(merged_df)} total\n"
        )