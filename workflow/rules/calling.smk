rule gatk_HaplotypeCaller:
    input:
        bam=str(ALIGN_DIR / "{sample}_recal.bam"),
        bai=str(ALIGN_DIR / "{sample}_recal.bam.bai"),
    output:
        vcf=str(VARIANT_DIR / "{sample}_variants.vcf")
    threads: 16
    log:
        str(LOG_DIR / "gatk_HaplotypeCaller_{sample}.log")
    container:
        CONTAINERS["gatk4"]
    shell:
        """
        mkdir -p {VARIANT_DIR} {LOG_DIR}
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
        bam=str(ALIGN_DIR / "{sample}_recal.bam"),
        bai=str(ALIGN_DIR / "{sample}_recal.bam.bai"),
    output:
        vcf=str(VARIANT_DIR / "{sample}_freebayes.vcf")
    threads: 8
    log:
        str(LOG_DIR / "freebayes_{sample}.log")
    container:
        CONTAINERS["freebayes"]
    shell:
        """
        mkdir -p {VARIANT_DIR} {LOG_DIR}
        freebayes \
            -f {REF} \
            --genotype-qualities \
            {input.bam} \
            > {output.vcf} 2> {log}
        """


rule filter_strand_bias:
    input:
        vcf=str(VARIANT_DIR / "{sample}_variants.vcf")
    output:
        vcf=str(VARIANT_DIR / "{sample}_no_strand_bias.vcf")
    log:
        str(LOG_DIR / "filter_strand_bias_{sample}.log")
    container:
        CONTAINERS["bcftools"]
    shell:
        r"""
        mkdir -p {LOG_DIR}
        set -o pipefail
        bcftools norm -f {REF} --check-ref w {input.vcf} 2>> {log} \
        | bcftools norm -f {REF} -d none 2>> {log} \
        | bcftools filter -i 'INFO/SOR > 0.5 & INFO/SOR < 2.0' -O v -o {output.vcf} - >> {log} 2>&1
        """


rule bcftools_filter:
    input:
        vcf=str(VARIANT_DIR / "{sample}_no_strand_bias.vcf")
    output:
        vcf=str(VARIANT_DIR / "{sample}_filtered.vcf.gz"),
        csi=str(VARIANT_DIR / "{sample}_filtered.vcf.gz.csi")
    params:
        DP=FILTERS["DP"],
        GQ=FILTERS["GQ"]
    threads: 8
    log:
        str(LOG_DIR / "bcftools_filter_{sample}.log")
    container:
        CONTAINERS["bcftools"]
    shell:
        """
        mkdir -p {LOG_DIR}
        bcftools filter -i '(FORMAT/DP > {params.DP} & FORMAT/GQ >= {params.GQ}) & \
        (GT="0/1" || GT="1/0" || GT="0/2" || GT="2/0" || GT="1/2" || GT="2/1" || GT="1/1" || GT="2/2")' \
        {input.vcf} -Oz -o {output.vcf} > {log} 2>&1
        bcftools index -f {output.vcf} >> {log} 2>&1
        """


rule bcftools_filter_fb:
    input:
        vcf=str(VARIANT_DIR / "{sample}_freebayes.vcf")
    output:
        vcf=str(VARIANT_DIR / "{sample}_freebayes_filtered.vcf.gz"),
        csi=str(VARIANT_DIR / "{sample}_freebayes_filtered.vcf.gz.csi")
    params:
        DP=FILTERS["DP"],
        GQ=FILTERS["GQ"]
    threads: 8
    log:
        str(LOG_DIR / "bcftools_filter_fb_{sample}.log")
    container:
        CONTAINERS["bcftools"]
    shell:
        r"""
        mkdir -p {LOG_DIR}
        bcftools norm -f {REF} --check-ref w {input.vcf} \
        | bcftools norm -f {REF} -d none \
        | bcftools filter -i '(FMT/DP >= {params.DP}) & (FMT/GQ >= {params.GQ}) & (GT!="0/0")' \
        -Oz -o {output.vcf} - > {log} 2>&1
        bcftools index -f {output.vcf} >> {log} 2>&1
        """


rule compare_vcf:
    input:
        fb=str(VARIANT_DIR / "{sample}_freebayes_filtered.vcf.gz"),
        hc=str(VARIANT_DIR / "{sample}_filtered.vcf.gz")
    output:
        unique=str(COMPARE_DIR / "{sample}" / "0000.vcf")
    log:
        str(LOG_DIR / "compare_vcf_{sample}.log")
    threads: 8
    container:
        CONTAINERS["bcftools"]
    shell:
        """
        mkdir -p {COMPARE_DIR} {LOG_DIR}
        bcftools isec -n-1 -c all {input.fb} {input.hc} -p {COMPARE_DIR}/{wildcards.sample} > {log} 2>&1
        """
