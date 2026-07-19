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
            -t {BED_PATH} \
            --genotype-qualities \
            {input.bam} \
            > {output.vcf} 2> {log}
        """


rule filter_strand_bias:
    input:
        vcf=str(VARIANT_DIR / "{sample}_variants.vcf")
    output:
        vcf=str(VARIANT_DIR / "{sample}_no_strand_bias.vcf")
    params:
        snp_qd_min=HARD_FILTERS["snps"]["QD_min"],
        snp_fs_max=HARD_FILTERS["snps"]["FS_max"],
        snp_sor_max=HARD_FILTERS["snps"]["SOR_max"],
        snp_mq_min=HARD_FILTERS["snps"]["MQ_min"],
        snp_mqranksum_min=HARD_FILTERS["snps"]["MQRankSum_min"],
        snp_readpos_min=HARD_FILTERS["snps"]["ReadPosRankSum_min"],
        indel_qd_min=HARD_FILTERS["indels"]["QD_min"],
        indel_fs_max=HARD_FILTERS["indels"]["FS_max"],
        indel_sor_max=HARD_FILTERS["indels"]["SOR_max"],
        indel_readpos_min=HARD_FILTERS["indels"]["ReadPosRankSum_min"],
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
        | bcftools filter -i '(
            TYPE="snp" &&
            (INFO/QD="." || INFO/QD >= {params.snp_qd_min}) &&
            (INFO/FS="." || INFO/FS <= {params.snp_fs_max}) &&
            (INFO/SOR="." || INFO/SOR <= {params.snp_sor_max}) &&
            (INFO/MQ="." || INFO/MQ >= {params.snp_mq_min}) &&
            (INFO/MQRankSum="." || INFO/MQRankSum >= {params.snp_mqranksum_min}) &&
            (INFO/ReadPosRankSum="." || INFO/ReadPosRankSum >= {params.snp_readpos_min})
          ) || (
            TYPE!="snp" &&
            (INFO/QD="." || INFO/QD >= {params.indel_qd_min}) &&
            (INFO/FS="." || INFO/FS <= {params.indel_fs_max}) &&
            (INFO/SOR="." || INFO/SOR <= {params.indel_sor_max}) &&
            (INFO/ReadPosRankSum="." || INFO/ReadPosRankSum >= {params.indel_readpos_min})
          )' -O v -o {output.vcf} - >> {log} 2>&1
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
        set -o pipefail
        bcftools norm -f {REF} --check-ref w {input.vcf} 2>> {log} \
        | bcftools norm -f {REF} -d none 2>> {log} \
        | bcftools filter -i '(FMT/DP >= {params.DP}) & (FMT/GQ >= {params.GQ}) & (GT!="0/0")' 2>> {log} \
        | bcftools sort -Oz -o {output.vcf} - >> {log} 2>&1
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
