rule fastp:
    input:
        r1=lambda wc: input_fastq(wc.sample, "R1"),
        r2=lambda wc: input_fastq(wc.sample, "R2"),
    output:
        fo1=str(FASTP_DIR / "{sample}_R1_fp.fastq.gz"),
        fo2=str(FASTP_DIR / "{sample}_R2_fp.fastq.gz"),
        html=str(FASTP_DIR / "{sample}_fastp.html"),
    threads: 8
    log:
        str(LOG_DIR / "fastp_{sample}.log")
    container:
        CONTAINERS["fastp"]
    shell:
        """
        mkdir -p {FASTP_DIR} {LOG_DIR}
        fastp \
            -w {threads} \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.fo1} \
            -O {output.fo2} \
            -h {output.html} \
            > {log} 2>&1
        """


rule bwa_mem:
    input:
        r1=str(FASTP_DIR / "{sample}_R1_fp.fastq.gz"),
        r2=str(FASTP_DIR / "{sample}_R2_fp.fastq.gz"),
    output:
        sam=str(ALIGN_DIR / "{sample}.sam")
    threads: 24
    log:
        str(LOG_DIR / "bwa_mem_{sample}.log")
    params:
        read_group=lambda wc: read_group_string(wc.sample)
    container:
        CONTAINERS["bwa"]
    shell:
        """
        mkdir -p {ALIGN_DIR} {LOG_DIR}
        bwa mem -t {threads} \
            -R "{params.read_group}" \
            {BWA_INDEX} {input.r1} {input.r2} \
            > {output.sam} 2> {log}
        """


rule samtools_sort:
    input:
        sam=str(ALIGN_DIR / "{sample}.sam")
    output:
        bam=str(ALIGN_DIR / "{sample}_grouped.bam")
    threads: 16
    log:
        str(LOG_DIR / "samtools_sort_{sample}.log")
    container:
        CONTAINERS["samtools"]
    shell:
        """
        mkdir -p {ALIGN_DIR} {LOG_DIR}
        samtools sort -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        """


rule samtools_index:
    input:
        bam=str(ALIGN_DIR / "{sample}_grouped.bam")
    output:
        bai=str(ALIGN_DIR / "{sample}_grouped.bam.bai")
    threads: 8
    log:
        str(LOG_DIR / "samtools_index_{sample}.log")
    container:
        CONTAINERS["samtools"]
    shell:
        """
        mkdir -p {LOG_DIR}
        samtools index {input.bam} -@ {threads} > {log} 2>&1
        """


rule gatk_MarkDuplicatesSpark:
    input:
        bam=str(ALIGN_DIR / "{sample}_grouped.bam"),
        bai=str(ALIGN_DIR / "{sample}_grouped.bam.bai")
    output:
        bam=str(ALIGN_DIR / "{sample}_marked.bam")
    threads: 8
    resources:
        mem_mb=lambda wc, attempt: RESOURCES["gatk_mem_mb"],
        tmpdir=lambda wc, attempt: RESOURCES["tmpdir"],
    log:
        str(LOG_DIR / "gatk_mark_duplicates_{sample}.log")
    params:
        remove_all_duplicates=lambda wc: str(
            DEDUP.get("remove_all_duplicates", False)
        ).lower()
    container:
        CONTAINERS["gatk4"]
    shell:
        """
        mkdir -p {ALIGN_DIR} {LOG_DIR}
        gatk --java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={resources.tmpdir}" \
            MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            --remove-all-duplicates {params.remove_all_duplicates} \
            --tmp-dir {resources.tmpdir} \
            > {log} 2>&1
        """


rule samtools_index_marked:
    input:
        bam=str(ALIGN_DIR / "{sample}_marked.bam")
    output:
        bai=str(ALIGN_DIR / "{sample}_marked.bam.bai")
    threads: 8
    log:
        str(LOG_DIR / "samtools_index_marked_{sample}.log")
    container:
        CONTAINERS["samtools"]
    shell:
        """
        mkdir -p {LOG_DIR}
        samtools index {input.bam} -@ {threads} > {log} 2>&1
        """


rule gatk_BaseRecalibratorSpark:
    input:
        bam=str(ALIGN_DIR / "{sample}_marked.bam"),
        bai=str(ALIGN_DIR / "{sample}_marked.bam.bai"),
    output:
        table=str(ALIGN_DIR / "{sample}_recal.table")
    threads: 16
    resources:
        mem_mb=lambda wc, attempt: RESOURCES["gatk_mem_mb"],
        tmpdir=lambda wc, attempt: RESOURCES["tmpdir"],
    log:
        str(LOG_DIR / "gatk_bqsr_{sample}.log")
    container:
        CONTAINERS["gatk4"]
    shell:
        """
        mkdir -p {LOG_DIR}
        gatk --java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={resources.tmpdir}" \
            BaseRecalibratorSpark \
            -I {input.bam} \
            -R {REF} \
            --known-sites {KNOWN_SITES} \
            -O {output.table} \
            --tmp-dir {resources.tmpdir} \
            --spark-master local[{threads}] \
            > {log} 2>&1
        """


rule gatk_ApplyBQSRSpark:
    input:
        bam=str(ALIGN_DIR / "{sample}_marked.bam"),
        bai=str(ALIGN_DIR / "{sample}_marked.bam.bai"),
        table=str(ALIGN_DIR / "{sample}_recal.table"),
    output:
        bam=str(ALIGN_DIR / "{sample}_recal.bam")
    threads: 16
    resources:
        mem_mb=lambda wc, attempt: RESOURCES["gatk_mem_mb"],
        tmpdir=lambda wc, attempt: RESOURCES["tmpdir"],
    log:
        str(LOG_DIR / "gatk_applybqsr_{sample}.log")
    container:
        CONTAINERS["gatk4"]
    shell:
        """
        mkdir -p {ALIGN_DIR} {LOG_DIR}
        gatk --java-options "-Xmx{resources.mem_mb}M -Djava.io.tmpdir={resources.tmpdir}" \
            ApplyBQSRSpark \
            -I {input.bam} \
            -R {REF} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} \
            --tmp-dir {resources.tmpdir} \
            --spark-master local[{threads}] \
            > {log} 2>&1
        """


rule samtools_index_recal:
    input:
        bam=str(ALIGN_DIR / "{sample}_recal.bam")
    output:
        bai=str(ALIGN_DIR / "{sample}_recal.bam.bai")
    threads: 8
    log:
        str(LOG_DIR / "samtools_index_recal_{sample}.log")
    container:
        CONTAINERS["samtools"]
    shell:
        """
        mkdir -p {LOG_DIR}
        samtools index {input.bam} -@ {threads} > {log} 2>&1
        """
