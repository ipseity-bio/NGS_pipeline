rule postprocess_reports:
    input:
        merged=expand(str(ANNOTATION_DIR / "{sample}_merged_vep.vcf"), sample=sample_names),
        filtered_vcf=expand(str(VARIANT_DIR / "{sample}_filtered.vcf.gz"), sample=sample_names),
        filtered_vcf_csi=expand(str(VARIANT_DIR / "{sample}_filtered.vcf.gz.csi"), sample=sample_names),
        freebayes_vcf=expand(str(VARIANT_DIR / "{sample}_freebayes.vcf"), sample=sample_names),
        freebayes_filtered_vcf=expand(str(VARIANT_DIR / "{sample}_freebayes_filtered.vcf.gz"), sample=sample_names),
        freebayes_filtered_vcf_csi=expand(str(VARIANT_DIR / "{sample}_freebayes_filtered.vcf.gz.csi"), sample=sample_names),
        bam=expand(str(ALIGN_DIR / "{sample}_recal.bam"), sample=sample_names),
        bai=expand(str(ALIGN_DIR / "{sample}_recal.bam.bai"), sample=sample_names),
        bed=str(BED_PATH),
        variant_summary=str(Path(POSTPROCESS["variant_summary"]).resolve()),
        pseudogene_bed=str(Path(REFERENCES["pseudogene_bed"]).resolve()),
    output:
        all_candidates=expand(str(REPORT_DIR / "{sample}_all_candidates.csv"), sample=sample_names),
        raw_output=expand(str(REPORT_DIR / "{sample}_raw_output.csv"), sample=sample_names),
        output_p=expand(str(REPORT_DIR / "{sample}_output_p.csv"), sample=sample_names),
        filtered_tp=expand(str(REPORT_DIR / "{sample}_raw_output_filtered_tp.csv"), sample=sample_names),
        filtered_tp_report=expand(str(REPORT_DIR / "{sample}_raw_output_filtered_tp_report.csv"), sample=sample_names),
        pseudogene_hits=expand(str(PSEUDOGENE_DIR / "{sample}_filtered_pseudogene_variants.txt"), sample=sample_names),
        coverage_summary=str(COVERAGE_SUMMARY),
        tp_zip=str(REPORT_DIR / "true_positive.zip"),
        tp_report_zip=str(REPORT_DIR / "true_positive_report.zip"),
        done=str(REPORT_DIR / ".postprocess.done"),
    log:
        str(LOG_DIR / "postprocess.log")
    threads: 1
    container:
        CONTAINERS["postprocess"]
    shell:
        """
        mkdir -p {REPORT_DIR} {PSEUDOGENE_DIR} {LOG_DIR}
        bash workflow/scripts/run_post_process.sh config/config.yaml > {log} 2>&1
        touch {output.done}
        """
