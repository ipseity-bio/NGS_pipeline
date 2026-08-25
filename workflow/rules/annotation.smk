ruleorder: ensembl_vep_hc > ensembl_vep_fb > merge_vep


rule ensembl_vep_hc:
    input:
        vcf=str(VARIANT_DIR / "{sample}_filtered.vcf.gz"),
        csi=str(VARIANT_DIR / "{sample}_filtered.vcf.gz.csi")
    output:
        vep_out=str(ANNOTATION_DIR / "{sample}_vep.vcf")
    log:
        str(LOG_DIR / "ensembl_vep_{sample}.log")
    threads: 8
    container:
        CONTAINERS["vep"]
    shell:
        """
        mkdir -p {ANNOTATION_DIR} {LOG_DIR}
        vep \
            --fork {threads} \
            --offline \
            --cache \
            --sift b \
            --polyphen b \
            --ccds \
            --symbol \
            --canonical \
            --protein \
            --biotype \
            --af \
            --pubmed \
            --uniprot \
            --variant_class \
            --gene_phenotype \
            --mirna \
            --check_existing \
            --hgvs \
            --fasta {REF} \
            --spdi \
            --hgvsg \
            --individual all \
            --allele_number \
            --dir_cache {VEP_CACHE} \
            --buffer_size 5000 \
            --force_overwrite \
            --assembly GRCh37 \
            --input_file {input.vcf} \
            --output_file {output.vep_out} \
            --tab \
            --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,IMPACT,Existing_variation,CLIN_SIG,SIFT,PolyPhen,CCDS,SYMBOL,ENSP,CANONICAL,BIOTYPE,AF,MAX_AF,PUBMED,VARIANT_CLASS,GENE_PHENO,SOMATIC,PHENO,HGVSg,HGVSc,HGVSp,HGVS_OFFSET,SPDI,ZYG,IND,ALLELE_NUM,REF_ALLELE" \
            --cache_version 115 \
            --refseq \
            --show_ref_allele \
            --max_af \
            > {log} 2>&1
        """


rule ensembl_vep_fb:
    input:
        vcf=str(COMPARE_DIR / "{sample}" / "0000.vcf")
    output:
        vep_out=str(ANNOTATION_DIR / "{sample}_fb_vep.vcf")
    log:
        str(LOG_DIR / "ensembl_vep_rsid_{sample}.log")
    threads: 8
    container:
        CONTAINERS["vep"]
    shell:
        """
        mkdir -p {ANNOTATION_DIR} {LOG_DIR}
        vep \
            --fork {threads} \
            --offline \
            --cache \
            --sift b \
            --polyphen b \
            --ccds \
            --symbol \
            --canonical \
            --protein \
            --biotype \
            --af \
            --pubmed \
            --uniprot \
            --variant_class \
            --gene_phenotype \
            --mirna \
            --check_existing \
            --hgvs \
            --fasta {REF} \
            --spdi \
            --hgvsg \
            --individual all \
            --allele_number \
            --dir_cache {VEP_CACHE} \
            --buffer_size 5000 \
            --force_overwrite \
            --assembly GRCh37 \
            --input_file {input.vcf} \
            --output_file {output.vep_out} \
            --tab \
            --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,IMPACT,Existing_variation,CLIN_SIG,SIFT,PolyPhen,CCDS,SYMBOL,ENSP,CANONICAL,BIOTYPE,AF,MAX_AF,PUBMED,VARIANT_CLASS,GENE_PHENO,SOMATIC,PHENO,HGVSg,HGVSc,HGVSp,HGVS_OFFSET,SPDI,ZYG,IND,ALLELE_NUM,REF_ALLELE" \
            --cache_version 115 \
            --refseq \
            --show_ref_allele \
            --max_af \
            > {log} 2>&1
        """


rule merge_vep:
    input:
        hc=str(ANNOTATION_DIR / "{sample}_vep.vcf"),
        fb=str(ANNOTATION_DIR / "{sample}_fb_vep.vcf")
    output:
        merged=str(ANNOTATION_DIR / "{sample}_merged_vep.vcf")
    log:
        str(LOG_DIR / "merge_vep_{sample}.log")
    container:
        CONTAINERS["python"]
    shell:
        """
        mkdir -p {ANNOTATION_DIR} {LOG_DIR}
        python3 workflow/scripts/merge_vep.py \
            --hc {input.hc} \
            --fb {input.fb} \
            --output {output.merged} \
            --log {log}
        """
