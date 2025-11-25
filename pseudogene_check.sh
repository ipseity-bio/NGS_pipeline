#!/bin/bash
set -euo pipefail

# Define the pseudogene BED file
PSEUDOGENE_BED="./pseudogenes_37.bed"

for VCF_FILE in *_filtered.vcf.gz; do
    [ -e "$VCF_FILE" ] || continue 

    BASE="${VCF_FILE%.vcf.gz}"

    bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' "$VCF_FILE" \
      | awk 'BEGIN{OFS="\t"}{
            id=$4;
            if(id=="." || id=="") id=$1":"$3;  # fallback ID if missing
            print $1,$2,$3,id
        }' > "${BASE}.bed"

    bedtools intersect \
        -a "${BASE}.bed" \
        -b "$PSEUDOGENE_BED" \
        -wa -wb > "${BASE}_pseudogene_variants.txt"
done
