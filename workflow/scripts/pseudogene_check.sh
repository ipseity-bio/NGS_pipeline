#!/usr/bin/env bash
set -euo pipefail

variant_dir="${1:?variant_dir is required}"
pseudogene_bed="${2:?pseudogene_bed is required}"
output_dir="${3:?output_dir is required}"

mkdir -p "$output_dir"

for vcf_file in "$variant_dir"/*_filtered.vcf.gz; do
    [ -e "$vcf_file" ] || continue

    base="$(basename "${vcf_file%.vcf.gz}")"
    bed_out="$output_dir/${base}.bed"
    hit_out="$output_dir/${base}_pseudogene_variants.txt"

    bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' "$vcf_file" \
      | awk 'BEGIN{OFS="\t"}{
            id=$4;
            if(id=="." || id=="") id=$1":"$3;
            print $1,$2,$3,id
        }' > "$bed_out"

    bedtools intersect -a "$bed_out" -b "$pseudogene_bed" -wa -wb > "$hit_out"
done
