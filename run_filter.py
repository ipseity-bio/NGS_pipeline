import filter

file_pattern = '*_merged_vep.vcf'
variant_summary_path = './variant_summary.txt'

filter.process_vcf(file_pattern, variant_summary_path)
