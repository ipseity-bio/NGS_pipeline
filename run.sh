snakemake --cores 32 --keep-going --snakefile snakemake.smk --configfile ./config.yaml

bash ./run_post_process.sh

zip true_positive *_tp.csv
zip true_positive_report *_tp_report.csv
