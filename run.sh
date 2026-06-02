#!/usr/bin/env bash
set -euo pipefail

config_path="${CONFIG_PATH:-config/config.yaml}"
snakefile_path="${SNAKEFILE_PATH:-workflow/Snakefile}"
cores="${CORES:-$(nproc)}"
sdm="${SNAKEMAKE_SDM:-apptainer}"
singularity_args="${SINGULARITY_ARGS:-}"

container_flag="--use-singularity"
case "$sdm" in
  apptainer|singularity|"")
    container_flag="--use-singularity"
    ;;
  *)
    echo "Unsupported SNAKEMAKE_SDM='$sdm' for this Snakemake version. Use 'apptainer' or 'singularity'." >&2
    exit 1
    ;;
esac

snakemake_args=(
  "$container_flag"
  --cores "$cores"
  --keep-going
  --snakefile "$snakefile_path"
  --configfile "$config_path"
)

if [[ -n "$singularity_args" ]]; then
  snakemake_args+=(--singularity-args "$singularity_args")
fi

snakemake "${snakemake_args[@]}"
