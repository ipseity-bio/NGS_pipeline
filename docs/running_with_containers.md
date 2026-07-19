# Running With Containers

## Standard Run

Run the full workflow:

```bash
CORES=16 bash run.sh
```

## Bind Mounts

If your input, reference, or cache paths live outside the project directory, bind them into the container runtime.

Typical example:

```bash
SINGULARITY_ARGS="--bind /home --bind /mnt" CORES=16 bash run.sh
```

If all required paths share a common parent, binding that parent is usually sufficient.

Common cases that require bind mounts:

- FASTQs in `/home/...` or `/mnt/...`
- reference FASTA and BWA index outside the repository
- VEP cache outside the repository
- output directory outside the repository

## Direct Snakemake Execution

```bash
snakemake \
  --use-singularity \
  --singularity-args "--bind /home --bind /mnt" \
  --cores 16 \
  --keep-going \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml
```

