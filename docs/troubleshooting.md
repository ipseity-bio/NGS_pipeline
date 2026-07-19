# Troubleshooting

## Container Bind Issues

If tools inside the container cannot see the reference, FASTQ, or cache paths, add bind mounts through `SINGULARITY_ARGS`.

Example:

```bash
SINGULARITY_ARGS="--bind /home --bind /mnt" CORES=16 bash run.sh
```

## FASTQs Not Detected

Check both:

- the files are located under `input_dir`
- the filenames follow `<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`

## VEP Cache Not Found

Confirm that:

- `references.vep_cache` points to the cache directory
- the cache path is visible inside the container via `SINGULARITY_ARGS`

## Known-Sites or Reference Path Errors

Confirm that:

- `references.ref`, `references.bwa_index`, and `references.known_sites` are correct in `config/config.yaml`
- those paths are bind-mounted if they are outside the project directory
