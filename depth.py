import os, glob, subprocess, shutil, sys
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor, as_completed

BED_FILE = "regions.bed"
DEPTHS = [1, 10, 20, 30, 40, 50, 100, 500]


def bed_total_bases(bed_path: str) -> int:
    intervals = {}
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid BED line (expected >=3 columns): {line}")
            chrom, start, end = parts[0], parts[1], parts[2]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                raise ValueError(f"Invalid BED coordinates: {line}")
            if end_i < start_i:
                raise ValueError(f"BED end < start: {line}")
            if end_i == start_i:
                continue 
            intervals.setdefault(chrom, []).append((start_i, end_i))

    total = 0
    for chrom, ivals in intervals.items():
        ivals.sort()
        cur_s, cur_e = ivals[0]
        for s, e in ivals[1:]:
            if s <= cur_e: 
                cur_e = max(cur_e, e)
            else:
                total += cur_e - cur_s
                cur_s, cur_e = s, e
        total += cur_e - cur_s
    return total


def check_prereqs(bed_path: str):
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    if shutil.which("samtools") is None:
        raise EnvironmentError("samtools not found on PATH. Please install samtools.")
    if os.path.getsize(bed_path) == 0:
        raise ValueError(f"BED file is empty: {bed_path}")


def process_bam(bam_file: str, bed_bases: int):
    prefix = bam_file.removesuffix("_recal.bam")
    out_txt = prefix + "_all_coverage.txt"

    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM not found: {bam_file}")
    if not (os.path.exists(bam_file + ".bai") or os.path.exists(prefix + "_aligned.bai")):
        print(f"Warning: BAM index not found for {bam_file}. Depth may be slow/fail.", file=sys.stderr)

    cmd = ["samtools", "depth", "-b", BED_FILE, "-a", bam_file]

    seen_bases = 0 
    sum_depth = 0    
    covered_counts = [0] * len(DEPTHS)

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    try:
        for line in p.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                depth = int(line.rsplit("\t", 1)[-1])
            except ValueError:
                continue

            seen_bases += 1
            sum_depth += depth

            if depth > 0:
                k = bisect_right(DEPTHS, depth)
                for i in range(k):
                    covered_counts[i] += 1
    finally:
        if p.stdout:
            p.stdout.close()

    stderr = p.stderr.read() if p.stderr else ""
    if p.stderr:
        p.stderr.close()
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(
            f"samtools failed on {bam_file} (exit {rc}).\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{stderr.strip()}"
        )

    denom = bed_bases

    if denom == 0:
        coverage_output = ["No non-zero depth bases found in BED regions."]
        mean_depth = 0.0
    else:
        coverage_output = []
        for depth_thr, c in zip(DEPTHS, covered_counts):
            percent = (c / denom) * 100.0
            coverage_output.append(
                f"Percentage of bases covered at {depth_thr}X: {percent:.2f}%"
            )
        mean_depth = sum_depth / denom
        coverage_output.append(f"Mean depth of coverage: {mean_depth:.2f}")

        if seen_bases != denom:
            print(
                f"Note: {bam_file}: samtools reported {seen_bases} positions, "
                f"BED length is {denom}. Denominator uses BED length.",
                file=sys.stderr,
            )

    with open(out_txt, "w") as f:
        f.write("\n".join(coverage_output))

    return bam_file, out_txt, mean_depth


def main():
    try:
        check_prereqs(BED_FILE)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        bed_bases = bed_total_bases(BED_FILE)
    except Exception as e:
        print(f"Error parsing BED: {e}", file=sys.stderr)
        sys.exit(1)

    bam_files = sorted(glob.glob("*_recal.bam"))
    if not bam_files:
        print("No *_recal.bam files found.")
        return

    workers = min(len(bam_files), os.cpu_count() or 1)
    print(f"Processing {len(bam_files)} BAMs with {workers} workers...")

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(process_bam, bam, bed_bases) for bam in bam_files]
        for fut in as_completed(futures):
            bam, out_txt, mean_depth = fut.result()
            print(f"Done: {bam} -> {out_txt} (mean depth {mean_depth:.2f})")


if __name__ == "__main__":
    main()