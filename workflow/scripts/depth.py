import argparse
import glob
import os
import shutil
import subprocess
import sys
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


DEPTHS = [1, 10, 20, 30, 40, 50, 100, 500]


def bed_total_bases(bed_path: str) -> int:
    intervals = {}
    with open(bed_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid BED line: {line}")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end < start:
                raise ValueError(f"BED end < start: {line}")
            if end == start:
                continue
            intervals.setdefault(chrom, []).append((start, end))

    total = 0
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
        cur_s, cur_e = chrom_intervals[0]
        for start, end in chrom_intervals[1:]:
            if start <= cur_e:
                cur_e = max(cur_e, end)
            else:
                total += cur_e - cur_s
                cur_s, cur_e = start, end
        total += cur_e - cur_s
    return total


def process_bam(bam_file: str, bed_file: str, bed_bases: int, output_dir: str):
    bam_path = Path(bam_file)
    sample = bam_path.name.removesuffix("_recal.bam")
    out_txt = Path(output_dir) / f"{sample}_all_coverage.txt"
    cmd = ["samtools", "depth", "-b", bed_file, "-a", bam_file]

    seen_bases = 0
    sum_depth = 0
    covered_counts = [0] * len(DEPTHS)

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    try:
        assert process.stdout is not None
        for line in process.stdout:
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
                for idx in range(bisect_right(DEPTHS, depth)):
                    covered_counts[idx] += 1
    finally:
        if process.stdout:
            process.stdout.close()

    stderr = process.stderr.read() if process.stderr else ""
    if process.stderr:
        process.stderr.close()
    rc = process.wait()
    if rc != 0:
        raise RuntimeError(f"samtools failed on {bam_file} (exit {rc}).\n{stderr.strip()}")

    coverage_output = []
    if bed_bases == 0:
        coverage_output.append("No non-zero depth bases found in BED regions.")
        mean_depth = 0.0
    else:
        for depth_thr, count in zip(DEPTHS, covered_counts):
            coverage_output.append(
                f"Percentage of bases covered at {depth_thr}X: {(count / bed_bases) * 100.0:.2f}%"
            )
        mean_depth = sum_depth / bed_bases
        coverage_output.append(f"Mean depth of coverage: {mean_depth:.2f}")
        if seen_bases != bed_bases:
            print(
                f"Note: {bam_file}: samtools reported {seen_bases} positions, BED length is {bed_bases}.",
                file=sys.stderr,
            )

    out_txt.write_text("\n".join(coverage_output), encoding="utf-8")
    return bam_file, str(out_txt), mean_depth


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed-file", required=True)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()

    if shutil.which("samtools") is None:
        raise EnvironmentError("samtools not found on PATH.")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    bed_bases = bed_total_bases(args.bed_file)
    bam_files = sorted(glob.glob(os.path.join(args.bam_dir, "*_recal.bam")))
    if not bam_files:
        print("No *_recal.bam files found.")
        return

    workers = min(len(bam_files), max(1, args.workers))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(process_bam, bam, args.bed_file, bed_bases, args.output_dir)
            for bam in bam_files
        ]
        for future in as_completed(futures):
            bam, out_txt, mean_depth = future.result()
            print(f"Done: {bam} -> {out_txt} (mean depth {mean_depth:.2f})")


if __name__ == "__main__":
    main()
