from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_vep(path: str):
    header_cols = None
    rows = []

    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            if raw.lstrip().startswith("##"):
                continue
            line = raw.rstrip("\n")
            if header_cols is None and line.lstrip().startswith("#"):
                header_cols = [
                    (column or "").lstrip("\ufeff").strip()
                    for column in re.split(r"\t+", line.lstrip()[1:].strip())
                ]
                continue
            if line.strip():
                row = line.split("\t")
                width = len(header_cols)
                if len(row) < width:
                    row = row + [""] * (width - len(row))
                elif len(row) > width:
                    row = row[:width]
                rows.append(row)

    if header_cols is None:
        raise ValueError(f"No single-# VEP header found in {path}")

    return header_cols, rows


def read_meta(path: str):
    meta = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                meta.append(line)
            elif line.lstrip().startswith("#"):
                break
    return meta


def main():
    parser = argparse.ArgumentParser(description="Merge HC and FB VEP tabular outputs.")
    parser.add_argument("--hc", required=True, help="HaplotypeCaller VEP output")
    parser.add_argument("--fb", required=True, help="FreeBayes VEP output")
    parser.add_argument("--output", required=True, help="Merged VEP output")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    header_hc, rows_hc = read_vep(args.hc)
    header_fb, rows_fb_raw = read_vep(args.fb)

    fb_index = {column: idx for idx, column in enumerate(header_fb)}
    rows_fb = []
    for row in rows_fb_raw:
        rows_fb.append([row[fb_index[column]] if column in fb_index else "" for column in header_hc])

    meta = read_meta(args.hc)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as out_handle:
        out_handle.writelines(meta)
        out_handle.write("#" + "\t".join(header_hc) + "\n")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows_hc)
        writer.writerows(rows_fb)

    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        f"Merged {len(rows_hc)} HC rows and {len(rows_fb)} FB rows into {len(rows_hc) + len(rows_fb)} records.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
