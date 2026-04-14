import argparse
from pathlib import Path

import pandas as pd


def sample_sort_key(sample: str):
    digits = "".join(ch for ch in sample if ch.isdigit())
    return (sample.rstrip(digits), int(digits) if digits else -1, sample)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-dir", required=True)
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()

    data = []
    for file_path in sorted(Path(args.coverage_dir).glob("*_all_coverage.txt")):
        content = file_path.read_text(encoding="utf-8").splitlines()
        sample = file_path.name.split("_all_coverage")[0]
        data.append(
            [
                sample,
                float(content[0].split(": ")[1].rstrip("%")),
                float(content[1].split(": ")[1].rstrip("%")),
                float(content[2].split(": ")[1].rstrip("%")),
                float(content[3].split(": ")[1].rstrip("%")),
                float(content[4].split(": ")[1].rstrip("%")),
                float(content[5].split(": ")[1].rstrip("%")),
                float(content[6].split(": ")[1].rstrip("%")),
                float(content[7].split(": ")[1].rstrip("%")),
                float(content[8].split(": ")[1]),
            ]
        )

    columns = [
        "Sample",
        "Percentage of bases covered at 1X",
        "Percentage of bases covered at 10X",
        "Percentage of bases covered at 20X",
        "Percentage of bases covered at 30X",
        "Percentage of bases covered at 40X",
        "Percentage of bases covered at 50X",
        "Percentage of bases covered at 100X",
        "Percentage of bases covered at 500X",
        "Mean depth of coverage",
    ]
    df = pd.DataFrame(data, columns=columns).fillna(0)
    df = df.sort_values(by="Sample", key=lambda col: col.map(sample_sort_key))
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
