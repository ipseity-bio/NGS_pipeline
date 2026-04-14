import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument("--report-dir", required=True)
    parser.add_argument("--qc-min-30x", type=float, default=30.0)
    parser.add_argument("--positive-control-min-30x", type=float, default=30.0)
    parser.add_argument("--ntc-max-30x", type=float, default=1.0)
    args = parser.parse_args()

    coverage_summary_path = Path(args.coverage_summary)
    report_dir = Path(args.report_dir)
    coverage_summary = pd.read_csv(coverage_summary_path)

    for index, row in coverage_summary.iterrows():
        sample_name = str(row["Sample"]).strip()
        raw_report_path = report_dir / f"{sample_name}_raw_output_filtered_tp_report.csv"

        if raw_report_path.exists():
            raw_report = pd.read_csv(raw_report_path)
            if raw_report["Gene"].isna().all() and raw_report["Transcript"].isna().all():
                coverage_summary.at[index, "Status"] = "Negative"
            else:
                coverage_summary.at[index, "Status"] = "Positive"

            if row["Percentage of bases covered at 30X"] < args.qc_min_30x:
                coverage_summary.at[index, "Status"] = "Qc Fail"
        else:
            if sample_name.startswith("NTC"):
                coverage_summary.at[index, "Status"] = (
                    "Pass"
                    if row["Percentage of bases covered at 30X"] < args.ntc_max_30x
                    else "Fail"
                )
            else:
                coverage_summary.at[index, "Status"] = "Qc Fail"

        if sample_name.startswith("PC"):
            coverage_summary.at[index, "Status"] = (
                "Pass"
                if row["Percentage of bases covered at 30X"] > args.positive_control_min_30x
                else "Fail"
            )

    coverage_summary.to_csv(coverage_summary_path, index=False)


if __name__ == "__main__":
    main()
