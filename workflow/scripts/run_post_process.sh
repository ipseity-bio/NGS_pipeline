#!/usr/bin/env bash
set -euo pipefail

config_path="${1:-config/config.yaml}"

read_config() {
  python3 - "$config_path" "$1" <<'PY'
import sys

config_path, key = sys.argv[1], sys.argv[2]


def parse_scalar(value: str):
    value = value.strip()
    if not value:
        return ""
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def load_simple_yaml(path: str):
    root = {}
    stack = [(-1, root)]

    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            indent = len(line) - len(line.lstrip(" "))
            stripped = line.strip()
            if ":" not in stripped:
                raise ValueError(f"Unsupported config line: {line}")

            field, value = stripped.split(":", 1)
            field = field.strip()
            value = value.strip()

            while stack and indent <= stack[-1][0]:
                stack.pop()

            parent = stack[-1][1]
            if value == "":
                parent[field] = {}
                stack.append((indent, parent[field]))
            else:
                parent[field] = parse_scalar(value)

    return root


config = load_simple_yaml(config_path)

parts = key.split(".")
value = config
for part in parts:
    value = value[part]
print(value)
PY
}

output_dir="$(read_config output_dir)"
bed_file="$(read_config references.bed_file)"
variant_summary="$(read_config postprocess.variant_summary)"
pseudogene_bed="$(read_config references.pseudogene_bed)"
coverage_subdir="$(read_config postprocess.coverage_dir)"
report_subdir="$(read_config postprocess.report_dir)"
qc_min_30x="$(read_config postprocess.qc_min_30x)"
pc_min_30x="$(read_config postprocess.positive_control_min_30x)"
ntc_max_30x="$(read_config postprocess.ntc_max_30x)"
require_protein_coding_for_rare_high_impact="$(read_config reporting.require_protein_coding_for_rare_high_impact)"

annotation_dir="$output_dir/annotation"
variant_dir="$output_dir/variants"
alignment_dir="$output_dir/alignment"
coverage_dir="$output_dir/$coverage_subdir"
report_dir="$output_dir/$report_subdir"
pseudogene_dir="$report_dir/pseudogene"
coverage_summary="$report_dir/$(basename "$output_dir")_coverage_summary.csv"
log="$report_dir/pipeline_timing_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$coverage_dir" "$report_dir" "$pseudogene_dir"

run_step() {
  local name="$1"; shift
  echo "==============================" | tee -a "$log"
  echo "[STEP] $name" | tee -a "$log"
  echo "[CMD ] $*" | tee -a "$log"
  echo "[START] $(date '+%F %T')" | tee -a "$log"
  if command -v /usr/bin/time >/dev/null 2>&1; then
    /usr/bin/time -f "[ELAPSED] %E  [CPU] user:%U sys:%S  [MAXRSS] %M KB" "$@" 2>>"$log"
  elif command -v time >/dev/null 2>&1; then
    time "$@" 2>>"$log"
  else
    echo "[WARN] time command not found; running step without timing metrics." | tee -a "$log"
    "$@" 2>>"$log"
  fi
  echo "[END] $(date '+%F %T')" | tee -a "$log"
  echo "" | tee -a "$log"
}

run_step "run_filter.py" \
  python3 workflow/scripts/run_filter.py \
  --annotation-dir "$annotation_dir" \
  --variant-dir "$variant_dir" \
  --variant-summary "$variant_summary" \
  --output-dir "$report_dir" \
  --require-protein-coding-for-rare-high-impact "$require_protein_coding_for_rare_high_impact"

run_step "depth.py" \
  python3 workflow/scripts/depth.py \
  --bed-file "$bed_file" \
  --bam-dir "$alignment_dir" \
  --output-dir "$coverage_dir"

run_step "pseudogene_check.sh" \
  bash workflow/scripts/pseudogene_check.sh "$variant_dir" "$pseudogene_bed" "$pseudogene_dir"

run_step "coverage_table.py" \
  python3 workflow/scripts/coverage_table.py \
  --coverage-dir "$coverage_dir" \
  --output-file "$coverage_summary"

run_step "True_positive_filter.py" \
  python3 workflow/scripts/True_positive_filter.py \
  --input-dir "$report_dir" \
  --output-dir "$report_dir" \
  --coverage-summary "$coverage_summary"

run_step "coverage_summary.py" \
  python3 workflow/scripts/coverage_summary.py \
  --coverage-summary "$coverage_summary" \
  --report-dir "$report_dir" \
  --qc-min-30x "$qc_min_30x" \
  --positive-control-min-30x "$pc_min_30x" \
  --ntc-max-30x "$ntc_max_30x"

zip -j "$report_dir/true_positive.zip" "$report_dir"/*_tp.csv >/dev/null
zip -j "$report_dir/true_positive_report.zip" "$report_dir"/*_tp_report.csv >/dev/null

echo "All steps done. Timing log: $log"
