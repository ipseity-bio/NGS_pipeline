#!/usr/bin/env bash
set -euo pipefail

LOG="pipeline_timing_$(date +%Y%m%d_%H%M%S).log"

run_step () {
  local name="$1"; shift
  echo "==============================" | tee -a "$LOG"
  echo "[STEP] $name" | tee -a "$LOG"
  echo "[CMD ] $*" | tee -a "$LOG"
  echo "[START] $(date '+%F %T')" | tee -a "$LOG"

  /usr/bin/time -f "[ELAPSED] %E  [CPU] user:%U sys:%S  [MAXRSS] %M KB" \
    "$@" 2>>"$LOG"

  echo "[END] $(date '+%F %T')" | tee -a "$LOG"
  echo "" | tee -a "$LOG"
}

run_step "run_filter.py" \
  python3 ./run_filter.py

run_step "depth.py" \
  python3 ./depth.py

run_step "pseudogene_check.sh" \
  bash ./pseudogene_check.sh

run_step "coverage_table.py" \
  python3 ./coverage_table.py

run_step "True_positive_filter.py" \
  python3 ./True_positive_filter.py

run_step "coverage_summary.py" \
  python3 ./coverage_summary.py

echo "All steps done. Timing log: $LOG"
