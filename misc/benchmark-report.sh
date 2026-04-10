#!/usr/bin/env bash
# =========================================================================
# Generate a Markdown report from benchmark results.
#
# Usage:
#   ./benchmark-report.sh <results-dir>              # single run
#   ./benchmark-report.sh <dir1> <dir2> ...           # compare runs
#   ./benchmark-report.sh benchmark-results/*/        # all runs
#
# Output: writes REPORT.md into the first results directory
#         (or into benchmark-results/ when comparing multiple)
# =========================================================================
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <results-dir> [<results-dir2> ...]"
  exit 1
fi

DIRS=("$@")
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"

# Where to write the report
if [[ ${#DIRS[@]} -eq 1 ]]; then
  REPORT="${DIRS[0]}/REPORT-${TIMESTAMP}.md"
else
  mkdir -p "${ROOT}/benchmark-results"
  REPORT="${ROOT}/benchmark-results/REPORT-${TIMESTAMP}.md"
fi

# ---------- Helpers ----------
read_field() {
  local file="$1" key="$2"
  grep "^${key}:" "$file" 2>/dev/null | sed "s/^${key}:[[:space:]]*//" || true
}

# ---------- Collect system info ----------
declare -a HOSTS=()
for dir in "${DIRS[@]}"; do
  dir="${dir%/}"
  sysfile="${dir}/system.txt"
  if [[ -f "$sysfile" ]]; then
    HOSTS+=("$sysfile")
  fi
done

# ---------- Start report ----------
{
  echo "# MDOODZ Benchmark Report"
  echo ""
  echo "_Generated: $(date -u +%Y-%m-%d\ %H:%M\ UTC)_"
  echo ""

  # --- System info ---
  echo "## System Information"
  echo ""
  if [[ ${#HOSTS[@]} -eq 1 ]]; then
    f="${HOSTS[0]}"
    echo "| Property | Value |"
    echo "|----------|-------|"
    echo "| Hostname | $(read_field "$f" hostname) |"
    echo "| OS | $(read_field "$f" os) $(read_field "$f" arch) |"
    echo "| CPUs | $(read_field "$f" cpus) |"
    echo "| RAM | $(read_field "$f" ram_mb) MB |"
    echo "| Compiler | $(read_field "$f" compiler) |"
    echo "| Date | $(read_field "$f" date) |"
    # Test results if available
    tp=$(read_field "$f" tests_passed)
    tf=$(read_field "$f" tests_failed)
    if [[ -n "$tp" || -n "$tf" ]]; then
      te=$(read_field "$f" tests_exit)
      status="PASS"
      [[ "$te" != "0" ]] && status="FAIL"
      echo "| Tests | ${tp:-0} passed, ${tf:-0} failed (${status}) |"
    fi
  else
    echo "| Property |"
    for f in "${HOSTS[@]}"; do
      printf " %s |" "$(read_field "$f" hostname)"
    done
    echo ""

    echo "|----------|"
    for f in "${HOSTS[@]}"; do printf "------|"; done
    echo ""

    for key in os arch cpus ram_mb compiler date; do
      printf "| %s |" "$key"
      for f in "${HOSTS[@]}"; do
        printf " %s |" "$(read_field "$f" "$key")"
      done
      echo ""
    done
  fi
  echo ""

  # --- Summary table ---
  echo "## Results Summary"
  echo ""
  echo "| Host | Grid | Threads | Steps | Avg Wall (s) | Avg Rheology (s) | Avg Assembly (s) | Avg Solve (s) | Avg Nit | Peak RSS (MB) |"
  echo "|------|------|---------|-------|-------------|-----------------|-----------------|--------------|---------|--------------|"

  for dir in "${DIRS[@]}"; do
    dir="${dir%/}"
    csv="${dir}/summary.csv"
    [[ -f "$csv" ]] || continue
    tail -n+2 "$csv" | while IFS=',' read -r host os arch grid nx nz threads steps avg_wall total_wall avg_rheo avg_asm avg_sol avg_nit peak_rss; do
      printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" \
        "$host" "$grid" "$threads" "$steps" "$avg_wall" "$avg_rheo" "$avg_asm" "$avg_sol" "$avg_nit" "$peak_rss"
    done
  done
  echo ""

  # --- Thread scaling ---
  echo "## Thread Scaling"
  echo ""
  echo "Speedup relative to single-thread baseline for same grid size."
  echo ""
  echo "| Host | Grid | Threads | Avg Wall (s) | Speedup | Efficiency |"
  echo "|------|------|---------|-------------|---------|------------|"

  for dir in "${DIRS[@]}"; do
    dir="${dir%/}"
    csv="${dir}/summary.csv"
    [[ -f "$csv" ]] || continue

    # Get unique grids
    grids_list=$(tail -n+2 "$csv" | cut -d',' -f4 | sort -u)
    for grid in $grids_list; do
      # Get baseline (1 thread)
      baseline=$(grep ",${grid}," "$csv" | grep ",1," | head -1 | cut -d',' -f9)
      if [[ -z "$baseline" || "$baseline" == "0" ]]; then
        # Use lowest thread count as baseline
        baseline=$(grep ",${grid}," "$csv" | sort -t',' -k7 -n | head -1 | cut -d',' -f9)
      fi
      [[ -z "$baseline" ]] && continue

      grep ",${grid}," "$csv" | sort -t',' -k7 -n | while IFS=',' read -r host os arch g nx nz threads steps avg_wall rest; do
        speedup=$(awk "BEGIN { printf \"%.2f\", ${baseline}/${avg_wall} }")
        efficiency=$(awk "BEGIN { printf \"%.0f%%\", 100.0*${baseline}/${avg_wall}/${threads} }")
        printf "| %s | %s | %s | %s | %sx | %s |\n" "$host" "$grid" "$threads" "$avg_wall" "$speedup" "$efficiency"
      done
    done
  done
  echo ""

  # --- Grid scaling ---
  echo "## Grid Scaling"
  echo ""
  echo "How wall time grows with problem size (at max thread count per host)."
  echo ""
  echo "| Host | Grid | Cells | Threads | Avg Wall (s) | Wall/Cell (µs) |"
  echo "|------|------|-------|---------|-------------|---------------|"

  for dir in "${DIRS[@]}"; do
    dir="${dir%/}"
    csv="${dir}/summary.csv"
    [[ -f "$csv" ]] || continue

    # Get max thread count
    max_thr=$(tail -n+2 "$csv" | cut -d',' -f7 | sort -rn | head -1)

    tail -n+2 "$csv" | grep ",${max_thr}," | sort -t',' -k5 -n | while IFS=',' read -r host os arch grid nx nz threads steps avg_wall rest; do
      cells=$((nx * nz))
      wall_per_cell=$(awk "BEGIN { printf \"%.2f\", ${avg_wall}/${cells}*1e6 }")
      printf "| %s | %s | %d | %s | %s | %s |\n" "$host" "$grid" "$cells" "$threads" "$avg_wall" "$wall_per_cell"
    done
  done
  echo ""

  # --- Memory ---
  echo "## Memory Usage"
  echo ""
  echo "| Host | Grid | Threads | Peak RSS (MB) | RSS/Cell (KB) |"
  echo "|------|------|---------|--------------|--------------|"

  for dir in "${DIRS[@]}"; do
    dir="${dir%/}"
    csv="${dir}/summary.csv"
    [[ -f "$csv" ]] || continue

    max_thr=$(tail -n+2 "$csv" | cut -d',' -f7 | sort -rn | head -1)
    tail -n+2 "$csv" | grep ",${max_thr}," | sort -t',' -k5 -n | while IFS=',' read -r host os arch grid nx nz threads steps avg_wall total_wall avg_rheo avg_asm avg_sol avg_nit peak_rss; do
      cells=$((nx * nz))
      rss_per_cell=$(awk "BEGIN { printf \"%.1f\", ${peak_rss}/${cells}*1024 }")
      printf "| %s | %s | %s | %s | %s |\n" "$host" "$grid" "$threads" "$peak_rss" "$rss_per_cell"
    done
  done
  echo ""

  # --- Per-timestep detail (first run only) ---
  echo "## Per-Timestep Detail (first run)"
  echo ""
  first_dir="${DIRS[0]%/}"
  first_perf=$(find "$first_dir" -name "perf.csv" -type f | head -1)
  if [[ -n "$first_perf" && -f "$first_perf" ]]; then
    tag="$(basename "$(dirname "$first_perf")")"
    echo "_Run: ${tag}_"
    echo ""
    echo "| Step | Wall (s) | Rheology (s) | Assembly (s) | Solve (s) | Nit | Particles | RSS (MB) |"
    echo "|------|---------|-------------|-------------|----------|-----|-----------|---------|"
    tail -n+2 "$first_perf" | while IFS=',' read -r step wall rheo asm sol nit npart neq_m neq_c rss ucpu scpu; do
      printf "| %s | %s | %s | %s | %s | %s | %s | %s |\n" \
        "$step" "$wall" "$rheo" "$asm" "$sol" "$nit" "$npart" "$rss"
    done
  fi
  echo ""

} > "$REPORT"

echo "Report written to: $REPORT"
