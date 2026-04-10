#!/usr/bin/env bash
# =========================================================================
# MDOODZ Performance Benchmark Script
#
# Runs RiftingBasic at multiple grid sizes and thread counts,
# collecting perf.csv from each run into a results directory.
#
# Usage:
#   ./benchmark.sh                    # default: all grids, all thread counts
#   ./benchmark.sh --grids "150x100 301x201"
#   ./benchmark.sh --threads "1 4 8"
#   ./benchmark.sh --steps 5
#   ./benchmark.sh --quick            # small grids, few steps
#
# Prerequisites:
#   make build SET=RiftingBasic       # or the script builds it for you
# =========================================================================
set -euo pipefail

# ---------- Configuration (override via flags) ----------
GRIDS="${BENCH_GRIDS:-150x100 301x201 501x501}"
THREADS="${BENCH_THREADS:-1 2 4 8 12 16}"
STEPS="${BENCH_STEPS:-10}"
SCENARIO="RiftingBasic"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_DIR="${ROOT}/benchmark-results/$(date +%Y%m%d-%H%M%S)"
QUICK=0
SKIP_BUILD=0
RUN_TESTS=0

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --grids)    GRIDS="$2"; shift 2 ;;
    --threads)  THREADS="$2"; shift 2 ;;
    --steps)    STEPS="$2"; shift 2 ;;
    --quick)    QUICK=1; shift ;;
    --skip-build) SKIP_BUILD=1; shift ;;
    --run-tests)  RUN_TESTS=1; shift ;;
    --help|-h)
      head -15 "$0" | grep '^#' | sed 's/^# *//'
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ $QUICK -eq 1 ]]; then
  GRIDS="150x100 301x201"
  THREADS="1 4"
  STEPS=3
fi

# ---------- Detect system ----------
OS="$(uname -s)"
ARCH="$(uname -m)"
NCPU=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
HOSTNAME_SHORT="$(hostname -s 2>/dev/null || hostname)"
RAM_MB=0
if [[ "$OS" == "Linux" ]]; then
  RAM_MB=$(awk '/MemTotal/ {printf "%d", $2/1024}' /proc/meminfo)
elif [[ "$OS" == "Darwin" ]]; then
  RAM_MB=$(( $(sysctl -n hw.memsize) / 1024 / 1024 ))
fi

echo "========================================"
echo "  MDOODZ Performance Benchmark"
echo "========================================"
echo "Host:     ${HOSTNAME_SHORT} (${OS} ${ARCH})"
echo "CPUs:     ${NCPU}"
echo "RAM:      ${RAM_MB} MB"
echo "Grids:    ${GRIDS}"
echo "Threads:  ${THREADS}"
echo "Steps:    ${STEPS}"
echo "Results:  ${RESULTS_DIR}"
echo "========================================"

# ---------- Build ----------
if [[ $SKIP_BUILD -eq 0 ]]; then
  echo ""
  echo ">>> Building ${SCENARIO} (OPT=ON, OMP=ON)..."
  cd "$ROOT"
  make build SET="${SCENARIO}" 2>&1 | tail -5
fi

EXEC_DIR="${ROOT}/cmake-exec/${SCENARIO}"
EXEC="${EXEC_DIR}/${SCENARIO}"
TXT_ORIG="${EXEC_DIR}/${SCENARIO}.txt"

if [[ ! -x "$EXEC" ]]; then
  echo "ERROR: Executable not found: $EXEC"
  echo "Run: make build SET=${SCENARIO}"
  exit 1
fi

# ---------- Create results directory ----------
mkdir -p "$RESULTS_DIR"

# Save system info
cat > "${RESULTS_DIR}/system.txt" <<EOF
hostname: ${HOSTNAME_SHORT}
os: ${OS}
arch: ${ARCH}
cpus: ${NCPU}
ram_mb: ${RAM_MB}
date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
uname: $(uname -a)
compiler: $(cc --version 2>/dev/null | head -1 || echo unknown)
EOF

# ---------- Helper: patch .txt file ----------
make_bench_txt() {
  local nx="$1" nz="$2" steps="$3" dst="$4"
  cp "$TXT_ORIG" "$dst"
  # Patch grid size
  sed -i.bak -E "s/^Nx[[:space:]]*=.*/Nx      = ${nx}/" "$dst"
  sed -i.bak -E "s/^Nz[[:space:]]*=.*/Nz      = ${nz}/" "$dst"
  # Patch timesteps
  sed -i.bak -E "s/^Nt[[:space:]]*=.*/Nt      = ${steps}/" "$dst"
  # Disable HDF5 output (pure compute benchmark)
  sed -i.bak -E "s/^writer[[:space:]]*=.*/writer = 0/" "$dst"
  # Ensure logging is minimal
  sed -i.bak -E "s/^log_timestamp[[:space:]]*=.*/log_timestamp = 0/" "$dst"
  rm -f "${dst}.bak"
}

# ---------- Run benchmarks ----------
echo ""
TOTAL=0
for grid in $GRIDS; do
  for thr in $THREADS; do
    TOTAL=$((TOTAL + 1))
  done
done

RUN=0
SUMMARY_CSV="${RESULTS_DIR}/summary.csv"
echo "host,os,arch,grid,nx,nz,threads,steps,avg_wall_s,total_wall_s,avg_rheology_s,avg_assembly_s,avg_solve_s,avg_nit,peak_rss_mb" > "$SUMMARY_CSV"

for grid in $GRIDS; do
  NX="${grid%x*}"
  NZ="${grid#*x}"

  for thr in $THREADS; do
    # Skip if thread count exceeds available CPUs
    if [[ $thr -gt $NCPU ]]; then
      echo "  [skip] ${grid} @ ${thr}t — only ${NCPU} CPUs available"
      continue
    fi

    RUN=$((RUN + 1))
    TAG="${grid}_${thr}t"
    RUN_DIR="${RESULTS_DIR}/${TAG}"
    mkdir -p "$RUN_DIR"

    echo ""
    echo ">>> [${RUN}/${TOTAL}] Grid ${NX}x${NZ}, ${thr} threads..."

    # Prepare .txt
    BENCH_TXT="${RUN_DIR}/bench.txt"
    make_bench_txt "$NX" "$NZ" "$STEPS" "$BENCH_TXT"

    # Run
    cd "$EXEC_DIR"
    rm -f perf.csv
    export OMP_NUM_THREADS="$thr"

    T0=$(date +%s)
    "./${SCENARIO}" "${RUN_DIR}/bench.txt" > "${RUN_DIR}/stdout.log" 2>&1 || true
    T1=$(date +%s)
    ELAPSED=$((T1 - T0))

    echo "    Completed in ${ELAPSED}s"

    # Collect results
    if [[ -f perf.csv ]]; then
      cp perf.csv "${RUN_DIR}/perf.csv"

      # Compute averages for summary
      AVG=$(awk -F',' 'NR>1 {
        n++; wall+=$2; rheo+=$3; asm+=$4; sol+=$5; nit+=$6; rss=$10
      } END {
        if(n>0) printf "%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.1f,%.1f", n, wall/n, wall, rheo/n, asm/n, sol/n, nit/n, rss
        else print "0,0,0,0,0,0,0,0"
      }' perf.csv)

      echo "${HOSTNAME_SHORT},${OS},${ARCH},${grid},${NX},${NZ},${thr},${AVG}" >> "$SUMMARY_CSV"
      echo "    $(head -1 perf.csv)"
      tail -1 perf.csv | awk -F',' '{printf "    Last step: wall=%.3fs, RSS=%.0fMB, nit=%s\n", $2, $10, $6}'
    else
      echo "    WARNING: No perf.csv produced"
    fi

    # Clean up HDF5 output
    rm -f "${EXEC_DIR}"/Output*.gzip.h5 "${EXEC_DIR}"/Particles*.gzip.h5 "${EXEC_DIR}"/Breakpoint*.dat
  done
done

# ---------- Print summary ----------
echo ""
echo "========================================"
echo "  Benchmark Complete"
echo "========================================"
echo ""
echo "Results in: ${RESULTS_DIR}/"
echo ""
echo "Summary:"
column -t -s',' "$SUMMARY_CSV" 2>/dev/null || cat "$SUMMARY_CSV"

# ---------- Run tests (optional) ----------
if [[ $RUN_TESTS -eq 1 ]]; then
  echo ""
  echo "========================================"
  echo "  Running Test Suite"
  echo "========================================"
  TEST_BUILD="${ROOT}/cmake-build-test"
  cmake -DTEST=ON -DOPT=ON -DOMP=ON -B "$TEST_BUILD" -S "$ROOT" > /dev/null 2>&1
  cmake --build "$TEST_BUILD" -j"${NCPU}" 2>&1 | tail -5
  TEST_LOG="${RESULTS_DIR}/tests.log"
  cd "$TEST_BUILD"
  ctest --extra-verbose --output-on-failure 2>&1 | tee "$TEST_LOG"
  TEST_EXIT=$?
  TESTS_PASSED=$(grep -c '\[  PASSED  \]' "$TEST_LOG" 2>/dev/null || echo 0)
  TESTS_FAILED=$(grep -c '\[  FAILED  \]' "$TEST_LOG" 2>/dev/null || echo 0)
  echo ""
  echo "Tests: ${TESTS_PASSED} passed, ${TESTS_FAILED} failed"
  # Save test summary for report
  echo "tests_passed: ${TESTS_PASSED}" >> "${RESULTS_DIR}/system.txt"
  echo "tests_failed: ${TESTS_FAILED}" >> "${RESULTS_DIR}/system.txt"
  echo "tests_exit: ${TEST_EXIT}" >> "${RESULTS_DIR}/system.txt"
fi

echo ""
echo "Files:"
echo "  ${RESULTS_DIR}/summary.csv      — one-line-per-run aggregates"
echo "  ${RESULTS_DIR}/system.txt       — hardware info"
echo "  ${RESULTS_DIR}/<grid>_<threads>t/perf.csv — per-timestep detail"
