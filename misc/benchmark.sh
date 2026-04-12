#!/usr/bin/env bash
# =========================================================================
# MDOODZ Performance Benchmark Script
#
# Runs a scenario at multiple grid sizes and thread counts,
# collecting perf.csv from each run into a results directory.
#
# Usage:
#   ./benchmark.sh                    # default: RiftingBasic, all grids, all threads
#   ./benchmark.sh --grids "150x100 301x201"
#   ./benchmark.sh --threads "1 4 8"
#   ./benchmark.sh --steps 5
#   ./benchmark.sh --quick            # small grids, few steps
#   ./benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default"
#   ./benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default" --validate
#
# Prerequisites:
#   make build SET=<scenario>         # or the script builds it for you
# =========================================================================
set -euo pipefail

# ---------- Configuration (override via flags) ----------
GRIDS="${BENCH_GRIDS:-150x100 301x201 501x501}"
THREADS="${BENCH_THREADS:-1 2 4 8 12 16}"
STEPS="${BENCH_STEPS:-10}"
SCENARIO="${BENCH_SCENARIO:-RiftingBasic}"
RESOLUTIONS=""
VALIDATE=0
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_DIR="${ROOT}/benchmark-results/$(date +%Y%m%d-%H%M%S)"
QUICK=0
SKIP_BUILD=0
RUN_TESTS=0

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --grids)       GRIDS="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --steps)       STEPS="$2"; shift 2 ;;
    --scenario)    SCENARIO="$2"; shift 2 ;;
    --resolutions) RESOLUTIONS="$2"; shift 2 ;;
    --validate)    VALIDATE=1; shift ;;
    --quick)       QUICK=1; shift ;;
    --skip-build)  SKIP_BUILD=1; shift ;;
    --run-tests)   RUN_TESTS=1; shift ;;
    --help|-h)
      head -17 "$0" | grep '^#' | sed 's/^# *//'
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

# ---------- Determine mode: resolution-based or grid-based ----------
MODE="grid"
if [[ -n "$RESOLUTIONS" ]]; then
  MODE="resolution"
fi

echo "========================================"
echo "  MDOODZ Performance Benchmark"
echo "========================================"
echo "Host:     ${HOSTNAME_SHORT} (${OS} ${ARCH})"
echo "CPUs:     ${NCPU}"
echo "RAM:      ${RAM_MB} MB"
echo "Scenario: ${SCENARIO}"
if [[ "$MODE" == "resolution" ]]; then
  echo "Mode:     resolution-based"
  echo "Resols:   ${RESOLUTIONS}"
else
  echo "Mode:     grid-based"
  echo "Grids:    ${GRIDS}"
fi
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

# ---------- Validation pass (--validate) ----------
if [[ $VALIDATE -eq 1 ]]; then
  echo ""
  echo "========================================"
  echo "  Validation Pass (3 steps, 1 thread)"
  echo "========================================"

  if [[ "$MODE" == "resolution" ]]; then
    for res in $RESOLUTIONS; do
      if [[ "$res" == "default" ]]; then
        SRC_TXT="${ROOT}/SETS/${SCENARIO}.txt"
      else
        SRC_TXT="${ROOT}/SETS/${SCENARIO}_${res}.txt"
      fi
      if [[ ! -f "$SRC_TXT" ]]; then
        echo "  ERROR: Resolution file not found: $SRC_TXT"
        exit 1
      fi

      echo "  Validating ${res}..."
      cd "$EXEC_DIR"
      cp "$SRC_TXT" "${SCENARIO}.txt"
      sed -i.bak -E "s/^Nt[[:space:]]*=.*/Nt      = 3/" "${SCENARIO}.txt"
      sed -i.bak -E "s/^t_end[[:space:]]*=.*/t_end   = 0/" "${SCENARIO}.txt"
      sed -i.bak -E "s/^writer[[:space:]]*=.*/writer = 0/" "${SCENARIO}.txt"
      rm -f "${SCENARIO}.txt.bak" perf.csv

      export OMP_NUM_THREADS=1
      if ! "./${SCENARIO}" > /tmp/mdoodz_validate.log 2>&1; then
        echo "  FAILED: ${res} crashed (exit code $?)"
        echo "  Last 20 lines of output:"
        tail -20 /tmp/mdoodz_validate.log
        exit 1
      fi

      # Check for NaN in perf.csv
      if [[ -f perf.csv ]] && grep -qi "nan\|inf" perf.csv; then
        echo "  FAILED: ${res} produced NaN/Inf in perf.csv"
        cat perf.csv
        exit 1
      fi

      NX=$(grep -E "^Nx[[:space:]]" "$SRC_TXT" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
      NZ=$(grep -E "^Nz[[:space:]]" "$SRC_TXT" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
      echo "  OK: ${res} (${NX}x${NZ}) — 3 steps completed"
    done
  else
    for grid in $GRIDS; do
      NX="${grid%x*}"
      NZ="${grid#*x}"
      echo "  Validating ${grid}..."
      cd "$EXEC_DIR"
      cp "${ROOT}/SETS/${SCENARIO}.txt" "${SCENARIO}.txt"
      sed -i.bak -E "s/^Nx[[:space:]]*=.*/Nx      = ${NX}/" "${SCENARIO}.txt"
      sed -i.bak -E "s/^Nz[[:space:]]*=.*/Nz      = ${NZ}/" "${SCENARIO}.txt"
      sed -i.bak -E "s/^Nt[[:space:]]*=.*/Nt      = 3/" "${SCENARIO}.txt"
      sed -i.bak -E "s/^writer[[:space:]]*=.*/writer = 0/" "${SCENARIO}.txt"
      rm -f "${SCENARIO}.txt.bak" perf.csv

      export OMP_NUM_THREADS=1
      if ! "./${SCENARIO}" > /tmp/mdoodz_validate.log 2>&1; then
        echo "  FAILED: ${grid} crashed (exit code $?)"
        tail -20 /tmp/mdoodz_validate.log
        exit 1
      fi
      echo "  OK: ${grid} — 3 steps completed"
    done
  fi

  echo ""
  echo "  All validations passed!"
  echo ""

  # Clean up: remove validation's patched .txt so it doesn't interfere with benchmark runs
  rm -f "${EXEC_DIR}/${SCENARIO}.txt"
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
scenario: ${SCENARIO}
mode: ${MODE}
date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
uname: $(uname -a)
compiler: $(cc --version 2>/dev/null | head -1 || echo unknown)
EOF

# ---------- Helper: patch .txt file (grid mode) ----------
make_bench_txt() {
  local nx="$1" nz="$2" steps="$3" dst="$4"
  cp "${ROOT}/SETS/${SCENARIO}.txt" "$dst"
  # Patch grid size
  sed -i.bak -E "s/^Nx[[:space:]]*=.*/Nx      = ${nx}/" "$dst"
  sed -i.bak -E "s/^Nz[[:space:]]*=.*/Nz      = ${nz}/" "$dst"
  # Patch timesteps
  sed -i.bak -E "s/^Nt[[:space:]]*=.*/Nt      = ${steps}/" "$dst"
  sed -i.bak -E "s/^t_end[[:space:]]*=.*/t_end   = 0/" "$dst"
  # Disable HDF5 output (pure compute benchmark)
  sed -i.bak -E "s/^writer[[:space:]]*=.*/writer = 0/" "$dst"
  # Ensure logging is minimal
  sed -i.bak -E "s/^log_timestamp[[:space:]]*=.*/log_timestamp = 0/" "$dst"
  rm -f "${dst}.bak"
}

# ---------- Helper: prepare resolution .txt file (resolution mode) ----------
make_res_txt() {
  local res="$1" steps="$2" dst="$3"
  if [[ "$res" == "default" ]]; then
    cp "${ROOT}/SETS/${SCENARIO}.txt" "$dst"
  else
    cp "${ROOT}/SETS/${SCENARIO}_${res}.txt" "$dst"
  fi
  # Patch timesteps
  sed -i.bak -E "s/^Nt[[:space:]]*=.*/Nt      = ${steps}/" "$dst"
  sed -i.bak -E "s/^t_end[[:space:]]*=.*/t_end   = 0/" "$dst"
  # Disable HDF5 output
  sed -i.bak -E "s/^writer[[:space:]]*=.*/writer = 0/" "$dst"
  sed -i.bak -E "s/^log_timestamp[[:space:]]*=.*/log_timestamp = 0/" "$dst"
  rm -f "${dst}.bak"
}

# ---------- Build run list ----------
declare -a RUNS_TAG=()
declare -a RUNS_GRID=()
declare -a RUNS_NX=()
declare -a RUNS_NZ=()
declare -a RUNS_THR=()
declare -a RUNS_RES=()

if [[ "$MODE" == "resolution" ]]; then
  for res in $RESOLUTIONS; do
    if [[ "$res" == "default" ]]; then
      SRC="${ROOT}/SETS/${SCENARIO}.txt"
    else
      SRC="${ROOT}/SETS/${SCENARIO}_${res}.txt"
    fi
    NX=$(grep -E "^Nx[[:space:]]" "$SRC" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
    NZ=$(grep -E "^Nz[[:space:]]" "$SRC" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
    for thr in $THREADS; do
      RUNS_TAG+=("${NX}x${NZ}_${thr}t")
      RUNS_GRID+=("${NX}x${NZ}")
      RUNS_NX+=("$NX")
      RUNS_NZ+=("$NZ")
      RUNS_THR+=("$thr")
      RUNS_RES+=("$res")
    done
  done
else
  for grid in $GRIDS; do
    NX="${grid%x*}"
    NZ="${grid#*x}"
    for thr in $THREADS; do
      RUNS_TAG+=("${grid}_${thr}t")
      RUNS_GRID+=("$grid")
      RUNS_NX+=("$NX")
      RUNS_NZ+=("$NZ")
      RUNS_THR+=("$thr")
      RUNS_RES+=("")
    done
  done
fi

TOTAL=${#RUNS_TAG[@]}

# ---------- Summary CSV header ----------
SUMMARY_CSV="${RESULTS_DIR}/summary.csv"
echo "host,os,arch,grid,nx,nz,threads,steps,avg_wall_s,total_wall_s,avg_rheology_s,avg_assembly_s,avg_solve_s,avg_thermal_s,avg_advection_s,avg_melting_s,avg_anisotropy_s,avg_gse_s,avg_output_s,avg_nit,peak_rss_mb" > "$SUMMARY_CSV"

# ---------- Run benchmarks ----------
for ((i=0; i<TOTAL; i++)); do
  TAG="${RUNS_TAG[$i]}"
  NX="${RUNS_NX[$i]}"
  NZ="${RUNS_NZ[$i]}"
  THR="${RUNS_THR[$i]}"
  RES="${RUNS_RES[$i]}"
  GRID="${RUNS_GRID[$i]}"

  # Skip if thread count exceeds available CPUs
  if [[ $THR -gt $NCPU ]]; then
    echo "  [skip] ${TAG} — only ${NCPU} CPUs available"
    continue
  fi

  RUN_DIR="${RESULTS_DIR}/${TAG}"
  mkdir -p "$RUN_DIR"

  echo ""
  echo ">>> [$((i+1))/${TOTAL}] Grid ${NX}x${NZ}, ${THR} threads..."

  # Prepare .txt — copy directly into exec dir (MDOODZ reads hardcoded ${SCENARIO}.txt)
  BENCH_TXT="${RUN_DIR}/bench.txt"
  if [[ -n "$RES" ]]; then
    make_res_txt "$RES" "$STEPS" "$BENCH_TXT"
  else
    make_bench_txt "$NX" "$NZ" "$STEPS" "$BENCH_TXT"
  fi

  # Run — place bench.txt as ${SCENARIO}.txt in exec dir (exe ignores CLI args)
  cd "$EXEC_DIR"
  cp "$BENCH_TXT" "${SCENARIO}.txt"
  rm -f perf.csv
  export OMP_NUM_THREADS="$THR"

  T0=$(date +%s)
  "./${SCENARIO}" > "${RUN_DIR}/stdout.log" 2>&1 || true
  T1=$(date +%s)
  ELAPSED=$((T1 - T0))

  echo "    Completed in ${ELAPSED}s"

  # Collect results
  if [[ -f perf.csv ]]; then
    cp perf.csv "${RUN_DIR}/perf.csv"

    # Compute averages for summary (21-column perf.csv)
    # Columns: 1=step 2=wall 3=time_ma 4=rheology 5=assembly 6=solve
    #          7=thermal 8=advection 9=free_surface 10=reseeding
    #          11=melting 12=anisotropy 13=gse 14=output
    #          15=nit 16=n_particles 17=neq_mom 18=neq_cont
    #          19=peak_rss 20=user_cpu 21=sys_cpu
    AVG=$(awk -F',' 'NR>1 {
      n++; wall+=$2; rheo+=$4; asm+=$5; sol+=$6;
      therm+=$7; adv+=$8; melt+=$11; aniso+=$12; gse+=$13; outp+=$14;
      nit+=$15; rss=$19
    } END {
      if(n>0) printf "%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.1f,%.1f",
        n, wall/n, wall, rheo/n, asm/n, sol/n, therm/n, adv/n, melt/n, aniso/n, gse/n, outp/n, nit/n, rss
      else print "0,0,0,0,0,0,0,0,0,0,0,0,0,0"
    }' perf.csv)

    echo "${HOSTNAME_SHORT},${OS},${ARCH},${GRID},${NX},${NZ},${THR},${AVG}" >> "$SUMMARY_CSV"
    echo "    $(head -1 perf.csv)"
    tail -1 perf.csv | awk -F',' '{printf "    Last step: wall=%.3fs, RSS=%.0fMB, nit=%s\n", $2, $19, $15}'
  else
    echo "    WARNING: No perf.csv produced"
  fi

  # Clean up HDF5 output
  rm -f "${EXEC_DIR}"/Output*.gzip.h5 "${EXEC_DIR}"/Particles*.gzip.h5 "${EXEC_DIR}"/Breakpoint*.dat
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
