#!/bin/bash
# add-gmg-upleg-fix §4.1 — gmg_levels bisection at 201² on SolVi-Dani 1000×.
#
# For each gmg_levels ∈ {2, 3, 4, 5, 6}, patch the base
# SolViGmgDefaultLevels201.txt fixture on the fly (writer_subfolder +
# gmg_levels + gmg_fgmres_max_restarts cap), run GmgLevelsSweepProbe,
# and dump the full MDOODZ log to lvl_N.log for analysis.
#
# Runtime estimate: ~90 inner FGMRES iters × ~2s per V-cycle = ~3 min/level
# × 5 levels ≈ 15 min on MacBook M1.
#
# After the sweep:
#   grep "GMG-FGMRES" lvl_*.log
# to see convergence vs stall at each level.

set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/cmake-build"
PROBE_BIN="${BUILD_DIR}/TESTS/GmgLevelsSweepProbe"
BASE_TXT="${REPO_ROOT}/TESTS/SolViBenchmark/SolViGmgDefaultLevels201.txt"
SET_DIR="${BUILD_DIR}/TESTS/SolViBenchmark"

if [[ ! -x "${PROBE_BIN}" ]]; then
  echo "probe binary not found at ${PROBE_BIN}" >&2
  echo "build with: cmake --build ${BUILD_DIR} --target GmgLevelsSweepProbe -j8" >&2
  exit 1
fi
if [[ ! -f "${BASE_TXT}" ]]; then
  echo "base fixture not found at ${BASE_TXT}" >&2
  exit 1
fi

# Cap wall time per restart: fewer restarts so stalled runs exit fast.
MAX_RESTARTS="${MAX_RESTARTS:-3}"
RESTART_LEN="${RESTART_LEN:-30}"

for LVL in 2 3 4 5 6; do
  LABEL="D1ProbeA_Lvl${LVL}"
  TMP_TXT="${SET_DIR}/${LABEL}.txt"
  LOG="${SCRIPT_DIR}/lvl_${LVL}.log"

  sed -e "s|writer_subfolder .*|writer_subfolder = ${LABEL}|" \
      -e "s|gmg_levels .*|gmg_levels              = ${LVL}|" \
      -e "s|gmg_fgmres_max_restarts .*|gmg_fgmres_max_restarts = ${MAX_RESTARTS}|" \
      -e "s|gmg_fgmres_restart .*|gmg_fgmres_restart      = ${RESTART_LEN}|" \
      "${BASE_TXT}" > "${TMP_TXT}"

  echo "===== gmg_levels = ${LVL} (log → ${LOG}) ====="
  cd "${BUILD_DIR}/TESTS"
  /usr/bin/time -l "${PROBE_BIN}" "SolViBenchmark/${LABEL}.txt" \
    > "${LOG}" 2>&1 || echo "  (probe exited non-zero at lvl=${LVL}; see ${LOG})"

  grep -E "GMG-FGMRES (config|converged|did not converge|restart)" "${LOG}" | head -15
  echo ""
done

echo "===== summary ====="
for LVL in 2 3 4 5 6; do
  LOG="${SCRIPT_DIR}/lvl_${LVL}.log"
  printf "lvl=%d | " "${LVL}"
  grep -E "GMG-FGMRES (converged|did not converge)" "${LOG}" | head -1 || echo "(no verdict)"
done
