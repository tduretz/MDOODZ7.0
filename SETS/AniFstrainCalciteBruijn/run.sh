#!/bin/bash
# Build, run, and plot the AniFstrainCalciteBruijn calibration in one shot.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC_DIR="$ROOT/cmake-exec/AniFstrainCalciteBruijn"
BINARY="$EXEC_DIR/AniFstrainCalciteBruijn"
LOG_DIR="$HERE/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/run_$(date +%Y%m%d_%H%M%S).log"

echo "[$(date '+%F %T')] building target AniFstrainCalciteBruijn" | tee -a "$LOG"
cmake --build "$ROOT/cmake-exec" --target AniFstrainCalciteBruijn >> "$LOG" 2>&1

echo "[$(date '+%F %T')] running $BINARY" | tee -a "$LOG"
cd "$EXEC_DIR"
mkdir -p output
"$BINARY" >> "$LOG" 2>&1
NH5=$(ls output/Output*.gzip.h5 2>/dev/null | wc -l | tr -d ' ')
echo "[$(date '+%F %T')] wrote $NH5 HDF5 outputs" | tee -a "$LOG"

echo "[$(date '+%F %T')] plotting MDOODZ vs Bruijn+11" | tee -a "$LOG"
cd "$HERE"
python3 mdoodz_vs_bruijn.py >> "$LOG" 2>&1
echo "[$(date '+%F %T')] done — see $HERE/img/mdoodz_vs_bruijn.png" | tee -a "$LOG"
