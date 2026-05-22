#!/bin/bash
# Build, run, and plot the AniFstrainCalcitePieri calibration in one shot.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC_DIR="$ROOT/cmake-exec/AniFstrainCalcitePieri"
BINARY="$EXEC_DIR/AniFstrainCalcitePieri"
LOG_DIR="$HERE/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/run_$(date +%Y%m%d_%H%M%S).log"

echo "[$(date '+%F %T')] building target AniFstrainCalcitePieri" | tee -a "$LOG"
cmake --build "$ROOT/cmake-exec" --target AniFstrainCalcitePieri >> "$LOG" 2>&1

echo "[$(date '+%F %T')] running $BINARY" | tee -a "$LOG"
cd "$EXEC_DIR"
# MDOODZ does not auto-create writer_subfolder; pre-create it.
mkdir -p output
"$BINARY" >> "$LOG" 2>&1
NH5=$(ls output/Output*.gzip.h5 2>/dev/null | wc -l | tr -d ' ')
echo "[$(date '+%F %T')] wrote $NH5 HDF5 outputs" | tee -a "$LOG"

echo "[$(date '+%F %T')] plotting MDOODZ vs Pieri+01" | tee -a "$LOG"
cd "$HERE"
python3 mdoodz_vs_pieri.py >> "$LOG" 2>&1
echo "[$(date '+%F %T')] done — see $HERE/img/mdoodz_vs_pieri.png" | tee -a "$LOG"
