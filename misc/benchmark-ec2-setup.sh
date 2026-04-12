#!/usr/bin/env bash
# =========================================================================
# MDOODZ EC2 Benchmark Setup & Run
#
# This script runs ON the EC2 instance (invoked by benchmark-aws.sh).
# It installs dependencies, builds MDOODZ, validates, and benchmarks.
#
# Usage (on EC2):
#   ./benchmark-ec2-setup.sh [--scenario NAME] [--resolutions "..."] [--threads "..."] [--steps N]
#
# Environment:
#   BENCH_S3_BUCKET  — S3 bucket for result uploads
#   BENCH_REGION     — AWS region (default: eu-central-1)
# =========================================================================
set -euo pipefail

SCENARIO="RiftingComprehensive"
RESOLUTIONS="lowres default medres highres"
THREADS="1 2 4 6 8 12 16"
STEPS="10"
BENCH_S3_BUCKET="${BENCH_S3_BUCKET:-mdoodz-bench-s3-bucket}"
BENCH_REGION="${BENCH_REGION:-eu-central-1}"
REPO_DIR="$HOME/MDOODZ7.0"
TIMEOUT_HOURS=6

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scenario)    SCENARIO="$2"; shift 2 ;;
    --resolutions) RESOLUTIONS="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --steps)       STEPS="$2"; shift 2 ;;
    --timeout)     TIMEOUT_HOURS="$2"; shift 2 ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

HOSTNAME_SHORT="$(hostname -s 2>/dev/null || hostname)"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
S3_PREFIX="s3://${BENCH_S3_BUCKET}/mdoodz-benchmarks/${HOSTNAME_SHORT}/${TIMESTAMP}"

echo "========================================"
echo "  EC2 Benchmark Setup"
echo "========================================"
echo "Scenario:    ${SCENARIO}"
echo "Resolutions: ${RESOLUTIONS}"
echo "Threads:     ${THREADS}"
echo "Steps:       ${STEPS}"
echo "Timeout:     ${TIMEOUT_HOURS}h"
echo "S3 dest:     ${S3_PREFIX}"
echo "========================================"

# ---------- Install dependencies ----------
echo ""
echo ">>> Installing dependencies..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
  build-essential cmake git \
  libsuitesparse-dev libhdf5-dev \
  libblas-dev liblapack-dev \
  > /dev/null 2>&1
echo "    Dependencies installed."

# ---------- Build ----------
echo ""
echo ">>> Building ${SCENARIO}..."
cd "$REPO_DIR"

# Create env.cmake for Ubuntu/GCC (no special paths needed — system packages)
cat > env.cmake <<'ENVEOF'
# EC2 Ubuntu — system packages, no special paths needed
set(OPT ON)
set(OMP ON)
ENVEOF

make clean 2>/dev/null || true
make build SET="${SCENARIO}" 2>&1 | tail -5
echo "    Build complete."

# ---------- Upload helper ----------
upload_results() {
  local results_dir="$1"
  if [[ -d "$results_dir" ]]; then
    echo "    Uploading results to S3..."
    aws s3 sync "$results_dir" "${S3_PREFIX}/" --region "$BENCH_REGION" --quiet 2>/dev/null || \
      echo "    WARNING: S3 upload failed"
  fi
}

# ---------- Run benchmark with timeout ----------
echo ""
echo ">>> Running benchmark (timeout: ${TIMEOUT_HOURS}h)..."

TIMEOUT_SECS=$((TIMEOUT_HOURS * 3600))
RESULTS_DIR=""

# Use a wrapper that captures the results dir path
BENCH_CMD="./misc/benchmark.sh \
  --scenario ${SCENARIO} \
  --resolutions \"${RESOLUTIONS}\" \
  --threads \"${THREADS}\" \
  --steps ${STEPS} \
  --skip-build \
  --validate"

# Run with timeout; capture output to find results dir
BENCH_LOG="/tmp/benchmark-run.log"
set +e
timeout "${TIMEOUT_SECS}" bash -c "$BENCH_CMD" 2>&1 | tee "$BENCH_LOG"
BENCH_EXIT=$?
set -e

# Extract results directory from output
RESULTS_DIR=$(grep "^Results in:" "$BENCH_LOG" | head -1 | sed 's/Results in: //' | tr -d '/' | xargs -I{} echo "{}")
# Fallback: find the latest benchmark-results dir
if [[ -z "$RESULTS_DIR" || ! -d "$RESULTS_DIR" ]]; then
  RESULTS_DIR=$(ls -td "${REPO_DIR}/benchmark-results/"*/ 2>/dev/null | head -1)
fi

if [[ $BENCH_EXIT -eq 124 ]]; then
  echo ""
  echo ">>> Benchmark timed out after ${TIMEOUT_HOURS}h — uploading partial results"
fi

# ---------- Generate report ----------
if [[ -n "$RESULTS_DIR" && -d "$RESULTS_DIR" ]]; then
  echo ""
  echo ">>> Generating report..."
  ./misc/benchmark-report.sh "$RESULTS_DIR" 2>/dev/null || true

  # Upload final results
  upload_results "$RESULTS_DIR"

  echo ""
  echo ">>> Results uploaded to: ${S3_PREFIX}/"
else
  echo ""
  echo ">>> WARNING: No results directory found"
fi

# ---------- Summary ----------
echo ""
echo "========================================"
echo "  EC2 Benchmark Complete"
echo "========================================"
echo "Exit code: ${BENCH_EXIT}"
echo "Results:   ${RESULTS_DIR:-none}"
echo "S3:        ${S3_PREFIX}/"
echo "========================================"
exit $BENCH_EXIT
