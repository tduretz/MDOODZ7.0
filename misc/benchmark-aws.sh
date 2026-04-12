#!/usr/bin/env bash
# =========================================================================
# MDOODZ AWS Benchmark Orchestrator
#
# Manages the full lifecycle: start EC2 → rsync → build & benchmark → download → stop
#
# Usage:
#   ./benchmark-aws.sh
#   ./benchmark-aws.sh --scenario RiftingComprehensive --resolutions "lowres default medres highres"
#   ./benchmark-aws.sh --instance-id i-xxx --ssh-key ~/.ssh/k.pem --s3-bucket mybucket
#
# Environment (override via flags):
#   BENCH_EC2_INSTANCE  — EC2 instance ID
#   BENCH_SSH_KEY       — Path to SSH private key
#   BENCH_S3_BUCKET     — S3 bucket name
#   BENCH_REGION        — AWS region
# =========================================================================
set -euo pipefail

# ---------- Defaults ----------
BENCH_EC2_INSTANCE="${BENCH_EC2_INSTANCE:-i-0df5754e400cb3ec1}"
BENCH_SSH_KEY="${BENCH_SSH_KEY:-$HOME/.ssh/key-mdoodz.pem}"
BENCH_S3_BUCKET="${BENCH_S3_BUCKET:-mdoodz-bench-s3-bucket}"
BENCH_SSH_USER="${BENCH_SSH_USER:-ubuntu}"
BENCH_REGION="${BENCH_REGION:-eu-central-1}"
SCENARIO="RiftingComprehensive"
RESOLUTIONS="lowres default medres highres"
THREADS="1 2 4 6 8 12 16"
STEPS="10"
TIMEOUT_HOURS="6"
SKIP_PREFLIGHT=0

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --instance-id)    BENCH_EC2_INSTANCE="$2"; shift 2 ;;
    --ssh-key)        BENCH_SSH_KEY="$2"; shift 2 ;;
    --s3-bucket)      BENCH_S3_BUCKET="$2"; shift 2 ;;
    --region)         BENCH_REGION="$2"; shift 2 ;;
    --scenario)       SCENARIO="$2"; shift 2 ;;
    --resolutions)    RESOLUTIONS="$2"; shift 2 ;;
    --threads)        THREADS="$2"; shift 2 ;;
    --steps)          STEPS="$2"; shift 2 ;;
    --timeout)        TIMEOUT_HOURS="$2"; shift 2 ;;
    --skip-preflight) SKIP_PREFLIGHT=1; shift ;;
    --help|-h)
      head -16 "$0" | grep '^#' | sed 's/^# *//'
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

export BENCH_EC2_INSTANCE BENCH_SSH_KEY BENCH_S3_BUCKET BENCH_SSH_USER BENCH_REGION

EXPANDED_KEY="${BENCH_SSH_KEY/#\~/$HOME}"
SSH_OPTS="-i ${EXPANDED_KEY} -o StrictHostKeyChecking=accept-new -o ConnectTimeout=15 -o ServerAliveInterval=30 -o ServerAliveCountMax=5"
INSTANCE_STARTED=0

echo "========================================"
echo "  MDOODZ AWS Benchmark"
echo "========================================"
echo "Instance: ${BENCH_EC2_INSTANCE}"
echo "Region:   ${BENCH_REGION}"
echo "Scenario: ${SCENARIO}"
echo "Resols:   ${RESOLUTIONS}"
echo "Threads:  ${THREADS}"
echo "Steps:    ${STEPS}"
echo "Timeout:  ${TIMEOUT_HOURS}h"
echo "S3:       s3://${BENCH_S3_BUCKET}/"
echo "========================================"

# ---------- Trap: stop instance only on explicit request ----------
STOP_ON_EXIT=0
cleanup() {
  if [[ $STOP_ON_EXIT -eq 1 && $INSTANCE_STARTED -eq 1 ]]; then
    echo ""
    echo ">>> Stopping instance ${BENCH_EC2_INSTANCE}..."
    aws ec2 stop-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
      --output text 2>/dev/null || echo "WARNING: Failed to stop instance"
    echo "    Instance stop requested."
  elif [[ $INSTANCE_STARTED -eq 1 ]]; then
    echo ""
    echo ">>> Instance ${BENCH_EC2_INSTANCE} left running (benchmark may still be active)."
    echo "    Stop manually: aws ec2 stop-instances --region ${BENCH_REGION} --instance-ids ${BENCH_EC2_INSTANCE}"
  fi
}
trap cleanup EXIT INT TERM

# ---------- Preflight ----------
if [[ $SKIP_PREFLIGHT -eq 0 ]]; then
  echo ""
  echo ">>> Running preflight checks..."
  if ! "${ROOT}/misc/benchmark-preflight.sh"; then
    echo ""
    echo "Preflight failed. Fix the issues above before running the benchmark."
    INSTANCE_STARTED=0  # don't try to stop in cleanup
    exit 1
  fi
fi

# ---------- Start instance ----------
echo ""
echo ">>> Starting instance..."
STATE=$(aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
  --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null)

if [[ "$STATE" == "running" ]]; then
  echo "    Instance already running."
  INSTANCE_STARTED=1
elif [[ "$STATE" == "stopped" ]]; then
  aws ec2 start-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" --output text >/dev/null
  INSTANCE_STARTED=1
  echo "    Start requested. Waiting for running state..."
  aws ec2 wait instance-running --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE"
  echo "    Instance running."
  # Short delay for SSH daemon to initialise
  sleep 10
else
  echo "ERROR: Instance is in state '${STATE}' — cannot start. Expected 'stopped' or 'running'."
  INSTANCE_STARTED=0
  exit 1
fi

# Get public IP (may change on restart)
PUB_IP=$(aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
  --query 'Reservations[0].Instances[0].PublicIpAddress' --output text 2>/dev/null)

if [[ -z "$PUB_IP" || "$PUB_IP" == "None" ]]; then
  echo "ERROR: No public IP assigned to instance."
  exit 1
fi
echo "    Public IP: ${PUB_IP}"

# ---------- Wait for SSH ----------
echo ""
echo ">>> Waiting for SSH..."
for i in $(seq 1 30); do
  if ssh $SSH_OPTS -o BatchMode=yes "${BENCH_SSH_USER}@${PUB_IP}" "echo ok" &>/dev/null; then
    echo "    SSH ready."
    break
  fi
  if [[ $i -eq 30 ]]; then
    echo "ERROR: SSH not available after 5 minutes. Check security group."
    exit 1
  fi
  echo "    Attempt $i/30..."
  sleep 10
done

# ---------- Rsync repo ----------
echo ""
echo ">>> Syncing repository to instance..."
rsync -az --delete \
  --exclude '.git' \
  --exclude 'cmake-build*' \
  --exclude 'cmake-exec' \
  --exclude 'benchmark-results' \
  --exclude 'visualtests-out' \
  --exclude 'env.cmake' \
  --exclude '*.h5' \
  --exclude '*.dat' \
  -e "ssh ${SSH_OPTS}" \
  "${ROOT}/" "${BENCH_SSH_USER}@${PUB_IP}:~/MDOODZ7.0/"
echo "    Sync complete."

# ---------- Run remote benchmark ----------
echo ""
echo ">>> Starting remote benchmark..."
echo "    This may take several hours. Instance will be stopped automatically on completion."
echo ""

# ---------- Run remote benchmark (detached — survives SSH drops) ----------
echo ""
echo ">>> Starting remote benchmark (detached)..."
echo "    Benchmark will continue even if this terminal disconnects."
echo ""

REMOTE_LOG="/tmp/mdoodz-bench-remote.log"
REMOTE_PID_FILE="/tmp/mdoodz-bench.pid"

ssh $SSH_OPTS "${BENCH_SSH_USER}@${PUB_IP}" \
  "export BENCH_S3_BUCKET='${BENCH_S3_BUCKET}' BENCH_REGION='${BENCH_REGION}'; \
   cd ~/MDOODZ7.0 && chmod +x misc/*.sh && \
   nohup ./misc/benchmark-ec2-setup.sh \
     --scenario '${SCENARIO}' \
     --resolutions '${RESOLUTIONS}' \
     --threads '${THREADS}' \
     --steps '${STEPS}' \
     --timeout '${TIMEOUT_HOURS}' \
     > ${REMOTE_LOG} 2>&1 & echo \$! > ${REMOTE_PID_FILE}"

echo "    Remote benchmark launched. Polling for completion..."

# ---------- Poll for completion ----------
POLL_INTERVAL=60
POLL_MAX=$((TIMEOUT_HOURS * 3600 / POLL_INTERVAL + 10))
for ((p=1; p<=POLL_MAX; p++)); do
  sleep "$POLL_INTERVAL"

  # Check if remote process is still running
  STILL_RUNNING=$(ssh $SSH_OPTS -o ConnectTimeout=10 "${BENCH_SSH_USER}@${PUB_IP}" \
    "cat ${REMOTE_PID_FILE} 2>/dev/null | xargs -I{} kill -0 {} 2>/dev/null && echo yes || echo no" 2>/dev/null) || {
    echo "    [poll $p] SSH failed — retrying next cycle (benchmark continues on EC2)"
    continue
  }

  # Show progress
  PROGRESS=$(ssh $SSH_OPTS -o ConnectTimeout=10 "${BENCH_SSH_USER}@${PUB_IP}" \
    "wc -l ~/MDOODZ7.0/benchmark-results/*/summary.csv 2>/dev/null | tail -1 | awk '{print \$1}'" 2>/dev/null) || PROGRESS="?"
  echo "    [poll $p] running=${STILL_RUNNING} summary_rows=${PROGRESS}"

  if [[ "$STILL_RUNNING" == "no" ]]; then
    echo "    Remote benchmark finished."
    # Grab exit code
    REMOTE_EXIT=$(ssh $SSH_OPTS "${BENCH_SSH_USER}@${PUB_IP}" \
      "wait \$(cat ${REMOTE_PID_FILE}) 2>/dev/null; echo \$?" 2>/dev/null) || REMOTE_EXIT=1
    break
  fi
done
REMOTE_EXIT=${REMOTE_EXIT:-1}

# ---------- Download results ----------
echo ""
echo ">>> Downloading results..."
mkdir -p "${ROOT}/benchmark-results"
REMOTE_LATEST=$(ssh $SSH_OPTS "${BENCH_SSH_USER}@${PUB_IP}" \
  "ls -td ~/MDOODZ7.0/benchmark-results/*/ 2>/dev/null | head -1" || true)

if [[ -n "$REMOTE_LATEST" ]]; then
  LOCAL_DIR="${ROOT}/benchmark-results/$(basename "$REMOTE_LATEST")"
  rsync -az -e "ssh ${SSH_OPTS}" \
    "${BENCH_SSH_USER}@${PUB_IP}:${REMOTE_LATEST}" "$LOCAL_DIR"
  echo "    Downloaded to: ${LOCAL_DIR}"

  # Generate local report
  if [[ -d "$LOCAL_DIR" ]]; then
    "${ROOT}/misc/benchmark-report.sh" "$LOCAL_DIR" 2>/dev/null || true
  fi
else
  echo "    No results to download — check S3 for partial results."
fi

# ---------- Done — stop instance now ----------
echo ""
echo "========================================"
echo "  AWS Benchmark Complete"
echo "========================================"
echo "Remote exit: ${REMOTE_EXIT}"
echo "Local results: ${LOCAL_DIR:-none}"
echo "S3: s3://${BENCH_S3_BUCKET}/mdoodz-benchmarks/"
echo "Instance will be stopped automatically."
echo "========================================"
STOP_ON_EXIT=1
