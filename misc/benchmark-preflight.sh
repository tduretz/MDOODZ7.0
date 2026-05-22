#!/usr/bin/env bash
# =========================================================================
# MDOODZ Benchmark Preflight — verify AWS environment before benchmarking
#
# Usage:
#   ./benchmark-preflight.sh
#   BENCH_EC2_INSTANCE=i-xxx BENCH_SSH_KEY=~/.ssh/k.pem BENCH_S3_BUCKET=b ./benchmark-preflight.sh
#
# Checks:
#   1. AWS CLI installed and configured
#   2. Required env vars set
#   3. EC2 instance exists
#   4. SSH key exists with correct permissions
#   5. SSH connectivity (instance must be running)
#   6. S3 bucket exists and is writable
# =========================================================================
set -euo pipefail

# ---------- Defaults (override via env) ----------
BENCH_EC2_INSTANCE="${BENCH_EC2_INSTANCE:-i-0df5754e400cb3ec1}"
BENCH_SSH_KEY="${BENCH_SSH_KEY:-$HOME/.ssh/key-mdoodz.pem}"
BENCH_S3_BUCKET="${BENCH_S3_BUCKET:-mdoodz-bench-s3-bucket}"
BENCH_SSH_USER="${BENCH_SSH_USER:-ubuntu}"
BENCH_REGION="${BENCH_REGION:-eu-central-1}"

PASS=0
FAIL=0
WARN=0

check_pass() { echo "  ✓ $1"; PASS=$((PASS+1)); }
check_fail() { echo "  ✗ $1"; FAIL=$((FAIL+1)); }
check_warn() { echo "  ⚠ $1"; WARN=$((WARN+1)); }

echo "========================================"
echo "  MDOODZ Benchmark Preflight"
echo "========================================"
echo ""

# --- 1. AWS CLI ---
echo "1. AWS CLI"
if command -v aws &>/dev/null; then
  check_pass "aws CLI found: $(aws --version 2>&1 | head -1)"
else
  check_fail "aws CLI not found — install: brew install awscli (macOS) or apt install awscli (Linux)"
fi

if aws sts get-caller-identity &>/dev/null; then
  CALLER=$(aws sts get-caller-identity --output text --query 'Account' 2>/dev/null)
  check_pass "AWS credentials configured (account: ${CALLER})"
else
  check_fail "AWS credentials not configured — run: aws configure"
fi
echo ""

# --- 2. Required env vars ---
echo "2. Environment variables"
echo "   BENCH_EC2_INSTANCE = ${BENCH_EC2_INSTANCE}"
echo "   BENCH_SSH_KEY      = ${BENCH_SSH_KEY}"
echo "   BENCH_S3_BUCKET    = ${BENCH_S3_BUCKET}"
echo "   BENCH_SSH_USER     = ${BENCH_SSH_USER}"
echo "   BENCH_REGION       = ${BENCH_REGION}"

[[ -n "$BENCH_EC2_INSTANCE" ]] && check_pass "BENCH_EC2_INSTANCE set" || check_fail "BENCH_EC2_INSTANCE not set"
[[ -n "$BENCH_SSH_KEY" ]]      && check_pass "BENCH_SSH_KEY set"      || check_fail "BENCH_SSH_KEY not set"
[[ -n "$BENCH_S3_BUCKET" ]]    && check_pass "BENCH_S3_BUCKET set"    || check_fail "BENCH_S3_BUCKET not set"
echo ""

# --- 3. EC2 instance ---
echo "3. EC2 instance"
if aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" &>/dev/null; then
  STATE=$(aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
    --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null)
  ITYPE=$(aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
    --query 'Reservations[0].Instances[0].InstanceType' --output text 2>/dev/null)
  check_pass "Instance found: ${BENCH_EC2_INSTANCE} (${ITYPE}, ${STATE})"
  if [[ "$STATE" != "running" ]]; then
    check_warn "Instance is ${STATE} — SSH check will be skipped. Start it for full preflight."
  fi
else
  check_fail "Instance not found: ${BENCH_EC2_INSTANCE} — check instance ID and region (${BENCH_REGION})"
fi
echo ""

# --- 4. SSH key ---
echo "4. SSH key"
EXPANDED_KEY="${BENCH_SSH_KEY/#\~/$HOME}"
if [[ -f "$EXPANDED_KEY" ]]; then
  check_pass "Key file exists: ${EXPANDED_KEY}"
  PERMS=$(stat -f '%Lp' "$EXPANDED_KEY" 2>/dev/null || stat -c '%a' "$EXPANDED_KEY" 2>/dev/null || echo "???")
  if [[ "$PERMS" == "600" || "$PERMS" == "400" ]]; then
    check_pass "Key permissions OK (${PERMS})"
  else
    check_warn "Key permissions are ${PERMS} — should be 600 or 400. Fix: chmod 600 ${EXPANDED_KEY}"
  fi
else
  check_fail "Key file not found: ${EXPANDED_KEY}"
fi
echo ""

# --- 5. SSH connectivity ---
echo "5. SSH connectivity"
if [[ "${STATE:-}" == "running" ]]; then
  PUB_IP=$(aws ec2 describe-instances --region "$BENCH_REGION" --instance-ids "$BENCH_EC2_INSTANCE" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' --output text 2>/dev/null)
  if [[ -n "$PUB_IP" && "$PUB_IP" != "None" ]]; then
    check_pass "Public IP: ${PUB_IP}"
    if ssh -i "$EXPANDED_KEY" -o ConnectTimeout=10 -o StrictHostKeyChecking=accept-new -o BatchMode=yes \
         "${BENCH_SSH_USER}@${PUB_IP}" "echo ok" &>/dev/null; then
      check_pass "SSH connection successful"
      REMOTE_CPUS=$(ssh -i "$EXPANDED_KEY" -o ConnectTimeout=10 "${BENCH_SSH_USER}@${PUB_IP}" "nproc" 2>/dev/null || echo "?")
      REMOTE_RAM=$(ssh -i "$EXPANDED_KEY" -o ConnectTimeout=10 "${BENCH_SSH_USER}@${PUB_IP}" \
        "awk '/MemTotal/ {printf \"%d\", \$2/1024}' /proc/meminfo" 2>/dev/null || echo "?")
      check_pass "Remote: ${REMOTE_CPUS} CPUs, ${REMOTE_RAM} MB RAM"
    else
      check_fail "SSH connection failed — check security group allows port 22 from your IP"
    fi
  else
    check_fail "No public IP — instance may be in a private subnet"
  fi
else
  check_warn "Skipping SSH check — instance not running"
fi
echo ""

# --- 6. S3 bucket ---
echo "6. S3 bucket"
if aws s3 ls "s3://${BENCH_S3_BUCKET}" --region "$BENCH_REGION" &>/dev/null; then
  check_pass "Bucket exists: s3://${BENCH_S3_BUCKET}"
  # Test write
  TEST_KEY="mdoodz-benchmarks/.preflight-test-$(date +%s)"
  if echo "preflight" | aws s3 cp - "s3://${BENCH_S3_BUCKET}/${TEST_KEY}" --region "$BENCH_REGION" &>/dev/null; then
    check_pass "Bucket is writable"
    aws s3 rm "s3://${BENCH_S3_BUCKET}/${TEST_KEY}" --region "$BENCH_REGION" &>/dev/null || true
  else
    check_fail "Bucket not writable — check IAM policy (s3:PutObject)"
  fi
else
  check_fail "Bucket not found: s3://${BENCH_S3_BUCKET} — create it or check region (${BENCH_REGION})"
fi
echo ""

# --- Summary ---
echo "========================================"
echo "  Preflight: ${PASS} passed, ${FAIL} failed, ${WARN} warnings"
echo "========================================"

if [[ $FAIL -gt 0 ]]; then
  echo ""
  echo "Fix the failures above before running benchmark-aws.sh"
  exit 1
fi

if [[ $WARN -gt 0 ]]; then
  echo ""
  echo "Warnings present — benchmark may still work but review them."
  exit 0
fi

echo ""
echo "All checks passed — ready to benchmark!"
exit 0
