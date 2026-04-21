#!/usr/bin/env bash
# Launch a pinned c7i.8xlarge in eu-west-1, tag it with the current change +
# git SHA, record the instance ID so teardown.sh can find it later. Refuses
# to launch if an instance with the same change tag is already alive (D8 in
# gmg-performance-harness spec: orphan detection).

source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

INSTANCE_TYPE="$(read_yaml instance_type || true)"
REGION="$(read_yaml region || true)"
INSTANCE_TYPE="${INSTANCE_TYPE:-c7i.8xlarge}"
REGION="${REGION:-eu-west-1}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --instance-type)
            shift
            INSTANCE_TYPE="$1"
            warn "instance-type override: ${INSTANCE_TYPE} (pinned default is c7i.8xlarge)"
            ;;
        --region)
            shift
            REGION="$1"
            ;;
        -h|--help)
            sed -n '1,8p' "$0"; exit 0
            ;;
        *)
            die "unknown argument: $1"
            ;;
    esac
    shift
done

require_aws_cli

log "provisioning ${INSTANCE_TYPE} in ${REGION} for change=${CHANGE_NAME}"

# 1. Orphan detection: refuse to launch if an instance with the same change
#    tag is already alive (states running, pending, stopping, stopped).
orphan_json="$(aws_cli ec2 describe-instances \
    --region "${REGION}" \
    --filters \
        "Name=tag:change,Values=${CHANGE_NAME}" \
        "Name=instance-state-name,Values=pending,running,stopping,stopped" \
    --query 'Reservations[].Instances[].InstanceId' \
    --output json 2>/dev/null || echo "[]")"

if [[ -n "${orphan_json}" ]] && ! [[ "${orphan_json}" =~ ^\[[[:space:]]*\]$ ]]; then
    warn "found orphan instance(s) tagged change=${CHANGE_NAME}:"
    echo "${orphan_json}" >&2
    warn "aborting launch. Run: ./teardown.sh  (or manually terminate the orphan)"
    exit 2
fi

# 2. Resolve the SSH key and AMI.
KEY_NAME="${MDOODZ_EC2_KEY_NAME:-mdoodz-ec2}"
# Amazon Linux 2023 x86_64 — resolved via SSM public parameter.
AMI_ID="$(aws_cli ssm get-parameter \
    --region "${REGION}" \
    --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
    --query 'Parameter.Value' \
    --output text 2>/dev/null || echo "ami-placeholder")"

GIT_SHA="$(current_git_sha)"

# 3. Launch.
launch_json="$(aws_cli ec2 run-instances \
    --region "${REGION}" \
    --image-id "${AMI_ID}" \
    --instance-type "${INSTANCE_TYPE}" \
    --key-name "${KEY_NAME}" \
    --count 1 \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=mdoodz-perf-${CHANGE_NAME}},{Key=project,Value=${PROJECT_TAG}},{Key=change,Value=${CHANGE_NAME}},{Key=git_sha,Value=${GIT_SHA}},{Key=instance_type,Value=${INSTANCE_TYPE}}]" \
    --output json)"

INSTANCE_ID="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["Instances"][0]["InstanceId"])' <<<"${launch_json}")"
[[ -n "${INSTANCE_ID}" && "${INSTANCE_ID}" != "None" ]] || die "failed to parse InstanceId from run-instances output"

log "launched ${INSTANCE_ID} (${INSTANCE_TYPE}, git_sha=${GIT_SHA})"

# 4. Persist state for the rest of the harness.
{
    echo "INSTANCE_ID=${INSTANCE_ID}"
    echo "INSTANCE_TYPE=${INSTANCE_TYPE}"
    echo "REGION=${REGION}"
    echo "GIT_SHA=${GIT_SHA}"
    echo "LAUNCHED_AT=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "PRICE_PER_HOUR_USD=$(instance_price_usd_per_hour "${INSTANCE_TYPE}")"
} >"${HARNESS_STATE_DIR}/active_instance.env"

# 5. Wait for it to be running + reachable, then print DNS.
aws_cli ec2 wait instance-running \
    --region "${REGION}" \
    --instance-ids "${INSTANCE_ID}" >/dev/null 2>&1 || warn "wait instance-running failed (may be a mock)"

PUBLIC_DNS="$(aws_cli ec2 describe-instances \
    --region "${REGION}" \
    --instance-ids "${INSTANCE_ID}" \
    --query 'Reservations[0].Instances[0].PublicDnsName' \
    --output text 2>/dev/null || echo "unknown")"

printf 'INSTANCE_ID=%s\n' "${INSTANCE_ID}"
printf 'PUBLIC_DNS=%s\n' "${PUBLIC_DNS}"
log "instance is up: ${INSTANCE_ID} @ ${PUBLIC_DNS}"
