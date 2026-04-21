#!/usr/bin/env bash
# Terminate every instance tagged with the current change. Idempotent: safe
# to call repeatedly. Registered as an EXIT/INT/TERM trap on every harness
# entry point.

source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

REGION="$(read_yaml region || true)"
REGION="${REGION:-eu-west-1}"

if [[ -f "${HARNESS_STATE_DIR}/active_instance.env" ]]; then
    # shellcheck disable=SC1091
    source "${HARNESS_STATE_DIR}/active_instance.env"
    REGION="${REGION:-eu-west-1}"
fi

require_aws_cli

# Collect every instance tagged for this change still in a non-terminated state.
instances_json="$(aws_cli ec2 describe-instances \
    --region "${REGION}" \
    --filters \
        "Name=tag:change,Values=${CHANGE_NAME}" \
        "Name=instance-state-name,Values=pending,running,stopping,stopped" \
    --query 'Reservations[].Instances[].InstanceId' \
    --output json 2>/dev/null || echo "[]")"

INSTANCE_IDS=()
while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    INSTANCE_IDS+=("${line}")
done < <(python3 -c 'import json,sys; [print(x) for x in json.loads(sys.stdin.read())]' <<<"${instances_json}")

if [[ ${#INSTANCE_IDS[@]} -eq 0 ]]; then
    log "no instances tagged change=${CHANGE_NAME} in ${REGION}; nothing to tear down"
    rm -f "${HARNESS_STATE_DIR}/active_instance.env"
    exit 0
fi

log "terminating instance(s): ${INSTANCE_IDS[*]}"
aws_cli ec2 terminate-instances \
    --region "${REGION}" \
    --instance-ids "${INSTANCE_IDS[@]}" >/dev/null

# Only block on the waiter in real mode; mocks don't implement `wait`.
if [[ "${MDOODZ_EC2_DRY_RUN}" != "1" ]]; then
    aws_cli ec2 wait instance-terminated \
        --region "${REGION}" \
        --instance-ids "${INSTANCE_IDS[@]}" >/dev/null \
        || warn "wait instance-terminated returned non-zero"
fi

rm -f "${HARNESS_STATE_DIR}/active_instance.env"
log "teardown complete"
