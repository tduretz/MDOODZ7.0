#!/usr/bin/env bash
# Shared helpers for the EC2 perf harness.
# Sourced by provision.sh, run_perf_sweep.sh, teardown.sh.

set -euo pipefail

# Allow a dry-run to route AWS calls through mocks/aws without hitting AWS.
: "${MDOODZ_EC2_DRY_RUN:=0}"

HARNESS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS_RESULTS_DIR="${HARNESS_DIR}/results"
HARNESS_STATE_DIR="${HARNESS_RESULTS_DIR}/.state"
HARNESS_RESULTS_CSV="${HARNESS_RESULTS_DIR}/results.csv"
HARNESS_GRIDS_FILE="${HARNESS_DIR}/grids.yaml"
HARNESS_MOCKS_DIR="${HARNESS_DIR}/mocks"

CHANGE_NAME="add-gmg-stokes-defence"
PROJECT_TAG="mdoodz-gmg"

mkdir -p "${HARNESS_RESULTS_DIR}" "${HARNESS_STATE_DIR}"

log()  { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die()  { log "FATAL: $*"; exit 1; }
warn() { log "WARN: $*"; }

# Prepend mocks/aws to PATH during dry-run so `aws ...` calls are intercepted.
aws_cli() {
    if [[ "${MDOODZ_EC2_DRY_RUN}" == "1" ]]; then
        PATH="${HARNESS_MOCKS_DIR}:${PATH}" aws "$@"
    else
        command aws "$@"
    fi
}

# Pull a scalar field from grids.yaml via python (avoids yq as a dependency).
read_yaml() {
    local field="$1"
    python3 - "$field" <<'PY'
import sys, os, re
field = sys.argv[1]
with open(os.environ["HARNESS_GRIDS_FILE"], "r") as f:
    txt = f.read()
# Tiny hand-rolled parser: first "key: value" match at any indent.
for line in txt.splitlines():
    m = re.match(r"\s*" + re.escape(field) + r":\s*(.+?)\s*$", line)
    if m:
        v = m.group(1)
        # strip inline comments + surrounding quotes
        v = re.sub(r"\s+#.*$", "", v).strip().strip('"').strip("'")
        print(v)
        break
PY
}

current_git_sha() {
    (cd "${HARNESS_DIR}/../.." && git rev-parse HEAD)
}

# On-demand instance price per hour for the known pinned types (USD).
# Prices are a static reference table used by the budget cap; they don't need
# to be exact to the cent. Update as needed.
instance_price_usd_per_hour() {
    case "$1" in
        c7i.8xlarge)  echo "1.4280" ;;
        c7g.8xlarge)  echo "1.1560" ;;
        m7i.8xlarge)  echo "1.6128" ;;
        *)            echo "2.0000" ;;  # conservative default
    esac
}

# Ensure the AWS CLI (real or mock) is available.
require_aws_cli() {
    if [[ "${MDOODZ_EC2_DRY_RUN}" == "1" ]]; then
        [[ -x "${HARNESS_MOCKS_DIR}/aws" ]] \
            || die "dry-run requested but mocks/aws is not executable"
    else
        command -v aws >/dev/null 2>&1 \
            || die "aws CLI not found on PATH; install it or set MDOODZ_EC2_DRY_RUN=1"
    fi
}

export HARNESS_GRIDS_FILE
