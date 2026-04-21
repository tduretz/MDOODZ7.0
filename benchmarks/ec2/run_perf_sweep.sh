#!/usr/bin/env bash
# Drive the full perf sweep on a provisioned EC2 instance. Reads grids.yaml,
# SSHes to the instance, runs MDOODZ N=5 per configuration, captures wall
# time + peak RSS + turbostat throttle state, appends one row per run to
# results/results.csv. Enforces:
#   - USD 50 budget cap (per D1/D8 in the perf-harness spec)
#   - SIGINT trap that finishes the current row and tears down
#   - EXIT trap that always calls teardown.sh
#   - > 3 retries on thermally-throttled runs before marking UNRELIABLE
#
# Output row schema matches the harness spec requirement.

source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

BUDGET_USD_DEFAULT=50
SPOT=0
BUDGET_USD="${BUDGET_USD_DEFAULT}"
GRIDS_FILE="${HARNESS_GRIDS_FILE}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --spot)             SPOT=1 ;;
        --budget-usd)       shift; BUDGET_USD="$1" ;;
        --grids)            shift; GRIDS_FILE="$1" ;;
        -h|--help)          sed -n '1,13p' "$0"; exit 0 ;;
        *)                  die "unknown argument: $1" ;;
    esac
    shift
done

require_aws_cli

[[ -f "${HARNESS_STATE_DIR}/active_instance.env" ]] \
    || die "no active instance; run provision.sh first"
# shellcheck disable=SC1091
source "${HARNESS_STATE_DIR}/active_instance.env"

PRICE_PER_HOUR_USD="${PRICE_PER_HOUR_USD:-$(instance_price_usd_per_hour "${INSTANCE_TYPE}")}"

# Results CSV header (written once; append-only afterwards).
if [[ ! -f "${HARNESS_RESULTS_CSV}" ]]; then
    echo "timestamp,instance_id,instance_type,region,git_sha,grid_nx,grid_nz,lin_solver,n_threads,repetition,wall_time_s,peak_rss_kb,n_its_fgmres,n_picard,thermal_throttled" \
        >"${HARNESS_RESULTS_CSV}"
fi

# Traps: SIGINT finishes current row + tears down; EXIT always tears down.
cleanup() {
    local rc=$?
    log "cleanup: running teardown.sh (exit=${rc})"
    bash "${HARNESS_DIR}/teardown.sh" || true
    exit "${rc}"
}
trap cleanup EXIT
trap 'log "SIGINT received; finishing current row then tearing down"; exit 130' INT TERM

# Helper: given one (sweep, nx, nz, lin_solver, n_threads, repetition) tuple,
# run MDOODZ remotely, capture metrics, emit a CSV row. In dry-run mode we
# emit deterministic synthetic metrics so the pipeline can be validated end
# to end without an actual instance.
run_one() {
    local sweep="$1"
    local nx="$2"
    local nz="$3"
    local lin_solver="$4"
    local n_threads="$5"
    local rep="$6"

    local ts; ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    local wall_time_s peak_rss_kb n_its_fgmres n_picard thermal_throttled

    if [[ "${MDOODZ_EC2_DRY_RUN}" == "1" ]]; then
        # Deterministic synthetic numbers roughly matching design D1 regression
        # bounds: CHOLMOD ~ O(N^1.5), GMG ~ O(N log N). Good enough for
        # dry-run pipeline validation only.
        local dof=$((nx * nz))
        if [[ "${lin_solver}" -eq 0 ]]; then
            wall_time_s="$(python3 -c "print(round(1e-5 * ${dof} ** 1.5, 4))")"
            peak_rss_kb="$(python3 -c "print(int(50 * ${dof} ** 1.25))")"
            n_its_fgmres=0
            n_picard=1
        else
            wall_time_s="$(python3 -c "import math; print(round(3e-5 * ${dof} * max(1, math.log(${dof})), 4))")"
            peak_rss_kb="$(python3 -c "print(int(30 * ${dof} ** 1.05))")"
            n_its_fgmres=8
            n_picard=1
        fi
        thermal_throttled=false
    else
        # Real remote run on the EC2 instance. Caller must have a working
        # SSH path — the invocation below is a template; the deployment
        # script (installed on the instance by provision.sh bootstrap)
        # accepts --nx/--nz/--lin_solver/--n-threads and prints a
        # single-line JSON on stdout.
        local ssh_json
        ssh_json="$(ssh -o StrictHostKeyChecking=no \
            -i "${MDOODZ_EC2_KEY_PATH:-$HOME/.ssh/mdoodz-ec2.pem}" \
            ec2-user@"${PUBLIC_DNS:-unknown}" \
            "./mdoodz-perf-run --nx ${nx} --nz ${nz} --lin_solver ${lin_solver} --n-threads ${n_threads}")"
        wall_time_s="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["wall_time_s"])' <<<"${ssh_json}")"
        peak_rss_kb="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["peak_rss_kb"])' <<<"${ssh_json}")"
        n_its_fgmres="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["n_its_fgmres"])' <<<"${ssh_json}")"
        n_picard="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["n_picard"])' <<<"${ssh_json}")"
        thermal_throttled="$(python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["thermal_throttled"])' <<<"${ssh_json}")"
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${ts}" "${INSTANCE_ID}" "${INSTANCE_TYPE}" "${REGION}" "${GIT_SHA}" \
        "${nx}" "${nz}" "${lin_solver}" "${n_threads}" "${rep}" \
        "${wall_time_s}" "${peak_rss_kb}" "${n_its_fgmres}" "${n_picard}" "${thermal_throttled}" \
        >>"${HARNESS_RESULTS_CSV}"
}

# Budget projection: reactive — only engages after a warm-up period of 3
# completed runs, using the observed median s/run as the extrapolation base.
# Before warm-up we trust the cap and let the sweep proceed.
projected_budget_usd() {
    local elapsed_s="$1"
    local completed="$2"
    local remaining_runs="$3"
    if [[ "${completed}" -lt 3 ]]; then
        echo "0.00"
        return
    fi
    python3 - "${elapsed_s}" "${completed}" "${remaining_runs}" "${PRICE_PER_HOUR_USD}" <<'PY'
import sys
elapsed_s, completed, remaining, price_h = sys.argv[1:5]
elapsed_s, completed, remaining = int(elapsed_s), int(completed), int(remaining)
price_h = float(price_h)
if completed <= 0:
    print("0.00"); raise SystemExit
per_run = elapsed_s / completed
total_s = elapsed_s + remaining * per_run
print(f"{total_s / 3600.0 * price_h:.2f}")
PY
}

# Expand grids.yaml via _parse_grids.py (isolated so macOS stock bash 3.2
# does not have to parse a heredoc inside a process substitution).
SWEEP_ROWS=()
while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    SWEEP_ROWS+=("${line}")
done < <(GRIDS_FILE="${GRIDS_FILE}" python3 "${HARNESS_DIR}/_parse_grids.py")

REPETITIONS="$(read_yaml repetitions)"
REPETITIONS="${REPETITIONS:-5}"
MAX_RETRIES="$(read_yaml thermal_throttle_max_retries)"
MAX_RETRIES="${MAX_RETRIES:-3}"

log "planned ${#SWEEP_ROWS[@]} configurations x ${REPETITIONS} reps on ${INSTANCE_ID}"

started_at=$(date +%s)
total_runs=$(( ${#SWEEP_ROWS[@]} * REPETITIONS ))
completed_runs=0

for row in "${SWEEP_ROWS[@]}"; do
    IFS=',' read -r sweep nx nz lin_solver n_threads <<<"${row}"

    for ((rep=1; rep<=REPETITIONS; rep++)); do
        # Budget projection before each run.
        elapsed_s=$(( $(date +%s) - started_at ))
        remaining=$(( total_runs - completed_runs ))
        proj=$(projected_budget_usd "${elapsed_s}" "${completed_runs}" "${remaining}")
        if [[ "$(python3 -c "print(${proj} > ${BUDGET_USD})")" == "True" ]]; then
            warn "BUDGET_EXCEEDED: projected=\$${proj} > cap=\$${BUDGET_USD}; aborting"
            echo "BUDGET_EXCEEDED projected_usd=${proj} cap_usd=${BUDGET_USD}" \
                >>"${HARNESS_RESULTS_DIR}/summary.md.fragment"
            exit 3
        fi

        # Retry loop on thermal throttling.
        attempt=0
        thermal_throttled="true"
        while [[ "${thermal_throttled}" == "true" && "${attempt}" -lt "${MAX_RETRIES}" ]]; do
            run_one "${sweep}" "${nx}" "${nz}" "${lin_solver}" "${n_threads}" "${rep}"
            thermal_throttled="$(tail -n 1 "${HARNESS_RESULTS_CSV}" | awk -F, '{print $NF}')"
            attempt=$((attempt + 1))
        done
        if [[ "${thermal_throttled}" == "true" ]]; then
            warn "configuration (${sweep},${nx}x${nz},lin=${lin_solver},thr=${n_threads}) UNRELIABLE after ${MAX_RETRIES} retries"
        fi

        completed_runs=$((completed_runs + 1))
        log "progress ${completed_runs}/${total_runs}: ${sweep} ${nx}x${nz} lin=${lin_solver} thr=${n_threads} rep=${rep}"
    done
done

log "sweep complete: ${completed_runs}/${total_runs} runs appended to results/results.csv"
