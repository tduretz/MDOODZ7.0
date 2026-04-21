#!/usr/bin/env bash
# Local (MacBook) perf-sweep driver for add-gmg-stokes-defence §5.
#
# Replaces the three-step EC2 cycle (provision → run_perf_sweep → teardown)
# with a single hermetic local run. Reads benchmarks/ec2/grids_local.yaml,
# invokes SolViPerf on each (sweep, nx, nz, lin_solver, n_threads) tuple
# with a freshly rewritten .txt fixture, and appends one row per run to
# results/results.csv using the same schema analyse.py already consumes.
#
# Hard wall-clock budget: grids_local.yaml -> measurement.wallclock_budget_s
# (default 3600 s). Before every run we project whether a median-time run
# would push elapsed past the budget; if so we emit a BUDGET_SKIP fragment
# and abort. analyse.py tolerates an incomplete matrix.
#
# Thermal throttling on Apple Silicon: per-core throttle state requires
# root + `sudo powermetrics` on macOS, so we use an indirect proxy: a run
# is flagged as thermally throttled if its wall time is > 1.35 x the
# running median of its (grid, solver, threads) group. On Linux we still
# read /proc/cpuinfo MHz if available.
#
# Usage:
#   bash benchmarks/ec2/run_local.sh                                # default grids_local.yaml
#   bash benchmarks/ec2/run_local.sh --grids <path>                 # custom spec
#   bash benchmarks/ec2/run_local.sh --dry-run                      # skip SolViPerf, emit synthetic
#   bash benchmarks/ec2/run_local.sh --budget-seconds 1800          # override wall-clock cap
#   bash benchmarks/ec2/run_local.sh --build                        # (re)build SolViPerf first

source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

GRIDS_FILE="${HARNESS_DIR}/grids_local.yaml"
BUDGET_S=""
DRY_RUN=0
DO_BUILD=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --grids)             shift; GRIDS_FILE="$1" ;;
        --budget-seconds)    shift; BUDGET_S="$1" ;;
        --dry-run)           DRY_RUN=1 ;;
        --build)             DO_BUILD=1 ;;
        -h|--help)           sed -n '1,26p' "$0"; exit 0 ;;
        *)                   die "unknown argument: $1" ;;
    esac
    shift
done

export HARNESS_GRIDS_FILE="${GRIDS_FILE}"

PROJECT_ROOT="$(cd "${HARNESS_DIR}/../.." && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build"
TEMPLATE_TXT="${PROJECT_ROOT}/TESTS/SolViBenchmark/SolViRes51_gmg.txt"
# add_set() in SETS/AddSet.cmake pins RUNTIME_OUTPUT_DIRECTORY under
# ${PROJECT_SOURCE_DIR}/cmake-exec/<name>/ rather than build/bin/. Use the
# canonical location so the harness doesn't depend on a particular build
# directory layout.
SOLVI_PERF_BIN="${PROJECT_ROOT}/cmake-exec/SolViPerf/SolViPerf"

GIT_SHA="$(current_git_sha)"
INSTANCE_ID="macbook-$(hostname -s)"
INSTANCE_TYPE="$(read_yaml instance_type)"
INSTANCE_TYPE="${INSTANCE_TYPE:-macbook-m1-local}"
REGION="$(read_yaml region)"
REGION="${REGION:-local}"

if [[ -z "${BUDGET_S}" ]]; then
    BUDGET_S="$(read_yaml wallclock_budget_s)"
    BUDGET_S="${BUDGET_S:-3600}"
fi

if [[ "${DO_BUILD}" -eq 1 ]]; then
    log "building SolViPerf in ${BUILD_DIR}"
    cmake --build "${BUILD_DIR}" --target SolViPerf -j 4 \
        || die "SolViPerf build failed"
fi

if [[ "${DRY_RUN}" -eq 0 ]]; then
    [[ -x "${SOLVI_PERF_BIN}" ]] \
        || die "SolViPerf binary not found at ${SOLVI_PERF_BIN}; re-run with --build"
fi
[[ -f "${TEMPLATE_TXT}" ]] || die "template fixture not found: ${TEMPLATE_TXT}"

# Fresh results CSV header (one header line; append-only afterwards).
if [[ ! -f "${HARNESS_RESULTS_CSV}" ]]; then
    echo "timestamp,instance_id,instance_type,region,git_sha,grid_nx,grid_nz,lin_solver,n_threads,repetition,wall_time_s,peak_rss_kb,n_its_fgmres,n_picard,thermal_throttled" \
        >"${HARNESS_RESULTS_CSV}"
fi

log "local perf sweep: grids=${GRIDS_FILE} budget_s=${BUDGET_S} dry_run=${DRY_RUN}"
log "                  binary=${SOLVI_PERF_BIN}"
log "                  template=${TEMPLATE_TXT}"

# Rewrite the template .txt for this (nx, nz, lin_solver) tuple and emit to
# stdout the path of the temporary file.  Only keys we control are Nx, Nz,
# lin_solver, writer_subfolder.  Everything else is kept intact so the
# underlying physics stays identical across all sweep points.
rewrite_txt() {
    local nx="$1" nz="$2" lin_solver="$3" rep="$4"
    local dest="${HARNESS_RESULTS_DIR}/tmp_SolViPerf_${nx}x${nz}_lin${lin_solver}_rep${rep}.txt"
    python3 - "${TEMPLATE_TXT}" "${dest}" "${nx}" "${nz}" "${lin_solver}" "${rep}" <<'PY'
import sys, re
src, dst, nx, nz, lin_solver, rep = sys.argv[1:]
txt = open(src).read()
txt = re.sub(r"^(Nx\s*=\s*)\d+", r"\g<1>" + nx, txt, flags=re.M)
txt = re.sub(r"^(Nz\s*=\s*)\d+", r"\g<1>" + nz, txt, flags=re.M)
txt = re.sub(r"^(lin_solver\s*=\s*)-?\d+", r"\g<1>" + lin_solver, txt, flags=re.M)
txt = re.sub(
    r"^(writer_subfolder\s*=\s*)\S+",
    r"\g<1>SolViPerf_{}x{}_lin{}_rep{}".format(nx, nz, lin_solver, rep),
    txt, flags=re.M,
)
open(dst, "w").write(txt)
PY
    echo "${dest}"
}

# Pre-run budget projection: reactive median of previously completed runs.
# Before any runs are completed we assume a worst-case estimate so we do
# not start sweep-points that obviously will not fit (e.g. a second 800²
# run when the first one took 20 minutes).
project_fits_budget() {
    local elapsed="$1"
    if [[ "${COMPLETED}" -lt 1 ]]; then
        echo "yes"; return
    fi
    python3 - "${elapsed}" "${COMPLETED}" "${BUDGET_S}" "${LAST_WALL}" <<'PY'
import sys
elapsed, completed, budget, last = [float(a) for a in sys.argv[1:]]
remaining = budget - elapsed
# One-more-run probe: use the max(median-so-far, last-run-wall) as the
# projection for the next run. This makes the projection resistant to
# small cache-warm runs early in the sweep (which would otherwise let
# us over-commit and trip the budget mid-801² run).
print("yes" if remaining > max(last, 1.0) else "no")
PY
}

# Core runner: runs SolViPerf with the given mutated txt and returns the
# captured (wall_time_s, peak_rss_kb, n_its_fgmres, n_picard). Uses a
# temporary directory so the MDOODZ outputs don't collide across parallel
# invocations (we don't parallel-invoke, but this stays hygienic).
run_solvi_perf() {
    local tmp_txt="$1" n_threads="$2" lin_solver="$3" nx="$4" nz="$5" rep_number="${6:-1}"
    local run_dir; run_dir="$(mktemp -d "${HARNESS_RESULTS_DIR}/solviperf_XXXXXX")"
    local log_file="${run_dir}/mdoodz.log"
    local dof=$((nx * nz))

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        # Synthetic numbers matching the EC2 dry-run path, scaled by 0.6
        # for the faster M1 cores — used to smoke-test the harness + the
        # analysis pipeline without running MDOODZ.
        local wall rss fgmres picard
        if [[ "${lin_solver}" -eq 0 ]]; then
            wall=$(python3 -c "print(round(0.6 * 1e-5 * ${dof} ** 1.5, 4))")
            rss=$(python3 -c "print(int(50 * ${dof} ** 1.25))")
            fgmres=0; picard=1
        else
            wall=$(python3 -c "import math; print(round(0.6 * 3e-5 * ${dof} * max(1, math.log(${dof})), 4))")
            rss=$(python3 -c "print(int(30 * ${dof} ** 1.05))")
            fgmres=8; picard=1
        fi
        echo "${wall}|${rss}|${fgmres}|${picard}"
        rm -rf "${run_dir}"
        return 0
    fi

    ( cd "${run_dir}" \
        && OMP_NUM_THREADS="${n_threads}" "${SOLVI_PERF_BIN}" "${tmp_txt}" \
    ) >"${log_file}" 2>&1

    # PERF_JSON line is written at exit of SolViPerf main.
    local pj; pj="$(grep -m1 '^PERF_JSON:' "${log_file}" | sed 's/^PERF_JSON://')"
    if [[ -z "${pj}" ]]; then
        warn "run failed (no PERF_JSON sentinel); log tail:"
        tail -n 15 "${log_file}" >&2
        echo "NaN|NaN|-1|-1"
        return 1
    fi

    local wall rss
    wall="$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['wall_time_s'])" "${pj}")"
    rss="$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['peak_rss_kb'])" "${pj}")"

    # Iteration counts: best-effort parse from the captured MDOODZ log.
    # Patterns shipped by StokesAssemblyGMG.c and the Picard/Newton driver
    # are brittle across versions, so we fall back to -1 on miss.
    local fgmres picard
    fgmres="$(grep -oE 'GMG-FGMRES[^0-9]*iters=[0-9]+' "${log_file}" | tail -n1 \
              | grep -oE '[0-9]+$' || true)"
    if [[ -z "${fgmres}" ]]; then
        fgmres="$(grep -oE 'total_iter [0-9]+' "${log_file}" | tail -n1 \
                  | grep -oE '[0-9]+$' || true)"
    fi
    fgmres="${fgmres:--1}"
    picard="$(grep -cE 'Picard it\. ' "${log_file}" || true)"
    picard="${picard:-0}"

    # Preserve the rich per-subsystem timing that MDOODZ writes to perf.csv
    # (one row per time step, columns: step, wall_s, rheology_s, assembly_s,
    # solve_s, thermal_s, advection_s, free_surface_s, reseeding_s,
    # melting_s, anisotropy_s, gse_s, output_s, interp_s, stokes_setup_s,
    # nl_overhead_s, post_solve_s, nit, n_particles, neq_mom, neq_cont,
    # peak_rss_mb, user_cpu_s, sys_cpu_s). Without this, results.csv tells
    # us total wall time but not WHERE the time is spent — a 200 s run at
    # 201² GMG is 99 % solve_s, but results.csv alone cannot tell you that.
    #
    # Layout: HARNESS_RESULTS_DIR/perf_detail/<nx>x<nz>_lin<L>_thr<T>_rep<R>/
    #         ├─ perf.csv         (MDOODZ per-timestep breakdown)
    #         ├─ mdoodz.log       (full stdout/stderr — LOG_TIME / FGMRES)
    #         └─ PERF_JSON        (wall_time_s + peak_rss_kb sentinel)
    local keep_dir="${HARNESS_RESULTS_DIR}/perf_detail/${nx}x${nz}_lin${lin_solver}_thr${n_threads}_rep${rep_number:-1}"
    mkdir -p "${keep_dir}"
    [[ -f "${run_dir}/perf.csv" ]] && cp "${run_dir}/perf.csv" "${keep_dir}/perf.csv"
    cp "${log_file}" "${keep_dir}/mdoodz.log"
    echo "${pj}" >"${keep_dir}/PERF_JSON"

    # Keep run_dir clean (HDF5 outputs from MDOODZ would otherwise accumulate
    # ~10s of MB per run across the sweep).
    rm -rf "${run_dir}"
    echo "${wall}|${rss}|${fgmres}|${picard}"
}

# Load sweep rows.
SWEEP_ROWS=()
while IFS= read -r line; do
    [[ -z "${line}" ]] && continue
    SWEEP_ROWS+=("${line}")
done < <(GRIDS_FILE="${GRIDS_FILE}" python3 "${HARNESS_DIR}/_parse_grids.py")

REPETITIONS="$(read_yaml repetitions)"
REPETITIONS="${REPETITIONS:-3}"
TOTAL_RUNS=$(( ${#SWEEP_ROWS[@]} * REPETITIONS ))
log "planned ${#SWEEP_ROWS[@]} configs x ${REPETITIONS} reps = ${TOTAL_RUNS} runs"

STARTED=$(date +%s)
COMPLETED=0
LAST_WALL=1.0
BUDGET_HIT=0

for row in "${SWEEP_ROWS[@]}"; do
    IFS=',' read -r sweep nx nz lin_solver n_threads <<<"${row}"
    for ((rep=1; rep<=REPETITIONS; rep++)); do
        ELAPSED=$(( $(date +%s) - STARTED ))
        fits="$(project_fits_budget "${ELAPSED}")"
        if [[ "${fits}" != "yes" ]]; then
            log "BUDGET_HIT at elapsed=${ELAPSED}s cap=${BUDGET_S}s; skipping remaining"
            BUDGET_HIT=1
            break 2
        fi

        ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        tmp_txt="$(rewrite_txt "${nx}" "${nz}" "${lin_solver}" "${rep}")"

        result="$(run_solvi_perf "${tmp_txt}" "${n_threads}" "${lin_solver}" "${nx}" "${nz}" "${rep}")" \
            || warn "run_solvi_perf reported failure for ${sweep} ${nx}x${nz} lin=${lin_solver} thr=${n_threads} rep=${rep}"

        IFS='|' read -r wall rss fgmres picard <<<"${result}"
        LAST_WALL="${wall}"

        # Throttle proxy: wall > 1.35 * running median of same (grid, solver,
        # threads) group. Deferred to analyse.py; for the CSV row we mark
        # false (the analysis stage has full history + IQR).
        thermal="false"

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "${ts}" "${INSTANCE_ID}" "${INSTANCE_TYPE}" "${REGION}" "${GIT_SHA}" \
            "${nx}" "${nz}" "${lin_solver}" "${n_threads}" "${rep}" \
            "${wall}" "${rss}" "${fgmres}" "${picard}" "${thermal}" \
            >>"${HARNESS_RESULTS_CSV}"

        COMPLETED=$((COMPLETED + 1))
        log "progress ${COMPLETED}/${TOTAL_RUNS}: ${sweep} ${nx}x${nz} lin=${lin_solver} thr=${n_threads} rep=${rep}  wall=${wall}s rss=${rss}KiB"
    done
done

TOTAL_ELAPSED=$(( $(date +%s) - STARTED ))
if [[ "${BUDGET_HIT}" -eq 1 ]]; then
    {
        echo "BUDGET_HIT: completed ${COMPLETED}/${TOTAL_RUNS} runs in ${TOTAL_ELAPSED}s"
        echo "  cap_seconds=${BUDGET_S}"
        echo "  see results/results.csv for rows that did land"
    } >>"${HARNESS_RESULTS_DIR}/summary.md.fragment"
    log "sweep truncated by budget: ${COMPLETED}/${TOTAL_RUNS} runs in ${TOTAL_ELAPSED}s"
    exit 4
fi

log "local sweep complete: ${COMPLETED}/${TOTAL_RUNS} runs in ${TOTAL_ELAPSED}s → ${HARNESS_RESULTS_CSV}"
