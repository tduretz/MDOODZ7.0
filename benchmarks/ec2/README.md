# EC2 Performance Harness

Reusable AWS EC2-backed perf-and-memory measurement infrastructure for MDOODZ.
First consumer is the `add-gmg-stokes-defence` change; any later change that
needs perf numbers can re-invoke this harness with a new `grids.yaml`.

## What it does

1. `provision.sh` — launches an `c7i.8xlarge` in `eu-west-1`, tags it with the
   current change name and git SHA, returns the instance ID and public DNS.
2. `run_perf_sweep.sh` — reads `grids.yaml`, SSHes to the instance, clones
   MDOODZ at the tagged SHA, builds it, runs each `(grid, lin_solver, threads)`
   configuration N=5 times, captures wall time + peak RSS + `turbostat`
   throttle state, appends one row per run to `results/results.csv`.
3. `teardown.sh` — terminates the tagged instance. Idempotent, safe to call
   multiple times, registered as a trap on every harness entry point so a
   signal interrupt still tears down.
4. `analyse.py` — aggregates `results.csv` into median + IQR tables, fits
   `log RSS = log a + α log DOF` for each solver, emits `summary.md` plus
   `results/figs/perf_*.png`.
5. `mocks/aws` — stub `aws` CLI for dry-run tests of the whole pipeline
   without hitting AWS.

## Workstation prerequisites

To run a real sweep (not a dry-run) you need:

- **AWS CLI v2** configured with credentials scoped to:
    - `ec2:RunInstances`, `ec2:TerminateInstances`, `ec2:DescribeInstances`,
      `ec2:CreateTags`, `ec2:DescribeSpotPriceHistory`
    - `budgets:ViewBudget` (optional, only if you want the budget cap to
      consult live AWS billing rather than a static price table)
- **An SSH key pair** registered in `eu-west-1` whose private half lives at
  `~/.ssh/mdoodz-ec2.pem` (override via `MDOODZ_EC2_KEY_PATH`).
- **A pinned MDOODZ commit SHA** on `origin` — the sweep checks out the exact
  SHA named in `grids.yaml` (falls back to `HEAD` if unset) so the measurements
  are reproducible and traceable per row in `results.csv`.

Nothing else. The MDOODZ build itself has no AWS dependency; the harness is
an external tool that happens to drive an AWS instance.

## Usage

### Real run

```bash
./provision.sh                    # returns INSTANCE_ID + PUBLIC_DNS
./run_perf_sweep.sh               # reads grids.yaml, appends to results.csv
./teardown.sh                     # clean shutdown
python analyse.py results/results.csv --out results/summary.md
```

`run_perf_sweep.sh` also registers `teardown.sh` as an `EXIT`/`INT` trap so a
Ctrl-C still tears down cleanly.

### Dry-run (no AWS)

Set `MDOODZ_EC2_DRY_RUN=1` to route all `aws` calls through `mocks/aws`:

```bash
MDOODZ_EC2_DRY_RUN=1 ./provision.sh
MDOODZ_EC2_DRY_RUN=1 ./run_perf_sweep.sh
MDOODZ_EC2_DRY_RUN=1 ./teardown.sh
```

The stub echoes every call so the end-to-end orchestration can be verified
without cost.

### Flags

- `./provision.sh --instance-type c7g.8xlarge` — launches a non-pinned type
  (exploratory Graviton runs). A warning is logged and the chosen type is
  recorded in every results row so defence aggregation can exclude non-pinned
  data.
- `./run_perf_sweep.sh --spot` — uses spot pricing for exploratory runs.
  Rejected for headline measurements per design D1.
- `./run_perf_sweep.sh --budget-usd 50` — override the USD 50 default cap.

## Output

Every run appends a row to `results/results.csv` with columns:

```
timestamp, instance_id, instance_type, region, git_sha,
grid_nx, grid_nz, lin_solver, n_threads, repetition,
wall_time_s, peak_rss_kb, n_its_fgmres, n_picard, thermal_throttled
```

The file is append-only across sweeps — you can re-run a month later on a
different SHA and both datasets coexist. `analyse.py` filters by `git_sha`
when you want a single-SHA view.

Generated artefacts:

- `results/results.csv` — raw one-row-per-run data
- `results/summary.md` — aggregated medians + IQRs + α fit + warnings
- `results/figs/perf_memory_scaling.png` — `log RSS` vs `log DOF`, two curves
- `results/figs/perf_walltime_comparison.png` — bar chart at 201×201
- `results/figs/perf_strong_scaling.png` — wall time vs `n_threads` at 201²

`results/` is gitignored — local, user-managed. A future change can add S3
sync if we want durable perf history.

## Safety

- **Orphan detection**: `provision.sh` refuses to launch if a
  `change = add-gmg-stokes-defence`-tagged instance is already alive. It
  prints the surviving instance ID and the terminate command.
- **Budget cap**: `run_perf_sweep.sh` tracks accumulated
  `instance-hours × instance price`. If the projection exceeds USD 50 the
  sweep writes a `BUDGET_EXCEEDED` marker to `summary.md` and calls
  `teardown.sh` before exiting.
- **SIGINT trap**: on Ctrl-C the current row finishes writing, `teardown.sh`
  fires, and the process exits cleanly.
- **Teardown-on-exit**: `trap teardown.sh EXIT INT TERM` is registered at the
  top of `run_perf_sweep.sh` so any error path still tears the instance down.

## Design references

- Decisions **D1–D7** in `openspec/changes/add-gmg-stokes-defence/design.md`
  pin the instance type, measurement protocol, output format, repository
  layout, and budget/teardown discipline.
- Contract in
  `openspec/changes/add-gmg-stokes-defence/specs/gmg-performance-harness/spec.md`.
