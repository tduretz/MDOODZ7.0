## ADDED Requirements

### Requirement: Harness lives under `benchmarks/ec2/` as reusable infrastructure

The performance harness SHALL live at repository-top-level `benchmarks/ec2/` and SHALL be structured as reusable infrastructure — not a one-off script tied to the GMG change. Future OpenSpec changes that need perf measurements SHALL be able to re-invoke it by supplying a grid-specification file without modifying the harness internals.

The directory SHALL contain at minimum: `README.md`, `provision.sh`, `run_perf_sweep.sh`, `analyse.py`, `teardown.sh`, and `grids.yaml`. A `results/` subdirectory is gitignored and user-managed.

#### Scenario: A later change re-uses the harness without modification

- **WHEN** a subsequent OpenSpec change (e.g. a two-phase-flow perf baseline) wants to run a different grid sweep
- **THEN** the change's tasks include a new `grids.yaml`-style spec file and re-invoke `run_perf_sweep.sh`, and no file under `benchmarks/ec2/` beyond `grids.yaml` and `results/` is modified

#### Scenario: Harness README documents the required workstation setup

- **WHEN** a new contributor reads `benchmarks/ec2/README.md` on a fresh machine
- **THEN** the README specifies: AWS CLI configured with credentials scoped to EC2 launch + termination + budget alerts, an SSH key pair registered in the target region, and a pinned MDOODZ commit SHA to ensure reproducible builds. No other workstation state is required.

### Requirement: EC2 instances are pinned to `c7i.8xlarge` in `eu-west-1`

The harness SHALL launch EC2 instances of type `c7i.8xlarge` in region `eu-west-1` for all perf measurements whose results will appear in the defence document. Instance type is fixed per design D1 for thermal and configuration consistency; deviations require a documented rationale in the results artefact.

#### Scenario: Default provisioning uses the pinned type

- **WHEN** `provision.sh` is invoked without arguments
- **THEN** it launches a `c7i.8xlarge` in `eu-west-1`, tags it with `project = mdoodz-gmg`, `change = add-gmg-stokes-defence`, and `git_sha = <current HEAD>`, and returns the instance ID and public DNS

#### Scenario: Instance-type override is allowed with explicit flag

- **WHEN** `provision.sh --instance-type c7g.8xlarge` is invoked for an exploratory Graviton run
- **THEN** it launches the specified type, and the resulting results CSV row carries the overridden type string so defence-aggregation scripts can exclude non-pinned runs

### Requirement: Every configuration is measured 5 times with median + IQR reporting and throttle rejection

For every `(grid_nx, grid_nz, lin_solver, n_threads)` combination in the sweep, the harness SHALL execute exactly 5 repetitions. Results SHALL be aggregated by reporting the median wall time and peak RSS plus the inter-quartile range as an uncertainty band. Runs in which thermal throttling is detected via `turbostat` or `/sys/class/powercap` during the measurement window SHALL be rejected and re-run.

#### Scenario: Non-throttled sweep produces statistically sufficient output

- **WHEN** the full grid sweep runs on a well-behaved instance without thermal events
- **THEN** `analyse.py` produces a `summary.md` with per-configuration median wall time, median peak RSS, IQR bands, and a scaling-law fit; every configuration has exactly 5 contributing repetitions

#### Scenario: Throttled run is rejected and re-run

- **WHEN** the harness detects thermal throttling during a measurement
- **THEN** that repetition is flagged with `thermal_throttled = true` in `results.csv`, excluded from aggregation, and re-run up to 3 times before the configuration is marked UNRELIABLE in `summary.md`

#### Scenario: Configuration with > 20% throttle rejections warns in summary

- **WHEN** more than 20% of repetitions for a given configuration were thermally rejected
- **THEN** `summary.md` includes a `WARNING` block naming the configuration and recommending either a different instance type or an off-peak reschedule in `eu-west-1`

### Requirement: Results format is one-row-per-run CSV with full provenance

The raw `results.csv` SHALL contain one row per individual run (not pre-aggregated) with at minimum the columns `timestamp, instance_id, instance_type, region, git_sha, grid_nx, grid_nz, lin_solver, n_threads, repetition, wall_time_s, peak_rss_kb, n_its_fgmres, n_picard, thermal_throttled`. The file SHALL be append-only across harness invocations so multiple sweeps accumulate in one canonical dataset.

#### Scenario: Row-level provenance enables future re-aggregation

- **WHEN** a reviewer asks a question after the defence has shipped (e.g. "what was the IQR at 81×81 with 4 threads?")
- **THEN** the answer is reproducible by re-running `analyse.py` on `results.csv` without relaunching any EC2 instances

#### Scenario: Multiple sweeps accumulate safely

- **WHEN** the harness is invoked a second time (e.g. a month later, on a different git SHA)
- **THEN** new rows are appended to `results.csv` without touching existing rows, and `analyse.py` can filter by `git_sha` for historical comparisons

### Requirement: Harness enforces a cost safeguard and guarantees teardown

The harness SHALL refuse to launch a new instance if a previous run's instance tagged `change = add-gmg-stokes-defence` is still alive. The harness SHALL enforce a hard budget cap of USD 50 per invocation: if the accumulated runtime × instance price exceeds this cap, the sweep aborts and `teardown.sh` is called. `teardown.sh` SHALL terminate instances on every harness exit path — success, failure, signal interruption, or timeout.

#### Scenario: Orphan instance blocks new run

- **WHEN** a previous harness run exited without calling `teardown.sh` (e.g. a network disconnect) and `provision.sh` is re-invoked
- **THEN** `provision.sh` detects the orphan instance, refuses to launch a new one, and prints the surviving instance ID plus the command to manually terminate it

#### Scenario: Budget cap triggers clean shutdown

- **WHEN** a long sweep is projected to exceed USD 50 mid-run
- **THEN** the harness stops issuing new runs, writes a `BUDGET_EXCEEDED` marker to `summary.md`, and invokes `teardown.sh` before exiting

#### Scenario: Signal interruption still tears down

- **WHEN** the operator sends SIGINT (Ctrl-C) during a running sweep
- **THEN** the harness traps the signal, finishes writing the current row, invokes `teardown.sh`, and exits cleanly with no instance left alive

### Requirement: Grid sweep covers the regression envelope declared in the previous change

The default `grids.yaml` SHALL specify a sweep that measures the two regression bounds declared in `add-gmg-stokes-solver`'s spec: memory-scaling exponent `α ≤ 1.2` and 201×201 wall-time ratio ≤ 1/3 vs CHOLMOD. Specifically it SHALL include:

- memory-scaling grids: 41×41, 81×81, 161×161, 321×321, 801×801 under `lin_solver ∈ {0, 3}` with `n_threads = 1`
- wall-time comparison: 201×201 under `lin_solver ∈ {0, 3}` with `n_threads = 1`
- strong-scaling sweep: 201×201 under `lin_solver = 3` with `n_threads ∈ {1, 2, 4, 8, 16, 32}`

#### Scenario: α measurement is derivable from the default sweep

- **WHEN** `analyse.py` runs on the default sweep's results
- **THEN** it produces a least-squares fit `log RSS = log a + α log DOF` from the five memory-scaling grids for each solver, and reports the fitted α with its 95% CI; for `lin_solver = 3` the central value SHALL be ≤ 1.2 (or the defence reports the shortfall honestly)

#### Scenario: Strong-scaling plot is derivable

- **WHEN** `analyse.py` runs on the default sweep's results
- **THEN** it produces a strong-scaling efficiency plot (wall time vs `n_threads` normalised to single-thread) for `lin_solver = 3` at 201×201, annotated with the observed plateau thread count
