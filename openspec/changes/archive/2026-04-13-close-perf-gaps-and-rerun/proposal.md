## Why

The current perf.csv instrumentation accounts for only ~50–60% of actual wall time per timestep. The "gap" — ~10s out of 20s at 1000×800/16t — is un-instrumented code inside the Newton/Picard loop (marker-to-grid interpolation, boundary conditions, residual evaluation, matrix scaling, sparse system alloc/free, grid-to-marker back-interpolation, stress rotation, strain updates) and between major phases (time-step setup, Courant dt selection, marker reseeding overhead). This makes the benchmark report misleading: the "Other/overhead" row dominates at 40–55%, hiding where time actually goes.

The benchmark was run on c5ad.4xlarge (AMD EPYC 7R32, 3.3 GHz, 16 vCPUs) but the results are unreliable because half the wall time is unaccounted for. Once instrumentation is complete, we need to re-run the same sweep on the same instance type to get an accurate breakdown.

## What Changes

- Add new timing variables to capture **all** currently un-instrumented phases in `Main_DOODZ.c`:
  - `dt_interp` — marker-to-grid and grid-to-marker interpolation
  - `dt_stokes_residual` — residual evaluation and convergence checks inside the nonlinear loop
  - `dt_bc` — boundary condition application
  - `dt_setup` — per-timestep setup (Courant dt, mesh updates, stress rotation, strain bookkeeping)
- Expand perf.csv with new columns for the above timers
- Add a `dt_gap` or validation check: `wall_s - sum(all_tracked)` should be < 5% of wall time
- Update `misc/benchmark.sh` summary.csv header to include new columns
- Update `misc/benchmark-report.sh` to display new columns
- Re-run the same 28-run sweep (4 resolutions × 7 thread counts × 10 steps) on c5ad.4xlarge with the improved instrumentation
- Produce an updated benchmark report with full phase coverage

## Capabilities

### New Capabilities
- `full-phase-coverage`: Ensure perf.csv tracked time accounts for ≥95% of wall time by instrumenting all un-timed code sections in the timestep loop
### Modified Capabilities
- `phase-timing`: Add new perf.csv columns (interp_s, stokes_residual_s, bc_s, setup_s) and gap validation

## Impact

- **Code**: `MDLIB/Main_DOODZ.c` — new timer variables, expanded fprintf for perf.csv header and data row
- **Benchmark scripts**: `misc/benchmark.sh` (summary.csv header), `misc/benchmark-report.sh` (new columns in tables)
- **AWS**: Re-run on existing c5ad.4xlarge instance (i-0df5754e400cb3ec1)
- **Skills**: `skill-logging` (perf.csv column table), `skill-code-glossary` (new timer variables)
- **perf.csv format**: **BREAKING** — additional columns mean older analysis scripts expecting 21 columns will need updating
