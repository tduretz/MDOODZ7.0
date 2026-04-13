## Context

The perf.csv currently tracks 11 per-step timers: `rheology_s`, `assembly_s`, `solve_s`, `thermal_s`, `advection_s`, `free_surface_s`, `reseeding_s`, `melting_s`, `anisotropy_s`, `gse_s`, `output_s`. Together they account for only 46–61% of wall time (measured across 1000×800 at 1–16 threads). The remaining 39–54% ("gap") is distributed across several un-instrumented code sections in `Main_DOODZ.c`.

Analysis of the timestep loop (lines ~455–1350) identifies these un-timed sections:

1. **Particle-to-grid interpolation** (~lines 519–660): ~20× `P2Mastah` calls, `CohesionFrictionDilationGrid`, `UpdateDensity`, `ComputeLithostaticPressure`, `SurfaceDensityCorrection`. Currently timed only as a single `LOG_TIME` message but **not** written to perf.csv.

2. **Stokes setup** (~lines 663–725): `RheologicalOperators` (initial), `ApplyBC`, `SetBCs`, `InitialiseSolutionFields`, `EvalNumberOfEquations`, `SAlloc` ×8–12, `cholmod_start`, `InitialiseSolutionVector`. Not timed at all.

3. **Non-linear loop overhead** (within the while loop, ~lines 727–910): Between the three timed phases (rheology, assembly, solve), there is: `MinMaxArrayTag` diagnostic output (when `noisy==1`), diagonal scaling (`ExtractDiagonalScale`, `ScaleMatrix`), `RheologicalOperators` calls for residual evaluation, `EvaluateStokesResidualDecoupled`, residual storage, Picard-to-Newton switching logic, `FreeSparseSystems`. None of these are individually timed.

4. **Post-solve particle updates** (~lines 920–1095): `UpdateParticleStress`, `UpdateParticlePressure`, `UpdateParticleGrainSize`, `MeltFractionGrid` + `UpdateAlphaCp` (second melting pass), `UpdateParticlePhi`, chemical solver, `UpdateParticleX`, second `UpdateDensity`, `UpdateParticleDensity`, `MassSourceTerm`, `AccumulatedStrainII`, `ComputeMeanQuantitesForTimeSeries`. Only melting and GSE have their own timers; the rest are unaccounted.

5. **Memory cleanup** (~lines 1260–1290): `SFree` ×8–12, `DoodzFree` ×3–6. Not timed.

## Goals / Non-Goals

**Goals:**
- Close the perf.csv instrumentation gap so that `sum(all_tracked_columns) / wall_s >= 0.95` for every timestep
- Add minimal new columns to perf.csv — group related operations into a few aggregate timers rather than one column per function call
- Maintain backward compatibility where possible (append new columns after existing ones)
- Re-run the full benchmark sweep on c5ad.4xlarge with the new instrumentation
- Produce an updated report where "Other/overhead" is < 5%

**Non-Goals:**
- Refactoring or optimising any of the timed code paths (this change is measurement only)
- Changing the Newton/Picard solver algorithm
- Switching to a different benchmark instance type
- Adding per-iteration perf.csv rows (keep one row per timestep)

## Decisions

### D1: Timer granularity — 4 new aggregate timers

Rather than adding a column per function call (which would make perf.csv unwieldy with 30+ columns), group the gap into 4 new timers:

| Timer variable | perf.csv column | What it measures |
|----------------|-----------------|------------------|
| `dt_interp` | `interp_s` | All P2Mastah calls + grid property setup (lines 519–660) |
| `dt_stokes_setup` | `stokes_setup_s` | Pre-loop Stokes initialization: BCs, alloc, RheologicalOperators (lines 663–725) |
| `dt_nl_overhead` | `nl_overhead_s` | Within-loop overhead not covered by rheology/assembly/solve: diagonal scaling, residual evaluation, convergence checks, Picard↔Newton logic, FreeSparseSystems (lines 727–910, excluding the 3 existing timers) |
| `dt_post_solve` | `post_solve_s` | Post-solve particle updates, density, strain, chemical solver, memory cleanup (lines 920–1290) |

**Rationale**: 4 columns is the minimum to cover all 5 identified gap regions (regions 4 and 5 are merged into `post_solve` since memory cleanup is trivial). Each maps to a coherent phase of the timestep.

**Alternative considered**: One catch-all `dt_other` column computed as `wall - sum(existing)`. Rejected because it provides no diagnostic value — you can't tell if interpolation or residual evaluation is the bottleneck.

### D2: Column ordering — append after existing columns

New perf.csv header:
```
step,wall_s,time_ma,rheology_s,assembly_s,solve_s,thermal_s,advection_s,free_surface_s,reseeding_s,melting_s,anisotropy_s,gse_s,output_s,interp_s,stokes_setup_s,nl_overhead_s,post_solve_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s
```

This preserves the first 14 data columns (positions 1–14) unchanged. The 4 new timing columns are inserted at positions 15–18, pushing `nit` and the metadata columns right. Total: 25 columns (was 21).

**Rationale**: Appending after `output_s` keeps existing analysis scripts working for columns 1–14 (most only use wall, rheology, assembly, solve, thermal). The `nit` and metadata columns shift but are rarely accessed by positional index.

### D3: `nl_overhead` implementation — accumulate per-iteration remainder

Inside the while loop, each iteration already times `dt_rheology`, `dt_assembly`, `dt_solve_iter`. The approach:

```c
double t_iter_start = omp_get_wtime();
// ... existing rheology timing ...
// ... existing assembly timing ...
// ... residual evaluation, scaling, convergence check ...
// ... existing solve timing ...
// ... FreeSparseSystems ...
double dt_iter_total = omp_get_wtime() - t_iter_start;
dt_nl_overhead_total += dt_iter_total - dt_rheology - dt_assembly - dt_solve_iter;
```

This captures everything in the iteration that isn't already timed without having to wrap each individual call. Simple and accurate.

### D4: Benchmark re-run — same parameters, same instance

Re-run on the existing c5ad.4xlarge (i-0df5754e400cb3ec1, 16 vCPU, 32 GB):
- Same 4 resolutions: lowres (100×80), default (200×160), medres (501×401), highres (1000×800)
- Same 7 thread counts: 1, 2, 4, 6, 8, 12, 16
- Same 10 steps per run
- Results stored in a new timestamped directory under `benchmark-results/`

### D5: benchmark.sh and benchmark-report.sh updates

- `benchmark.sh`: Update `summary.csv` header and awk extraction to include the 4 new columns
- `benchmark-report.sh`: Add new columns to the Results Summary table; replace "Other/overhead" with exact tracked percentages; add a "Coverage" column showing `tracked/wall × 100%`

## Risks / Trade-offs

**[Timing overhead]** → Each new `omp_get_wtime()` call adds ~50ns. With 4 new timers × ~10 calls each, that's < 1µs per step — negligible vs. the 17–37s step wall time.

**[perf.csv column shift]** → Scripts that access `nit`, `n_particles`, `peak_rss_mb` by column number (15–21) will break. → Mitigation: The benchmark scripts (`benchmark.sh`, `benchmark-report.sh`) are the only consumers and will be updated simultaneously. Document the column change in skill-logging.

**[nl_overhead accuracy]** → The subtraction approach (`iter_total - rheo - asm - solve`) may include tiny timing imprecision. → Acceptable: at ~4s overhead per step, sub-millisecond error is irrelevant.

**[Stokes-disabled runs]** → When `mechanical == 0`, `stokes_setup_s`, `nl_overhead_s` will be 0.0. → Fine: those phases simply don't execute.
