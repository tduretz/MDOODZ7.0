## ADDED Requirements

### Requirement: Interpolation timer in perf.csv
The system SHALL record the cumulative wall-clock time for all particle-to-grid interpolation operations (all `P2Mastah` calls, `CohesionFrictionDilationGrid`, `UpdateDensity`, `ComputeLithostaticPressure`, `SurfaceDensityCorrection`, and related grid property setup) into variable `dt_interp` and write it to perf.csv as column `interp_s`.

#### Scenario: Interpolation time recorded per step
- **WHEN** the particle-to-grid interpolation block completes in a timestep
- **THEN** the elapsed wall-clock time SHALL be written to the `interp_s` column of perf.csv

#### Scenario: Interpolation timer covers all P2Mastah calls
- **WHEN** a timestep executes with `mechanical == 1`
- **THEN** `interp_s` SHALL include time for all `P2Mastah` calls, `CohesionFrictionDilationGrid`, `UpdateDensity`, `ComputeLithostaticPressure`, and `SurfaceDensityCorrection`

### Requirement: Stokes setup timer in perf.csv
The system SHALL record the cumulative wall-clock time for the Stokes pre-loop initialization (boundary conditions, `RheologicalOperators`, `ApplyBC`, `SetBCs`, `InitialiseSolutionFields`, `EvalNumberOfEquations`, sparse system allocation, `cholmod_start`, `InitialiseSolutionVector`) into variable `dt_stokes_setup` and write it to perf.csv as column `stokes_setup_s`.

#### Scenario: Stokes setup time recorded per step
- **WHEN** the Stokes initialization block completes before the nonlinear iteration loop
- **THEN** the elapsed wall-clock time SHALL be written to the `stokes_setup_s` column of perf.csv

#### Scenario: Stokes setup is zero when mechanical is disabled
- **WHEN** a timestep executes with `mechanical == 0`
- **THEN** `stokes_setup_s` SHALL be `0.0000`

### Requirement: Nonlinear loop overhead timer in perf.csv
The system SHALL record the cumulative wall-clock time for all operations inside the Newton/Picard iteration loop that are NOT already captured by `rheology_s`, `assembly_s`, or `solve_s`. This includes diagonal scaling, residual evaluation, convergence checks, Picard-to-Newton switching, and `FreeSparseSystems`. The value SHALL be stored in `dt_nl_overhead` and written to perf.csv as column `nl_overhead_s`.

#### Scenario: NL overhead computed as iteration remainder
- **WHEN** a nonlinear iteration completes
- **THEN** `dt_nl_overhead` SHALL accumulate `(iteration_total_time - dt_rheology - dt_assembly - dt_solve_iter)` for that iteration

#### Scenario: NL overhead accumulated across all iterations
- **WHEN** a timestep completes N nonlinear iterations
- **THEN** `nl_overhead_s` in perf.csv SHALL equal the sum of per-iteration overhead across all N iterations

#### Scenario: NL overhead is zero when mechanical is disabled
- **WHEN** a timestep executes with `mechanical == 0`
- **THEN** `nl_overhead_s` SHALL be `0.0000`

### Requirement: Post-solve timer in perf.csv
The system SHALL record the cumulative wall-clock time for all post-solve operations: particle stress/pressure/grain-size updates, second melting pass, chemical solver, density updates, strain accumulation, free surface cleanup, sparse memory deallocation, and timestep bookkeeping. The value SHALL be stored in `dt_post_solve` and written to perf.csv as column `post_solve_s`.

#### Scenario: Post-solve time recorded per step
- **WHEN** the post-solve block completes (from end of nonlinear loop to start of thermal solver)
- **THEN** the elapsed wall-clock time SHALL be written to the `post_solve_s` column of perf.csv

#### Scenario: Post-solve includes particle updates and memory cleanup
- **WHEN** a timestep executes with `mechanical == 1`
- **THEN** `post_solve_s` SHALL include `UpdateParticleStress`, `UpdateParticlePressure`, `UpdateParticleGrainSize`, second `MeltFractionGrid`/`UpdateAlphaCp`, `UpdateParticlePhi`, chemical solver, `UpdateParticleX`, `UpdateDensity`, `UpdateParticleDensity`, `AccumulatedStrainII`, `ComputeMeanQuantitesForTimeSeries`, and all `SFree`/`DoodzFree` calls

### Requirement: perf.csv 25-column format
The perf.csv header SHALL be:
```
step,wall_s,time_ma,rheology_s,assembly_s,solve_s,thermal_s,advection_s,free_surface_s,reseeding_s,melting_s,anisotropy_s,gse_s,output_s,interp_s,stokes_setup_s,nl_overhead_s,post_solve_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s
```
Total: 25 columns. The first 14 data columns (positions 1–14) SHALL remain unchanged from the 21-column format. New columns SHALL be at positions 15–18. Metadata columns (`nit` through `sys_cpu_s`) SHALL shift to positions 19–25.

#### Scenario: Header format
- **WHEN** perf.csv is created at the start of a simulation
- **THEN** the header line SHALL contain exactly 25 comma-separated column names as specified above

#### Scenario: Data row format
- **WHEN** a timestep completes and a row is written to perf.csv
- **THEN** the row SHALL contain exactly 25 comma-separated values in the order matching the header

### Requirement: Coverage validation target
The sum of all tracked timing columns (`rheology_s + assembly_s + solve_s + thermal_s + advection_s + free_surface_s + reseeding_s + melting_s + anisotropy_s + gse_s + output_s + interp_s + stokes_setup_s + nl_overhead_s + post_solve_s`) SHALL account for at least 95% of `wall_s` for every timestep.

#### Scenario: Coverage above 95% at highres
- **WHEN** a 1000×800 RiftingComprehensive simulation runs for 10 steps at any thread count (1–16)
- **THEN** for every step, `sum(tracked) / wall_s >= 0.95`

#### Scenario: Coverage above 95% at lowres
- **WHEN** a 100×80 RiftingComprehensive simulation runs for 10 steps at any thread count (1–16)
- **THEN** for every step, `sum(tracked) / wall_s >= 0.95`

### Requirement: benchmark.sh summary.csv updated
The `misc/benchmark.sh` summary.csv header SHALL include new average columns for the 4 new timers (`avg_interp_s`, `avg_stokes_setup_s`, `avg_nl_overhead_s`, `avg_post_solve_s`) and the awk extraction logic SHALL compute their averages from the 25-column perf.csv.

#### Scenario: Summary CSV includes new columns
- **WHEN** `benchmark.sh` produces a summary.csv
- **THEN** the header SHALL include `avg_interp_s,avg_stokes_setup_s,avg_nl_overhead_s,avg_post_solve_s`

### Requirement: benchmark-report.sh displays new columns
The `misc/benchmark-report.sh` output SHALL include the 4 new timing columns in the Results Summary table and SHALL display a Coverage percentage column (`tracked/wall × 100%`) instead of an "Other/overhead" row.

#### Scenario: Report table shows new timers
- **WHEN** `benchmark-report.sh` generates a report
- **THEN** the Results Summary table SHALL include columns for Interp, Stokes Setup, NL Overhead, and Post-Solve

#### Scenario: Report shows coverage percentage
- **WHEN** `benchmark-report.sh` generates a report
- **THEN** each row SHALL include a Coverage column showing tracked time as a percentage of wall time
