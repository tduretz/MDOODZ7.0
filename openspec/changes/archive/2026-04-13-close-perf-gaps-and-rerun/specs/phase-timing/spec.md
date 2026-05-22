## ADDED Requirements

### Requirement: Timing instrumentation for Stokes setup phase
The system SHALL measure wall-clock time for the Stokes pre-loop initialization block (BCs, RheologicalOperators, sparse allocation, solution field initialization) using `omp_get_wtime()` and emit the elapsed time via `LOG_TIME`.

#### Scenario: Stokes setup timed
- **WHEN** the Stokes setup block executes before the nonlinear iteration loop
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Stokes setup: <elapsed> sec"

### Requirement: Timing instrumentation for nonlinear loop overhead
The system SHALL measure the wall-clock time for each nonlinear iteration that is not covered by the existing rheology, assembly, and solve timers (diagonal scaling, residual evaluation, convergence checks, Picard↔Newton logic, FreeSparseSystems) and emit the per-iteration overhead via `LOG_TIME`.

#### Scenario: NL overhead timed per iteration
- **WHEN** a nonlinear iteration completes
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Iteration NN: rheology=X assembly=Y solve=Z overhead=W sec"

### Requirement: Timing instrumentation for post-solve phase
The system SHALL measure wall-clock time for the post-solve particle updates and memory cleanup block using `omp_get_wtime()` and emit the elapsed time via `LOG_TIME`.

#### Scenario: Post-solve timed
- **WHEN** the post-solve block completes after the nonlinear loop
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Post-solve updates: <elapsed> sec"

## MODIFIED Requirements

### Requirement: Per-timestep timing summary
The system SHALL emit a summary at the end of each timestep showing the total wall-clock time for the timestep and a breakdown by ALL major phases (interpolation, stokes setup, nonlinear solve total including overhead, thermal, advection, post-solve, free surface, reseeding, melting, anisotropy, GSE, I/O) plus the coverage percentage.

#### Scenario: Timestep summary with full coverage
- **WHEN** a timestep completes
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Step <N> total: <elapsed> sec" followed by a breakdown listing ALL timed phases and a coverage line showing `tracked/wall × 100%`

### Requirement: Per-iteration timing summary
The system SHALL emit a summary of all timed phases at the end of each nonlinear iteration, showing the elapsed time for each phase within that iteration including the overhead time.

#### Scenario: Iteration summary
- **WHEN** a nonlinear iteration completes
- **THEN** a `LOG_TIME` message SHALL be emitted listing the elapsed time for rheology, assembly, solver, and overhead phases for that iteration
