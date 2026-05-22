## ADDED Requirements

### Requirement: Timing instrumentation for particle interpolation phase
The system SHALL measure wall-clock time for the particle-to-grid interpolation phase (all `P2Mastah` calls and related operations) using `omp_get_wtime()` and emit the elapsed time via `LOG_TIME`.

#### Scenario: Particle interpolation timed
- **WHEN** the particle interpolation phase executes in the main timestep loop
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Particle interpolation: <elapsed> sec"

### Requirement: Timing instrumentation for rheology evaluation
The system SHALL measure wall-clock time for `UpdateNonLinearity` and `RheologicalOperators` calls within the nonlinear iteration loop.

#### Scenario: Rheology evaluation timed per iteration
- **WHEN** a nonlinear iteration executes rheology evaluation
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Rheology evaluation: <elapsed> sec"

### Requirement: Timing instrumentation for Stokes assembly
The system SHALL measure wall-clock time for `BuildStokesOperatorDecoupled` and `BuildJacobianOperatorDecoupled` calls.

#### Scenario: Stokes assembly timed
- **WHEN** Stokes matrix assembly executes within a nonlinear iteration
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Stokes assembly: <elapsed> sec"

### Requirement: Timing instrumentation for direct solver
The system SHALL measure wall-clock time for the CHOLMOD/UMFPACK direct solver phase (Cholesky analysis, factorization, and solve).

#### Scenario: Direct solver timed
- **WHEN** the direct solver executes
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Direct solver: <elapsed> sec"

### Requirement: Timing instrumentation for thermal solver
The system SHALL measure wall-clock time for the thermal solver phase.

#### Scenario: Thermal solver timed
- **WHEN** the thermal solver executes in the main timestep loop
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Thermal solver: <elapsed> sec"

### Requirement: Timing instrumentation for advection
The system SHALL measure wall-clock time for the advection phase (marker advection, `RogerGunther` and related calls).

#### Scenario: Advection timed
- **WHEN** the advection phase executes in the main timestep loop
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Advection: <elapsed> sec"

### Requirement: Timing instrumentation for I/O
The system SHALL measure wall-clock time for HDF5 output operations (`WriteOutputHDF5`, `WriteOutputHDF5Particles`).

#### Scenario: I/O timed
- **WHEN** HDF5 output is written
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "HDF5 output: <elapsed> sec"

### Requirement: Per-iteration timing summary
The system SHALL emit a summary of all timed phases at the end of each nonlinear iteration, showing the elapsed time for each phase within that iteration.

#### Scenario: Iteration summary
- **WHEN** a nonlinear iteration completes
- **THEN** a `LOG_TIME` message SHALL be emitted listing the elapsed time for rheology, assembly, and solver phases for that iteration

### Requirement: Per-timestep timing summary
The system SHALL emit a summary at the end of each timestep showing the total wall-clock time for the timestep and a breakdown by major phase (interpolation, nonlinear solve total, thermal, advection, I/O).

#### Scenario: Timestep summary
- **WHEN** a timestep completes
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Step <N> total: <elapsed> sec" followed by a breakdown of particle interpolation, nonlinear solve, thermal, advection, and I/O times

### Requirement: Timing uses omp_get_wtime
All timing measurements SHALL use `omp_get_wtime()` for wall-clock accuracy. The existing `#ifdef _OMP_` fallback to `clock()/CLOCKS_PER_SEC` SHALL be preserved for non-OpenMP builds.

#### Scenario: Timing with OpenMP
- **WHEN** the code is compiled with `_OMP_` defined
- **THEN** timing measurements SHALL use the real `omp_get_wtime()` function

#### Scenario: Timing without OpenMP
- **WHEN** the code is compiled without `_OMP_`
- **THEN** timing measurements SHALL fall back to `clock()/CLOCKS_PER_SEC`
