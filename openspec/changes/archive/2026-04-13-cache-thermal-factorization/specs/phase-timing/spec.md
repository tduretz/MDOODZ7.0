## MODIFIED Requirements

### Requirement: Per-timestep timing summary
The system SHALL emit a summary at the end of each timestep showing the total wall-clock time for the timestep and a breakdown by ALL major phases (interpolation, stokes setup, nonlinear solve total including overhead, thermal, advection, post-solve, free surface, reseeding, melting, anisotropy, GSE, I/O) plus the coverage percentage.

The thermal phase timing SHALL correctly measure the wall-clock time of `EnergyDirectSolve` regardless of whether the CHOLMOD context is created locally or provided externally via a persistent `DirectSolver`.

#### Scenario: Timestep summary with full coverage
- **WHEN** a timestep completes
- **THEN** a `LOG_TIME` message SHALL be emitted with the format "Step <N> total: <elapsed> sec" followed by a breakdown listing ALL timed phases and a coverage line showing `tracked/wall × 100%`

#### Scenario: Thermal timing unaffected by factorization caching
- **WHEN** the thermal solver uses a cached `cholmod_factor` (skipping `cholmod_analyze`)
- **THEN** the `dt_thermal` measurement SHALL still capture the full wall-clock time of `EnergyDirectSolve` including matrix assembly, numeric factorization, and solve
