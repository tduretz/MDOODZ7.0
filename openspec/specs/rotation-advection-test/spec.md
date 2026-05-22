## ADDED Requirements

### Requirement: Rigid body rotation advection test
The system SHALL provide a GTest (`TESTS/RotationAdvectionTests.cpp`) that prescribes a solid-body rotation velocity field, advects markers carrying a compositional anomaly for one full revolution, and measures the L2 error of the composition field against the initial condition. The test SHALL compare L2 errors across `reseed_mode` 0, 1, and 2.

#### Scenario: Velocity field is prescribed, Stokes is skipped
- **GIVEN** a domain centered at the origin with a prescribed angular velocity ω
- **WHEN** the simulation advances one time step
- **THEN** the velocity field SHALL be `Vx = -ω·z`, `Vz = ω·x` (rigid body rotation about the origin), and the Stokes solver SHALL NOT be called (`mechanical=0`)

#### Scenario: Compositional anomaly is a circular disk
- **GIVEN** a simulation domain [-0.5, 0.5] × [-0.5, 0.5] (non-dimensional, unit scaling)
- **WHEN** particles are initialized
- **THEN** particles with `(x - x0)² + (z - z0)² < R²` SHALL be assigned phase 1 (anomaly), and all others phase 0 (background), where `x0`, `z0`, `R` are configurable via `user0`, `user1`, `user2`

#### Scenario: One full revolution returns composition to initial state
- **WHEN** the simulation runs for exactly `T = 2π/ω` time (one revolution)
- **THEN** the analytical composition field SHALL be identical to the initial condition (every marker returns to its starting position)

#### Scenario: L2 error is computed on cell centres
- **GIVEN** the simulation has completed one revolution
- **WHEN** the L2 error is computed
- **THEN** it SHALL be the relative L2 norm of the cell-centre phase proportion field: `L2 = sqrt(sum((C_num - C_init)²) / sum(C_init²))`, where `C_init` and `C_num` are the phase-1 volume fraction on cell centres at initial and final time

### Requirement: L2 comparison across reseed_modes
The test SHALL run three simulations (reseed_mode 0, 1, 2) with identical parameters and compare L2 errors to demonstrate the relative quality of each reseeding strategy.

#### Scenario: Three modes produce comparable or improving L2
- **WHEN** all three reseed_mode simulations complete
- **THEN** the test SHALL report L2(mode 0), L2(mode 1), L2(mode 2) and assert that all three are below 0.5 (sanity bound)

#### Scenario: Mode 2 does not degrade accuracy
- **WHEN** the L2 errors are compared
- **THEN** L2(mode 2) SHALL be less than or equal to 1.5 × L2(mode 1) (mode 2 must not be significantly worse than mode 1)

### Requirement: Test parameter file configuration
The test SHALL use a dedicated parameter file `TESTS/RotationAdvection/RotationAdvection.txt` with unit scaling (`eta=1, L=1, V=1, T=1`), `mechanical=0`, `thermal=0`, `advection=1`, and `Courant=0.25`.

#### Scenario: Grid resolution is suitable for CI
- **GIVEN** the test must complete in CI within a few minutes
- **THEN** the grid SHALL be 51×51 with 4×4 markers per cell

#### Scenario: Time stepping covers one revolution
- **GIVEN** angular velocity ω and total time T = 2π/ω
- **WHEN** the simulation runs
- **THEN** the number of time steps SHALL be sufficient to complete one full revolution with `Courant ≤ 0.25`

#### Scenario: Output written at initial and final step only
- **THEN** `writer_step` SHALL be set so that HDF5 output is written at the first and last step (to compare initial vs final composition)
