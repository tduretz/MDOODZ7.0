## ADDED Requirements

### Requirement: Picard iteration convergence test
The CI test suite SHALL include a GTest case that runs a nonlinear problem with `Newton=0` (Picard iteration only) and verifies convergence.

#### Scenario: Picard converges for power-law rheology
- **WHEN** a simulation runs with `Newton=0`, `pwlv=40` (dry olivine), and a small grid for 1 time step with sufficient nonlinear iterations (`nit=50`)
- **THEN** the final momentum residual (`/Iterations/rx_abs`) SHALL be below the specified tolerance (`lin_abs_mom`)

#### Scenario: Picard requires more iterations than Newton
- **WHEN** the same scenario is run with Picard (`Newton=0`, `nit=50`) and the iteration count is compared to the Newton variant
- **THEN** Picard SHALL require more iterations than Newton for the same convergence tolerance

### Requirement: Adaptive time stepping test
The CI test suite SHALL include a GTest case that verifies adaptive time stepping adjusts dt based on Courant constraints.

#### Scenario: Time step adapts to velocity field
- **WHEN** a simulation runs with `constant_dt=0`, `Courant=0.5`, for 3 time steps with a velocity field that varies between steps
- **THEN** the time step sizes (from `/Model/Params[7]` across output files) SHALL NOT all be identical — the solver SHALL have adapted dt
