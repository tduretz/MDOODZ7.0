## ADDED Requirements

### Requirement: Diffusion creep convergence test
The CI test suite SHALL include a GTest case that runs a ShearTemplate-like scenario with `linv` as the only active flow law (all other mechanisms disabled) and verifies Newton solver convergence.

#### Scenario: Linear diffusion creep converges in one iteration
- **WHEN** a simulation runs with `cstv=0`, `pwlv=0`, `linv=1`, `expv=0`, `gbsv=0` on a small grid (≤51×51) for 1 time step
- **THEN** the Newton solver SHALL converge in exactly 1 iteration (linear problem)

#### Scenario: Diffusion creep with grain size produces correct viscosity range
- **WHEN** a simulation runs with `linv=1` and a grain size `gs_ref` of 1e-3 m
- **THEN** the maximum viscosity in the output (`/Centers/eta_n`) SHALL fall within physically reasonable bounds (1e18–1e24 Pa·s after scaling)

### Requirement: Peierls creep convergence test
The CI test suite SHALL include a GTest case that runs a scenario with `expv` as the primary flow law and verifies nonlinear convergence.

#### Scenario: Peierls creep requires multiple Newton iterations
- **WHEN** a simulation runs with `expv=1` (Peierls mechanism) on a small grid for 1 time step with a background strain rate of 1e-15 s⁻¹
- **THEN** the Newton solver SHALL converge in more than 1 and fewer than 20 iterations

### Requirement: Grain boundary sliding convergence test
The CI test suite SHALL include a GTest case that runs a scenario with `gbsv` as the primary flow law.

#### Scenario: GBS creep convergence
- **WHEN** a simulation runs with `gbsv=1` on a small grid for 1 time step
- **THEN** the Newton solver SHALL converge in fewer than 20 iterations

### Requirement: Composite creep convergence test
The CI test suite SHALL include a GTest case that activates both `pwlv` and `linv` simultaneously to verify composite rheology.

#### Scenario: Composite pwl+lin convergence
- **WHEN** a simulation runs with `pwlv=40` (dry olivine) and `linv=1` on a small grid for 1 time step
- **THEN** the Newton solver SHALL converge in more than 1 and fewer than 15 iterations

#### Scenario: Composite viscosity is bounded by end-member viscosities
- **WHEN** a composite creep simulation completes
- **THEN** the effective viscosity at each cell centre SHALL be less than or equal to the viscosity that would result from the weaker mechanism alone
