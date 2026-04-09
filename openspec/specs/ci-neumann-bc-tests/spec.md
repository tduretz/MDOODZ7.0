## MODIFIED Requirements

### Requirement: Stress boundary condition enforcement
The test suite SHALL verify that Neumann (stress/traction) boundary conditions (BC type 13) produce a stress field consistent with the applied traction.

#### Scenario: Boundary stress matches applied traction
- **WHEN** a simulation applies a known normal stress on one boundary using BC type 13
- **THEN** the deviatoric stress `sxxd` near that boundary SHALL be non-zero and consistent with the applied value (within 20%)

#### Scenario: Solver converges with stress BCs
- **WHEN** stress BCs are applied instead of velocity BCs
- **THEN** the solver SHALL converge (iteration count > 0 and residual below tolerance)

#### Scenario: Stress BC pure shear L2 error
- **WHEN** a pure-shear simulation runs with stress BCs applying σ'_xx = 2η·ε̇ on E/W boundaries with known constant viscosity η and strain rate ε̇
- **THEN** the relative L2 error of the full sxxd field against the analytical σ'_xx = 2η·ε̇ SHALL be less than 1e-2
- **AND** the absolute L2 error of the sxz field SHALL be less than 1e-6 (shear stress is identically zero in pure shear)
