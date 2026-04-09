## ADDED Requirements

### Requirement: Stress boundary condition enforcement
The test suite SHALL verify that Neumann (stress/traction) boundary conditions (BC type 13) produce a stress field consistent with the applied traction.

#### Scenario: Boundary stress matches applied traction
- **WHEN** a simulation applies a known normal stress on one boundary using BC type 13
- **THEN** the deviatoric stress `sxxd` near that boundary SHALL be non-zero and consistent with the applied value (within 20%)

#### Scenario: Solver converges with stress BCs
- **WHEN** stress BCs are applied instead of velocity BCs
- **THEN** the solver SHALL converge (iteration count > 0 and residual below tolerance)
