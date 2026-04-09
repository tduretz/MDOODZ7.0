## ADDED Requirements

### Requirement: Periodic boundary condition test
The CI test suite SHALL include a GTest case that runs a simulation with `periodic_x=1` and verifies field continuity across the lateral boundaries.

#### Scenario: Velocity is periodic across x boundaries
- **WHEN** a simulation runs with `periodic_x=1`, `shear_style=1` (simple shear), and a circular inclusion
- **THEN** the horizontal velocity (Vx) at the left boundary SHALL match the horizontal velocity at the right boundary within machine precision

#### Scenario: Pressure is periodic across x boundaries
- **WHEN** the same periodic simulation completes
- **THEN** the pressure field at the leftmost cell column SHALL match the pressure at the rightmost cell column within solver tolerance

### Requirement: Free-slip boundary verification
The CI test suite SHALL include a GTest case that verifies free-slip boundary conditions produce zero shear stress at the domain boundaries.

#### Scenario: Zero shear stress at free-slip walls
- **WHEN** a simulation runs with free-slip BCs on all four boundaries and a pure shear background strain rate
- **THEN** the shear stress (`/Vertices/sxz`) at all boundary vertices SHALL be zero (< 1e-10 × stress scale)

### Requirement: No-slip boundary verification
The CI test suite SHALL include a GTest case that verifies no-slip (Dirichlet zero velocity) boundary conditions.

#### Scenario: Zero tangential velocity at no-slip walls
- **WHEN** a simulation runs with no-slip BCs on the top and bottom boundaries
- **THEN** the tangential velocity component at those boundaries SHALL be exactly zero
