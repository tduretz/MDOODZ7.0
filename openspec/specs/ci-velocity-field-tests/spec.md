## ADDED Requirements

### Requirement: Pure shear velocity field validation
The test suite SHALL verify that velocity fields Vx and Vz match the imposed pure-shear boundary conditions: Vx = ε̇·x on boundaries.

#### Scenario: Vx at boundary matches imposed strain rate
- **WHEN** a pure-shear simulation runs with known ε̇ on a non-dimensional domain [-0.5, 0.5]
- **THEN** `max(Vx)` SHALL be approximately ε̇ × 0.5 (within 20%)

#### Scenario: Vz at boundary matches imposed strain rate
- **WHEN** a pure-shear simulation runs with known ε̇
- **THEN** `min(Vz)` SHALL be approximately -ε̇ × 0.5 (within 20%)

#### Scenario: Pure shear Vx L2 error
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates
- **THEN** the relative L2 error of the Vx field against Vx_ana(x) = −ε̇·x SHALL be less than 5e-2 (with 10:1 viscosity inclusion)
- **NOTE** positive bkg_strain_rate = shortening in x, so Vx = −ε̇·x, Vz = +ε̇·z

#### Scenario: Pure shear Vz L2 error
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates
- **THEN** the relative L2 error of the Vz field against Vz_ana(z) = +ε̇·z SHALL be less than 5e-2 (with 10:1 viscosity inclusion)

### Requirement: No-slip velocity verification
The test suite SHALL verify that velocity is approximately zero at the boundaries when no-slip conditions are applied.

#### Scenario: Boundary velocity is near zero
- **WHEN** a no-slip simulation runs
- **THEN** `max(|Vx|)` at boundary nodes SHALL be less than 1e-10 times the maximum interior velocity
