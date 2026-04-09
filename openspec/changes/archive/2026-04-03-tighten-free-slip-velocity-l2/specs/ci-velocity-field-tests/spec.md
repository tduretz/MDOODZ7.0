## MODIFIED Requirements

### Requirement: Pure shear velocity field validation
The test suite SHALL verify that velocity fields Vx and Vz match the imposed pure-shear boundary conditions: Vx = ε̇·x on boundaries.

#### Scenario: Vx at boundary matches imposed strain rate
- **WHEN** a pure-shear simulation runs with known ε̇ on a non-dimensional domain [-0.5, 0.5]
- **THEN** `max(Vx)` SHALL be approximately ε̇ × 0.5 (within 20%)

#### Scenario: Vz at boundary matches imposed strain rate
- **WHEN** a pure-shear simulation runs with known ε̇
- **THEN** `min(Vz)` SHALL be approximately -ε̇ × 0.5 (within 20%)

#### Scenario: Pure shear Vx L2 error
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates on a 21×21 grid with a 10:1 viscosity inclusion
- **THEN** the relative L2 error of the Vx field against Vx_ana(x) = −ε̇·x SHALL be less than 3e-2
- **NOTE** Tightened from 5e-2 to 3e-2 based on measured L2 ≈ 0.024 (1.25× margin)

#### Scenario: Pure shear Vz L2 error
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates on a 21×21 grid with a 10:1 viscosity inclusion
- **THEN** the relative L2 error of the Vz field against Vz_ana(z) = +ε̇·z SHALL be less than 3e-2
- **NOTE** Tightened from 5e-2 to 3e-2 based on measured L2 ≈ 0.024 (1.25× margin)
