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
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates on a 21×21 grid with a 10:1 viscosity inclusion
- **THEN** the relative L2 error of the Vx field against Vx_ana(x) = −ε̇·x SHALL be less than 3e-2
- **NOTE** Tightened from 5e-2 to 3e-2 based on measured L2 ≈ 0.024 (1.25× margin)

#### Scenario: Pure shear Vz L2 error
- **WHEN** a pure-shear simulation runs with known ε̇ and domain coordinates on a 21×21 grid with a 10:1 viscosity inclusion
- **THEN** the relative L2 error of the Vz field against Vz_ana(z) = +ε̇·z SHALL be less than 3e-2
- **NOTE** Tightened from 5e-2 to 3e-2 based on measured L2 ≈ 0.024 (1.25× margin)

### Requirement: No-slip velocity verification
The test suite SHALL verify that velocity is approximately zero at the boundaries when no-slip conditions are applied.

#### Scenario: Boundary velocity is near zero
- **WHEN** a no-slip simulation runs
- **THEN** `max(|Vx|)` at boundary nodes SHALL be less than 1e-10 times the maximum interior velocity

### Requirement: Simple shear velocity field validation
The test suite SHALL verify that the velocity field under periodic simple shear with homogeneous viscosity matches the analytical linear profile Vx(z) = γ̇·z, Vz = 0.

#### Scenario: Simple shear Vx L2 error
- **WHEN** a simple shear simulation runs with `shear_style=1`, `periodic_x=1`, `bkg_strain_rate=1`, homogeneous viscosity (η=1 for all phases, no inclusion), on a 21×21 non-dimensional domain [−0.5, 0.5]²
- **THEN** the relative L2 error of the Vx field against Vx_ana(z) = γ̇·z SHALL be less than 1e-6
- **NOTE** A linear profile is reproduced exactly by the FD stencil; the residual is solver/floating-point noise only

#### Scenario: Simple shear Vz near zero
- **WHEN** a simple shear simulation runs with the same homogeneous setup
- **THEN** max(|Vz|) SHALL be less than 1e-6

#### Scenario: Simple shear boundary Vx values
- **WHEN** a simple shear simulation runs with `bkg_strain_rate=1` on domain [−0.5, 0.5]
- **THEN** Vx at the bottom boundary SHALL be approximately −0.5 (within 1%) and Vx at the top boundary SHALL be approximately +0.5 (within 1%)
