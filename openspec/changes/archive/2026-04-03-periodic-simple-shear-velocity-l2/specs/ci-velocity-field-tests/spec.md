## ADDED Requirements

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
