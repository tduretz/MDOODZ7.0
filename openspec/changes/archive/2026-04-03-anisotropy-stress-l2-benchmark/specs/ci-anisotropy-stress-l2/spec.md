## ADDED Requirements

### Requirement: Spatial L2 error norm for anisotropic deviatoric stress sxxd
The CI test suite SHALL include a GTest case that computes the spatial L2 error norm of the numerical sxxd field against the closed-form analytical value at every cell-centre node, for a homogeneous pure-shear setup with a fixed director angle.

#### Scenario: sxxd L2 error is below threshold
- **WHEN** a single-step pure-shear simulation runs with director angle θ=30°, η=1, δ=6, Exx=−1 on an 11×11 grid (StressAngle.txt) and the sxxd field is read from HDF5 output
- **THEN** the L2 error between the numerical sxxd array (size (Nx-1)×(Nz-1)) and a constant analytical vector computed via the rotation formula SHALL be less than 1e-6

### Requirement: Spatial L2 error norm for anisotropic deviatoric stress szzd
The CI test suite SHALL include a GTest case that computes the spatial L2 error norm of the numerical szzd field against the closed-form analytical value at every cell-centre node.

#### Scenario: szzd L2 error is below threshold
- **WHEN** the same single-step pure-shear simulation output is used and the szzd field is read from HDF5 output
- **THEN** the L2 error between the numerical szzd array (size (Nx-1)×(Nz-1)) and a constant analytical vector (szzd = −sxxd by incompressibility) SHALL be less than 1e-6

### Requirement: Spatial L2 error norm for anisotropic shear stress sxz
The CI test suite SHALL include a GTest case that computes the spatial L2 error norm of the numerical sxz field against the closed-form analytical value at every vertex node.

#### Scenario: sxz L2 error is below threshold
- **WHEN** the same single-step pure-shear simulation output is used and the sxz field is read from HDF5 output
- **THEN** the L2 error between the numerical sxz array (size Nx×Nz) and a constant analytical vector computed via the rotation formula SHALL be less than 1e-6

### Requirement: L2 test reuses existing simulation and analytical helpers
The StressAnisotropyL2 test SHALL reuse the existing `StressAngle.txt` parameter file and the existing `analyticalStressSxxd` / `analyticalStressSxz` helper functions. No new parameter files or analytical derivations SHALL be added.

#### Scenario: No new parameter files created
- **WHEN** the StressAnisotropyL2 test is added to `AnisotropyBenchmarkTests.cpp`
- **THEN** it SHALL reference only `StressAngle.txt` and no additional `.txt` files SHALL be created in `TESTS/AnisotropyBenchmark/`

#### Scenario: Analytical values come from existing helpers
- **WHEN** the test constructs analytical reference vectors
- **THEN** it SHALL call `analyticalStressSxxd` and `analyticalStressSxz` with the same parameters (η=1, θ=30°, δ=6, Exx=−1) to obtain the constant analytical values
