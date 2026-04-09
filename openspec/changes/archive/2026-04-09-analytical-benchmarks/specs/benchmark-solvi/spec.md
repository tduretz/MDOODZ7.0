## ADDED Requirements

### Requirement: SolVi L2 error benchmark
The CI test suite SHALL include a GTest case that runs the SolVi benchmark (circular viscous inclusion under pure shear, Schmid & Podladchikov 2003) and measures L2 error norms for velocity and pressure against the analytical solution `eval_anal_Dani()`.

#### Scenario: L2 error for Vx is below threshold
- **WHEN** the SolVi benchmark runs at 51×51 resolution with viscosity contrast mc/mm = 1000, inclusion radius = 0.2, and pure shear boundary conditions derived from the analytical solution
- **THEN** the relative L2 error for Vx SHALL be less than 1e-2

#### Scenario: L2 error for Vz is below threshold
- **WHEN** the SolVi benchmark runs at 51×51 resolution under the same configuration
- **THEN** the relative L2 error for Vz SHALL be less than 1e-2

#### Scenario: L2 error for pressure is below threshold
- **WHEN** the SolVi benchmark runs at 51×51 resolution under the same configuration
- **THEN** the relative L2 error for P SHALL be less than 1e-1

#### Scenario: L2 error for deviatoric stress is below threshold
- **WHEN** the SolVi benchmark runs at 51×51 resolution under the same configuration
- **THEN** the relative L2 error for σ'_xx SHALL be less than 1e-1

### Requirement: SolVi grid-convergence order
The CI test suite SHALL include a GTest case that runs the SolVi benchmark at 3 resolutions (21×21, 41×41, 81×81) and verifies that the velocity L2 error converges at approximately second order in the grid spacing.

#### Scenario: Velocity converges at second order
- **WHEN** the SolVi benchmark runs at resolutions 21×21, 41×41, and 81×81 under the same configuration
- **THEN** the convergence order $p = \log(L_2^{coarse} / L_2^{fine}) / \log(h_{coarse} / h_{fine})$ for Vx SHALL be ≥ 1.8

#### Scenario: Pressure converges at near-second order
- **WHEN** the SolVi benchmark runs at resolutions 21×21, 41×41, and 81×81
- **THEN** the convergence order for P SHALL be ≥ 1.5

### Requirement: SolVi analytical boundary conditions
The SolVi test `.cpp` file SHALL apply the analytical solution as boundary conditions: Dirichlet (type=0) on E/W faces and analytical velocity (type=11) on N/S faces, matching the pattern in `SETS/AnisotropyDabrowski.c`.

#### Scenario: Boundary conditions use analytical solution values
- **WHEN** the SolVi test setup callbacks are invoked for boundary positions
- **THEN** SetBCVx and SetBCVz SHALL call the ported `eval_anal_Dani()` function to compute the boundary velocity values at each coordinate

### Requirement: SolVi parameter files
The SolVi test suite SHALL include `.txt` parameter files for each resolution (21×21, 41×41, 51×51, 81×81) in a `TESTS/SolViBenchmark/` directory. Each file SHALL configure: single time step (`Nt=1`), no thermal coupling, viscosity contrast mc=1e3, inclusion radius via `user1`.

#### Scenario: Parameter files exist for all resolutions
- **WHEN** the test build system copies the `SolViBenchmark/` directory
- **THEN** there SHALL be separate `.txt` files for resolutions 21, 41, 51, and 81

### Requirement: SolVi reads coordinates from HDF5
The L2 error computation SHALL use the coordinate arrays stored in the HDF5 `Model/` group (`xc_coord`, `zc_coord` for centers; `xg_coord`, `zg_coord` for vertices) to reconstruct each grid point's physical position for the analytical evaluation.

#### Scenario: Coordinates come from HDF5 output
- **WHEN** L2 error is computed for any field
- **THEN** the test SHALL read coordinate arrays from the HDF5 file rather than recomputing them from domain parameters
