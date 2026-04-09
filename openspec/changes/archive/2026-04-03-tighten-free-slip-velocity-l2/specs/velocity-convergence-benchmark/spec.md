## ADDED Requirements

### Requirement: Pure shear velocity grid-convergence study
The test suite SHALL verify that the Stokes solver produces the expected convergence rate for pure shear with a 10:1 viscosity inclusion by running at multiple grid resolutions and computing pair-wise convergence orders.

#### Scenario: Multi-resolution L2 error computation
- **WHEN** the pure shear simulation (ε̇ = 1, η_matrix = 1, η_inclusion = 10, r = 0.05) runs at resolutions 21×21, 41×41, and 81×81
- **THEN** the relative L2 error for Vx and Vz against analytical Vx = −ε̇·x, Vz = +ε̇·z SHALL be computed at each resolution

#### Scenario: Monotonic error decrease with resolution
- **WHEN** L2 errors are computed at all three resolutions
- **THEN** L2(Vx) and L2(Vz) SHALL decrease monotonically from coarsest to finest grid

#### Scenario: Convergence order for Vx
- **WHEN** pair-wise convergence orders are computed as log(L2_coarse / L2_fine) / log(h_coarse / h_fine)
- **THEN** the convergence order for Vx SHALL be at least 0.8 for each resolution pair

#### Scenario: Convergence order for Vz
- **WHEN** pair-wise convergence orders are computed
- **THEN** the convergence order for Vz SHALL be at least 0.8 for each resolution pair
