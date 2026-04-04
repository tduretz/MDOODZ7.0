## MODIFIED Requirements

### Requirement: Anisotropic stress invariant vs orientation angle
The CI test suite SHALL include GTest cases that run single-step pure-shear simulations with prescribed director angles and verify the computed stress against the analytical formula: rotate strain rate to the principal anisotropy plane, apply D = 2η_N·diag(1, 1, 1/δ), rotate stress back. Verification SHALL include both mean-value checks and spatial L2 error norms.

#### Scenario: Stress invariant matches analytical at reference angle
- **WHEN** a single-step pure-shear simulation runs with director angle θ=30°, η_N=1, δ=6, and known strain rate
- **THEN** the computed stress invariant τ_II SHALL match the analytical value within 0.01%

#### Scenario: Stress invariant is angle-dependent
- **WHEN** a single-step pure-shear simulation runs at director angle θ=0° and again at θ=45°
- **THEN** the stress invariants SHALL differ, confirming orientation dependence of the anisotropic constitutive law

#### Scenario: Spatial L2 norms validate full stress field
- **WHEN** the single-step pure-shear simulation with θ=30° completes and stress fields are read from HDF5
- **THEN** the spatial L2 error norms for sxxd, szzd, and sxz SHALL each be less than 1e-6
