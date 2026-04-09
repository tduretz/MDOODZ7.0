## ADDED Requirements

### Requirement: Director evolution under simple shear
The CI test suite SHALL include a GTest case that runs a homogeneous simple-shear simulation with anisotropy enabled and verifies the director field against the analytical Mühlhaus director equation: ṅ = Wn − Dn(nᵀn) + (nᵀDn)n.

#### Scenario: Director L2 error is below threshold after multiple steps
- **WHEN** a homogeneous simple-shear simulation runs for ~50 time steps with constant anisotropy factor δ=10, initial director at 45°, and simple shear rate γ̇=1
- **THEN** the L2 error between the numerical director field and the analytically integrated director SHALL be less than 1e-2

#### Scenario: Director rotates toward the shear plane
- **WHEN** the simulation completes all time steps under simple shear
- **THEN** the director Ny component SHALL have decreased from its initial value (director rotates toward the shear plane)

### Requirement: Anisotropic stress invariant vs orientation angle
The CI test suite SHALL include a GTest case that runs single-step pure-shear simulations with prescribed director angles and verifies the computed stress invariant against the analytical formula: rotate strain rate to the principal anisotropy plane, apply D = 2η_N·diag(1, 1, 1/δ), rotate stress back.

#### Scenario: Stress invariant matches analytical at reference angle
- **WHEN** a single-step pure-shear simulation runs with director angle θ=30°, η_N=1, δ=3, and known strain rate
- **THEN** the computed stress invariant τ_II SHALL match the analytical value within 5%

#### Scenario: Stress invariant is angle-dependent
- **WHEN** a single-step pure-shear simulation runs at director angle θ=0° and again at θ=45°
- **THEN** the stress invariants SHALL differ, confirming orientation dependence of the anisotropic constitutive law

### Requirement: Anisotropy benchmark parameter files
The anisotropy benchmark SHALL include `.txt` parameter files in `TESTS/AnisotropyBenchmark/` configuring: anisotropy enabled (`anisotropy=1`), appropriate anisotropy factor, homogeneous single-phase domain, and the correct boundary conditions for simple shear or pure shear.

#### Scenario: Simple shear parameter file exists
- **WHEN** the test build system copies the `AnisotropyBenchmark/` directory
- **THEN** there SHALL be a `.txt` file configuring homogeneous simple shear with anisotropy enabled

#### Scenario: Pure shear parameter file exists for stress-vs-angle test
- **WHEN** the test build system copies the `AnisotropyBenchmark/` directory
- **THEN** there SHALL be a `.txt` file configuring single-step pure shear with anisotropy enabled
