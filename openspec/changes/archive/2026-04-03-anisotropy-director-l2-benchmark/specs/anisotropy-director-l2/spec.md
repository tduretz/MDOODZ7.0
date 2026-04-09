## ADDED Requirements

### Requirement: Director evolution L2 error against analytical ODE
The `DirectorEvolution` test in `AnisotropyBenchmarkTests.cpp` SHALL compute the L2 error of the numerical director angle field against the analytical closed-form solution θ(t) = arctan(tan(θ₀) − γ̇·t), where θ₀ is the initial director angle, γ̇ is the simple shear rate, and t is the total elapsed time.

#### Scenario: L2 error below threshold at baseline dt
- **WHEN** a homogeneous simple-shear simulation runs with initial director angle θ₀ = 45°, shear rate γ̇ = 1.0 (`bkg_strain_rate=0.5`, `shear_style=1`), dt = 0.05, and Nt = 10 steps (total strain γ = 0.5)
- **THEN** the L2 error between the numerical director angle field (computed from HDF5 `nx`, `nz` via atan2) and the analytical θ_final = arctan(tan(45°) − 1.0·0.5) SHALL be less than 0.01 radians

#### Scenario: L2(Vx) = L2(Vz) symmetry not applicable
- **WHEN** the L2 error is computed over all interior cells
- **THEN** the error SHALL be spatially uniform (all cells have the same analytical value since the field is homogeneous)

### Requirement: Director dt-convergence order
The test suite SHALL include a convergence-order test that runs the director evolution at 3 different time step sizes for the same total time and verifies first-order convergence of the forward Euler integration scheme.

#### Scenario: First-order convergence in dt
- **WHEN** the director evolution runs at dt = 0.1 (Nt=5), dt = 0.05 (Nt=10), and dt = 0.025 (Nt=20) for total time t = 0.5
- **THEN** the convergence order computed from the finest pair (dt=0.05 → dt=0.025) SHALL be ≥ 0.8

#### Scenario: Error decreases with smaller dt
- **WHEN** the L2 errors are computed at all three dt values
- **THEN** L2(dt=0.025) < L2(dt=0.05) < L2(dt=0.1)

### Requirement: Anisotropic stress L2 error against analytical formula
The `StressAnisotropy` test SHALL compute the analytical stress for homogeneous pure shear with transverse isotropy and compare against the numerical stress field. The analytical formula: rotate strain rate to anisotropy frame using director angle, apply σ_nn = 2η·ε_nn and σ_ns = 2η·ε_ns/δ, rotate stress back to lab frame.

#### Scenario: Stress invariant relative error below 2%
- **WHEN** a single-step pure-shear simulation runs with director angle θ=30°, constant viscosity η=1, anisotropy factor δ=6, and strain rate ε̇=1.0
- **THEN** the mean stress invariant τ_II SHALL match the analytical value within 2% relative error

#### Scenario: Mean deviatoric stress matches analytical
- **WHEN** the same simulation completes
- **THEN** the mean sxxd SHALL match the analytical sxxd (computed from the rotation formula) within 2% relative error

### Requirement: Director evolution parameter files for dt convergence
The test suite SHALL include parameter files in `TESTS/AnisotropyBenchmark/` for the dt-convergence test.

#### Scenario: Coarse dt parameter file
- **WHEN** the test build system copies `AnisotropyBenchmark/`
- **THEN** there SHALL be `DirectorEvolution_dt01.txt` with dt=0.1, Nt=5, total time=0.5

#### Scenario: Fine dt parameter file
- **WHEN** the test build system copies `AnisotropyBenchmark/`
- **THEN** there SHALL be `DirectorEvolution_dt025.txt` with dt=0.025, Nt=20, total time=0.5

### Requirement: Documentation of analytical solutions
The analytical solutions for the anisotropy benchmarks SHALL be documented in `TESTS/AnalyticalSolutions.md` and referenced from the relevant skills.

#### Scenario: AnalyticalSolutions.md contains anisotropy section
- **WHEN** reviewing `TESTS/AnalyticalSolutions.md`
- **THEN** there SHALL be a section documenting: the director evolution ODE (dθ/dt = −γ̇·cos²θ), the closed-form solution (θ(t) = arctan(tan(θ₀) − γ̇·t)), the transverse isotropy stress formula, measured L2 errors, and convergence orders

#### Scenario: skill-testing-guide references anisotropy benchmarks
- **WHEN** reviewing `.github/skills/skill-testing-guide/SKILL.md`
- **THEN** the anisotropy benchmark tests SHALL be listed in the test coverage section with their L2 error thresholds

#### Scenario: skill-anisotropy references benchmark accuracy
- **WHEN** reviewing `.github/skills/skill-anisotropy/SKILL.md`
- **THEN** there SHALL be a reference to the quantitative benchmark tests and measured director evolution accuracy
