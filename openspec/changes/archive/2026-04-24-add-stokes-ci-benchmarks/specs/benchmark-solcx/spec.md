## ADDED Requirements

### Requirement: SolCx L2 error benchmark
The CI test suite SHALL include a GTest case that runs the SolCx benchmark (step-viscosity Stokes flow, Zhong 1996) at 51×51 resolution and measures L2 velocity errors and L1 pressure error against the Velic analytical solution.

#### Scenario: L2 error for Vx is below threshold
- **WHEN** SolCx runs at 51×51 with η = 1 for x < 0.5, η = 10⁶ for x ≥ 0.5, density ρ = sin(πz)·cos(πx), gravity (0, 1), and analytical Dirichlet BCs on all four walls
- **THEN** the relative L2 error for Vx SHALL be below a threshold calibrated post-first-run at 2–5× the measured error

#### Scenario: L2 error for Vz is below threshold
- **WHEN** SolCx runs at 51×51 under the same configuration
- **THEN** the relative L2 error for Vz SHALL be below a threshold calibrated post-first-run at 2–5× the measured error

#### Scenario: L1 error for pressure is below threshold
- **WHEN** SolCx runs at 51×51 under the same configuration
- **THEN** the absolute L1 error for P SHALL be below a threshold calibrated post-first-run at 2–5× the measured error

### Requirement: SolCx grid-convergence order
The CI test suite SHALL include a GTest case that runs SolCx at three resolutions (21×21, 41×41, 81×81) and verifies that velocity and pressure errors decrease at specified convergence rates.

#### Scenario: Velocity converges at better than first order
- **WHEN** SolCx runs at 21×21, 41×41, and 81×81 under the same configuration
- **THEN** the convergence order `p = log(L2_coarse / L2_fine) / log(h_coarse / h_fine)` SHALL be ≥ 1.5 for both Vx and Vz across the 41→81 pair

#### Scenario: Pressure converges at a measurable rate despite the η jump
- **WHEN** SolCx runs at 21×21, 41×41, and 81×81 under the same configuration
- **THEN** the convergence order of L1(P) SHALL be ≥ 0.5 across the 41→81 pair — permissive floor reflecting the known P-order reduction across discontinuous-η interfaces for staggered FD methods

### Requirement: SolCx analytical Dirichlet boundary conditions
Each SolCx test parameter file SHALL apply the Velic analytical velocity as Dirichlet boundary conditions on all four walls, following the precedent of the SolVi benchmark.

#### Scenario: Boundary conditions evaluate the analytical solution
- **WHEN** the SolCx `SetBCVx` and `SetBCVz` callbacks are invoked for any boundary coordinate
- **THEN** they SHALL call the ported Velic `eval_anal_solcx()` function to compute the exact boundary velocity at that coordinate and return it as a type-0 (Dirichlet) boundary condition

### Requirement: SolCx uses the SetGridDensity callback
The SolCx test SHALL use the `SetGridDensity` callback (defined by the `set-grid-density-callback` capability) to prescribe the spatially-varying body-force density `ρ(x, z) = sin(π·z)·cos(π·x)`. The corresponding `.txt` parameter files SHALL declare `density_model = 0` on every phase so that the startup-validation rule in `set-grid-density-callback` is satisfied.

#### Scenario: SetGridDensity is registered and returns sin(πz)·cos(πx)
- **WHEN** the SolCx `MdoodzSetup` is constructed
- **THEN** `SetParticles_ff.SetGridDensity` SHALL be set to a function that returns `sin(π · coordinates.z) · cos(π · coordinates.x)`

#### Scenario: density_model = 0 on every phase
- **WHEN** the SolCx `.txt` files are parsed
- **THEN** `density_model` SHALL be 0 on both phases, satisfying the validation block in `set-grid-density-callback`

#### Scenario: Grid density matches the analytical formula
- **WHEN** the SolCx benchmark's HDF5 output is read after the single Stokes solve completes
- **THEN** for every cell centre `(x_c, z_c)`, `Centers/rho_n[c]` SHALL equal `sin(π · z_c) · cos(π · x_c)` up to scaling-factor precision

### Requirement: SolCx parameter files
The SolCx test suite SHALL include `.txt` parameter files in `TESTS/SolCx/` for resolutions 21×21, 41×41, 51×51, and 81×81. Each file SHALL configure: `Nt = 1`, `Nb_phases = 2`, `eta0[0] = 1`, `eta0[1] = 1e6`, all flow-law indices zero, `density_model = 0` on every phase, `gx = 0`, `gz = 1`, no thermal coupling (`thermal = 0`), no elasticity (`elastic = 0`), no plasticity (`plast = 0`), `writer_subfolder = SolCx<NNN>`, `coh_soft = 0`.

#### Scenario: Parameter files exist for all tested resolutions
- **WHEN** the test build copies `TESTS/SolCx/` into the CMake binary directory
- **THEN** `.txt` files SHALL exist for `Nx ∈ {21, 41, 51, 81}`

#### Scenario: Phase assignment produces the η jump at x = 0.5
- **WHEN** the SolCx `SetPhase` callback is invoked at marker coordinates (x, z)
- **THEN** it SHALL return `0` for `x < 0.5` and `1` for `x ≥ 0.5`, so that the marker inherits `eta0[0] = 1` or `eta0[1] = 1e6` respectively

### Requirement: SolCx analytical solution port
A C++ port of the Velic SolCx reference solution (Chebyshev-series evaluator) SHALL be included at `TESTS/SolCx/SolCxAnalytical.cpp` with a scalar function that returns exact `(Vx, Vz, P)` at any `(x, z) ∈ [0, 1]²`.

#### Scenario: Analytical solution matches tabulated reference values
- **WHEN** the ported analytical function is evaluated at three tabulated reference points from ASPECT's SolCx benchmark test
- **THEN** the returned velocities and pressure SHALL match the ASPECT values to within 1e-10 relative error

### Requirement: SolCx convergence figure
The SolCx benchmark SHALL ship a PNG convergence figure produced by a committed gnuplot script and embedded in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md). The figure SHALL be a log-log L2-error-vs-grid-spacing plot across the three convergence-sweep resolutions, with reference slopes for visual interpretation.

#### Scenario: Gnuplot script is committed and produces the figure
- **WHEN** a developer runs `gnuplot SolCx/plot_solcx.gp` from the CMake build TESTS/ directory after `SolCxBenchmarkTests` has run
- **THEN** the script SHALL produce `SolCx/solcx_convergence.png` from the `solcx_convergence.dat` file that the test writes during `GridConvergence`

#### Scenario: Figure shows log-log convergence with reference slopes
- **WHEN** the generated figure is viewed
- **THEN** it SHALL plot `L2(Vx)`, `L2(Vz)`, and `L1(P)` versus `h = 1/Nx` on log-log axes for the three convergence-sweep resolutions, with dashed reference lines at slopes 1 and 2, and a legend identifying each series

#### Scenario: Figure is embedded in AnalyticalSolutions.md
- **WHEN** `TESTS/AnalyticalSolutions.md` is read
- **THEN** the SolCx section SHALL embed the PNG via `![SolCx grid-convergence plot: ...](SolCx/solcx_convergence.png)` positioned after the "Code Assertions" subsection, with accompanying prose explaining how to regenerate it

#### Scenario: Test writes convergence data file
- **WHEN** `SolCxBenchmarkTests.GridConvergence` completes successfully
- **THEN** it SHALL write `solcx_convergence.dat` to the current directory with four columns: `h L2(Vx) L2(Vz) L1(P)`, one row per resolution tested

### Requirement: SolCx benchmark documentation
The SolCx benchmark SHALL be documented as a new numbered section in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) following the existing SolVi (§1) structure.

#### Scenario: Documentation section exists with the canonical subsections
- **WHEN** `TESTS/AnalyticalSolutions.md` is read
- **THEN** it SHALL contain a "SolCx Benchmark (Step Viscosity)" section with subsections "Problem Statement", "Analytical Solution", "Implementation Details", "Measured L2 Errors", "Convergence Orders", and "Code Assertions"

#### Scenario: Measured values and code assertions are recorded
- **WHEN** the SolCx documentation is finalised after the first successful CI run
- **THEN** "Measured L2 Errors" SHALL contain measured `L2(Vx)`, `L2(Vz)`, `L1(P)` values at each tested resolution (21, 41, 51, 81); "Convergence Orders" SHALL contain measured orders for the 41→81 pair; and "Code Assertions" SHALL reproduce the exact `EXPECT_LT(...)` / `EXPECT_GE(...)` lines from the GTest source with inline margin comments (e.g., `// measured 1.2e-2 (2.1× margin)`)
