## Requirements

### Requirement: Analytical L2 benchmark for paleowattmeter steady-state grain size

The CI test suite SHALL include a GoogleTest case that runs a homogeneous pure-shear simulation with calcite paleowattmeter grain-size evolution (`gs = 10`) and L2-compares the resulting grid grain-size field against the closed-form analytical steady state $d_{ss} = (B_g\,\dot\varepsilon_{II}\,\tau_{II}\,p / A_g)^{-1/(p+1)}$ derived from the local equilibrium of grain growth and stress-driven reduction.

#### Scenario: L2 error of grain size is below threshold

- **WHEN** a single-step (`Nt = 1`) homogeneous pure-shear simulation runs with a single phase (no inclusion), `elastic = 0`, `thermal = 0`, dislocation creep + paleowattmeter `gs = 10`, fixed temperature, and `bkg_strain_rate` chosen so that $d_{ss}$ falls in $[10^{-6}, 10^{-2}]$ m, and `mesh.d_n` is read from `Output00001.gzip.h5`
- **THEN** the relative L2 error between the numerical `mesh.d_n` array (size $(N_x-1)\times(N_z-1)$) and a constant analytical vector at $d_{ss}$ SHALL be less than `1e-4`

#### Scenario: Mean grain size matches analytical to 1 percent

- **WHEN** the same simulation output is used
- **THEN** the mean of `mesh.d_n` over the interior cells SHALL be within 1% of the analytical $d_{ss}$ in SI units

### Requirement: Benchmark detects NaN propagation from the local Newton iteration

The CI test SHALL fail loudly if any cell of `mesh.d_n` is non-positive or non-finite, so that a regression in `LocalIterationViscoElasticGrainSize` (NaN, Inf, negative grain size) is caught at the test layer rather than via a downstream "Cell went empty" exit.

#### Scenario: Grain size is positive and finite everywhere

- **WHEN** the benchmark simulation completes and `mesh.d_n` is read from `Output00001.gzip.h5`
- **THEN** the minimum of `mesh.d_n` over interior cells SHALL be strictly greater than 0 and finite

### Requirement: Benchmark detects accidental loss of homogeneity

The CI test SHALL fail if `mesh.d_n` deviates from spatial uniformity by more than 0.1% (max/min ratio), so that interpolation bugs or boundary artifacts that introduce spurious spatial structure into a homogeneous problem are caught.

#### Scenario: Grain size field is uniform within tolerance

- **WHEN** the benchmark simulation completes and `mesh.d_n` is read from `Output00001.gzip.h5`
- **THEN** `(max(d_n) / min(d_n)) - 1` SHALL be less than `1e-3`

### Requirement: Benchmark uses the established L2 helper and HDF5 reader

The benchmark SHALL use the existing `computeL2Error` helper from `TESTS/TestHelpers.h` and the existing `getMin/MaxFieldValue` HDF5-reader helpers, with no new analytical-error or HDF5 utility code added at the test-helper layer.

#### Scenario: No new test-helper utilities added

- **WHEN** the benchmark is added to `TESTS/RheologyCreepTests.cpp` (or a new `GrainSizeTests.cpp`)
- **THEN** the implementation SHALL call only existing helpers from `TestHelpers.h` and SHALL NOT introduce new files under `TESTS/` named `*Helpers.h`, `*Helpers.cpp`, or similar

### Requirement: Benchmark documented in `TESTS/AnalyticalSolutions.md`

The new test SHALL be documented as a numbered section in `TESTS/AnalyticalSolutions.md` (matching the format of §3.4 Stress L2 Error Norm), with the analytical formula, parameter values, expected L2 magnitude, and code-assertion block, and SHALL appear as a row in the summary table at the bottom of that document.

#### Scenario: AnalyticalSolutions.md documents the benchmark

- **WHEN** `TESTS/AnalyticalSolutions.md` is read
- **THEN** it SHALL contain a section titled with "Grain Size" that quotes the closed-form $d_{ss}$ formula, lists the constants $K_g$, $Q_g$, $\gamma$, $\lambda$, $c_g$, $p$ used by `gs = 10`, states the parameter file path, gives the expected L2 magnitude, and includes the GTest assertion block

#### Scenario: Summary table contains a row for the benchmark

- **WHEN** the summary table at the bottom of `TESTS/AnalyticalSolutions.md` is read
- **THEN** it SHALL contain a row for the grain-size L2 test with columns Test, Source, Analytical formula, Assertion, Threshold

### Requirement: Benchmark generates a gnuplot visualisation

The benchmark SHALL emit a `.dat` file as a test side-effect (containing the MDOODZ measurement and the analytical reference) and a checked-in gnuplot script that renders a PNG comparing them, following the pattern established by §3.2 director (`director_convergence.dat` + `plot_director_convergence.gp`) and §7 SolCx (`solcx_convergence.dat` + `plot_solcx.gp`). The PNG SHALL be embedded in the new `TESTS/AnalyticalSolutions.md` section via a markdown image link.

The plot SHALL show the paleowattmeter $d_{ss}(\tau_{II})$ analytical curve over a stress range that brackets the test point, with the single MDOODZ measurement overlaid as a marker. This is the canonical Austin & Evans (2002) paleowattmeter plot and visually validates that the MDOODZ point lands on the analytical curve.

#### Scenario: Test writes a `.dat` file with MDOODZ and analytical values

- **WHEN** the grain-size L2 GTest case completes
- **THEN** a file SHALL exist in the test working directory containing at minimum: the prescribed $\dot\varepsilon_{II}$, prescribed $T$, computed $\tau_{II}$, analytical $d_{ss}$, and measured mean of `mesh.d_n` — in a gnuplot-readable column format

#### Scenario: A gnuplot script renders the comparison PNG

- **WHEN** `gnuplot TESTS/RheologyCreep/plot_grain_size.gp` is run from the build directory after the test
- **THEN** a PNG SHALL be produced (e.g. `grain_size_benchmark.png`) showing the analytical $d_{ss}(\tau_{II})$ curve and the overlaid MDOODZ measurement on log-log axes (stress vs grain size)

#### Scenario: PNG is embedded in `TESTS/AnalyticalSolutions.md`

- **WHEN** the new "Grain Size" section in `TESTS/AnalyticalSolutions.md` is rendered
- **THEN** it SHALL contain an image link of the form `![<alt>](<path-to-png>)` pointing at the generated PNG, sibling to the existing `director_benchmark.png`, `thermoelastic_benchmark.png`, and `solcx_convergence.png` references

### Requirement: PinchSwellGSE scenario stops crashing at initialization

The `PinchSwellGSE` scenario SHALL run to completion at its configured `Nt` without exiting via `LOG_ERR("Cell went empty!!! Exiting...")` or any other NaN-driven failure. Physical correctness of the scenario itself is validated by the analytical L2 benchmark above, not by Pinch&Swell — this requirement only states that the scenario is no longer broken.

#### Scenario: PinchSwellGSE completes its time loop

- **WHEN** `cmake-exec/PinchSwellGSE/PinchSwellGSE` is run with `cmake-exec/PinchSwellGSE/PinchSwellGSE.txt`
- **THEN** the process SHALL exit with status 0 and `Output*.gzip.h5` files SHALL exist for at least one time step beyond initialization

### Requirement: Dead `gsel` switch removed from input parser

The per-phase `gsel` switch SHALL be removed from `MDLIB/InputOutput.c` so that the misleading "Grain size evolution activated" log line is no longer emitted. `gsel = N` lines in existing `.txt` files SHALL be silently ignored (no parser warning, no functional effect — `ReadMatProps` only warns on *expected-but-missing* parameters, not on unknown lines present in the file). The actual grain-size-evolution control SHALL remain the per-phase flow-law index `gs` (used at `RheologyDensity.c:522`).

#### Scenario: Input file no longer logs "Grain size evolution activated"

- **WHEN** any `.txt` parameter file is loaded by MDOODZ (with or without a `gsel = N` line)
- **THEN** the log SHALL NOT contain the line "Grain size evolution activated"

#### Scenario: GSE remains armed by `gs` flow-law index

- **WHEN** any `.txt` parameter file sets `gs = 10` for some phase (and no `gsel` line)
- **THEN** the local Newton iteration in `LocalIterationViscoElasticGrainSize` SHALL be invoked for cells of that phase

### Requirement: NaN in the local Newton iteration fails fast with a diagnostic

If the local Newton iteration in `LocalIterationViscoElasticGrainSize` produces a non-finite `eta_ve`, the simulation SHALL exit immediately via `LOG_ERR` with a message that includes the phase, temperature, strain-rate invariant, and last `eta_ve`. The simulation SHALL NOT propagate the NaN into `mesh.d_n` and rely on the downstream "Cell went empty" check.

#### Scenario: NaN-producing iteration exits immediately

- **WHEN** the local Newton iteration in `LocalIterationViscoElasticGrainSize` produces an `eta_ve` for which `isnan(eta_ve)` or `isinf(eta_ve)` is true
- **THEN** MDOODZ SHALL emit a `LOG_ERR` containing the phase index, temperature, `Eii`, and last `eta_ve`, and SHALL exit with non-zero status before the harmonic average at `RheologyDensity.c:1294` runs

### Requirement: Analytical L2 benchmark covers the coupled dislocation+diffusion regime

The CI test suite SHALL include a 5-strain-rate sweep test (`RheologyCreep.GrainSizeSweepCoupled`) that runs the calcite paleowattmeter scenario with **both** dislocation creep (`pwlv = 15`, Renner et al. 2002) and diffusion creep (`linv = 15`, Calcite Herwegh 2003) active simultaneously, in addition to the wattmeter (`gs = 10`). The test SHALL L2-compare the resulting grid grain-size field at each of five strain rates against an **analytical coupled steady state** computed inside the test by solving the 1D fixed-point on $\tau_{II}$ that satisfies $\dot\varepsilon_{II}^{tot} = \dot\varepsilon_{II}^{pwl}(\tau_{II}) + \dot\varepsilon_{II}^{lin}(\tau_{II}, d_{ss})$ with $d_{ss} = (B_g\dot\varepsilon_{II}^{pwl}\tau_{II}p/A_g)^{-1/(p+1)}$. The test's solver SHALL use a different iteration scheme than MDOODZ's (e.g. secant on $\tau_{II}$ rather than bisection-then-Newton on $\eta_{ve}$) so that a self-consistent-but-wrong MDOODZ result is exposed by L2 disagreement rather than masked.

#### Scenario: Coupled sweep converges to analytical d_ss across 4 decades

- **WHEN** five single-step (`Nt = 1`) homogeneous pure-shear simulations run with `pwlv = 15`, `linv = 15`, `gs = 10`, `elastic = 0`, `thermal = 0`, fixed temperature, single phase (`user1 = 0`), and `bkg_strain_rate ∈ {10⁻¹⁶, 10⁻¹⁵, 10⁻¹⁴, 10⁻¹³, 10⁻¹²}` s⁻¹
- **THEN** for each strain rate, the relative L2 error between `mesh.d_n` (from `Output00001.gzip.h5`) and a constant-vector at the analytical $d_{ss}$ (solved from the 1D coupled fixed-point at that Eii_total) SHALL be less than `5e-3`
- **AND** the mean of `mesh.d_n` SHALL be within 1% of the analytical $d_{ss}$
- **AND** `min(mesh.d_n) > 0` and finite
- **AND** `max(mesh.d_n) / min(mesh.d_n) − 1 < 1e-3` (homogeneity)

#### Scenario: Test's analytical solver converges

- **WHEN** the test's coupled fixed-point solver runs at any of the five strain rates
- **THEN** it SHALL converge to a $\tau_{II}$ satisfying $|Eii_{total} - Eii_{pwl}(\tau_{II}) - Eii_{lin}(\tau_{II}, d_{ss})| / Eii_{total} < 10^{-12}$ within 50 iterations
- **AND** the resulting $d_{ss}$ SHALL be larger than the dislocation-only $d_{ss}$ at the same Eii_total (because diffusion takes a fraction of the strain rate, leaving less for dislocation, which the wattmeter rewards with a larger steady-state grain size)

#### Scenario: Single base fixture file is used for the coupled sweep

- **WHEN** the test source tree is inspected
- **THEN** exactly one fixture SHALL exist for the coupled sweep, [TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt](TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt), with `pwlv = 15`, `linv = 15`, `gs = 10`
- **AND** the test SHALL run this base file 5 times with `bkg_strain_rate` and `writer_subfolder` injected per iteration via MDLIB's `MutateInput` callback hook
- **AND** the dislocation-only `RheologyCreep.GrainSizeSweep` test SHALL similarly share its base ([TESTS/RheologyCreep/GrainSizeSteadyState.txt](TESTS/RheologyCreep/GrainSizeSteadyState.txt)) — no per-strain-rate fixture files SHALL exist

#### Scenario: Coupled-sweep test failure points to the coupled-Jacobian path

- **WHEN** `RheologyCreep.GrainSizeSweepCoupled` fails while `RheologyCreep.GrainSizeSweep` (dislocation only) still passes
- **THEN** the failure SHALL be diagnostic of a regression in the coupled-mechanism Jacobian path of `LocalIterationViscoElasticGrainSize` (specifically the `dfdeta += params.m_lin * Eii_lin * dddeta / d_ve` term and the diffusion strain-rate evaluation), not in the dislocation-only path

### Requirement: 2D smoke test gates the integrated PinchSwellGSE code path

The CI test suite SHALL include a 2D smoke test (`RheologyCreep.PinchSwellGSESmoke`) that runs a downsized version of the Schmalholz & Duretz (2017) pinch-and-swell scenario through MDOODZ's full integrated pipeline (Stokes solve, advection, P2G interpolation, harmonic averaging, particle reseeding) and asserts that grain size remains physically valid and develops the expected spatial heterogeneity.

#### Scenario: Downsized PinchSwellGSE completes with healthy d_n field

- **WHEN** [TESTS/RheologyCreep/PinchSwellGSESmoke.txt](TESTS/RheologyCreep/PinchSwellGSESmoke.txt) — a copy of [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt) with `Nx = Nz = 51` and `Nt = 10` — is run by `RheologyCreep.PinchSwellGSESmoke`
- **THEN** the simulation SHALL exit with status 0 and `Output00010.gzip.h5` SHALL exist
- **AND** all values in `mesh.d_n` SHALL be strictly positive and finite (catches NaN propagation in the integrated code path)

#### Scenario: Heterogeneity emerges in the layer

- **WHEN** the same simulation output is read
- **THEN** `max(mesh.d_n) / min(mesh.d_n)` SHALL exceed 5 (catches "scenario silently produced uniform field" regressions — the wattmeter MUST drive heterogeneity given the cosine-perturbed pinch-swell stress field)

#### Scenario: Grain reduction occurs in the layer phase

- **WHEN** the same simulation output is read
- **THEN** `min(mesh.d_n)` SHALL be less than the layer phase's `gs_ref` (1×10⁻³ m) — confirming the wattmeter actively reduced grain size where stress is highest, not bypassed silently

#### Scenario: Smoke test runtime is bounded

- **WHEN** `RheologyCreep.PinchSwellGSESmoke` runs on a typical CI machine
- **THEN** total wall time SHALL be under 30 seconds (51×51 × 10 steps; expected ~5 s)

### Requirement: Coverage extensions documented in `TESTS/AnalyticalSolutions.md`

The new coupled sweep and 2D smoke tests SHALL be documented in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) as a "Coverage extensions" subsection appended to §8 (rather than as a new top-level section), since they validate the same paleowattmeter physics on different code paths.

#### Scenario: §8 contains a Coverage extensions subsection

- **WHEN** [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** §8 SHALL contain a subsection titled "Coverage extensions" (or equivalent) describing both `RheologyCreep.GrainSizeSweepCoupled` and `RheologyCreep.PinchSwellGSESmoke`, what each catches beyond the headline sweep, the parameter file paths, the analytical or invariant references, and the measured numerical results

#### Scenario: Summary table includes the new tests

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain at least two new rows: one for the coupled sweep (analytical reference: paleowattmeter steady state), and one for the 2D smoke (reference: invariants on `d_n` finite/positive/heterogeneous)
