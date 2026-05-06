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
