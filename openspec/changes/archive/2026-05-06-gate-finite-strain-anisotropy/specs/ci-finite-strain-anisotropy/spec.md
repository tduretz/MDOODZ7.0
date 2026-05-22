## ADDED Requirements

### Requirement: Analytical L2 benchmark for finite-strain anisotropy under simple shear

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainSimpleShear`) that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1` and per-phase `ani_fstrain = 1`, and asserts that MDOODZ's measured per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the closed-form analytical solution `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4 + 1)` while unsaturated. The reference is derived from the polar decomposition of the deformation gradient `F(t) = [[1, γ(t)], [0, 1]]` (with `γ = γ_dot · t`), giving `C = F^T·F = [[1, γ], [γ, 1+γ²]]`, eigenvalues `(1 + γ²/2) ± γ·sqrt(γ²/4 + 1)`, and `FS_AR = e1/e2 = λ_max` (since `λ_min · λ_max = 1`).

The test SHALL also write a `.dat` file with per-step `(t, γ, FS_AR_analytical, δ_clamped, δ_mdoodz_mean)` for the gnuplot CI visualisation.

**Note on output indexing**: MDOODZ writes `Output<k>.gzip.h5` at the END of step `k`, but the `mesh.aniso_factor_n` field in that file reflects `FS_AR` computed at the START of step `k` (from the deformation gradient at the END of step `k-1`). So the analytical reference for `Output<k>` is evaluated at `t = (k-1)·dt`.

**Pure shear is NOT exercised in CI**: MDOODZ's `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` configuration segfaults at the start of step 2 (after step 1 writes its output successfully) under multiple `pure_shear_ALE` settings (0, 1, -1). The pure-shear analytical formula `FS_AR(t) = exp(2·ε_dot·t)` is documented in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) for completeness, but evolved pure-shear is out of scope for this gating change. Investigation of the pure-shear segfault is left to a follow-up change. The simple-shear formula tested here exercises the same `FiniteStrainAspectRatio` → `mesh.FS_AR_n` → `AnisoFactorEvolv` → `mesh.aniso_factor_n` code path.

#### Scenario: Simple-shear δ(γ) matches the closed-form FS_AR expression to FP precision

- **WHEN** a multi-step homogeneous simple-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max` set high enough that no step saturates (e.g. `ani_fac_max = 1e6`), `bkg_strain_rate = 1.0` so that `γ_dot = 2·bkg_strain_rate = 2`, `Nt = 5`, `dt = 0.5`, single phase, and `Centers/ani_fac` is read at each output step
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and `1 + γ²/2 + γ·sqrt(γ²/4 + 1)` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector at the same value SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved through P2G interpolation

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4` — homogeneous setup must remain homogeneous after particle-to-grid interpolation

### Requirement: Analytical L2 benchmark for the saturation clamp

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainSaturation`) that explicitly exercises the `min(FS_AR, ani_fac_max)` clamp in [AnisotropyRoutines.c:46-54](MDLIB/AnisotropyRoutines.c#L46) by setting a small `ani_fac_max` and running a simple-shear simulation past the saturation strain `γ_sat` (the strain at which the simple-shear `FS_AR(γ)` formula reaches `ani_fac_max`). Steps before `γ_sat` SHALL match the unbounded FS_AR formula and steps after SHALL be clamped at `ani_fac_max`. This is the entire physical content of the `min()` form and must be gated. Simple shear is used (rather than pure shear) for the same reason as the previous requirement.

#### Scenario: Saturation clamp pins δ at ani_fac_max for γ > γ_sat

- **WHEN** a homogeneous simple-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max = 4`, `bkg_strain_rate = 1.0` (so `γ_dot = 2`), `Nt = 5`, `dt = 0.5` (so γ at output step `k` is `γ_dot · (k-1) · dt ∈ {0, 1, 2, 3, 4}`; saturation strain γ_sat ≈ 1.94 falls between step 2 and step 3)
- **THEN** at output steps 1–2 (γ ∈ {0, 1}, unsaturated): mean of `Centers/ani_fac` SHALL match the analytical FS_AR with relative error less than `1e-6`
- **AND** at output steps 3–5 (γ ∈ {2, 3, 4}, saturated): mean of `Centers/ani_fac` SHALL equal `ani_fac_max = 4.0` with absolute error less than `1e-9`
- **AND** spatial homogeneity (`(max/min) - 1 < 1e-4`) SHALL hold at every step

### Requirement: Single base fixture file with MutateInput per-test overrides

The benchmark SHALL use exactly one shared base fixture file [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) (with `Nb_phases = 1`, `anisotropy = 1`, `ani_fstrain = 1`, sensible defaults for the rest), and SHALL NOT introduce per-test fixture files. Each test case SHALL inject its differentiating parameters (`bkg_strain_rate`, `Nt`, `dt`, `dt0`, `dt_start`, `ani_fac_max`, simple-shear vs pure-shear BC mode, `writer_subfolder`) via MDLIB's `MutateInput` callback hook, matching the established convention in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp) (`mutateSweepInput`) and [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp) (`mutateDirOverride`).

#### Scenario: Exactly one base fixture exists

- **WHEN** the test source tree under [TESTS/AnisotropyBenchmark/](TESTS/AnisotropyBenchmark/) is inspected
- **THEN** exactly one new fixture SHALL exist for the finite-strain anisotropy tests, named `AniFstrainEvolution.txt`
- **AND** no per-test or per-strain-rate variants (e.g. `AniFstrainEvolution_PureShear.txt`, `AniFstrainEvolution_Saturation.txt`) SHALL be added

#### Scenario: Both test cases share the base via MutateInput

- **WHEN** the GoogleTest cases `AniFstrainSimpleShear` and `AniFstrainSaturation` are inspected
- **THEN** each SHALL call `RunMDOODZ` on the single shared `AniFstrainEvolution.txt` and rely on a `MutateInput` callback to override its differentiating parameters
- **AND** the override callback SHALL support per-test injection of `bkg_strain_rate`, `Nt`, `dt`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, and `writer_subfolder`
- **AND** the `dt` override path SHALL set `model.dt`, `model.dt0`, AND `model.dt_start` (to avoid the dt-reset bug documented during the `extend-gse-ci-coverage` debug — `AdvectionRoutines.c:685` resets `dt` from `dt_start`)
- **AND** the `writer_subfolder` override path SHALL `free()` the existing heap pointer and `strdup()` the new value (to avoid the heap-ownership crash documented during the same debug)

### Requirement: Benchmark detects NaN propagation from finite-strain integration

The CI test SHALL fail loudly if any cell of `mesh.aniso_factor_n` (or upstream `mesh.FS_AR_n`) is non-positive or non-finite at any output step, so that a regression in `FiniteStrainAspectRatio` ([RheologyParticles.c:634-678](MDLIB/RheologyParticles.c#L634)) or in the deformation-gradient integration ([Main_DOODZ.c:1275](MDLIB/Main_DOODZ.c#L1275)) is caught at the test layer rather than via a downstream solver failure or silent `δ = 1` everywhere.

#### Scenario: aniso_factor_n is positive and finite everywhere

- **WHEN** any of the three benchmark simulations completes and `mesh.aniso_factor_n` is read from each `Output*.gzip.h5`
- **THEN** at every output step, the minimum of `mesh.aniso_factor_n` over interior cells SHALL be strictly greater than or equal to 1.0 and finite (FS_AR is bounded below by 1 by construction; δ inherits this lower bound)

### Requirement: Benchmark documented in `TESTS/AnalyticalSolutions.md`

The new tests SHALL be documented as a numbered section in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) (matching the format of §8 "Grain Size Evolution"), with the closed-form references for `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear), the derivation from the polar decomposition of `F`, the `min(·, ani_fac_max)` clamp behaviour, the parameter values used by the tests, expected L2 magnitudes, the GTest assertion blocks, and a row in the summary table at the bottom of the document.

#### Scenario: AnalyticalSolutions.md documents the benchmark

- **WHEN** [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain a new numbered section titled with "Finite-Strain Anisotropy" that quotes both closed-form FS_AR formulae, derives them from `F` via `C = F^T·F` and the eigenvalues of `U = √C`, states the `min(FS_AR, ani_fac_max)` saturation form with the [AnisotropyRoutines.c:46-54](MDLIB/AnisotropyRoutines.c#L46) reference, lists the parameter file path, gives expected L2 magnitudes, and includes the GTest assertion blocks for all three cases

#### Scenario: Summary table contains a row for the benchmark

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain at least one new row for the finite-strain anisotropy benchmark with columns Test, Source, Analytical formula, Assertion, Threshold (or whatever columns the existing table uses)

### Requirement: Benchmark generates a CI visualisation `aniso_factor_evolution.png`

The benchmark SHALL emit a `.dat` file as a test side-effect (containing the analytical `FS_AR`, the clamped `δ = min(FS_AR, ani_fac_max)`, and the MDOODZ measurements per step) and a checked-in gnuplot script that renders a PNG comparing them, following the pattern established by §3.2 director (`director_convergence.dat` + `plot_director_convergence.gp`) and §8 grain size (`grain_size_benchmark.dat` + `plot_grain_size.gp`). The PNG SHALL be embedded in the new `TESTS/AnalyticalSolutions.md` section via a markdown image link.

The plot SHALL show, on a single panel, the unbounded `FS_AR(γ)` analytical curve (solid), the clamped MDOODZ form `min(FS_AR, ani_fac_max)` (dashed), the MDOODZ measurements as overlaid markers, and grey reference lines at `y = ani_fac_max` (saturation level) and `x = γ_sat` (saturation strain). This visually validates the knee shape of the `min()` clamp.

#### Scenario: Test writes a `.dat` file with MDOODZ and analytical values

- **WHEN** any of the three GTest cases completes
- **THEN** a file SHALL exist in the test working directory containing at minimum: per output step the time `t` (or strain `γ`), the analytical `FS_AR`, the clamped `δ_ana = min(FS_AR, ani_fac_max)`, and the measured mean of `mesh.aniso_factor_n` — in a gnuplot-readable column format

#### Scenario: A gnuplot script renders the comparison PNG

- **WHEN** `gnuplot TESTS/AnisotropyBenchmark/plot_aniso_factor.gp` (or equivalent path) is run from the build directory after the tests
- **THEN** a PNG `aniso_factor_evolution.png` SHALL be produced showing `FS_AR(γ)`, the clamped form, the MDOODZ measurements, and the saturation-level / saturation-strain reference lines

#### Scenario: PNG is embedded in `TESTS/AnalyticalSolutions.md`

- **WHEN** the new "Finite-Strain Anisotropy" section in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is rendered
- **THEN** it SHALL contain an image link of the form `![<alt>](<path-to-png>)` pointing at `aniso_factor_evolution.png`, sibling to the existing `director_benchmark.png`, `grain_size_benchmark.png`, and `solcx_convergence.png` references

### Requirement: Research note `SETS/AnisoFstrainResearch/hansen2012_comparison.md`

The change SHALL deliver a first-class research-grade markdown record at `SETS/AnisoFstrainResearch/hansen2012_comparison.md` documenting the form mismatch between MDOODZ's `δ = min(FS_AR, ani_fac_max)` and the empirical fit from Hansen, Zimmerman, and Kohlstedt (2012) *Nature* 492, 415–418, doi:10.1038/nature11671 (their Eq. 7: `δ(γ) = 13·tanh(0.25γ) + 1`). The note SHALL co-locate with a runnable scenario file `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` that reproduces the comparison conditions. The research note SHALL be flagged as future-work motivation, NOT used as a CI gate (digitisation noise of Hansen+12 lab data is comparable to the form mismatch we'd be flagging — using it as a CI threshold would conflate "MDOODZ has a bug" with "the model is approximate").

#### Scenario: Research note exists with required content

- **WHEN** `SETS/AnisoFstrainResearch/hansen2012_comparison.md` is read
- **THEN** it SHALL contain (a) the live MDOODZ formula `δ = min(FS_AR, ani_fac_max)` with code references to [AnisotropyRoutines.c:46-54](MDLIB/AnisotropyRoutines.c#L46) and [RheologyParticles.c:634-678](MDLIB/RheologyParticles.c#L634), (b) the Hansen+12 fit `δ(γ) = 13·tanh(0.25γ) + 1` with the asymptote at 14 and half-saturation at γ ≈ 2.2, (c) a quantitative comparison at several γ values (table or short prose) showing both forms diverge in shape but converge in asymptote, (d) a discussion noting the `min()` form is too sharp at the knee and rises too fast at low strain, (e) a future-work section pointing at the commented-out `erfc` blend in `AnisoFactorEvolv` and the option of replacing with a tanh form, explicitly marked **out of scope for this change**, and (f) a full citation to Hansen, Zimmerman, Kohlstedt (2012)
- **AND** it SHALL NOT be referenced from any `TESTS/` file or CI script (it is research, not regression gate)

#### Scenario: Companion scenario exists and runs

- **WHEN** the research directory is inspected
- **THEN** `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` SHALL exist with `anisotropy = 1`, `ani_fstrain = 1`, simple-shear BCs, and parameters chosen to reach γ ≈ 4 (sufficient to span the Hansen+12 fit's saturation region)
- **AND** running `cmake-exec/AniFstrainSimpleShear/AniFstrainSimpleShear` against this `.txt` SHALL exit with status 0 — the scenario is for manual research runs, but it must not be silently broken

### Requirement: Benchmark reuses existing test-helper utilities

The benchmark SHALL use the existing `computeL2Error` helper from [TESTS/TestHelpers.h](TESTS/TestHelpers.h) and the existing `getMin/MaxFieldValue` HDF5-reader helpers, with no new analytical-error or HDF5 utility code added at the test-helper layer. Test-local state structures (e.g. `AniFstrainOverride`) MAY live in file scope inside `AnisotropyBenchmarkTests.cpp` or be promoted to `TestHelpers.h` if and only if the same shape is shared with another test file (decision deferable to implementation).

#### Scenario: No new test-helper utilities added

- **WHEN** the benchmark is added to [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp)
- **THEN** the implementation SHALL call only existing helpers from [TESTS/TestHelpers.h](TESTS/TestHelpers.h) and SHALL NOT introduce new files under `TESTS/` named `*Helpers.h`, `*Helpers.cpp`, or similar
