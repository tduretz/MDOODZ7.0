# ci-finite-strain-anisotropy Specification

## Purpose

Analytical-benchmark CI coverage and the `ani_fstrain` / `aniso_db` enumeration contract for MDOODZ's finite-strain viscous-anisotropy module: the closed-form δ(FS_AR) validation under simple and pure shear, the saturation clamp, the per-mineral calibrated δ-forms (Hansen olivine, calcite, quartz, and the BRACKETED polyphase estimates), the grain-size-coupled δ modifier, and the `ani_fstrain` dispatch values — including `ani_fstrain = 3` (temperature-dependent δ-relaxation, see the `aniso-delta-relaxation` capability).
## Requirements
### Requirement: Analytical L2 benchmark for finite-strain anisotropy under simple shear

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainSimpleShear`) that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1` and per-phase `ani_fstrain = 1`, and asserts that MDOODZ's measured per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the closed-form analytical solution `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4 + 1)` while unsaturated. The reference is derived from the polar decomposition of the deformation gradient `F(t) = [[1, γ(t)], [0, 1]]` (with `γ = γ_dot · t`), giving `C = F^T·F = [[1, γ], [γ, 1+γ²]]`, eigenvalues `(1 + γ²/2) ± γ·sqrt(γ²/4 + 1)`, and `FS_AR = e1/e2 = λ_max` (since `λ_min · λ_max = 1`).

**Note on output indexing**: MDOODZ writes `Output<k>.gzip.h5` at the END of step `k`, but the `mesh.aniso_factor_n` field in that file reflects `FS_AR` computed at the START of step `k` (from the deformation gradient at the END of step `k-1`). So the analytical reference for `Output<k>` is evaluated at `t = (k-1)·dt`.

The pure-shear analytical reference (`FS_AR(t) = exp(2·ε_dot·t)`) is also exercised in CI by the sibling requirement *Analytical L2 benchmark for finite-strain anisotropy under pure shear*. Both deformation modes share the same `FiniteStrainAspectRatio` → `mesh.FS_AR_n` → `AnisoFactorEvolv` → `mesh.aniso_factor_n` code path; running both gives independent confirmation of the path under the two BC types most commonly used in research scenarios.

#### Scenario: Simple-shear δ(γ) matches the closed-form FS_AR expression to FP precision

- **WHEN** a multi-step homogeneous simple-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max` set high enough that no step saturates (e.g. `ani_fac_max = 1e6`), `bkg_strain_rate = 1.0` so that `γ_dot = 2·bkg_strain_rate = 2`, `Nt = 5`, `dt = 0.5`, single phase, and `Centers/ani_fac` is read at each output step
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and `1 + γ²/2 + γ·sqrt(γ²/4 + 1)` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector at the same value SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved through P2G interpolation

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4` — homogeneous setup must remain homogeneous after particle-to-grid interpolation

### Requirement: Analytical L2 benchmark for finite-strain anisotropy under pure shear

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainPureShear`) that runs a homogeneous prescribed pure-shear simulation with `anisotropy = 1` and per-phase `ani_fstrain = 1`, and asserts that MDOODZ's measured per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the closed-form analytical solution `FS_AR(t) = exp(2·ε̇·t)` while unsaturated. The reference is derived from the polar decomposition of the deformation gradient `F(t) = diag(exp(-ε̇·t), exp(ε̇·t))` (with prescribed pure-shear velocity field `Vx = -x·ε̇`, `Vz = z·ε̇`), giving principal stretches `(e1, e2) = (exp(ε̇·t), exp(-ε̇·t))` and `FS_AR = e1/e2 = exp(2·ε̇·t)`.

The same output-indexing offset documented for the simple-shear test applies: `Output<k>.gzip.h5` reflects `FS_AR` at `t = (k-1)·dt`, because `mesh.aniso_factor_n` in step `k` is computed at the START of the step from `particles.F` at the END of step `k-1`.

#### Scenario: Pure-shear δ(t) matches exp(2·ε̇·t) within Stokes-solve integration tolerance

- **WHEN** a multi-step homogeneous pure-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max` set high enough that no step saturates (e.g. `ani_fac_max = 1e6`), `bkg_strain_rate` small enough for stable ALE mesh evolution (e.g. `Eii = 0.1`), single phase, `Nt = 5`, `dt = 0.5`, and `Centers/ani_fac` is read at each output step
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and `exp(2·ε̇·(k-1)·dt)` SHALL be less than `1e-2`
- **AND** the relative L2 error between the spatial field and a constant analytical vector at the same value SHALL be less than `1e-2`

**Note on threshold**: pure-shear F is integrated through MDOODZ's Stokes solve + advection chain (forward-Euler-class), so residuals grow ~linearly with step count (observed: `~6e-4` at step 2 → `~3e-3` at step 5). This is significantly looser than the simple-shear test's `1e-6` because simple-shear's `F = [[1, γ̇·t], [0, 1]]` is exact-by-construction under periodic BCs while pure-shear's F evolves via Stokes solve + ALE remeshing.

#### Scenario: Spatial homogeneity is preserved through P2G interpolation under pure shear

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4` — homogeneous setup must remain homogeneous after particle-to-grid interpolation, including across the per-step ALE mesh reshape if `pure_shear_ALE > 0`

#### Scenario: Pure-shear simulation runs to completion through Nt steps

- **WHEN** the test executes `RunMDOODZ` on the shared `AniFstrainEvolution.txt` fixture with the pure-shear override (`shear_style = 0`, `periodic_x = 0`, configured `pure_shear_ALE`)
- **THEN** the process SHALL exit with status 0 (no SIGSEGV at start of step 2, no NaN-driven exit, no `LOG_ERR`-triggered abort)
- **AND** `Output00001.gzip.h5` through `Output<Nt>.gzip.h5` SHALL all exist with valid `Centers/ani_fac` datasets

### Requirement: Analytical L2 benchmark for the saturation clamp

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainSaturation`) that explicitly exercises the `min(FS_AR, ani_fac_max)` clamp in [AnisotropyRoutines.c:46-54](../../MDLIB/AnisotropyRoutines.c#L46) by setting a small `ani_fac_max` and running a simple-shear simulation past the saturation strain `γ_sat` (the strain at which the simple-shear `FS_AR(γ)` formula reaches `ani_fac_max`). Steps before `γ_sat` SHALL match the unbounded FS_AR formula and steps after SHALL be clamped at `ani_fac_max`. This is the entire physical content of the `min()` form and must be gated. Simple shear is used (rather than pure shear) for the same reason as the previous requirement.

The saturation test SHALL also write the `.dat` file `aniso_factor_evolution.dat` (columns: step, t, γ, FS_AR_analytical, δ_clamped_analytical, δ_mdoodz_mean) used by the gnuplot CI visualisation. The saturation test owns the `.dat` write because its data exhibits the `min()` clamp visually (which is the entire point of the plot).

#### Scenario: Saturation clamp pins δ at ani_fac_max for γ > γ_sat

- **WHEN** a homogeneous simple-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max = 4`, `bkg_strain_rate = 1.0` (so `γ_dot = 2`), `Nt = 5`, `dt = 0.5` (so γ at output step `k` is `γ_dot · (k-1) · dt ∈ {0, 1, 2, 3, 4}`; saturation strain γ_sat ≈ 1.94 falls between step 2 and step 3)
- **THEN** at output steps 1–2 (γ ∈ {0, 1}, unsaturated): mean of `Centers/ani_fac` SHALL match the analytical FS_AR with relative error less than `1e-6`
- **AND** at output steps 3–5 (γ ∈ {2, 3, 4}, saturated): mean of `Centers/ani_fac` SHALL equal `ani_fac_max = 4.0` with absolute error less than `1e-9`
- **AND** spatial homogeneity (`(max/min) - 1 < 1e-4`) SHALL hold at every step

### Requirement: Single base fixture file with MutateInput per-test overrides

The benchmark SHALL use exactly one shared base fixture file [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](../../TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) (with `Nb_phases = 1`, `anisotropy = 1`, `ani_fstrain = 1`, sensible defaults for the rest), and SHALL NOT introduce per-test fixture files. Each test case SHALL inject its differentiating parameters (`bkg_strain_rate`, `Nt`, `dt`, `dt0`, `dt_start`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, `writer_subfolder`) via MDLIB's `MutateInput` callback hook, matching the established convention in [TESTS/RheologyCreepTests.cpp](../../TESTS/RheologyCreepTests.cpp) (`mutateSweepInput`) and [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp) (`mutateDirectorInput`).

#### Scenario: Exactly one base fixture exists

- **WHEN** the test source tree under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** exactly one new fixture SHALL exist for the finite-strain anisotropy tests, named `AniFstrainEvolution.txt`
- **AND** no per-test or per-strain-rate variants (e.g. `AniFstrainEvolution_PureShear.txt`, `AniFstrainEvolution_Saturation.txt`) SHALL be added

#### Scenario: All test cases share the base via MutateInput

- **WHEN** the GoogleTest cases `AniFstrainSimpleShear`, `AniFstrainSaturation`, and `AniFstrainPureShear` are inspected
- **THEN** each SHALL call `RunMDOODZ` on the single shared `AniFstrainEvolution.txt` and rely on a `MutateInput` callback to override its differentiating parameters
- **AND** the override callback SHALL support per-test injection of `bkg_strain_rate`, `Nt`, `dt`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, and `writer_subfolder`
- **AND** the `dt` override path SHALL set `model.dt`, `model.dt0`, AND `model.dt_start` (to avoid the dt-reset bug — `AdvectionRoutines.c:685` resets `dt` from `dt_start`)
- **AND** the `writer_subfolder` override path SHALL `free()` the existing heap pointer and `strdup()` the new value (to avoid the heap-ownership crash on shutdown when the original heap-owned pointer is overwritten with a string literal)

### Requirement: Benchmark detects NaN propagation from finite-strain integration

The CI test SHALL fail loudly if any cell of `mesh.aniso_factor_n` (or upstream `mesh.FS_AR_n`) is non-positive or non-finite at any output step, so that a regression in `FiniteStrainAspectRatio` ([RheologyParticles.c:634-678](../../MDLIB/RheologyParticles.c#L634)) or in the deformation-gradient integration ([Main_DOODZ.c:1275](../../MDLIB/Main_DOODZ.c#L1275)) is caught at the test layer rather than via a downstream solver failure or silent `δ = 1` everywhere.

#### Scenario: aniso_factor_n is positive and finite everywhere

- **WHEN** any of the benchmark simulations completes and `mesh.aniso_factor_n` is read from each `Output*.gzip.h5`
- **THEN** at every output step, the minimum of `mesh.aniso_factor_n` over interior cells SHALL be greater than or equal to `0.999` (FS_AR is bounded below by 1.0 by construction; the 0.001 slack accommodates the float32 round-trip through HDF5) and finite

### Requirement: Benchmark documented in `TESTS/AnalyticalSolutions.md`

The tests SHALL be documented as a numbered section in [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) (matching the format of §8 "Grain Size Evolution"), with the closed-form references for `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear), the derivation from the polar decomposition of `F`, the `min(·, ani_fac_max)` clamp behaviour, the parameter values used by the tests, expected L2 magnitudes, the GTest assertion blocks, the off-by-one output-indexing note, and the embedded PNG. The section SHALL contain measured-accuracy tables for **both** active CI tests (simple-shear and pure-shear), with parallel structure. The previous "(NOT in CI)" footnote for pure shear (§9.7 in the prior revision) SHALL be replaced by a measured-accuracy table plus a "Historical note (resolved)" subsection describing the two MDLIB bugs that originally caused the deferral and how they were fixed.

#### Scenario: AnalyticalSolutions.md documents both deformation modes as active CI

- **WHEN** [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) is read
- **THEN** §9 SHALL contain measured-accuracy tables for **both** the simple-shear test and the pure-shear test, each with per-step `(γ or t, FS_AR_ana, δ_mean MDOODZ, spatial L2 or rel. error)` columns
- **AND** the prior §9.7 footnote acknowledging the pure-shear segfault as a known MDOODZ bug SHALL be removed (or rewritten as a historical note pointing at the change that resolved it, including code refs to the two MDLIB fixes)
- **AND** the closed-form derivations for both `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear) SHALL remain present

#### Scenario: Summary table contains rows for the benchmark

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain at least one row per active CI test (simple-shear with threshold `1e-6`, saturation with threshold `1e-9`, pure-shear with threshold `1e-2`) with columns Test, Source, Analytical formula, Assertion, Threshold

### Requirement: Benchmark generates a CI visualisation `aniso_factor_evolution.png`

The benchmark SHALL emit a `.dat` file as a test side-effect (containing the analytical `FS_AR`, the clamped `δ = min(FS_AR, ani_fac_max)`, and the MDOODZ measurements per step) and a checked-in gnuplot script that renders a PNG comparing them, following the pattern established by §3.2 director (`director_convergence.dat` + `plot_director_convergence.gp`) and §8 grain size (`grain_size_benchmark.dat` + `plot_grain_size.gp`). The PNG SHALL be embedded in the new `TESTS/AnalyticalSolutions.md` section via a markdown image link.

The plot SHALL show, on a single panel, the unbounded `FS_AR(γ)` analytical curve (solid), the clamped MDOODZ form `min(FS_AR, ani_fac_max)` (dashed), the MDOODZ measurements as overlaid markers, and grey reference lines at `y = ani_fac_max` (saturation level) and `x = γ_sat` (saturation strain). This visually validates the knee shape of the `min()` clamp.

#### Scenario: Test writes a `.dat` file with MDOODZ and analytical values

- **WHEN** the `AniFstrainSaturation` GTest case completes
- **THEN** a file `aniso_factor_evolution.dat` SHALL exist in the test working directory containing at minimum: per output step the time `t`, the strain `γ`, the analytical `FS_AR`, the clamped `δ_ana = min(FS_AR, ani_fac_max)`, and the measured mean of `mesh.aniso_factor_n` — in a gnuplot-readable column format

#### Scenario: A gnuplot script renders the comparison PNG

- **WHEN** `gnuplot TESTS/AnisotropyBenchmark/plot_aniso_factor.gp` (or equivalent path) is run from the build directory after the tests
- **THEN** a PNG `aniso_factor_evolution.png` SHALL be produced showing `FS_AR(γ)`, the clamped form, the MDOODZ measurements, and the saturation-level / saturation-strain reference lines

#### Scenario: PNG is embedded in `TESTS/AnalyticalSolutions.md`

- **WHEN** the "Finite-Strain Anisotropy" section in [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) is rendered
- **THEN** it SHALL contain an image link of the form `![<alt>](<path-to-png>)` pointing at `aniso_factor_evolution.png`, sibling to the existing `director_benchmark.png`, `grain_size_benchmark.png`, and `solcx_convergence.png` references

### Requirement: Hansen-calibrated olivine δ-form is available as `ani_fstrain = 2`

[`AnisoFactorEvolv`](../../MDLIB/AnisotropyRoutines.c#L46) SHALL provide a third form of the saturation map, dispatched when the per-phase parameter `ani_fstrain == 2`:

```
δ(FS_AR) = min( 24.5 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR`  (analytical inversion of the simple-shear identity `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4 + 1)`),
- `M(γ) = M_inf · (1 − exp(−γ / γ_e))`  (sigmoidal saturation of the olivine fabric-strength index),
- The mineral-specific δ(FS_AR) closed form is dispatched via a **per-phase function pointer** `materials->aniso_delta_fn[k]` populated by `ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c). Case 1 = Hansen olivine assigns the static helper `anisoDelta_HansenOlivine`, which encapsulates the calibration constants `M_inf = 0.536`, `γ_e = 3.96` and the linear-fit slope `24.5` directly inside the function body. Future minerals (calcite, quartz, …) add new cases that assign their own static helpers — `AnisotropyRoutines.c` stays mineral-agnostic and only invokes the pointer.
- The case selector is the per-phase input parameter `aniso_db` (default 0; auto-set to 1 when `ani_fstrain == 2` is requested without explicit selection, both in [MDLIB/InputOutput.c](../../MDLIB/InputOutput.c) at parse time and post-`MutateInput` in [MDLIB/Main_DOODZ.c](../../MDLIB/Main_DOODZ.c)).
- Values come from a least-squares fit to the **combined 38-point Hansen-group olivine torsion dataset**: 26 datapoints from Hansen, Zhao, Zimmerman & Kohlstedt (2014) *EPSL* 387, 157 Table 1 — itself a unified reanalysis of Bystricky+00, Zhang & Karato 2005, and Hansen+12a/b/c samples — plus 12 datapoints from Hansen, Warren, Zimmerman & Kohlstedt (2016) *EPSL* 445, 92 Part 1 Fig. 6 + Table 1, from 5 new samples not in Hansen+14.
- the linear coefficient `24.5` and the intercept `+1` are taken **verbatim** from Hansen, Zimmerman & Kohlstedt (2012) *Nature* 492, 415, Fig. 3b: `δ = (24.5 ± 5.5) · M + 1`.

The closed-form δ-asymptote (before the `ani_fac_max` cap) is `δ_∞ = 24.5 · M_inf + 1 = 14.13`, agreeing with Hansen+16 Part 1 Eq. 3-4's stress-exponent-aware direct calculation `δ = (F_s/F_w)^n = (1.39/0.73)^4.1 = 14.02` to **100.8%** — independent two-method validation that the Hansen-olivine asymptote is `δ_∞ ≈ 14`. The `ani_fac_max` cap remains the user-facing safety bound, identical in role and effect to its role under `ani_fstrain == 1`.

The form is **opt-in**. Existing `ani_fstrain ∈ {0, 1}` dispatch arms remain byte-identical — no existing scenario, fixture, or CI test is affected unless the operator deliberately sets `ani_fstrain = 2`.

#### Scenario: ani_fstrain = 2 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called with `FS_AR` set to the analytical simple-shear FS_AR at γ = 5.0 (`FS_AR ≈ 27.299`) and `ani_fstrain = 2`, with `ani_fac_max = 100` (uncapped)
- **THEN** the returned δ SHALL be within `1e-6` of `24.5 · 0.536 · (1 − exp(−5.0 / 3.96)) + 1 ≈ 10.40` (the calibrated formula evaluated symbolically)
- **AND** when called with `FS_AR` at γ = 10.0 (`FS_AR ≈ 102.005`), the returned δ SHALL be within `1e-6` of `24.5 · 0.536 · (1 − exp(−10.0 / 3.96)) + 1 ≈ 13.05`

#### Scenario: ani_fstrain = 2 saturates at δ_∞ ≈ 14.13 for large strain

- **WHEN** `AnisoFactorEvolv` is called with `FS_AR = 1e6` (γ_eff ≈ 1000) and `ani_fstrain = 2`, with `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `24.5 · M_inf + 1 = 14.13`

#### Scenario: ani_fstrain = 2 reduces to δ ≈ 1 at zero strain

- **WHEN** `AnisoFactorEvolv` is called with `FS_AR = 1.0` (γ_eff = 0, isotropic) and `ani_fstrain = 2`, with `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-9` of `1.0` — `M(0) = 0`, so δ = 24.5·0 + 1 = 1

#### Scenario: ani_fac_max cap is honoured under ani_fstrain = 2

- **WHEN** `AnisoFactorEvolv` is called with `FS_AR` chosen such that the uncapped Hansen+12 expression yields δ > `ani_fac_max` (e.g. `FS_AR = 1e6` with `ani_fac_max = 4`)
- **THEN** the returned δ SHALL equal `ani_fac_max` exactly (within `1e-12` absolute), matching the cap semantics of the existing `min` form

#### Scenario: ani_fstrain = 1 dispatch arm is byte-identical to its prior behaviour

- **WHEN** `AnisoFactorEvolv` is called with `ani_fstrain = 1` for any value of `FS_AR ∈ [1, 1e6]` and any `ani_fac_max ∈ [1, 1e6]`
- **THEN** the returned value SHALL equal `min(FS_AR, ani_fac_max)` to bitwise precision — the new dispatch arm SHALL NOT alter the existing arm's IEEE-754 output

### Requirement: CI test `AniFstrainHansenOlivine` exercises the new form

The CI test suite SHALL include a GoogleTest case `AnisotropyBenchmark.AniFstrainHansenOlivine` that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1` and `ani_fstrain = 2` and verifies that the per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the analytical reference (the calibrated formula evaluated symbolically) to FP precision over the simulated strain range.

The test SHALL share the existing base fixture [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](../../TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) via `MutateInput`, consistent with the existing simple-shear, pure-shear, and saturation tests in [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp). No new fixture file SHALL be added.

The reference δ at output step `k` is evaluated from the analytical FS_AR at `γ = γ_dot · (k-1) · dt`, then mapped through `γ_eff → M → δ` using the calibrated constants `M_inf = 0.536`, `γ_e = 3.96` from [MDLIB/AnisotropyRoutines.c](../../MDLIB/AnisotropyRoutines.c). The same off-by-one output-indexing rule documented for the existing simple-shear test applies (`Output<k>.gzip.h5` reflects FS_AR computed at the START of step `k`, i.e. `t = (k-1)·dt`).

#### Scenario: AniFstrainHansenOlivine matches the calibrated formula per step

- **WHEN** the test runs simple-shear with `anisotropy = 1`, `ani_fstrain = 2`, `ani_fac_max = 100` (uncapped), `bkg_strain_rate = 1.0` so `γ_dot = 2`, `Nt = 5`, `dt = 0.5` (γ ∈ {0, 1, 2, 3, 4} at output steps 1..5), single phase
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and the analytical reference `δ_ana(γ) = 24.5 · 0.536 · (1 − exp(−γ_eff(FS_AR(γ)) / 3.96)) + 1` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved under ani_fstrain = 2

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4`

#### Scenario: AniFstrainHansenOlivine reuses the shared base fixture

- **WHEN** the test source under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** no new `.txt` fixture file SHALL be added for this test
- **AND** the test SHALL call `RunMDOODZ` on `AniFstrainEvolution.txt` and rely on a `MutateInput` callback that sets `ani_fstrain = 2`, `ani_fac_max = 100`, and any other per-test parameters

### Requirement: Research note `misc/aniso_fstrain/notes/hansen2012_comparison.md`

The change SHALL deliver a first-class research-grade markdown record at `misc/aniso_fstrain/notes/hansen2012_comparison.md` documenting the form mismatch between MDOODZ's `δ = min(FS_AR, ani_fac_max)` and the empirical fit from Hansen, Zimmerman, and Kohlstedt (2012) *Nature* 492, 415–418, doi:10.1038/nature11671. The note SHALL co-locate with a runnable scenario file `misc/aniso_fstrain/notes/AniFstrainSimpleShear.txt` that reproduces the comparison conditions. The research note SHALL be flagged as future-work motivation for the `ani_fstrain = 1` (`min`) form, NOT used as a CI gate (digitisation noise of Hansen+12 lab data is comparable to the form mismatch we'd be flagging — using it as a CI threshold would conflate "MDOODZ has a bug" with "the model is approximate").

The note SHALL additionally contain a section (§7 or later) documenting the calibrated `ani_fstrain = 2` form delivered by the `add-hansen-olivine-delta-form` change, with: (a) the calibration procedure used to fit `M_inf` and `γ_e` from the **combined 38-point Hansen-group olivine torsion dataset** (26 from Hansen+14 Table 1 + 12 from Hansen+16 Part 1 Fig. 6 + Table 1), (b) the closed-form expressions for `γ_eff(FS_AR)` and `M(γ)`, (c) the resulting δ-asymptote `δ_∞ = 14.13` (linear formula) confirmed by Hansen+16 Eq. 3-4's stress-aware direct calc `δ = (1.39/0.73)^4.1 = 14.02` to **100.8% agreement** (independent two-method validation), and how it compares to typical user-side `ani_fac_max` values (4, 6, 8), (d) a quantitative comparison plot (`misc/aniso_fstrain/img/olivine_hansen_calibration.png`) showing the candidate curves on two panels — Panel (a): M(γ) with all 38 + 3 datapoints and per-dataset fits; Panel (b): δ(γ) with `min`, `erfc`, Hansen+12 3-pt, Hansen+16 12-pt, and combined 38-pt fits overlaid with lab data, (e) the RMS error of each form against the combined 38-point reference (γ ∈ [0, 19]): `min` ≈ 2.87, `erfc` ≈ 3.11, Hansen+12 3-pt ≈ 2.07, Hansen+16 12-pt ≈ 1.93, **combined 38-pt ≈ 1.85**, (f) a note on the **structural caveat** that `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` is a simple-shear inversion and reinterpreting it under other deformation modes is interpretation-questionable, (g) a note on the **stress-orientation caveat** from Hansen+16 Part 2 §4.1.1 / Fig. 10: real δ depends on alignment between principal stress and texture's strong axis (~10× spread), and MDOODZ's scalar δ corresponds to the worst-case (texture-aligned) value, (h) a note on the **ε_c ≈ 1 co-evolution** finding from Hansen+16 Part 2 Eq. 10: texture and grain size evolve on the same strain timescale, consistent with MDOODZ's wattmeter coupling.

The accompanying calibration script SHALL be checked in at `misc/aniso_fstrain/calibrate/calibrate_olivine_hansen.py` (Python, no scipy required — uses a 2D coarse-then-fine grid search to minimise SSE on `M_inf`, `γ_e`), reproducible from a clean Python environment.

#### Scenario: Research note exists with required content

- **WHEN** `misc/aniso_fstrain/notes/hansen2012_comparison.md` is read
- **THEN** it SHALL contain (a) the live MDOODZ formula `δ = min(FS_AR, ani_fac_max)` with code references to [AnisotropyRoutines.c:46-54](../../MDLIB/AnisotropyRoutines.c#L46) and [RheologyParticles.c:634-678](../../MDLIB/RheologyParticles.c#L634), (b) the Hansen+12 fit `δ(γ) = 13·tanh(0.25γ) + 1` with the asymptote at 14 and half-saturation at γ ≈ 2.2, (c) a quantitative comparison at several γ values (table or short prose) showing both forms diverge in shape but converge in asymptote, (d) a discussion noting the `min()` form is too sharp at the knee and rises too fast at low strain, (e) a future-work section pointing at the commented-out `erfc` blend in `AnisoFactorEvolv` and the option of replacing with a tanh form, and (f) a full citation to Hansen, Zimmerman, Kohlstedt (2012)
- **AND** it SHALL contain a section (§7 or later) documenting the calibrated `ani_fstrain = 2` form: the calibration procedure (38-point combined least-squares fit to Hansen+14 Table 1 + Hansen+16 Part 1), the closed-form `M(γ_eff(FS_AR))` mapping, the values `M_inf = 0.536` and `γ_e = 3.96`, the δ-asymptote `δ_∞ = 14.13` (validated to 100.8% by Hansen+16 Eq. 3-4 stress-aware direct calc 14.02), the comparison plot at `misc/aniso_fstrain/img/olivine_hansen_calibration.png`, the RMS-error comparison of all candidate forms (`min` ≈ 2.87, `erfc` ≈ 3.11, Hansen+12 3-pt ≈ 2.07, Hansen+16 12-pt ≈ 1.93, combined 38-pt ≈ 1.85), the simple-shear caveat for `γ_eff(FS_AR)`, the Hansen+16 Part 2 stress-orientation caveat, and the Hansen+16 Part 2 ε_c ≈ 1 co-evolution finding
- **AND** it SHALL reference the calibration script `misc/aniso_fstrain/calibrate/calibrate_olivine_hansen.py` as the reproducible artifact behind the calibration values
- **AND** it SHALL NOT be referenced from any `TESTS/` file or CI script (it remains research, not regression gate)

#### Scenario: Companion scenario exists

- **WHEN** the research directory is inspected
- **THEN** `misc/aniso_fstrain/notes/AniFstrainSimpleShear.txt` SHALL exist with `anisotropy = 1`, `ani_fstrain = 1`, simple-shear BCs, and parameters chosen to reach γ ≈ 4 (sufficient to span the Hansen+12 fit's saturation region)

#### Scenario: Calibration script and plot are checked in

- **WHEN** `misc/aniso_fstrain/` is inspected
- **THEN** `calibrate_hansen_olivine.py` SHALL exist as an executable Python script that, when run, prints the per-dataset fit summary table (Hansen+12, Hansen+14, Hansen+16 Part 1, combined), both candidate δ-asymptotes (linear formula and Hansen+16 Eq. 3-4 stress-aware) with their percentage agreement, and the RMS errors of all candidate δ-forms vs the combined 38-point dataset, AND writes `hansen_olivine_calibration.png` to the same directory
- **AND** the script SHALL NOT depend on scipy (numpy + matplotlib only)

### Requirement: Calcite-calibrated δ-form is available as `aniso_db = 2`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a calcite-calibrated entry, dispatched when the per-phase parameter `aniso_db == 2` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( a · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` (analytical inversion of the simple-shear identity, mineral-independent),
- `M(γ) = M_inf · (1 − exp(−γ / γ_e))` (saturating exponential),
- `M_inf = 0.41` and `γ_e = 3.71` (free LSQ fit to the combined 18-point Pieri+01 + Bruijn+11 + Barnhoorn+04 calcite torsion dataset, J→M-converted via Skemer+05; Barnhoorn high-J samples filtered for J ≤ 15 where the linear conversion is valid),
- `a = 6.0` (calcite δ-vs-M slope; bracketed by Tommasi+09-style VPSC bounds for calcite saturated CPO ≈ 2-5).

The δ-vs-M slope `a = 6.0` is mineral-specific (calcite-only), distinct from the olivine `a = 24.5` taken from Hansen+12 Fig 3b. The slope reflects calcite's lower viscous-anisotropy-per-fabric-strength ratio (Barnhoorn+04: CPO contributes ~1/3 of calcite weakening, the rest from grain-size reduction — already captured by MDOODZ's wattmeter).

The closed-form δ-asymptote (before the `ani_fac_max` cap) is `δ_∞ = a · M_inf + 1 = 3.46`, significantly lower than olivine's 14.13.

The form is **opt-in** via explicit `aniso_db = 2` per phase. The auto-default (`ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1`) remains **olivine-specific** — flipping it on existing scenarios would be a silent semantic change. Calcite scenarios MUST set `aniso_db = 2` explicitly. Existing `ani_fstrain ∈ {0, 1, 3}` and `aniso_db = 1` (olivine) dispatch arms remain byte-identical.

#### Scenario: aniso_db = 2 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer (set by `ReadDataAnisotropy(case 2)`) with `FS_AR` set to the analytical simple-shear FS_AR at γ = 5.0 (`FS_AR ≈ 27.299`), `ani_fstrain = 2`, and `ani_fac_max = 100` (uncapped)
- **THEN** the returned δ SHALL be within `1e-6` of `6.0 · 0.41 · (1 − exp(−5.0 / 3.71)) + 1 ≈ 2.82`
- **AND** when called with `FS_AR` at γ = 10.0 (`FS_AR ≈ 102.005`), the returned δ SHALL be within `1e-6` of `6.0 · 0.41 · (1 − exp(−10.0 / 3.71)) + 1 ≈ 3.29`

#### Scenario: aniso_db = 2 saturates at δ_∞ ≈ 3.46 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR = 1e6` (γ_eff ≈ 1000), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `6.0 · M_inf + 1 = 3.46`

#### Scenario: aniso_db = 2 reduces to δ ≈ 1 at zero strain

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR = 1.0` (γ_eff = 0, isotropic), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-9` of `1.0` — `M(0) = 0`, so δ = 6·0 + 1 = 1

#### Scenario: ani_fac_max cap is honoured under aniso_db = 2

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR` chosen such that the uncapped calcite formula yields δ > `ani_fac_max` (e.g. `FS_AR = 1e6` with `ani_fac_max = 2`)
- **THEN** the returned δ SHALL equal `ani_fac_max` exactly (within `1e-12` absolute), matching the cap semantics of all other dispatch arms

#### Scenario: aniso_db = 1 (olivine) dispatch arm is unaffected by the calcite addition

- **WHEN** `AnisoFactorEvolv` is called via the olivine function pointer (set by `ReadDataAnisotropy(case 1)`) for any value of `FS_AR ∈ [1, 1e6]`, `ani_fstrain = 2`, `ani_fac_max ∈ [1, 1e6]`
- **THEN** the returned δ SHALL equal the value the previous (olivine-only) implementation would have returned, to bitwise precision — adding case 2 SHALL NOT alter case 1's IEEE-754 output

#### Scenario: Auto-default for ani_fstrain = 2 without aniso_db remains olivine

- **WHEN** an input file (or a `MutateInput` callback) sets `ani_fstrain[k] = 2` without setting `aniso_db[k]` (left at default 0)
- **THEN** the auto-default in [MDLIB/InputOutput.c](../../MDLIB/InputOutput.c) at parse time AND in [MDLIB/Main_DOODZ.c](../../MDLIB/Main_DOODZ.c) post-`MutateInput` SHALL set `aniso_db[k] = 1` (Hansen olivine), NOT `aniso_db[k] = 2`
- **AND** the calcite calibration SHALL only be selected when the user OR the `MutateInput` callback explicitly sets `aniso_db[k] = 2`

### Requirement: CI test `AniFstrainCalcite` exercises the calcite form

The CI test suite SHALL include a GoogleTest case `AnisotropyBenchmark.AniFstrainCalcite` that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 2` and verifies that the per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the analytical reference (the calibrated calcite formula evaluated symbolically) to FP precision over the simulated strain range.

The test SHALL share the existing base fixture [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](../../TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) via `MutateInput`, consistent with the existing `AniFstrainSimpleShear`, `AniFstrainPureShear`, `AniFstrainSaturation`, and `AniFstrainHansenOlivine` tests in [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp). No new fixture file SHALL be added.

The reference δ at output step `k` is evaluated from the analytical FS_AR at `γ = γ_dot · (k-1) · dt`, then mapped through `γ_eff → M → δ` using the calibrated constants `M_inf = 0.41`, `γ_e = 3.71`, `slope = 6.0` from the static helper `anisoDelta_Calcite` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c). The same off-by-one output-indexing rule documented for the existing simple-shear test applies (`Output<k>.gzip.h5` reflects FS_AR computed at the START of step `k`, i.e. `t = (k-1)·dt`).

The test runs to γ = 10 (Nt = 11, dt = 0.5, γ_dot = 2) — sufficient to reach 95% of the asymptote (δ ≈ 3.29 at γ = 10 vs δ_∞ = 3.46). Running further provides no additional verification.

#### Scenario: AniFstrainCalcite matches the calibrated formula per step

- **WHEN** the test runs simple-shear with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 2`, `ani_fac_max = 100` (uncapped), `bkg_strain_rate = 1.0` so `γ_dot = 2`, `Nt = 11`, `dt = 0.5` (γ ∈ {0, 1, 2, ..., 10} at output steps 1..11), single phase
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and the analytical reference `δ_ana(γ) = 6.0 · 0.41 · (1 − exp(−γ_eff(FS_AR(γ)) / 3.71)) + 1` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved under calcite calibration

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4`

#### Scenario: AniFstrainCalcite reuses the shared base fixture and explicitly opts in to aniso_db = 2

- **WHEN** the test source under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** no new `.txt` fixture file SHALL be added for this test
- **AND** the test SHALL call `RunMDOODZ` on `AniFstrainEvolution.txt` and rely on a `MutateInput` callback that sets `ani_fstrain[0] = 2` AND `aniso_db[0] = 2` (the latter is required because the auto-default routes to olivine)

### Requirement: Research note `misc/aniso_fstrain/notes/calcite_calibration.md`

The change SHALL deliver a research-grade markdown record at `misc/aniso_fstrain/notes/calcite_calibration.md` documenting the calibration of the `aniso_db = 2` calcite form. The structure SHALL mirror §7 of [misc/aniso_fstrain/notes/hansen2012_comparison.md](../../misc/aniso_fstrain/notes/hansen2012_comparison.md) (the Hansen olivine research note) for consistency.

The note SHALL contain:

- (a) the functional form (Eq + constants `M_inf = 0.50`, `γ_e = 3.50`, `a = 6.0`),
- (b) the calibration procedure (least-squares fit on Pieri+01 + Bruijn+11 J-index data, J→M conversion via Skemer+05),
- (c) the **asymptote uncertainty**: `δ_∞ = a · M_inf + 1 = 4.0` is calibrated to match the empirical `ani_fac_max = 4` PinchSwell calibration; bracketed by Tommasi+09-style VPSC bounds (calcite δ_saturated ≈ 2-5). NOT a clean Hansen+12 Fig 3b-equivalent direct lab fit. This is the weakest link in the calibration.
- (d) the comparison plot at `misc/aniso_fstrain/img/calcite_calibration.png`, showing the new calibrated curve overlaid with Pieri+01 + Bruijn+11 J-index lab data (J→M-converted) and the existing min(FS_AR, 4) form for comparison,
- (e) the per-datapoint verification table (γ, J_obs, M_obs, M_fit, δ_fit) listing both Pieri+01 (γ ∈ {0, 1, 2, 5, 11}) and Bruijn+11 (γ ∈ {1, 2.6, 5}) datapoints,
- (f) a discussion of the **Pieri+01 vs Bruijn+11 disagreement at γ=5** (P087 J=9.76 vs Bruijn J≈4) and how the fit smooths through it,
- (g) a discussion of the **Pieri+01 γ=11 non-monotonic outlier** (J=8.92 < J=9.76 at γ=5; recrystallization randomization) and the structural decision to smooth through it,
- (h) the simple-shear `γ_eff(FS_AR)` caveat, identical to the olivine note,
- (i) explicit non-applicability to **HPT calcite (Schuster+18)** at confining pressures > 1.6 GPa where the calcite → CaCO₃-II phase transition kicks in,
- (j) full citations to Pieri+01, Bruijn+11, Schuster+18, Skemer+05, Barnhoorn+04, Tommasi+09 (or whichever VPSC-bound references underpin the slope choice).

The accompanying calibration script SHALL be checked in at `misc/aniso_fstrain/calibrate/calibrate_calcite.py` (Python, no scipy required — uses 2D grid search to fit `M_inf` and `γ_e` with `slope` fixed at 6.0).

#### Scenario: Calcite research note exists with required content

- **WHEN** `misc/aniso_fstrain/notes/calcite_calibration.md` is read
- **THEN** it SHALL contain sections covering each of (a)-(j) above
- **AND** it SHALL reference the calibration script `calibrate_calcite.py` and the comparison plot `mdoodz_vs_calcite_data.png` via relative paths
- **AND** it SHALL NOT be referenced from any `TESTS/` `.cpp` test file or CI script (research-grade, not regression gate)

#### Scenario: Calcite calibration script and plot are checked in

- **WHEN** `misc/aniso_fstrain/data/calcite/` is inspected
- **THEN** `calibrate_calcite.py` SHALL exist as an executable Python script that, when run, prints the fit constants `M_inf, γ_e`, the asymptote `δ_∞`, the verification table at all Pieri+01 + Bruijn+11 datapoints, and writes `mdoodz_vs_calcite_data.png` to the same directory
- **AND** the script SHALL NOT depend on scipy (numpy + matplotlib only)

### Requirement: Quartz Pennacchioni+10 δ-form is available as `aniso_db = 3`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a quartz-calibrated entry from the Pennacchioni+10 ODF-J dataset, dispatched when the per-phase parameter `aniso_db == 3` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( 3.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` (analytical inversion of the simple-shear identity, mineral-independent),
- `M(γ) = M_inf · (1 − exp(−γ / γ_e))` (saturating exponential),
- `M_inf = 0.75` and `γ_e = 3.50` (free LSQ fit to the 17-point Pennacchioni+10 natural-shear-zone dataset, J→M-converted via Skemer+05; LSQ centroid (0.736, 3.59) rounded to (0.75, 3.50) for clean numerics),
- `slope = 3.0` (quartz δ-vs-M slope; calibrated from VPSC-style Heilbronner-Tullis bounds for quartz aggregate viscous anisotropy at saturated prism-a CPO).

The closed-form δ-asymptote (before the `ani_fac_max` cap) is `δ_∞ = slope · M_inf + 1 = 3.25`, smaller than calcite's 4.0 and dramatically smaller than olivine's 14.13 — reflecting quartz's lower viscous-anisotropy magnitude even at saturated CPO.

**Calibrated for T ≥ 500°C, prism-a slip regime.** Below ~500°C basal-a slip becomes dominant and the fabric kinematics change; the calibration may not apply. For T < 500°C scenarios, falling back to `ani_fstrain = 1` with `ani_fac_max = 3` per phase is a reasonable approximation.

The form is **opt-in** via explicit `aniso_db = 3` per phase. The auto-default (`ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1`) remains **olivine-specific**. Existing `ani_fstrain ∈ {0, 1, 3}` and `aniso_db ∈ {1, 2}` dispatch arms remain byte-identical.

#### Scenario: aniso_db = 3 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the quartz function pointer (set by `ReadDataAnisotropy(case 3)`) with `FS_AR` set to the analytical simple-shear FS_AR at γ = 4.0 (`FS_AR ≈ 17.944`), `ani_fstrain = 2`, and `ani_fac_max = 100` (uncapped)
- **THEN** the returned δ SHALL be within `1e-6` of `3.0 · 0.75 · (1 − exp(−4.0 / 3.50)) + 1 ≈ 2.53`
- **AND** when called with `FS_AR` at γ = 8.0 (`FS_AR ≈ 66.0`), the returned δ SHALL be within `1e-6` of `3.0 · 0.75 · (1 − exp(−8.0 / 3.50)) + 1 ≈ 3.04`

#### Scenario: aniso_db = 3 saturates at δ_∞ ≈ 3.25 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the quartz function pointer with `FS_AR = 1e6`, `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `3.0 · M_inf + 1 = 3.25`

#### Scenario: aniso_db = 3 reduces to δ ≈ 1 at zero strain

- **WHEN** `AnisoFactorEvolv` is called via the quartz function pointer with `FS_AR = 1.0`, `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-9` of `1.0` — `M(0) = 0`, so δ = 3·0 + 1 = 1

#### Scenario: aniso_db ∈ {1, 2} dispatch arms are unaffected by the quartz addition

- **WHEN** `AnisoFactorEvolv` is called via the olivine function pointer (`aniso_db = 1`) OR the calcite function pointer (`aniso_db = 2`) for any value of `FS_AR ∈ [1, 1e6]`
- **THEN** the returned δ SHALL equal the value the previous (pre-quartz) implementation would have returned, to bitwise precision — adding case 3 SHALL NOT alter case 1 or case 2's IEEE-754 output

### Requirement: CI test `AniFstrainQuartz` exercises the quartz form

The CI test suite SHALL include a GoogleTest case `AnisotropyBenchmark.AniFstrainQuartz` that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 3` and verifies that the per-cell `aniso_factor_n` matches the analytical reference (the calibrated quartz formula evaluated symbolically) to FP precision over the simulated strain range.

The test SHALL share the existing base fixture `AniFstrainEvolution.txt` via `MutateInput`, consistent with the existing tests. No new fixture file SHALL be added.

The test runs to γ = 8 (Nt = 9, dt = 0.5, γ_dot = 2) — sufficient to reach 99% of the asymptote (δ ≈ 3.06 at γ = 8 vs δ_∞ = 3.1). Quartz saturates fast enough that γ = 8 is plenty.

#### Scenario: AniFstrainQuartz matches the calibrated formula per step

- **WHEN** the test runs simple-shear with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 3`, `ani_fac_max = 100` (uncapped), `bkg_strain_rate = 1.0` so `γ_dot = 2`, `Nt = 9`, `dt = 0.5` (γ ∈ {0, 1, 2, ..., 8} at output steps 1..9), single phase
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` and the analytical reference `δ_ana(γ) = 3.0 · 0.75 · (1 − exp(−γ_eff(FS_AR(γ)) / 3.50)) + 1` SHALL be less than `1e-6`
- **AND** the relative L2 error SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved under quartz calibration

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4`

#### Scenario: AniFstrainQuartz reuses the shared base fixture and explicitly opts in to aniso_db = 3

- **WHEN** the test source under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** no new `.txt` fixture file SHALL be added for this test
- **AND** the test SHALL call `RunMDOODZ` on `AniFstrainEvolution.txt` and rely on a `MutateInput` callback that sets `ani_fstrain[0] = 2` AND `aniso_db[0] = 3` (the latter is required because the auto-default routes to olivine)

### Requirement: Research note `misc/aniso_fstrain/notes/quartz_calibration.md`

The change SHALL deliver a research-grade markdown record at `misc/aniso_fstrain/notes/quartz_calibration.md` documenting the calibration of the `aniso_db = 3` quartz form. The structure SHALL mirror §7 of [misc/aniso_fstrain/notes/hansen2012_comparison.md](../../misc/aniso_fstrain/notes/hansen2012_comparison.md) and [misc/aniso_fstrain/notes/calcite_calibration.md](../../misc/aniso_fstrain/notes/calcite_calibration.md) for consistency.

The note SHALL contain:

- (a) the functional form (Eq + constants for case 3 `(M_inf = 0.75, γ_e = 3.50, slope = 3.0)` and case 6 `(M_inf = 0.20, γ_e = 4.00, slope = 11.0)`),
- (b) the calibration procedure (free LSQ fit on Pennacchioni+10 ODF J 17 points and Blackford+24 pfJ 38 points, J→M conversion via Skemer+05),
- (c) the **asymptote uncertainty**: `δ_∞ ≈ 3.2` is calibrated from VPSC-style Heilbronner-Tullis bounds (quartz δ_saturated ≈ 2-3); NOT a clean direct lab fit. Documented as the weakest link.
- (d) the comparison plot at `misc/aniso_fstrain/img/quartz_calibration.png`, showing both calibrated curves overlaid with Pennacchioni+10 + Blackford+24 lab data and the previous min-form alternatives.
- (e) the per-datapoint verification table covering all 17 Pennacchioni+10 datapoints with regime labels (WDV/MDV/SDV).
- (f) the **explicit T-range caveat**: ~500°C, prism-a slip dominant.
- (g) the **strain-measure conversion** for Blackford+24 (γ_oct → γ_simple_eq via 2·sinh(ε/√2) plane-strain inversion).
- (h) discussion of the **ODF J vs pfJ metric mismatch** that motivates two separate modes.
- (i) full citation to Pennacchioni+10, Blackford+24, Heilbronner & Tullis 2006 (VPSC bounds), and Skemer+05.

The accompanying calibration script SHALL be at `misc/aniso_fstrain/calibrate/calibrate_quartz.py` and the HDF5-overlay script at `misc/aniso_fstrain/mdoodz_compare/mdoodz_vs_quartz.py`. Both pure stdlib + numpy + matplotlib, no scipy.

#### Scenario: Quartz research note exists with required content

- **WHEN** `misc/aniso_fstrain/notes/quartz_calibration.md` is read
- **THEN** it SHALL contain sections covering each of (a)-(i) above
- **AND** it SHALL document both calibrated forms (case 3: M_inf=0.75, γ_e=3.50, slope=3.0, δ_∞=3.25; case 6: M_inf=0.20, γ_e=4.00, slope=11.0, δ_∞=3.20)
- **AND** it SHALL NOT be referenced from any `TESTS/` `.cpp` test file or CI script (research-grade, not regression gate)

#### Scenario: Quartz calibration scripts and plots are checked in

- **WHEN** `misc/aniso_fstrain/data/quartz/` is inspected
- **THEN** `calibrate_quartz.py` SHALL exist as an executable Python script that, when run, loads BOTH the Pennacchioni+10 17-point and Blackford+24 38-point datasets, runs free LSQ fits per dataset and combined, prints the verification tables, prints RMS errors of all candidate δ-forms (min-form alternatives + cases 3 and 6) against each dataset, and writes `mdoodz_vs_quartz_data.png` (3-panel plot showing both calibrated curves) to the same directory
- **AND** `mdoodz_vs_quartz_comparison.py` SHALL exist for HDF5-overlay rendering after CI tests run, handling both `AniFstrainQuartz` and `AniFstrainQuartzBlackford` build directories
- **AND** the scripts SHALL NOT depend on scipy (numpy + matplotlib only)

### Requirement: Calcite mid-low-T regime δ-form is available as `aniso_db = 4`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a calcite mid-low-T entry, dispatched when `aniso_db == 4`:

```
δ(FS_AR) = min( 6.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with `M_inf = 0.22`, `γ_e = 1.37`, `slope = 6.0`. Asymptote `δ_∞ = 2.33`.

Calibrated by free LSQ fit to the 7-point Barnhoorn+04 dataset at 500°C (4 pts) + 600°C (3 pts), J→M-converted via Skemer+05. Fit RMS = 0.073 — much tighter than the T-averaged calibration because subgrain-rotation recrystallization dominates this regime and produces a consistent peak-then-stable behaviour the saturating-exponential form captures well.

**Recommended for T < 650°C calcite scenarios.** SR-recryst randomizes CPO at high γ, capping fabric at modest strength.

#### Scenario: aniso_db = 4 returns the calibrated formula at γ = 4 and γ = 6

- **WHEN** `AnisoFactorEvolv` is called via the calcite-low-T function pointer with `FS_AR ≈ 17.944` (γ=4), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `6.0 · 0.22 · (1 − exp(−4.0 / 1.37)) + 1 ≈ 2.26`
- **AND** at `FS_AR ≈ 39.4` (γ=6), the returned δ SHALL be within `1e-6` of `6.0 · 0.22 · (1 − exp(−6.0 / 1.37)) + 1 ≈ 2.31`

#### Scenario: aniso_db = 4 saturates at δ_∞ ≈ 2.33

- **WHEN** `AnisoFactorEvolv` is called via the calcite-low-T function pointer with `FS_AR = 1e6`
- **THEN** the returned δ SHALL be within `1e-3` of `2.33`

### Requirement: Calcite high-T regime δ-form is available as `aniso_db = 5`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a calcite high-T entry, dispatched when `aniso_db == 5`:

```
δ(FS_AR) = min( 6.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with `M_inf = 0.83`, `γ_e = 8.24`, `slope = 6.0`. Asymptote `δ_∞ = 5.99`.

Calibrated by free LSQ fit to the 11-point combined Pieri+01 (5 pts) + Bruijn+11 (3 pts) + Barnhoorn+04(727°C, 3 pts after J ≤ 15 filter) dataset, J→M-converted via Skemer+05. Fit RMS = 0.169 — looser than mid-low-T because the published lab data has internal disagreement (Pieri+01 saturates at γ ≈ 5 while Barnhoorn+04 727°C continues strengthening past γ = 46).

**Recommended for T > 650°C calcite scenarios.** Grain-boundary-migration recrystallization dominates this regime and preferentially eliminates unfavorable orientations, strengthening CPO continuously even at very high strain.

The large γ_e = 8.24 reflects this slow, continuous strengthening; for very-high-γ scenarios (γ > 20) the calibration extrapolates beyond the linear-conversion-valid Barnhoorn data.

#### Scenario: aniso_db = 5 returns the calibrated formula at γ = 5 and γ = 14

- **WHEN** `AnisoFactorEvolv` is called via the calcite-high-T function pointer with `FS_AR ≈ 27.3` (γ=5), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `6.0 · 0.83 · (1 − exp(−5.0 / 8.24)) + 1 ≈ 3.27`
- **AND** at `γ=14` (`FS_AR ≈ 196`), the returned δ SHALL be within `1e-6` of `6.0 · 0.83 · (1 − exp(−14.0 / 8.24)) + 1 ≈ 5.07`

#### Scenario: aniso_db = 5 saturates at δ_∞ ≈ 5.99

- **WHEN** `AnisoFactorEvolv` is called via the calcite-high-T function pointer with `FS_AR = 1e6`
- **THEN** the returned δ SHALL be within `1e-3` of `5.99`

### Requirement: CI tests `AniFstrainCalciteLowT` and `AniFstrainCalciteHighT` exercise the 2-regime calcite forms

The CI test suite SHALL include `AnisotropyBenchmark.AniFstrainCalciteLowT` (γ ∈ [0, 6], `aniso_db = 4`) and `AnisotropyBenchmark.AniFstrainCalciteHighT` (γ ∈ [0, 14], `aniso_db = 5`), each asserting per-step δ matches the analytical reference to `1e-6` relative tolerance.

#### Scenario: AniFstrainCalciteLowT matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 7`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 4`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold

#### Scenario: AniFstrainCalciteHighT matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 15`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 5`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold

### Requirement: Quartz Blackford+24 δ-form is available as `aniso_db = 6`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a second quartz-calibrated entry from the Blackford+24 pole-figure-J dataset, dispatched when `aniso_db == 6` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( 11.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `M(γ) = 0.20 · (1 − exp(−γ / 4.00))` (Blackford+24-derived saturating exponential),
- `slope = 11.0` (rescaled from case 3's 3.0 to compensate for pfJ < ODF J on the same physical texture: pfJ-based M is ~5× smaller, so slope must be ~5× larger to give the same physical δ_∞).

Asymptote `δ_∞ = 11.0 · 0.20 + 1 = 3.20` matches case 3's δ_∞ = 3.25 to within rounding. The two modes describe the same physical mineral measured on different texture metrics; case 6 is appropriate when calibrating against modern EBSD-derived pole-figure J values, case 3 against texture-goniometry ODF J values.

Calibrated by free LSQ fit to the 38-point Blackford+24 dataset (Tectonics 43, e2023TC008166). Strain reported in the paper as Nadai natural strain ε; converted to simple-shear-equivalent γ via the plane-strain inversion `γ_eq = 2·sinh(ε/√2)` (samples described as plane-strain to slightly constrictional).

The form is **opt-in** via explicit `aniso_db = 6` per phase. The auto-default routes to olivine, so quartz-Blackford users MUST set `aniso_db = 6` explicitly. Cases 1–5 dispatch arms remain byte-identical.

#### Scenario: aniso_db = 6 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the quartz-Blackford function pointer (set by `ReadDataAnisotropy(case 6)`) with `FS_AR ≈ 17.944` (γ = 4), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `11.0 · 0.20 · (1 − exp(−4.0 / 4.00)) + 1 ≈ 2.39`
- **AND** at `FS_AR ≈ 66.0` (γ = 8), the returned δ SHALL be within `1e-6` of `11.0 · 0.20 · (1 − exp(−8.0 / 4.00)) + 1 ≈ 2.90`

#### Scenario: aniso_db = 6 saturates at δ_∞ ≈ 3.20 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the quartz-Blackford function pointer with `FS_AR = 1e6`, `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `11.0 · 0.20 + 1 = 3.20`

### Requirement: CI test `AniFstrainQuartzBlackford` exercises the Blackford form

The CI test suite SHALL include `AnisotropyBenchmark.AniFstrainQuartzBlackford` (γ ∈ [0, 8], `aniso_db = 6`), asserting per-step δ matches the analytical reference to `1e-6` relative tolerance.

#### Scenario: AniFstrainQuartzBlackford matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 9`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 6`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold for `δ_ana(γ) = 11.0 · 0.20 · (1 − exp(−γ_eff(FS_AR(γ)) / 4.00)) + 1`

### Requirement: Damped olivine δ-form is available as `aniso_db = 7`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a damped-olivine entry, dispatched when `aniso_db == 7` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( 24.5 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `M(γ) = 0.20 · (1 − exp(−γ / 1.30))` (free LSQ fit on combined non-Hansen 93-sample dataset),
- `slope = 24.5` — same as case 1, because Hansen+12 Fig 3b's M↔δ linear relation is universal across olivine fabric strength.

Asymptote `δ_∞ = 24.5 · 0.20 + 1 = 5.90` — about half of case 1's 14.13, reflecting damped fabric in natural / wet / pre-CPO olivine.

The 93-sample non-Hansen calibration combines four independent datasets:
- Tasaka 2016 (JGR-SE 121:92): 12 wet Fo50 lab torsion samples
- Boneh & Skemer 2014 (EPSL 406:213): 16 Aheim dunite samples with pre-existing CPO
- Kumamoto+19 (JGR-SE 124:12763): 30 natural Josephine Peridotite shear-zone samples
- Bernard+19 (G3 20:3469): 35 natural xenoliths with X/Z aspect-ratio strain proxies (Bernard's "octahedral" strain mapped to γ_eq via plane-strain inversion)

Each population's strain is converted to a simple-shear-equivalent γ before the LSQ fit; conversions and rationale are documented in `misc/aniso_fstrain/olivine_extended/olivine_calibration.md`.

The form is **opt-in** via explicit `aniso_db = 7` per phase. The auto-default `ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1` remains olivine-fresh-CPO; users wanting the damped mode (mantle xenoliths, natural shear-zone mylonites, wet/Fe-rich olivine) MUST set `aniso_db = 7` explicitly. Cases 1–6 dispatch arms remain byte-identical.

#### Scenario: aniso_db = 7 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the damped-olivine function pointer (set by `ReadDataAnisotropy(case 7)`) with `FS_AR` set to the analytical simple-shear FS_AR at γ = 2.0 (`FS_AR ≈ 5.83`), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `24.5 · 0.20 · (1 − exp(−2.0 / 1.30)) + 1 ≈ 4.85`
- **AND** at `FS_AR ≈ 17.944` (γ = 4), the returned δ SHALL be within `1e-6` of `24.5 · 0.20 · (1 − exp(−4.0 / 1.30)) + 1 ≈ 5.67`

#### Scenario: aniso_db = 7 saturates at δ_∞ ≈ 5.90 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the damped-olivine function pointer with `FS_AR = 1e6`, `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `24.5 · 0.20 + 1 = 5.90`

### Requirement: CI test `AniFstrainHansenOlivineDamped` exercises the damped olivine form

The CI test suite SHALL include `AnisotropyBenchmark.AniFstrainHansenOlivineDamped` (γ ∈ [0, 8], `aniso_db = 7`), asserting per-step δ matches the analytical reference to `1e-6` relative tolerance.

#### Scenario: AniFstrainHansenOlivineDamped matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 9`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 7`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold for `δ_ana(γ) = 24.5 · 0.20 · (1 − exp(−γ_eff(FS_AR(γ)) / 1.30)) + 1`

### Requirement: Olivine calibration note documents both case 1 and case 7

The change SHALL deliver a research-grade markdown record at `misc/aniso_fstrain/olivine_extended/olivine_calibration.md` documenting both calibrated olivine forms, the strain-measure conversions for the four Population B datasets, and the rationale for using a single damped mode rather than splitting into 3 (wet / pre-CPO / natural).

#### Scenario: Olivine calibration note exists

- **WHEN** `misc/aniso_fstrain/olivine_extended/olivine_calibration.md` is read
- **THEN** it SHALL document both calibrated forms (case 1: M_inf=0.536, γ_e=3.96, δ_∞=14.13; case 7: M_inf=0.20, γ_e=1.30, δ_∞=5.90)
- **AND** it SHALL document strain-measure conversions for Boneh axial strain (γ_eq = 2·sinh(ε/2)) and Bernard X/Z aspect ratio (γ_eq = √(X/Z) − 1/√(X/Z))
- **AND** it SHALL include the comparison plot `olivine_all_datasets.png`
- **AND** it SHALL contain the per-scenario recommendation table for case 1 vs case 7

### Requirement: Quartz-mica polyphase BRACKETED δ-form is available as `aniso_db = 8`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a quartz-mica polyphase entry, dispatched when `aniso_db == 8` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( 15.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with `M(γ) = 0.40 · (1 − exp(−γ / 2.50))`, `slope = 15.0`, asymptote `δ_∞ = 7.0`.

**This is a BRACKETED estimate, not a direct LSQ fit.** Constants are bracketed from disparate literature sources (Tokle+23 saturation γ; Dempsey+11 mica J at saturation; Holyoke & Tullis 2006 interconnection threshold; Bos & Spiers 2002 micromechanical mixing) because no co-located (γ, J_bulk) lab dataset exists for quartz-mica aggregates. Quality tier matches `pwlv = 16` Ranalli felsic granulite in the power-law database — a defensible point estimate from a literature compilation, not a direct lab fit.

Suitable for **typical 10–25% mica content** crustal mylonites (granitoid mylonite, quartzite mylonite with muscovite films, two-mica granite, mica schist). Outside this range users should fall back to `ani_fstrain = 1` with chosen `ani_fac_max`.

The form is **opt-in** via explicit `aniso_db = 8`. Cases 1–7 remain byte-identical.

#### Scenario: aniso_db = 8 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the quartz-mica function pointer with `FS_AR ≈ 5.83` (γ = 2), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `15.0 · 0.40 · (1 − exp(−2.0 / 2.50)) + 1 ≈ 4.31`
- **AND** at `FS_AR ≈ 17.944` (γ = 4), within `1e-6` of `15.0 · 0.40 · (1 − exp(−4.0 / 2.50)) + 1 ≈ 5.79`

#### Scenario: aniso_db = 8 saturates at δ_∞ ≈ 7.0 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the quartz-mica function pointer with `FS_AR = 1e6`
- **THEN** the returned δ SHALL be within `1e-3` of `15.0 · 0.40 + 1 = 7.0`

### Requirement: Plagioclase polyphase BRACKETED δ-form is available as `aniso_db = 9`

`ReadDataAnisotropy()` SHALL provide a plagioclase polyphase entry, dispatched when `aniso_db == 9`:

```
δ(FS_AR) = min( 8.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with `M(γ) = 0.12 · (1 − exp(−γ / 2.00))`, `slope = 8.0`, asymptote `δ_∞ = 1.96 ≈ 2.0`.

**This is a BRACKETED estimate, not a direct LSQ fit.** Constants are bracketed from Mehl & Hirth 2008 (saturated monophase plag M ≈ 0.12), Marti et al. 2018 (qualitative weak CPO via DPC + GBS), Gómez Barreiro et al. 2015 (Morin anorthosite max 2.9 m.r.d.), and Stünitz et al. 2003 (slip systems). The saturating-exponential form **over-predicts at γ > 8** because plagioclase CPO can WEAKEN at high strain due to DisGBS / diffusion creep (Mehl & Hirth polyphase plag layers M ≈ 0.04). Quality tier matches `pwlv = 17` Ranalli mafic granulite.

Suitable for **monophase plagioclase mylonite at saturation** (gabbro mylonite, anorthosite, mafic granulite, basalt, diabase). For polyphase plag-cpx layers at γ > 8, users should fall back to `ani_fstrain = 1, ani_fac_max ≈ 1.2–1.5`.

The form is **opt-in** via explicit `aniso_db = 9`. Cases 1–8 remain byte-identical.

#### Scenario: aniso_db = 9 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the plagioclase function pointer with `FS_AR ≈ 5.83` (γ = 2), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `8.0 · 0.12 · (1 − exp(−2.0 / 2.00)) + 1 ≈ 1.61`
- **AND** at `FS_AR ≈ 66.0` (γ = 8), within `1e-6` of `8.0 · 0.12 · (1 − exp(−8.0 / 2.00)) + 1 ≈ 1.94`

#### Scenario: aniso_db = 9 saturates at δ_∞ ≈ 1.96 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the plagioclase function pointer with `FS_AR = 1e6`
- **THEN** the returned δ SHALL be within `1e-3` of `8.0 · 0.12 + 1 = 1.96`

### Requirement: BRACKETED modes are flagged in LOG_INFO output

The LOG_INFO message emitted by `ReadDataAnisotropy` for cases 8 and 9 SHALL contain the literal string `"BRACKETED"` to flag the looser provenance. The LOG_INFO messages for cases 1–7 SHALL NOT contain that string.

#### Scenario: Case 8 LOG_INFO contains "BRACKETED"

- **WHEN** `ReadDataAnisotropy` is called with `aniso_db = 8`
- **THEN** the emitted log message SHALL contain the substring `"BRACKETED"`
- **AND** SHALL contain `"comparable rigor to pwlv=16 Ranalli granulite"` to flag the parallel with the power-law database convention

#### Scenario: Case 9 LOG_INFO contains "BRACKETED"

- **WHEN** `ReadDataAnisotropy` is called with `aniso_db = 9`
- **THEN** the emitted log message SHALL contain the substring `"BRACKETED"`
- **AND** SHALL contain `"comparable rigor to pwlv=17 Ranalli mafic granulite"`

### Requirement: CI tests `AniFstrainQuartzMicaBracketed` and `AniFstrainPlagioclaseBracketed`

The CI test suite SHALL include `AnisotropyBenchmark.AniFstrainQuartzMicaBracketed` (γ ∈ [0, 8], `aniso_db = 8`) and `AnisotropyBenchmark.AniFstrainPlagioclaseBracketed` (γ ∈ [0, 8], `aniso_db = 9`), each asserting per-step δ matches the analytical reference to `1e-6` relative tolerance.

The tests gate the IMPLEMENTATION (MDOODZ matches the formula) — they do NOT validate the science. Users are responsible for choosing whether bracketed-grade rigor is appropriate for their scenario.

#### Scenario: AniFstrainQuartzMicaBracketed matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 9`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 8`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold for `δ_ana(γ) = 15.0 · 0.40 · (1 − exp(−γ_eff(FS_AR(γ)) / 2.50)) + 1`

#### Scenario: AniFstrainPlagioclaseBracketed matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 9`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 9`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold for `δ_ana(γ) = 8.0 · 0.12 · (1 − exp(−γ_eff(FS_AR(γ)) / 2.00)) + 1`

### Requirement: `quartz_calibration.md` documents BRACKETED modes alongside direct-fit cases

The research note `misc/aniso_fstrain/notes/quartz_calibration.md` SHALL document both bracketed cases (8 and 9) in dedicated sections (§11 quartz-mica, §12 plagioclase) with explicit "BRACKETED" labelling, comparison to the `pwlv = 16` / `pwlv = 17` Ranalli convention, mode selection table by rock type, and the three structural reasons each case is bracketed rather than direct-fit.

#### Scenario: Bracketed sections are clearly distinguishable from direct-fit sections

- **WHEN** `quartz_calibration.md` §11 and §12 are read
- **THEN** they SHALL contain the literal word "BRACKETED" in the section heading or first paragraph
- **AND** SHALL contrast quality tier explicitly with cases 1–7 ("comparable rigor to pwlv=16/17", "not Hansen+14-grade")
- **AND** SHALL NOT report an RMS-of-fit metric (since no regression was performed)

### Requirement: Grain-size-coupled δ modifier (d-limiter)

`AnisoFactorEvolv` SHALL apply a grain-size-coupled modifier to the calibrated δ value when two new per-phase parameters are both set: `aniso_d_threshold > 0` and the local grain size `d > 0`. The modifier captures CPO randomization at small grain sizes (DisGBS / GBM regime) — physical observation, e.g., Mehl & Hirth (2008) for plagioclase, Précigout & Hirth (2014) for olivine.

**Functional form** (applied AFTER `aniso_delta_fn(FS_AR)` returns the base δ, BEFORE the `aniso_fac_max` cap):

```
if d_threshold > 0 AND d > 0 AND d < d_threshold:
    f      = (d / d_threshold) ^ d_decay      ∈ (0, 1)
    δ_eff  = 1 + f · (δ − 1)
else:
    δ_eff  = δ unchanged
```

with `d_decay = 1.0` giving a linear blend. `d_decay > 1` makes the modifier kick in only near `d = 0`; `d_decay < 1` gives a gentler decay starting earlier.

**Limit cases that the formula SHALL satisfy** (testable algebraically):

- `d → 0`: `f → 0`, `δ_eff → 1.0` (complete CPO randomization, isotropic limit)
- `d → d_threshold` from below: `f → 1`, `δ_eff → δ` (continuous transition)
- `d ≥ d_threshold`: modifier inactive, `δ_eff = δ` exactly
- `d_threshold = 0` (default): modifier inactive regardless of `d`, `δ_eff = δ` exactly
- `d_decay = 0`: `f = 1` for any `d > 0` below threshold (modifier gives full base δ; degenerate but well-defined)

**Backward compatibility requirement**: when `aniso_d_threshold[k] = 0` for all phases `k` (the default for cases 1–8), `AnisoFactorEvolv` SHALL produce δ values byte-identical to the pre-d-limiter implementation. This is testable: the existing 18 anisotropy CI tests SHALL continue to pass to FP precision (1e-6 relative tolerance) without any test fixture changes.

**Per-phase storage**: two new fields are added to `mat_prop` in [`MDLIB/include/mdoodz.h`](../../MDLIB/include/mdoodz.h):

```c
double aniso_d_threshold[20];   // [m]; 0 = inactive
double aniso_d_decay[20];       // dimensionless; default 1.0
```

**Parser**: [`MDLIB/InputOutput.c`](../../MDLIB/InputOutput.c) `ReadParameters` parses both fields per phase via `ReadMatProps` with defaults `(threshold=0.0, decay=1.0)`. The threshold is divided by `scaling.L` so it lives in scaled units alongside `gs_ref` and `mesh->d_n`.

**Per-case auto-default**: in `ReadDataAnisotropy` (FlowLaws.c), `case 9` SHALL set `aniso_d_threshold = 5e-6 m` and `aniso_d_decay = 1.0` automatically IF the user has not overridden them in the input file. This anchors plagioclase BRACKETED to Mehl & Hirth (2008) DisGBS observations. Cases 1–8 SHALL leave both fields at their parsed defaults (typically 0.0). The auto-default SHALL be detectable: the LOG_INFO message for case 9 SHALL announce the d-coupling defaults when applied, and SHALL NOT announce them when the user has set non-zero values explicitly.

**Function signature**: `AnisoFactorEvolv` SHALL take three additional arguments compared to its pre-d-limiter signature:

```c
double AnisoFactorEvolv( double FS_AR, double aniso_fac_max, int ani_fstrain,
                         double (*aniso_delta_fn)(double FS_AR),
                         double grain_size,           // NEW
                         double aniso_d_threshold,    // NEW
                         double aniso_d_decay );      // NEW
```

The function pointer signature `(*aniso_delta_fn)(double FS_AR)` SHALL remain unchanged — d-coupling is a uniform post-processing step in `AnisoFactorEvolv`, not embedded in the per-mineral helpers. This keeps cases 1–9 byte-identical at the helper level and concentrates the d-coupling logic in one place.

**Vertex (shear-node) grain size**: the existing `mesh->d_n` field is per-cell-center only; `mesh->d` does not have a vertex-grid version. For the 6 call sites in `AnisotropyRoutines.c::UpdateAnisoFactor`, the 3 center calls use `mesh->d_n[c0]` directly, and the 3 vertex calls compute `d_at_vertex` as a 2×2 average of the surrounding cell-center values:

```c
for ( int dj = -1; dj <= 0; dj++ ) for ( int di = -1; di <= 0; di++ ) {
    int ii = k + di, jj = l + dj;
    if ( ii >= 0 && ii < Ncx && jj >= 0 && jj < Ncz
         && mesh->d_n[ii + jj*Ncx] > 0.0 ) {
        d_at_vertex += mesh->d_n[ii + jj*Ncx];
        count++;
    }
}
if ( count > 0 ) d_at_vertex /= count;
```

This averaging is sufficient for the d-modifier (which is a smooth power-law decay, not a sharp threshold), but the implementation SHOULD note that vertex-level d resolution is one half-cell coarser than center-level.

#### Scenario: d-modifier inactive when threshold is zero (default for cases 1–8)

- **WHEN** `AnisoFactorEvolv` is called with `aniso_d_threshold = 0.0` for any combination of `(FS_AR, ani_fstrain, aniso_delta_fn, grain_size, aniso_d_decay)`
- **THEN** the returned δ SHALL equal the pre-d-limiter result EXACTLY (byte-identity to FP precision)
- **AND** the existing 18 anisotropy CI tests (cases 1–9 across olivine, calcite, quartz, plagioclase) SHALL continue to pass without any test fixture modifications

#### Scenario: d-modifier inactive when grain size exceeds threshold

- **WHEN** `AnisoFactorEvolv` is called with `aniso_d_threshold = 5e-6`, `grain_size = 1e-5` (10 µm > 5 µm threshold), `aniso_d_decay = 1.0`, `ani_fstrain = 2`, and any base δ from `aniso_delta_fn(FS_AR)`
- **THEN** the returned δ SHALL equal the base δ unchanged (capped by `aniso_fac_max`)

#### Scenario: d-modifier reduces δ continuously below threshold

- **WHEN** `AnisoFactorEvolv` is called with `aniso_d_threshold = 5e-6`, `aniso_d_decay = 1.0`, base δ = 1.96 (case 9 asymptote), `aniso_fac_max = 100`, and grain_size values 0, 1e-6, 2.5e-6, 5e-6
- **THEN** the returned δ values SHALL be 1.0, 1.192, 1.480, 1.96 respectively (within 1e-6 numerical tolerance)
- **AND** the function SHALL be monotonic non-decreasing in `grain_size` over `[0, ∞)`

#### Scenario: case 9 sets d-coupling defaults automatically

- **WHEN** `ReadDataAnisotropy(case 9)` is called via a phase where the user has not set `aniso_d_threshold` (parsed value is 0.0)
- **THEN** `mat->aniso_d_threshold[k]` SHALL be set to `5.0e-6 / scaling->L` (5 µm in scaled units)
- **AND** `mat->aniso_d_decay[k]` SHALL be set to 1.0
- **AND** the LOG_INFO output SHALL contain "case 9 default d-coupling" or equivalent identifying string

#### Scenario: case 9 respects user override

- **WHEN** the input file specifies a non-zero `aniso_d_threshold` for a phase using `aniso_db = 9`
- **THEN** `ReadDataAnisotropy(case 9)` SHALL leave the user-set value unchanged
- **AND** the LOG_INFO SHALL NOT announce the default (because the default branch is bypassed)

### Requirement: Benchmark reuses existing test-helper utilities

The benchmark SHALL use the existing `computeL2Error` helper from [TESTS/TestHelpers.h](../../TESTS/TestHelpers.h) and the existing `getMin/MaxFieldValue` HDF5-reader helpers, with no new analytical-error or HDF5 utility code added at the test-helper layer. Test-local state structures (e.g. `AniFstrainOverride`) MAY live in file scope inside `AnisotropyBenchmarkTests.cpp` or be promoted to `TestHelpers.h` if and only if the same shape is shared with another test file (in this change: kept file-scope).

#### Scenario: No new test-helper utilities added

- **WHEN** the benchmark is added to [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp)
- **THEN** the implementation SHALL call only existing helpers from [TESTS/TestHelpers.h](../../TESTS/TestHelpers.h) and SHALL NOT introduce new files under `TESTS/` named `*Helpers.h`, `*Helpers.cpp`, or similar

### Requirement: `ani_fstrain = 3` enables temperature-dependent δ-relaxation

The `ani_fstrain` enumeration SHALL define value `3` as: the `ani_fstrain == 2` δ-dispatch (mineral-calibrated `aniso_delta_fn(FS_AR)` selected by `aniso_db`) **plus** a temperature-dependent kinetic relaxation of the anisotropy factor δ toward the isotropic limit. The relaxation physics, state, integration scheme, inputs, and CI coverage are specified by the `aniso-delta-relaxation` capability. The previous meaning of `ani_fstrain == 3` (an upstream `T < ani_T_threshold` deformation-gradient freeze) and the associated `ani_T_threshold` per-phase parameter SHALL no longer exist.

#### Scenario: ani_fstrain enumeration documents value 3

- **WHEN** the `ani_fstrain` enumeration is documented (code comments and `TESTS/AnalyticalSolutions.md`)
- **THEN** value `3` SHALL be described as "`ani_fstrain == 2` dispatch + temperature-dependent δ-relaxation"
- **AND** `ani_T_threshold` SHALL NOT appear as a per-phase parameter

#### Scenario: ani_fstrain values 0, 1, 2 are unaffected

- **WHEN** this change is applied
- **THEN** the dispatch arms and CI tests for `ani_fstrain` in {0, 1, 2} SHALL be byte-identical to their prior behaviour

#### Scenario: T-threshold-freeze tests are removed

- **WHEN** the anisotropy CI suite is built after this change
- **THEN** the `AniFstrainT_Threshold_Freeze_Sentinel` and `AniFstrainT_Threshold_Freeze_Active` GoogleTest cases SHALL no longer exist
- **AND** the suite SHALL instead contain the `aniso-delta-relaxation` analytical-unit-test cases

