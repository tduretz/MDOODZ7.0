## ADDED Requirements

### Requirement: Analytical L2 benchmark for finite-strain anisotropy under pure shear

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainPureShear`) that runs a homogeneous prescribed pure-shear simulation with `anisotropy = 1` and per-phase `ani_fstrain = 1`, and asserts that MDOODZ's measured per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the closed-form analytical solution `FS_AR(t) = exp(2·ε̇·t)` while unsaturated. The reference is derived from the polar decomposition of the deformation gradient `F(t) = diag(exp(-ε̇·t), exp(ε̇·t))` (with prescribed pure-shear velocity field `Vx = -x·ε̇`, `Vz = z·ε̇`), giving principal stretches `(e1, e2) = (exp(ε̇·t), exp(-ε̇·t))` and `FS_AR = e1/e2 = exp(2·ε̇·t)`.

This requirement is added once the step-2 segfault under `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` is fixed in MDLIB. Before the fix, this test case was deferred (see archived [`gate-finite-strain-anisotropy/tasks.md` §5.1](../../changes/archive/2026-05-06-gate-finite-strain-anisotropy/tasks.md)).

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

## MODIFIED Requirements

### Requirement: Analytical L2 benchmark for finite-strain anisotropy under simple shear

The CI test suite SHALL include a GoogleTest case (`AnisotropyBenchmark.AniFstrainSimpleShear`) that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1` and per-phase `ani_fstrain = 1`, and asserts that MDOODZ's measured per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the closed-form analytical solution `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4 + 1)` while unsaturated. The reference is derived from the polar decomposition of the deformation gradient `F(t) = [[1, γ(t)], [0, 1]]` (with `γ = γ_dot · t`), giving `C = F^T·F = [[1, γ], [γ, 1+γ²]]`, eigenvalues `(1 + γ²/2) ± γ·sqrt(γ²/4 + 1)`, and `FS_AR = e1/e2 = λ_max` (since `λ_min · λ_max = 1`).

**Note on output indexing**: MDOODZ writes `Output<k>.gzip.h5` at the END of step `k`, but the `mesh.aniso_factor_n` field in that file reflects `FS_AR` computed at the START of step `k` (from the deformation gradient at the END of step `k-1`). So the analytical reference for `Output<k>` is evaluated at `t = (k-1)·dt`.

The pure-shear analytical reference (`FS_AR(t) = exp(2·ε̇·t)`) is also exercised in CI by the sibling requirement *Analytical L2 benchmark for finite-strain anisotropy under pure shear*. Both deformation modes share the same `FiniteStrainAspectRatio` → `mesh.FS_AR_n` → `AnisoFactorEvolv` → `mesh.aniso_factor_n` code path; running both gives independent confirmation of the path under the two BC types most commonly used in research scenarios.

#### Scenario: Simple-shear δ(γ) matches the closed-form FS_AR expression to FP precision

- **WHEN** a multi-step homogeneous simple-shear simulation runs with `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max` set high enough that no step saturates (e.g. `ani_fac_max = 1e6`), `bkg_strain_rate = 1.0` so that `γ_dot = 2·bkg_strain_rate = 2`, `Nt = 5`, `dt = 0.5`, single phase, and `Centers/ani_fac` is read at each output step
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and `1 + γ²/2 + γ·sqrt(γ²/4 + 1)` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector at the same value SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved through P2G interpolation

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4` — homogeneous setup must remain homogeneous after particle-to-grid interpolation

### Requirement: Single base fixture file with MutateInput per-test overrides

The benchmark SHALL use exactly one shared base fixture file [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](../../TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) (with `Nb_phases = 1`, `anisotropy = 1`, `ani_fstrain = 1`, sensible defaults for the rest), and SHALL NOT introduce per-test fixture files. Each test case SHALL inject its differentiating parameters (`bkg_strain_rate`, `Nt`, `dt`, `dt0`, `dt_start`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, `writer_subfolder`) via MDLIB's `MutateInput` callback hook, matching the established convention in [TESTS/RheologyCreepTests.cpp](../../TESTS/RheologyCreepTests.cpp) (`mutateSweepInput`) and [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp) (`mutateDirectorInput`).

#### Scenario: Exactly one base fixture exists

- **WHEN** the test source tree under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** exactly one fixture SHALL exist for the finite-strain anisotropy tests, named `AniFstrainEvolution.txt`
- **AND** no per-test or per-strain-rate variants (e.g. `AniFstrainEvolution_PureShear.txt`, `AniFstrainEvolution_Saturation.txt`) SHALL be added

#### Scenario: All test cases share the base via MutateInput

- **WHEN** the GoogleTest cases `AniFstrainSimpleShear`, `AniFstrainSaturation`, and `AniFstrainPureShear` are inspected
- **THEN** each SHALL call `RunMDOODZ` on the single shared `AniFstrainEvolution.txt` and rely on a `MutateInput` callback to override its differentiating parameters
- **AND** the override callback SHALL support per-test injection of `bkg_strain_rate`, `Nt`, `dt`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, and `writer_subfolder`
- **AND** the `dt` override path SHALL set `model.dt`, `model.dt0`, AND `model.dt_start` (to avoid the dt-reset bug — `AdvectionRoutines.c:685` resets `dt` from `dt_start`)
- **AND** the `writer_subfolder` override path SHALL `free()` the existing heap pointer and `strdup()` the new value (to avoid the heap-ownership crash on shutdown when the original heap-owned pointer is overwritten with a string literal)

### Requirement: Benchmark documented in `TESTS/AnalyticalSolutions.md`

The new tests SHALL be documented as a numbered section in [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) (matching the format of §8 "Grain Size Evolution"), with the closed-form references for `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear), the derivation from the polar decomposition of `F`, the `min(·, ani_fac_max)` clamp behaviour, the parameter values used by the tests, expected L2 magnitudes, the GTest assertion blocks, the off-by-one output-indexing note, and the embedded PNG. The previous "pure shear is NOT exercised in CI" footnote (§9.7 in the prior revision) SHALL be replaced by a measured-accuracy table for the pure-shear test, parallel to the existing simple-shear table.

#### Scenario: AnalyticalSolutions.md documents both deformation modes as active CI

- **WHEN** [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) is read
- **THEN** §9 SHALL contain measured-accuracy tables for **both** the simple-shear test and the pure-shear test, each with per-step `(γ or t, FS_AR_ana, δ_mean MDOODZ, spatial L2)` columns
- **AND** the prior §9.7 footnote acknowledging the pure-shear segfault as a known MDOODZ bug SHALL be removed (or rewritten as a historical note pointing at this change as the resolution)
- **AND** the closed-form derivations for both `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear) SHALL remain present

#### Scenario: Summary table contains rows for the benchmark

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](../../TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain at least one row per active CI test (simple-shear, saturation, pure-shear) with columns Test, Source, Analytical formula, Assertion, Threshold (or whatever columns the existing table uses)
