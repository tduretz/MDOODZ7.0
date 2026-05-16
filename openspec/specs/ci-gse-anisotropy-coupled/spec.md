## Requirements

### Requirement: Grain-size evolution active under finite-strain anisotropy

When the per-phase `gs` flow-law index is non-zero (i.e. the same gate as in [RheologyDensity.c:477](MDLIB/RheologyDensity.c#L477): `if ( materials->gs[phase] != 0 )`, which selects any of the GSE laws dispatched by [`ReadDataGSE`](MDLIB/FlowLaws.c#L844): `gs ∈ {9, 10, 11, 12, 40, …}`) AND `model.anisotropy = 1`, the per-particle grain size SHALL evolve according to the same law as it would under `model.anisotropy = 0`. Specifically, [`ViscosityConciseAniso`](MDLIB/AnisotropyRoutines.c#L82) SHALL update `*d` (the grain-size output) using the same wattmeter/piezometer-aware local iteration that [`ViscosityConcise`](MDLIB/RheologyDensity.c#L385) routes through, NOT leave `*d = d0` unchanged as the current inline iteration does. The natural implementation is to call `LocalIterationViscoElasticGrainSize` from `ViscosityConciseAniso` when `gs != 0`, but any equivalent in-function update that produces the same `*d` to FP precision is acceptable.

#### Scenario: GSE evolves under anisotropy at the same rate as without anisotropy in homogeneous flow

- **WHEN** a homogeneous calcite-paleowattmeter simulation runs with `gs = 10` and identical strain-rate-controlling parameters in two configurations: (a) `anisotropy = 0`, and (b) `anisotropy = 1` with `aniso_factor = 1.0` and `ani_fstrain = 0` (effectively isotropic)
- **THEN** the steady-state grain size `mesh.d_n` SHALL match between (a) and (b) to within `5e-3` relative error, demonstrating that anisotropy doesn't disable GSE
- **AND** when (b) is repeated with `aniso_factor > 1` or `ani_fstrain = 1` active, `mesh.d_n` SHALL still evolve (i.e., NOT remain equal to the input `d0_ref`); the absolute value may differ from (a) but it MUST evolve

#### Scenario: GSE under anisotropy uses the rotation-invariant strain-rate

- **WHEN** the GSE iteration runs from inside `ViscosityConciseAniso`
- **THEN** the iteration's `Eii` input SHALL be `sqrt(I2(E_rot))`, where `E_rot` is the rotated-frame strain rate already computed at [AnisotropyRoutines.c:176-179](MDLIB/AnisotropyRoutines.c#L176)
- **AND** because `I2` is rotation-invariant, this value SHALL be numerically equal (to FP precision) to `sqrt(I2(E))` of the un-rotated strain rate

### Requirement: CI smoke test for the coupled GSE + anisotropy regime

The CI test suite SHALL include a GoogleTest case (e.g. `RheologyCreep.PinchSwellGSEAniso`) that runs a downsized version of `SETS/PinchSwellGSEAniso.txt` (51×51 grid, `Nt = 10`, both `gs = 10` and `ani_fstrain = 1` active in the calcite layer phase) through MDOODZ's full integrated pipeline and asserts physical-invariant sentinels for BOTH subsystems' state. The test SHALL be sentinel-class, not analytical-L2, because the coupled physics has no closed-form or published reference.

The assertion block SHALL fail if either subsystem is silently disabled (the regression failure mode the merge fixes): if GSE is silently no-op'd under anisotropy (the bug as it stands today), `mesh.d_n` stays uniform → `max(d)/min(d) ≈ 1` violates the heterogeneity assertion; if anisotropy is silently disabled, `mesh.aniso_factor_n` stays at the input `aniso_factor` everywhere → `max(δ)` doesn't exceed 1 violates the evolution assertion.

#### Scenario: Coupled scenario completes with healthy d_n and aniso_factor_n fields

- **WHEN** `TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt` is run by `RheologyCreep.PinchSwellGSEAniso` (51×51 × Nt=10 downsized variant of `SETS/PinchSwellGSEAniso.txt`)
- **THEN** the simulation SHALL exit with status 0 and `Output00010.gzip.h5` SHALL exist
- **AND** all values in `mesh.d_n` SHALL be strictly positive and finite
- **AND** all values in `mesh.aniso_factor_n` (`Centers/ani_fac` in HDF5) SHALL be greater than or equal to `0.999` and finite (FS_AR is bounded below by 1.0)

#### Scenario: GSE is actually evolving under anisotropy (catches the silent no-op)

- **WHEN** the same simulation output is read
- **THEN** `max(mesh.d_n) / min(mesh.d_n)` SHALL exceed `5.0` — the wattmeter is driving heterogeneity in the layer, AS IT SHOULD when GSE is active. If GSE is silently no-op'd (the pre-fix state), grain size stays uniform per-phase and this ratio stays near 1.

#### Scenario: Anisotropy is actually developing (catches if aniso silently disabled)

- **WHEN** the same simulation output is read
- **THEN** `max(mesh.aniso_factor_n)` SHALL exceed `1.0` — the layer's strain-ellipse aspect ratio has grown above the IC of 1, demonstrating that `ani_fstrain = 1` is active.
- **AND** `max(mesh.aniso_factor_n)` SHALL be less than or equal to `ani_fac_max + 1e-9` — the saturation clamp respected.

#### Scenario: Smoke test runtime is bounded

- **WHEN** `RheologyCreep.PinchSwellGSEAniso` runs on a typical CI machine
- **THEN** total wall time SHALL be under 30 seconds (similar in shape and size to the existing `PinchSwellGSESmoke` which runs in ~3 s).

### Requirement: Coupled scenario file `SETS/PinchSwellGSEAniso.txt`

The change SHALL deliver a runnable scenario file at `SETS/PinchSwellGSEAniso.txt` that exercises GSE + finite-strain anisotropy together at full research-scale dimensions. It SHALL be derived from `SETS/PinchSwellGSE.txt` with the calcite layer phase modified to enable anisotropy: `anisotropy = 1` (in the model switches), and on the layer phase `aniso_angle`, `aniso_factor = 1.0`, `ani_fstrain = 1`, `ani_fac_max` (modest value, e.g. 4) added. The matrix phase SHALL keep `ani_fstrain = 0` and `aniso_factor = 1.0` (no anisotropy in the matrix), mirroring the `RiftingComprehensive.txt` pattern.

#### Scenario: Scenario file exists and runs

- **WHEN** the scenario `SETS/PinchSwellGSEAniso.txt` is run by its built executable (`cmake-exec/PinchSwellGSEAniso/PinchSwellGSEAniso`)
- **THEN** the simulation SHALL execute at least one full timestep without exiting non-zero
- **AND** `Output00001.gzip.h5` SHALL contain valid `Centers/d` (grain size) AND `Centers/ani_fac` (anisotropy factor) datasets
- **AND** the scenario file's matrix phase SHALL have `ani_fstrain = 0` (not 1), keeping the matrix isotropic

### Requirement: Existing isolated regimes unaffected

The MDLIB change SHALL NOT regress the isolated GSE-only or anisotropy-only test suites. Specifically:

- `RheologyCreep.GrainSize*` (4 tests: `GrainSizeSteadyState`, `GrainSizeSweep`, `GrainSizeSweepCoupled`, `PinchSwellGSESmoke`) MUST continue to pass with the same numerical results (within FP / float32-HDF5 noise) as before.
- `AnisotropyBenchmark.AniFstrain*` (3 tests: `AniFstrainSimpleShear`, `AniFstrainSaturation`, `AniFstrainPureShear`) MUST continue to pass with the same numerical results as before.
- `AnisotropyBenchmark.{DirectorEvolution, DirectorDtConvergence, StressAnisotropy, StressAnisotropyL2, DirectorEvolutionInterpMode1, DirectorEvolutionInterpMode2}` MUST continue to pass.

The MDLIB change SHALL be gated to take effect ONLY when `materials->gs[phase] != 0` (the same gate as the existing GSE-active check at [RheologyDensity.c:477](MDLIB/RheologyDensity.c#L477)). For phases with `gs == 0` (no GSE) under `anisotropy = 1` — which includes ALL existing live aniso scenarios (`SETS/RiftingComprehensive*.txt`, `SETS/CollisionPolarCartesianAniso.txt`) — the existing inline iteration MUST run unchanged.

#### Scenario: Aniso-only scenarios produce byte-identical output post-fix

- **WHEN** `AnisotropyBenchmark.AniFstrainSimpleShear` is run before and after the MDLIB change (with `gs = 0` for the only phase, the standard config for all `AnisotropyBenchmark` tests)
- **THEN** the per-step `mesh.aniso_factor_n` mean SHALL match between the two runs to within `1e-12` (i.e., effectively unchanged — this scenario doesn't enable GSE, so the new GSE branch should not execute)

#### Scenario: GSE-only scenarios produce byte-identical output post-fix

- **WHEN** `RheologyCreep.GrainSizeSteadyState` is run before and after the MDLIB change (with `anisotropy = 0`, the standard config for all `RheologyCreep` tests except `PinchSwellGSEAniso`)
- **THEN** the resulting `mesh.d_n` mean SHALL match between the two runs to within `1e-12` (this scenario doesn't enable anisotropy, so the new aniso path's GSE branch is not entered)

### Requirement: Documentation in `TESTS/AnalyticalSolutions.md`

The coupled regime SHALL be documented in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) as a new subsection — either §8.3 (extending the GSE section's "Coverage extensions") or §10 (new top-level section). The documentation SHALL state the silent-no-op bug that this change resolves, the fix approach (calling `LocalIterationViscoElasticGrainSize` from `ViscosityConciseAniso`), the assertion block for the smoke test, and the rationale for sentinel-class (no analytical L2) thresholds.

#### Scenario: AnalyticalSolutions.md documents the coupled regime

- **WHEN** [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain a subsection (or section) describing `RheologyCreep.PinchSwellGSEAniso` with: the parameter file paths, the assertion block, an explanation of WHY each sentinel exists (specifically that `max(d)/min(d) > 5` catches GSE silent no-op and `max(δ) > 1` catches aniso silent disable), and a reference to the MDLIB code change at `ViscosityConciseAniso`.

#### Scenario: Summary table contains a row for the coupled smoke test

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain a row for the coupled smoke test with columns Test, Source, Analytical formula (= "sentinel — no closed form"), Assertion, Threshold.

### Requirement: Reuse existing test-helper utilities

The smoke test SHALL use the existing `getMin/MaxFieldValue` HDF5-reader helpers from [TESTS/TestHelpers.h](TESTS/TestHelpers.h) and SHALL NOT introduce new files under `TESTS/` named `*Helpers.h`, `*Helpers.cpp`, or similar. No new utility functions are needed — reading two scalar mesh fields (`Centers/d` and `Centers/ani_fac`) and computing min/max/mean is fully covered by the existing helpers.

#### Scenario: No new test-helper utilities added

- **WHEN** the smoke test is added to [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp)
- **THEN** the implementation SHALL call only existing helpers from `TestHelpers.h` and SHALL NOT introduce new helper files
