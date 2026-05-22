## Why

Grain Size Evolution (GSE) and finite-strain anisotropy (`ani_fstrain = 1`) are now both fully CI-gated in isolation:

- GSE: `ci-grain-size-evolution` capability, `RheologyCreep.GrainSizeSteadyState/Sweep/SweepCoupled/PinchSwellGSESmoke` (4 tests, both dislocation-only and coupled regimes).
- Anisotropy: `ci-finite-strain-anisotropy` capability, `AnisotropyBenchmark.AniFstrainSimpleShear/Saturation/PureShear` (3 tests, including the recently-fixed pure-shear path).

But the **coupled regime — both `gs = 10` and `ani_fstrain = 1` active in the same phase** — has zero CI coverage and zero production usage. No scenario in `SETS/` currently exercises both at once: `PinchSwellGSE.txt` has `gs = 10` but no anisotropy; `RiftingComprehensive.txt` and `CollisionPolarCartesianAniso.txt` have `ani_fstrain = 1` but no `gs` flow law (the `gsel = 1` lines in those files were a dead switch removed in `fix-grain-size-evolution`).

Two pieces of MDLIB code write to overlapping state under the coupled regime:

1. **[AnisotropyRoutines.c:531-641](MDLIB/AnisotropyRoutines.c#L531)** writes `mesh->d_n[c0]` from per-phase `phase_perc_n` weights. This runs every step inside the anisotropy P2G pass.
2. **[RheologyDensity.c:117 `LocalIterationViscoElasticGrainSize`](MDLIB/RheologyDensity.c#L117)** updates the per-cell grain size via the local Newton iteration. Independently writes to `mesh->d_n`.

If both run in the coupled regime, write ordering and which one is the source of truth is unverified. Additionally, [`EffectiveStrainRate`](MDLIB/RheologyDensity.c#L116) rotates the strain rate into the anisotropy frame (`ani, d1, d2` parameters) — but `LocalIterationViscoElasticGrainSize` consumes a strain rate that may or may not be the rotated one. If the GSE iteration sees the un-rotated strain rate while the rest of the rheology sees the rotated one, the wattmeter gets the wrong driver under anisotropy.

A coupled scenario `PinchSwellGSEAniso` (calcite layer with both GSE and finite-strain anisotropy) would reveal latent bugs and demonstrate the "anisotropy enhances necking" interaction physically.

**Why now:** the next research direction (smoother δ-saturation per Hansen+12) is gated on having the coupled scenario actually run, because parameter calibration needs concrete trajectories to fit against. Doing this merge first → a working `PinchSwellGSEAniso` baseline → informed decisions about δ-form refinements later.

## What Changes

- **Investigate** the code-level interactions between GSE and anisotropy in MDLIB:
  - Write ordering on `mesh->d_n` (anisotropy P2G vs GSE local iteration).
  - Strain-rate input to `LocalIterationViscoElasticGrainSize` when `anisotropy = 1` (rotated vs un-rotated frame).
  - Any other latent assumptions in either path that break when both are on.
- **Resolve any conflicts** with minimum-scope fixes consistent with the [`fix-pure-shear-fstrain-segfault` design D3](openspec/changes/archive/2026-05-06-fix-pure-shear-fstrain-segfault/design.md) scope wall (one-line / few-line edits; pause and escalate if architectural).
- **Add `SETS/PinchSwellGSEAniso.txt`** — a calcite layer scenario that turns on GSE (`gs = 10`, calcite paleowattmeter) AND finite-strain anisotropy (`ani_fstrain = 1`) in the same phase. Based on `SETS/PinchSwellGSE.txt`, with the anisotropy parameters added.
- **Add a CI smoke test** `RheologyCreep.PinchSwellGSEAniso` (or similar) that runs a downsized version of the new scenario through MDOODZ's full integrated pipeline and asserts physical-invariant sentinels (no NaN, both `mesh.d_n` and `mesh.aniso_factor_n` finite, layer develops both grain-size heterogeneity and anisotropy strength).
- **Document the coupled regime** as a new subsection of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8 or §9 (whichever fits — leaning toward §8.3 since the test is a smoke/sentinel rather than an analytical L2 gate).

## Capabilities

### New Capabilities

- `ci-gse-anisotropy-coupled`: CI smoke gate for the coupled GSE + finite-strain anisotropy regime via `PinchSwellGSEAniso`. Sentinel-class test (finite, positive, heterogeneity emerging, both subsystems' state evolving) — NOT an analytical L2 gate, because the coupled physics has no closed form.

### Modified Capabilities

(potentially `ci-grain-size-evolution` and/or `ci-finite-strain-anisotropy` if the investigation reveals that one of those specs needs to acknowledge the coupled write to `mesh.d_n` or the anisotropy-aware strain-rate path. Decision deferred to design.md after investigation.)

## Impact

- **Code under change**: bounded by the investigation. Two suspected hot zones: `AnisotropyRoutines.c:531-641` (mesh.d_n write inside P2G) and `LocalIterationViscoElasticGrainSize` (strain-rate input). Either or both may need a small change. If the investigation finds the architecture is fine and just needs a scenario + smoke test, code under change is ZERO.
- **Tests**: 1 new GTest case `RheologyCreep.PinchSwellGSEAniso` in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp). New fixture `TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt` (downsized 51×51 × Nt=10 variant of the SETS scenario, mirroring the `PinchSwellGSESmoke` pattern).
- **Documentation**: new subsection in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) describing the coupled scenario, the assertion block, and any new physics gating.
- **CI runtime**: ~3 s estimated (51×51 × 10 steps, similar to `PinchSwellGSESmoke`).
- **Risk**: **medium-high** — touches the interaction surface between two non-trivial subsystems. Both subsystems are now well-tested in isolation, but their interaction is, by definition, unproven. Mitigations:
  - Full regression suites (`AnisotropyBenchmarkTests` 9, `RheologyCreepTests` 12) before declaring done.
  - Manual smoke run of `RiftingComprehensive.txt` for ≥ 2 timesteps (uses anisotropy without GSE — should be unaffected).
  - Manual smoke run of `PinchSwellGSE.txt` for ≥ 2 timesteps (uses GSE without anisotropy — should be unaffected).
- **Out of scope**:
  - Smoother δ-saturation form (e.g., tanh-fit to Hansen+12 or activating the commented-out `erfc` blend in `AnisoFactorEvolv`). The `min(FS_AR, ani_fac_max)` form stays. Future work — see `SETS/AnisoFstrainResearch/hansen2012_comparison.md`.
  - Quantitative L2 comparison against published GSE+aniso results (Schmalholz & Duretz 2017's pinch-swell paper had GSE only, no anisotropy; no equivalent published reference for this combination).
  - Per-mineral parameterisation review of the anisotropy `ani_fac_max` and `aniso_angle` defaults in calcite scenarios. The values in `PinchSwellGSEAniso.txt` will be set to plausible defaults, not calibrated.
  - Restructuring the `mesh.d_n` write authority. If the investigation finds that the dual-writer pattern is actually broken (e.g. anisotropy P2G clobbers GSE's value mid-iteration), the fix is in scope; redesigning the write authority itself (e.g. introducing a single source of truth) is a bigger change and would be its own follow-up.
