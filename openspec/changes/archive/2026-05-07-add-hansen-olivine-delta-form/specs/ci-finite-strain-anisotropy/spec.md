## ADDED Requirements

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

## MODIFIED Requirements

### Requirement: Research note `SETS/AnisoFstrainResearch/hansen2012_comparison.md`

The change SHALL deliver a first-class research-grade markdown record at `SETS/AnisoFstrainResearch/hansen2012_comparison.md` documenting the form mismatch between MDOODZ's `δ = min(FS_AR, ani_fac_max)` and the empirical fit from Hansen, Zimmerman, and Kohlstedt (2012) *Nature* 492, 415–418, doi:10.1038/nature11671. The note SHALL co-locate with a runnable scenario file `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` that reproduces the comparison conditions. The research note SHALL be flagged as future-work motivation for the `ani_fstrain = 1` (`min`) form, NOT used as a CI gate (digitisation noise of Hansen+12 lab data is comparable to the form mismatch we'd be flagging — using it as a CI threshold would conflate "MDOODZ has a bug" with "the model is approximate").

The note SHALL additionally contain a section (§7 or later) documenting the calibrated `ani_fstrain = 2` form delivered by this change, with: (a) the calibration procedure used to fit `M_inf` and `γ_e` from the **combined 38-point Hansen-group olivine torsion dataset** (26 from Hansen+14 Table 1 + 12 from Hansen+16 Part 1 Fig. 6 + Table 1), (b) the closed-form expressions for `γ_eff(FS_AR)` and `M(γ)`, (c) the resulting δ-asymptote `δ_∞ = 14.13` (linear formula) confirmed by Hansen+16 Eq. 3-4's stress-aware direct calc `δ = (1.39/0.73)^4.1 = 14.02` to **100.8% agreement** (independent two-method validation), and how it compares to typical user-side `ani_fac_max` values (4, 6, 8), (d) a quantitative comparison plot (`SETS/AnisoFstrainResearch/hansen_olivine_calibration.png`) showing the candidate curves on two panels — Panel (a): M(γ) with all 38 + 3 datapoints and per-dataset fits; Panel (b): δ(γ) with `min`, `erfc`, Hansen+12 3-pt, Hansen+16 12-pt, and combined 38-pt fits overlaid with lab data, (e) the RMS error of each form against the combined 38-point reference (γ ∈ [0, 19]): `min` ≈ 2.87, `erfc` ≈ 3.11, Hansen+12 3-pt ≈ 2.07, Hansen+16 12-pt ≈ 1.93, **combined 38-pt ≈ 1.85**, (f) a note on the **structural caveat** that `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` is a simple-shear inversion and reinterpreting it under other deformation modes is interpretation-questionable, (g) a note on the **stress-orientation caveat** from Hansen+16 Part 2 §4.1.1 / Fig. 10: real δ depends on alignment between principal stress and texture's strong axis (~10× spread), and MDOODZ's scalar δ corresponds to the worst-case (texture-aligned) value, (h) a note on the **ε_c ≈ 1 co-evolution** finding from Hansen+16 Part 2 Eq. 10: texture and grain size evolve on the same strain timescale, consistent with MDOODZ's wattmeter coupling.

The accompanying calibration script SHALL be checked in at `SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py` (Python, no scipy required — uses a 2D coarse-then-fine grid search to minimise SSE on `M_inf`, `γ_e`), reproducible from a clean Python environment.

#### Scenario: Research note exists with required content

- **WHEN** `SETS/AnisoFstrainResearch/hansen2012_comparison.md` is read
- **THEN** it SHALL contain (a) the live MDOODZ formula `δ = min(FS_AR, ani_fac_max)` with code references to [AnisotropyRoutines.c:46-54](../../MDLIB/AnisotropyRoutines.c#L46) and [RheologyParticles.c:634-678](../../MDLIB/RheologyParticles.c#L634), (b) the Hansen+12 fit `δ(γ) = 13·tanh(0.25γ) + 1` with the asymptote at 14 and half-saturation at γ ≈ 2.2, (c) a quantitative comparison at several γ values (table or short prose) showing both forms diverge in shape but converge in asymptote, (d) a discussion noting the `min()` form is too sharp at the knee and rises too fast at low strain, (e) a future-work section pointing at the commented-out `erfc` blend in `AnisoFactorEvolv` and the option of replacing with a tanh form, and (f) a full citation to Hansen, Zimmerman, Kohlstedt (2012)
- **AND** it SHALL contain a section (§7 or later) documenting the calibrated `ani_fstrain = 2` form: the calibration procedure (38-point combined least-squares fit to Hansen+14 Table 1 + Hansen+16 Part 1), the closed-form `M(γ_eff(FS_AR))` mapping, the values `M_inf = 0.536` and `γ_e = 3.96`, the δ-asymptote `δ_∞ = 14.13` (validated to 100.8% by Hansen+16 Eq. 3-4 stress-aware direct calc 14.02), the comparison plot at `SETS/AnisoFstrainResearch/hansen_olivine_calibration.png`, the RMS-error comparison of all candidate forms (`min` ≈ 2.87, `erfc` ≈ 3.11, Hansen+12 3-pt ≈ 2.07, Hansen+16 12-pt ≈ 1.93, combined 38-pt ≈ 1.85), the simple-shear caveat for `γ_eff(FS_AR)`, the Hansen+16 Part 2 stress-orientation caveat, and the Hansen+16 Part 2 ε_c ≈ 1 co-evolution finding
- **AND** it SHALL reference the calibration script `SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py` as the reproducible artifact behind the calibration values
- **AND** it SHALL NOT be referenced from any `TESTS/` file or CI script (it remains research, not regression gate)

#### Scenario: Companion scenario exists

- **WHEN** the research directory is inspected
- **THEN** `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` SHALL exist with `anisotropy = 1`, `ani_fstrain = 1`, simple-shear BCs, and parameters chosen to reach γ ≈ 4 (sufficient to span the Hansen+12 fit's saturation region)

#### Scenario: Calibration script and plot are checked in

- **WHEN** `SETS/AnisoFstrainResearch/` is inspected
- **THEN** `calibrate_hansen_olivine.py` SHALL exist as an executable Python script that, when run, prints the per-dataset fit summary table (Hansen+12, Hansen+14, Hansen+16 Part 1, combined), both candidate δ-asymptotes (linear formula and Hansen+16 Eq. 3-4 stress-aware) with their percentage agreement, and the RMS errors of all candidate δ-forms vs the combined 38-point dataset, AND writes `hansen_olivine_calibration.png` to the same directory
- **AND** the script SHALL NOT depend on scipy (numpy + matplotlib only)
