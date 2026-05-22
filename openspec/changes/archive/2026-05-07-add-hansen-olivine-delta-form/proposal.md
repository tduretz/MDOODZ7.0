## Why

The current finite-strain anisotropy saturation forms in MDOODZ — `δ = min(FS_AR, ani_fac_max)` (default, `ani_fstrain = 1`) and the `erfc`-blend form (build-time `-DANISO_ERFC_BLEND`, prototyped during the overnight visualisation queue) — both **fail to reproduce empirical olivine δ(γ) data** from Hansen, Zimmerman & Kohlstedt (2012) *Nature* 492, 415. An analytical comparison performed in this session (script + plot in [SETS/AnisoFstrainResearch/hansen2012_erfc_comparison.png](SETS/AnisoFstrainResearch/hansen2012_erfc_comparison.png)) showed:

- Both forms underestimate δ at low strain (γ ∈ [0, 1.5]) by a factor of ~3 because they grow as `~γ` near origin while Hansen+12's empirical δ(γ) grows as `~3.25γ`.
- The `erfc` form's RMS error against Hansen+12 (over γ ∈ [0, 5], the lab range) is **2.70**, vs **2.03** for the `min` form — i.e. `erfc` is 33% **worse**, not better. The smoother knee did not help because the dominant deviation is at low strain, not at the saturation transition.
- The structural cause: both MDOODZ forms use `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4+1)` (strain-ellipse aspect ratio) as the input variable. Hansen+12 instead reports δ as a function of the **fabric-strength index M** — a CPO measure unrelated to FS_AR.

Hansen+12's explicit linear fit (Fig 3b) is

```
δ = (24.5 ± 5.5) · M + 1
```

with M the fabric-strength index (0 for random fabric, 1 for single crystal). Two follow-up papers materially expand the calibration data:

- **Hansen, Zhao, Zimmerman & Kohlstedt (2014)** *EPSL* 387, 157–168 ("Protracted fabric evolution in olivine"): 26 (γ, M) datapoints over γ ∈ [0, 18.7], unified using consistent M-index methodology (Skemer+05). Hansen+14 reanalyzed previously published samples from Bystricky+00, Zhang & Karato 2005, and Hansen+12a/b/c, so its Table 1 is one harmonised dataset.
- **Hansen, Warren, Zimmerman & Kohlstedt (2016)** *EPSL* 445, 92–103 Part 1: 12 datapoints from 5 NEW samples (PT0716, PT0718, PT0751, PT0752, PT0754) not in Hansen+14, plus a stress-exponent-aware direct calculation `δ = (F_s/F_w)^n` with `F_w = 0.73`, `F_s = 1.39`, `n = 4.1` → **`δ_∞ ≈ 14`**, materially below the 16.22 a 3-point Hansen+12 fit's `24.5·M_inf+1` extrapolation predicts.

Hansen+16 Part 2 (JGR 121, 7137) develops a full pseudo-Taylor + director-method micromechanical model whose textural-evolution panel (Fig. 7) reproduces both Hansen+14 and Hansen+16 Part 1 data and confirms `M(γ)` saturates near 0.5–0.6 by γ ≈ 10. The model is too heavy for MDOODZ's per-cell rheology pass (per-grain integration), but Part 1's empirical formula is appropriate for a scalar δ-form.

M is not tracked in MDOODZ (would require a CPO model). A calibrated mapping `M(FS_AR)` lets us evaluate the Hansen olivine δ-formula using MDOODZ's existing deformation gradient. Calibration uses the **combined Hansen-group olivine torsion dataset (38 datapoints: 26 from Hansen+14 + 12 from Hansen+16 Part 1)**. Best-fit constants: `M_inf = 0.536`, `γ_e = 3.96`, linear-formula asymptote `24.5·0.536 + 1 = 14.13` — agreeing with Hansen+16 Eq. 3-4's stress-aware direct calc (14.02) to **100.8%**, an independent two-method validation that the Hansen-olivine δ-asymptote is `δ_∞ ≈ 14`. RMS error against the 38-point dataset: 1.85, vs 2.07 for the Hansen+12 3-point fit and 2.87 for the existing `min` form.

This change adds that translation as a third δ-form, gated by `ani_fstrain = 2` (additive — the existing `ani_fstrain = 1` form is preserved unchanged for backward compatibility).

## What Changes

- **Add a calibrated mapping** `M_from_FS_AR(FS_AR)` in [MDLIB/AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) — a sigmoidal fit through the three Hansen+12 (FS_AR, M) datapoints (FS_AR computed from γ via the analytical simple-shear formula). Form to be determined during implementation; candidates: `M = a · tanh(b · log(FS_AR))` (one calibrated parameter beyond asymptote a) or `M = M_max · (1 - exp(-c · (FS_AR-1)))`.
- **Add a new `ani_fstrain = 2` dispatch arm** in [`AnisoFactorEvolv`](MDLIB/AnisotropyRoutines.c#L46) returning `δ = 24.5·M(FS_AR) + 1`, capped at `ani_fac_max` for safety.
- **Backward compatibility (hard requirement)**: `ani_fstrain = 0` (off) and `ani_fstrain = 1` (current `min` form) unchanged byte-for-byte. All existing CI tests pass without modification. All existing `SETS/RiftingComprehensive*.txt`, `SETS/CollisionPolarCartesianAniso.txt`, `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt`, etc. behave identically. The new form is **opt-in only**.
- **Remove or deprecate the `-DANISO_ERFC_BLEND` build-time toggle** — it was an in-source prototype for our overnight queue; its functionality is now absorbed (and deemed empirically inferior) by the new dispatch. Decision deferable to design — keep as `ani_fstrain = 3` for completeness, or remove entirely.
- **Add 1 new CI test** `AnisotropyBenchmark.AniFstrainHansenOlivine` validating the new form against the analytical reference (the calibrated `δ(γ_eff(FS_AR))` formula evaluated symbolically) over the lab strain range γ ∈ [0, 4].
- **Update [SETS/AnisoFstrainResearch/hansen2012_comparison.md](SETS/AnisoFstrainResearch/hansen2012_comparison.md)** with: the calibration procedure, the comparison plot showing `min` / `erfc` / Hansen-calibrated curves vs the Hansen+16 12-point dataset, the RMS-error comparison, and the FS_AR-vs-M structural caveat (Hansen+16 Part 2 §4.1.1 stress-orientation dependence + ε_c ≈ 1 co-evolution of texture and grain size).
- **Add the calibration script** [SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py](SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py) (Python, no scipy required — already produced this session, as a standalone reproducible artifact). Calibration plot at [SETS/AnisoFstrainResearch/hansen_olivine_calibration.png](SETS/AnisoFstrainResearch/hansen_olivine_calibration.png).

## Capabilities

### New Capabilities

(none — we extend the existing capability rather than add a new one)

### Modified Capabilities

- `ci-finite-strain-anisotropy`: extend the `ani_fstrain` enumeration documentation to include the new value `= 2` for the Hansen-calibrated olivine form (Hansen+12 Fig 3b linear formula composed with M(γ) fit to the Hansen+16 Part 1 12-point dataset). Add a requirement documenting the analytical reference and the new test case.

## Impact

- **Code under change**: MDLIB additions confined to [`AnisotropyRoutines.c`](MDLIB/AnisotropyRoutines.c) — 1 new helper (`M_from_FS_AR`), 1 new dispatch arm in `AnisoFactorEvolv`, ~30 lines total. No changes to existing arms.
- **Tests**: 1 new GTest case `AnisotropyBenchmark.AniFstrainHansenOlivine` in [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp). Reuses the existing `AniFstrainEvolution.txt` base fixture via `MutateInput` (sets `ani_fstrain = 2` per phase) — no new fixture file.
- **Documentation**: updated [SETS/AnisoFstrainResearch/hansen2012_comparison.md](SETS/AnisoFstrainResearch/hansen2012_comparison.md) with a new section on the calibrated form. Calibration plot at [SETS/AnisoFstrainResearch/hansen_olivine_calibration.png](SETS/AnisoFstrainResearch/hansen_olivine_calibration.png). Calibration script at [SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py](SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py).
- **CI runtime**: +1 test × short timesteps × small homogeneous box → ~0.2 s. Within budget.
- **Risk**: **low**. Purely additive — the new form is gated by `ani_fstrain == 2` which is currently a no-op default value. No existing scenario uses this value. Backward compatibility is structural (different dispatch arm), not requiring careful preservation logic.
- **Out of scope**:
  - Modeling M directly via a CPO model (VPSC, D-Rex, Tommasi+2000-style, Hansen+16 Part 2 director method). That remains future work — `M(FS_AR)` is an empirical translation, not a physics model of fabric evolution.
  - Calibrating the form against minerals other than olivine. The Hansen-group lineage only publishes data for olivine; other minerals (calcite Barnhoorn+04, quartz, mica) need their own calibrated `M(FS_AR)`. Per-mineral parameterisation is future work — for now, the new form ships with olivine-only constants.
  - Quantitative validation of MDOODZ scenarios against Hansen+16 lab data (digitisation noise + sample-to-sample scatter is comparable to the form mismatch, same reasoning as in `gate-finite-strain-anisotropy`).
  - **Stress-orientation-aware viscous anisotropy** (Hansen+16 Part 1 Eq. 2 / Part 2 Fig. 10): the real δ depends on the angle between principal stress and the texture's strong axis (~10× spread). MDOODZ's δ remains a scalar magnitude — i.e. an upper-bound, texture-aligned δ. A future change could replace the scalar with a Hill-tensor parameterisation (6 components per cell). Out of scope here.
  - **Replace `FS_AR → M` shortcut with a director-method texture tracker** (Hansen+16 Part 2 Eqs. 5–9). Would require advecting per-cell ODFs / per-particle slip-system rotations. Out of scope.
  - Phase-3 unified `η(d, δ, F-orientation)` constitutive coupling. Still a separate research direction.
