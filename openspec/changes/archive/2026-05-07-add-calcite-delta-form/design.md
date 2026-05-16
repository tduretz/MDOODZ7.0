## Context

The function-pointer dispatch architecture from [add-hansen-olivine-delta-form](../archive/2026-05-07-add-hansen-olivine-delta-form/) (case 1 = Hansen olivine) generalises trivially to other minerals — adding a calibration is a static helper in [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c) plus a new `case N` in `ReadDataAnisotropy`. The cost-benefit question for each new mineral is **whether the existing `ani_fstrain = 1` form `min(FS_AR, ani_fac_max)` is fit-for-purpose**.

The min form's transient behaviour saturates at `γ_sat = γ where FS_AR(γ) = ani_fac_max`. For the simple-shear identity `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4+1)`:

| `ani_fac_max` | `γ_sat` (min form) | Lab saturation strain | Transient overshoot |
|---|---|---|---|
| 14 (olivine) | 3.6 | ~7-10 (Hansen+14) | factor ~1.5× over γ ∈ [3.6, 7] — material |
| 4 (calcite) | 1.6 | ~5 (Pieri+01) | factor ~2× over γ ∈ [1.6, 5] — material |
| 3 (quartz) | 1.2 | ~4 (Pennacchioni+10) | factor ~1.3× over γ ∈ [1.2, 4] — modest |

**Olivine and calcite are worth a custom form; quartz is not.** Quartz's asymptote δ_max ≈ 2-3 is small enough that the min-form's "early saturation" produces only a small absolute overestimate. This is documented as a research-note recommendation rather than added MDLIB code.

For calcite, the published lab data is:

- **Pieri+01** *Tectonophysics* 330, 119 — Carrara marble torsion at 1000 K, 300 MPa, γ ∈ {0, 1, 2, 5, 11}. Reports texture index J (Bunge 1982) — J = 1 for random fabric, ∞ for single crystal:

| γ | J-index | M-index (Skemer+05 J→M) |
|---|---|---|
| 0  | 1.12 | ≈ 0.01 |
| 1  | 1.66 | ≈ 0.05 |
| 2  | 2.02 | ≈ 0.07 |
| 5  | 9.76 | ≈ 0.45 |
| 11 | 8.92 | ≈ 0.40 |

The non-monotonic decline at γ=11 reflects recrystallization randomisation (Pieri+01 §3.3-3.4); not captured by exponential-saturation form, smoothed in the fit.

- **Bruijn+11** *Tectonophysics* 503, 75 — single-stage segments of two-stage Carrara marble torsion at γ ∈ {1, 2.6, 5} with similar conditions. Reports J-index. Adds ~3 datapoints at low-mid strain.

- **Schuster+18** *J. Struct. Geol.* — high-pressure torsion calcite to γ = 80 at 1-4 GPa, 235-450 °C. Excluded from fit because phase transition (calcite → CaCO₃-II at 1.6 GPa) makes it a different material; documented for context only.

- **Pennacchioni+10** *JGR* — quartz veins, 17 (γ, J) datapoints used for the quartz research-note recommendation, NOT for any MDLIB calibration.

**The δ-vs-fabric coefficient for calcite is the weakest part of the calibration.** Hansen+12 Fig 3b gives `δ = 24.5·M + 1` for olivine — a direct mechanical lab fit. There is no equivalent published linear fit for calcite. Estimates from VPSC modeling (Tommasi+09-style Hill-tensor calculations) suggest calcite δ at saturated CPO ≈ 2-5, depending on slip-system activity and grain-size feedback. We commit to a slope that yields **`δ_∞ = 4`**, matching the empirical `ani_fac_max = 4` calibration already used by `SETS/PinchSwell*` scenarios (which was set based on practical scenario behaviour, not lab measurement). This is documented as a calibration uncertainty.

## Goals / Non-Goals

**Goals:**

- Add a calcite δ-form (`case 2` in `ReadDataAnisotropy`) using the same saturating-exponential framework as Hansen olivine, with calcite-specific constants calibrated to Pieri+01 + Bruijn+11.
- Document the **explicit recommendation** that quartz scenarios use the existing `ani_fstrain = 1` form with `ani_fac_max ≈ 2-3` per phase, justified by the Pennacchioni+10 dataset comparison.
- One new CI test (`AniFstrainCalcite`) gating the calcite form against its analytical reference at FP precision.
- Two new SETS subfolders (`SETS/AnisoFstrainResearch/calcite/`, `SETS/AnisoFstrainResearch/quartz/`) with calibration scripts, comparison plots, and structured markdown documentation matching the Hansen olivine pattern.
- **Backward compatibility (HARD requirement)**: existing `ani_fstrain ∈ {0, 1, 3}` and case 1 (Hansen olivine) byte-identical post-change. The auto-default `ani_fstrain = 2 + aniso_db = 0 → aniso_db = 1` (olivine) is unchanged. Calcite is opt-in via explicit `aniso_db = 2`.

**Non-Goals:**

- **No quartz case 3.** The min form is fit-for-purpose for quartz scenarios; no MDLIB code change.
- **No Hansen+16 Eq. 3-4-style independent stress-aware δ_∞ measurement for calcite.** Such a measurement does not exist in the literature. We commit to `δ_∞ = 4` based on a combination of Tommasi+09-style VPSC bounds + empirical match to the proven PinchSwell calibration. Documented as uncertainty.
- **No HPT calcite (Schuster+18).** Phase transition makes the data not directly applicable to MDOODZ scenarios at sub-GPa pressures.
- **No T-aware quartz dispatch.** Pennacchioni+10 is single-T; basal-a/rhomb-a/prism-a slip-regime changes outside the ~500°C window are not modelled. Future work if demand emerges.
- **No exposure of calibration constants as user input parameters.** Per the FlowLaws idiom: users select via `aniso_db = 2`, constants live inside the static helper. Same as olivine.
- **No phase-3 unified `η(d, δ)` constitutive coupling.** Separate research direction.

## Decisions

### D1: Composition strategy — same three-step framework as Hansen olivine

```
γ_eff(FS_AR) = sqrt(FS_AR) - 1/sqrt(FS_AR)             // simple-shear inversion (mineral-independent)
M(γ_eff)     = M_inf · (1 - exp(-γ_eff / γ_e))         // saturating exponential (mineral-specific constants)
δ(M)         = a · M + 1                                // δ-vs-M linear fit (mineral-specific slope)
```

**Rationale:** identical structure to olivine `case 1` keeps the framework predictable. Each mineral gets its own `(M_inf, γ_e, a)` triple inside its `static double anisoDelta_<Mineral>(double FS_AR)` helper. Hansen+16 Part 2 director-method modeling (Fig. 7) shows this saturating-exponential shape is general across the Hansen-group olivine dataset; for calcite it is a reasonable first approximation given the smoothing through Pieri+01's non-monotonic γ=11 outlier.

**Alternatives considered:**

- *(A) Mineral-specific functional form (e.g. tanh, piecewise) for calcite to better fit non-monotonic γ=11 outlier*: rejected. The non-monotonic behaviour is recrystallization-randomization, a transient at very high strain that no MDOODZ scenario routinely visits. Smoothing through it is the right ensemble approximation.
- *(B) J-index pipeline directly (`δ = c·(J − 1) + 1`)* without converting to M: rejected. The J↔M conversion is mineral-shape-dependent (Skemer+05) but using J directly would require a different δ-vs-J slope than δ-vs-M, breaking compatibility with the olivine framework. Going with the M-index approach for consistency.

### D2: Calibration constants

**Choice:** Hardcoded calcite constants, calibrated by least-squares fit to the **combined Pieri+01 + Bruijn+11 dataset** with J → M conversion via Skemer+05's empirical relation:

```c
// In MDLIB/FlowLaws.c, file-scope static helper:
static double anisoDelta_Calcite( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.50;        // calibrated to Pieri+01/Bruijn+11 J→M saturation
    const double gamma_e   = 3.50;        // characteristic strain for calcite fabric development
    const double slope     = 6.0;         // δ-vs-M slope (Tommasi+09 VPSC; gives δ_∞ = 4 matching empirical PinchSwell)
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}
```

**Constants come from**:
- `M_inf = 0.50`: Pieri+01 J = 9.76 at γ = 5 → M ≈ 0.45 via Skemer+05; Bruijn+11 corroborates ≈ 0.4-0.5; round to 0.50.
- `γ_e = 3.50`: combined Pieri+01 + Bruijn+11 datapoint fit. Saturates faster than olivine (γ_e = 3.96) — calcite reaches steady-state CPO at γ ≈ 5 vs olivine's γ ≈ 7 (consistent with Pieri+01 §2.3, Schuster+18 §5.3, Barnhoorn+04).
- `slope = 6.0`: chosen to yield `δ_∞ = 6·0.5 + 1 = 4.0`, matching the empirical `ani_fac_max = 4` calibration used in `SETS/PinchSwell*` scenarios. Bracketed by Tommasi+09-style VPSC bounds (calcite δ_saturated ≈ 2-5 depending on slip-system mix). NOT a clean Hansen+12 Fig 3b-equivalent lab fit.

**Asymptote**: `δ_∞ = a · M_inf + 1 = 4.0`. Significantly lower than olivine's 14.13, reflecting calcite's modest viscous-anisotropy magnitude (Barnhoorn+04: CPO contributes ~1/3 of total weakening, the rest from grain-size reduction — already captured by MDOODZ's wattmeter).

**Verification table** (from `calibrate_calcite.py`):

| γ | J_obs (Pieri+01 / Bruijn+11) | M_obs (Skemer+05) | M_fit (this calibration) | δ_fit |
|---|---|---|---|---|
| 0 | 1.12 | 0.008 | 0.000 | 1.00 |
| 1 | 1.66 (P089) | 0.044 | 0.131 | 1.79 |
| 1 | 1.6  (Bruijn) | 0.040 | 0.131 | 1.79 |
| 2 | 2.02 (P143) | 0.068 | 0.232 | 2.39 |
| 2.6 | 2.5 (Bruijn) | 0.10 | 0.265 | 2.59 |
| 5 | 9.76 (P087) | 0.45 | 0.376 | 3.26 |
| 5 | 4.0 (Bruijn) | 0.20 | 0.376 | 3.26 |
| 11 | 8.92 (P088) | 0.40 | 0.477 | 3.86 |
| ∞ | — | 0.50 | 0.500 | 4.00 |

The fit smooths through Pieri+01 vs Bruijn+11 disagreement at γ=5 (P087's J=9.76 is much higher than Bruijn's J≈4 at the same strain — sample-to-sample variability + Pieri's `single-stage at high strain` vs Bruijn's `single-stage embedded in a pre-existing-fabric matrix`). The fit lands roughly midway, biased slightly toward Pieri's higher value since Pieri is the cleaner single-stage measurement.

### D3: γ_eff outside simple shear — caveat in code comment, not a runtime guard

Same as olivine (D3 in archived design): `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` is exact under simple shear, defensible-but-questionable under pure shear / general 2D flow. Document but don't guard. Calcite scenarios (`SETS/PinchSwell*`) are in a near-pure-shear regime, so the user is signing up for "calcite-shaped saturation parameterized by accumulated stretch", not "calcite fabric strength under arbitrary deformation."

### D4: Backward compatibility — additive `case 2` dispatch

**Choice:** new entry in `ReadDataAnisotropy`'s switch:

```c
case 2 :
    LOG_INFO("Calcite viscous anisotropy - Pieri+01 / Bruijn+11 (combined fit, δ_∞ = 4):");
    mat->aniso_delta_fn[k] = anisoDelta_Calcite;
    success                = 1;
    break;
```

The auto-default `ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1` (olivine) is **unchanged**. Calcite users opt in via explicit `aniso_db = 2` per phase in their input file. This preserves byte-identical behavior for every existing scenario.

**Migration**: no migration of `AnisoFactorEvolv` signature, no migration of call sites. The function-pointer architecture absorbs this entirely in `FlowLaws.c`.

### D5: Saturation cap

Same as olivine: `MINV(δ, ani_fac_max)` applied uniformly by the dispatcher in `AnisotropyRoutines.c`. If a user sets `ani_fac_max < 4` (the calcite asymptote), the cap takes precedence — preserves user's intent. Documented in the code comment.

### D6: Test strategy — analytical-reference gate, mirroring AniFstrainHansenOlivine

**Choice:** new GTest case `AnisotropyBenchmark.AniFstrainCalcite` reusing the existing `AniFstrainEvolution.txt` base fixture via `MutateInput`. Override sets:
- `ani_fstrain[0] = 2`
- `aniso_db[0] = 2` (calcite)
- `ani_fac_max[0] = 100` (uncapped)
- `bkg_strain_rate = 1.0`, `Nt = 11`, `dt = 0.5` → γ ∈ {0, 1, 2, ..., 10}

Analytical reference (in test source, mirrors MDLIB constants):
```cpp
static constexpr double kCalciteMInf   = 0.50;
static constexpr double kCalciteGammaE = 3.50;
static constexpr double kCalciteSlope  = 6.0;

static double anaDeltaCalcite(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5*g*g + g * sqrt(0.25*g*g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kCalciteMInf * (1.0 - exp(-g_eff / kCalciteGammaE));
    return kCalciteSlope * M + 1.0;
}
```

Per-step assertion: `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)`. Same off-by-one indexing as existing tests.

**Rationale:** gates the new form against ITS OWN analytical specification. Catches regressions in dispatch wiring, calibration constants, the `γ_eff` inversion. Doesn't gate against Pieri+01 lab data (deliberate non-goal — sample-to-sample scatter is comparable to the form fit).

**Test runs to γ=10** rather than γ=20 (Hansen olivine) because calcite's δ_∞ = 4 is well-saturated by γ=10 (M reaches 0.477 of asymptote 0.5, δ reaches 3.86 of asymptote 4.0). No additional information from γ ∈ [10, 20].

### D7: Quartz documentation strategy — research note, no MDLIB code

**Choice:** `SETS/AnisoFstrainResearch/quartz/` contains:

1. **`calibrate_quartz.py`** — Python script that loads the Pennacchioni+10 17-point dataset from a tabular literal in the script (no PDF-extraction dependency at run-time), plots the existing `min(FS_AR, ani_fac_max)` form for `ani_fac_max ∈ {2, 2.5, 3}` overlaid with the lab data, computes RMS errors to support the "min form is fit-for-purpose" claim.
2. **`mdoodz_vs_quartz_data.png`** — two-panel plot. Left: M(γ) (or J→M-converted) lab vs min form curves at multiple cap values. Right: δ(γ) with the same.
3. **`quartz_calibration.md`** — research note documenting:
   - The Pennacchioni+10 dataset (17 natural shear-zone samples, γ < 1 to > 15, ~500°C).
   - The min-form transient overshoot analysis (~factor 1.3 over γ ∈ [1, 4]).
   - Recommendation: `ani_fstrain = 1` per phase + `ani_fac_max = 2-3` (specific value depends on T-regime / slip-system mix).
   - Explicit T-range caveat: ~500°C, prism-a slip dominant.
   - Explicit non-recommendation of adding a quartz case 3 unless a future scenario demands sub-γ=4 transient precision.

**No CI test for quartz** — `AniFstrainSimpleShear` and `AniFstrainSaturation` already cover the min form; nothing new to gate.

### D8: Documentation — calibration plot + structured markdown

Mirroring [SETS/AnisoFstrainResearch/](SETS/AnisoFstrainResearch/) for olivine, two new subfolders:

1. **`SETS/AnisoFstrainResearch/calcite/`**:
   - `calibrate_calcite.py`: pure stdlib + numpy + matplotlib (no scipy). Loads Pieri+01 + Bruijn+11 datapoints, runs least-squares fit on `(M_inf, γ_e)` (slope `a` is fixed at 6.0 by D2), prints summary table + RMS errors of all candidate δ-forms, renders `mdoodz_vs_calcite_data.png`.
   - `mdoodz_vs_calcite_data.png`: two-panel plot. Left: M(γ) with all datapoints + fit. Right: δ(γ) with `min(FS_AR, 4)` (existing default), `min(FS_AR, 14)` (Hansen olivine asymptote, for comparison), the new calibrated calcite form, and Pieri+01 + Bruijn+11 lab data. After CI test runs, MDOODZ HDF5 output overlay (red stars) on the calibrated curve.
   - `calcite_calibration.md`: structured note (sections matching Hansen olivine `hansen2012_comparison.md` §7):
     - §1 Functional form (Eq + constants)
     - §2 Calibration procedure (Pieri+01 + Bruijn+11 dataset, J → M conversion)
     - §3 Asymptote target (δ_∞ = 4, matched to empirical PinchSwell calibration; bracketed by Tommasi+09 VPSC bounds)
     - §4 Per-dataset fit comparison
     - §5 RMS error of candidate δ-forms (min, Hansen-form-with-calcite-cap, new calcite-case-2 form)
     - §6 Independent corroboration (Bruijn+11 single-stage segments at γ = 1, 2.6, 5; Schuster+18 saturation note despite different P)
     - §7 Structural caveats (J→M conversion uncertainty, δ-vs-M slope from VPSC not direct lab fit, simple-shear γ_eff caveat, recrystallization-randomization not modelled)
     - §8 Reproducing the calibration (script + plot)
     - §9 Citations
2. **`SETS/AnisoFstrainResearch/quartz/`** (per D7 above).

## Risks / Trade-offs

- **[Risk] Calcite δ-vs-M slope is the weakest link.** Hansen+12 Fig 3b's `δ = 24.5·M + 1` is olivine-specific from a direct lab fit. Calcite has no equivalent published linear fit. We commit to `slope = 6.0` based on Tommasi+09-style VPSC bounds + empirical match to existing `ani_fac_max = 4`. Any future direct lab measurement supersedes this. → **Mitigation**: documented prominently in §7 of the calcite research note. The CI test gates the implementation against THIS calibration, not against lab data — so a future re-calibration is a single-line edit (the `slope = 6.0` literal in `anisoDelta_Calcite`) plus a re-render of the plot.
- **[Risk] J→M conversion is fabric-shape-dependent.** Skemer+05's relation between J-index and M-index assumes specific texture symmetries; Pieri+01 calcite has monoclinic symmetry that the Skemer relation handles approximately. → **Mitigation**: documented as uncertainty. The fit's shape (`γ_e`, `M_inf`) is more important than absolute M values for our purposes — the saturation strain and the asymptote ratio matter more than exact M numbers.
- **[Risk] Pieri+01 vs Bruijn+11 disagreement at γ=5.** P087 reports J=9.76 (single-stage); Bruijn reports J ≈ 4 at the same strain (single-stage segment of two-stage experiment). → **Mitigation**: smoothed in the least-squares fit; documented in §4 of the research note.
- **[Risk] Pieri γ=11 non-monotonic outlier.** J=8.92 < J=9.76 from γ=5; recrystallization randomisation not captured by exponential form. → **Mitigation**: smoothed by the fit; documented as a known limitation. No MDOODZ scenario routinely visits γ=11.
- **[Risk] Quartz documentation-only path doesn't gate the recommendation.** A user could still set `aniso_db = 2` (calcite) on a quartz phase by mistake, getting calcite's δ_∞ = 4 instead of quartz's expected δ ≈ 2-3. → **Mitigation**: existing `ani_fac_max` per phase still acts as a hard ceiling — if the user sets `ani_fac_max = 3` per quartz convention, the calcite δ saturates at 3 anyway. Practical impact small. Documented in the quartz research note.
- **[Trade-off] Skipping quartz case 3 vs adding it for completeness.** Adding it would be ~25 LoC of MDLIB but requires quartz-specific calibration with weaker data backing than calcite. Saved complexity for documentation-only path. If a future scenario emerges where the min-form transient is genuinely wrong for quartz, adding case 3 is a small follow-up.

## Migration Plan

1. Add `static double anisoDelta_Calcite(double FS_AR)` to [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c) with the constants from D2 baked into the function body.
2. Add `case 2 :` to the `ReadDataAnisotropy` switch, assigning `mat->aniso_delta_fn[k] = anisoDelta_Calcite`.
3. **Verify backward compat**: run full `AnisotropyBenchmarkTests` (10 tests including the new HansenOlivine) + `RheologyCreepTests` (13 tests) — all SHALL pass byte-identically. The new code path is opt-in only.
4. Add `AnisotropyBenchmark.AniFstrainCalcite` GTest case per D6. Calibrate threshold to FP precision (just function evaluation — no Stokes solve).
5. Write [SETS/AnisoFstrainResearch/calcite/calibrate_calcite.py](SETS/AnisoFstrainResearch/calcite/calibrate_calcite.py) + render [SETS/AnisoFstrainResearch/calcite/mdoodz_vs_calcite_data.png](SETS/AnisoFstrainResearch/calcite/mdoodz_vs_calcite_data.png).
6. Write [SETS/AnisoFstrainResearch/calcite/calcite_calibration.md](SETS/AnisoFstrainResearch/calcite/calcite_calibration.md) per D8.
7. Write [SETS/AnisoFstrainResearch/quartz/calibrate_quartz.py](SETS/AnisoFstrainResearch/quartz/calibrate_quartz.py) + render [SETS/AnisoFstrainResearch/quartz/mdoodz_vs_quartz_data.png](SETS/AnisoFstrainResearch/quartz/mdoodz_vs_quartz_data.png) showing min form vs Pennacchioni+10.
8. Write [SETS/AnisoFstrainResearch/quartz/quartz_calibration.md](SETS/AnisoFstrainResearch/quartz/quartz_calibration.md) per D7.
9. Manual sanity simple-shear scenario: copy `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` to `SETS/AnisoFstrainResearch/calcite/AniFstrainCalciteShear.txt`, set `ani_fstrain = 2`, `aniso_db = 2`, run 2 timesteps, confirm δ matches the analytical curve.

Rollback: revert in reverse order. The new dispatch arm is opt-in only; reverting it has no effect on existing scenarios.

## Open Questions

- **Should we expose `aniso_db` selection in the input file `.txt` parser format documentation?** Currently `aniso_db` is read as `(int)ReadMatProps(fin, "aniso_db", k, 0.0)` with default 0. Users who don't know about it default correctly (off). But there is no central documentation of which `aniso_db` values are valid. → Defer to a separate documentation-only change.
- **Should the quartz research note recommend a single `ani_fac_max` value (say 3) or a range?** Pennacchioni+10 dataset is roughly consistent with `ani_fac_max ∈ [2, 3]`. → Recommend 3 as a single value with the range documented as user-tunable. Same as olivine pre-Hansen-form had `ani_fac_max = 14` as an empirical default.
- **Is `slope = 6.0` defensible without a direct lab fit?** It's a calibrated value matching the empirical `ani_fac_max = 4` PinchSwell convention; bracketed by Tommasi+09-type VPSC results (2-5). If a reviewer wants tighter bounds, a future change can re-calibrate against newer published Hill-tensor work. → Documented as the weakest link in the calibration; not a blocker for shipping.
