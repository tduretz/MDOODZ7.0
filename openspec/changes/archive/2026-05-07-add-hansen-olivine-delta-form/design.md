## Context

Two facts established in the conversation that motivated this change:

1. **The `min(FS_AR, ani_fac_max)` form (current default, `ani_fstrain = 1`) and the `erfc`-blend form (build-time `-DANISO_ERFC_BLEND` prototype) both fail to fit Hansen+12 olivine data.** RMS error over γ ∈ [0, 5]: `min` = 2.03, `erfc` = 2.70. Both forms grow as `~γ` near the origin while Hansen+12 reports δ growing as `~3.25γ` — a structural ~3× shortfall in low-strain slope that no smoothing of the saturation knee can fix. See [SETS/AnisoFstrainResearch/hansen2012_erfc_comparison.png](SETS/AnisoFstrainResearch/hansen2012_erfc_comparison.png).
2. **The structural cause** is the input variable: both forms use `FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4+1)` (strain-ellipse aspect ratio from polar decomposition of `F`), but Hansen+12 reports δ as a function of the **fabric-strength index M** — a CPO measure not tracked by MDOODZ.

**The fundamental definition of δ** (Hansen+16 Part 1 Eq. 2): for a textured aggregate, the viscous-anisotropy magnitude is the ratio of strain rates measured under stress applied in the **weak** vs **strong** orientation, raised to the inverse stress exponent:

```
δ = (ε̇_e,w / ε̇_e,s) · (σ_e,s / σ_e,w)^n
```

(subscript `e` = von Mises equivalent, `w`/`s` = weak/strong orientation). All published δ formulas (Hansen+12 Fig 3b, Hansen+16 Eq. 3-4) are *derived from this definition*. MDOODZ's scalar `aniso_factor_n` field is this δ — the magnitude of viscous anisotropy at a given local texture state.

The Hansen group has published **three** explicit closed-form expressions for δ derived from this definition:

1. **Hansen+12 Fig 3b — δ-vs-fabric-strength linear fit:**
   ```
   δ = (24.5 ± 5.5) · M + 1
   ```
   M is the Skemer+05 fabric-strength index (0 for random fabric, 1 for single crystal). This is the only Hansen+12 expression that maps a CPO scalar to δ; the δ(γ) curve in Fig 3a is computed from supplementary-info evolution equations and is not given as a closed form.

2. **Hansen+16 Part 1 Eq. 3-4 — stress-exponent-aware steady-state δ:**
   ```
   δ = (F_s / F_w)^n      with F_w = 0.73, F_s = 1.39, n = 4.1
   ```
   `F_w`, `F_s` parameterise the texture's effect on strain rate when stress is in the weak vs strong orientation; `n` is the stress exponent. Evaluates to `δ ≈ 14` for a fully developed olivine texture. This gives the **asymptote** but not the δ(γ) ramp-up.

3. **Hansen+16 Part 2 Eq. 10 — coupled grain-size evolution:**
   ```
   ḋ = ε̇ · (d_ss − d) / ε_c        with ε_c = 1
   ```
   `d_ss` is the steady-state wattmeter grain size, `ε_c` is the critical strain for grain-size evolution. Hansen+16 Part 2 demonstrates that texture and grain size **co-evolve on the same strain timescale `ε_c ≈ 1`**. MDOODZ already implements this form in `LocalIterationViscoElasticGrainSize`; the implication for our change is that no separate textural timescale is needed.

Hansen+14 (EPSL 387, 157) does **not** publish a closed-form δ-formula but provides the **largest tabulated (γ, M) dataset** (Table 1, 26 datapoints unifying samples from Bystricky+00, Zhang & Karato 2005, and Hansen+12a/b/c under consistent M-index methodology). It also defines the **shape parameter K** (Eq. 4) which is informative for fabric-evolution diagnostics but does not enter the δ-formula directly.

Hansen+16 Part 2 also publishes a **micromechanical pseudo-Taylor + director-method model** with calibrated slip-system strengths `τ_0 = {0.30, 0.27, 1.29}` and rotation parameters `f^α = {1.35, 1.45, 0.07, 0.14}` — this is a per-grain integration scheme too heavy for MDOODZ's per-cell rheology pass and is documented here as future-work direction, not as part of this change.

**Independent corroboration of the M(γ) saturation shape**: Hansen+16 Part 2 Fig. 7 row 2 overlays the Hansen+14 26-point torsion dataset with the director-method model curve and shows the model reaches `M ≈ 0.6` by γ ≈ 15 along an exponential-like saturation. **The shape of our 38-point empirical fit `M(γ) = 0.536·(1 − exp(−γ/3.96))` is independently confirmed by their physically-motivated micromechanical model** — we are not pattern-matching the data with an arbitrary curve but with the same shape Hansen+16 Part 2's slip-system rotation kinematics produce. This is the meta-validation: Hansen+16 Eq. 3-4 validates the *asymptote*, Hansen+16 Part 2 Fig. 7 validates the *shape*, and the 38-point combined dataset constrains the *constants*.

This change uses formula (1) — Hansen+12 Fig 3b's linear `δ = 24.5·M + 1` — as the **δ(M) mapping**, with M evaluated as `M(γ_eff(FS_AR))` via a calibrated exponential-saturation. Formula (2) — Hansen+16 Part 1 Eq. 3-4 — is used as the **independent asymptote target** (`δ_∞ ≈ 14`) that validates the fit's limiting behaviour. Formula (3) — Hansen+16 Part 2 Eq. 10 — is already in MDOODZ via the wattmeter; we cite it as supporting evidence that scalar δ(γ) is the right level of approximation given the existing grain-size coupling.

## Goals / Non-Goals

**Goals:**

- Add a third δ-saturation form (`ani_fstrain = 2`) that takes Hansen+12's `δ = 24.5·M + 1` literally, with M computed from MDOODZ's existing `FS_AR` field via a calibrated mapping `M(FS_AR)`.
- The new form SHALL match the two non-trivial Hansen+12 datapoints (γ=4, γ=10) **exactly** by construction, and the noise-floor γ=0.03 point to ~1% (with the residual recognised as M-index sampling noise on a near-random fabric, not real fabric strength).
- **Backward compatibility (HARD requirement)**: existing `ani_fstrain ∈ {0, 1}` byte-identical post-change; all existing CI tests and live `SETS/` scenarios unaffected. The new form is opt-in only via the new dispatch value.
- All calibration documentation, comparison plot, and the calibration script live in [SETS/AnisoFstrainResearch/](SETS/AnisoFstrainResearch/) — not in CI test fixtures or MDLIB.
- One new CI test gating the analytical reference at the calibration datapoints.

**Non-Goals:**

- **No CPO model.** M is not tracked dynamically. The new form approximates M from FS_AR via a calibrated empirical relation. The structural simplification (M depends on grain orientations, FS_AR on bulk strain — they are different quantities) is acknowledged and called out.
- **No per-mineral calibration.** Constants are hardcoded for olivine (Hansen+12). Calcite (Barnhoorn+04, Pieri+01) and quartz/mica use cases would need separate calibrations with different constants — out of scope for this change.
- **No quantitative CI gating against Hansen+12 lab data.** Same reasoning as `gate-finite-strain-anisotropy`: digitisation noise in lab data is comparable to model approximations; using lab data as a CI threshold would conflate "MDOODZ has a bug" with "the model is approximate." The new CI test gates against the **analytical Hansen+12 fit at the calibration points**, not against the underlying lab measurements.
- **No removal of the existing `-DANISO_ERFC_BLEND` build-time toggle.** Decision deferred to D6 below.
- **No phase-3 unified `η(d, δ)` constitutive coupling.** Still future work.
- **No user-tunable calibration constants from the input file.** Per-mineral constants (`M_inf`, `γ_e`, slope) are baked into static helper functions in `FlowLaws.c`, not exposed as `ReadMatProps` parameters. Users select calibrations via `aniso_db = 1, 2, ...` (an enumeration of vetted calibrations), not by typing constants. This matches the FlowLaws idiom — `pwlv = 40` selects "Hirth & Kohlstedt 2003 dry olivine dislocation creep" with all its constants, not free numbers.

## Decisions

### D1: Composition strategy — three-step `δ(FS_AR) = 24.5·M(γ_eff(FS_AR)) + 1`

**Choice:** Build the new δ-form as a composition:

```
γ_eff(FS_AR) = sqrt(FS_AR) - 1/sqrt(FS_AR)        // step 1: invert simple-shear FS_AR(γ)
M(γ_eff)     = M_inf · (1 - exp(-γ_eff / γ_e))    // step 2: Hansen+12 fabric-evolution shape
δ(M)         = 24.5 · M + 1                       // step 3: Hansen+12 Fig 3b linear fit (verbatim)
```

Step 1 is the analytical inversion of the simple-shear relation `FS_AR = 1 + γ²/2 + γ·√(γ²/4+1)` — derivation: `FS_AR = (s² + s·√(s²+4))/2` where `s = γ`, giving `γ = √FS_AR - 1/√FS_AR`. Verifies trivially: `γ=0 → FS_AR=1 → γ_eff=0`; `γ=4 → FS_AR=17.94 → γ_eff=4`; `γ=10 → FS_AR=102 → γ_eff=10`.

**Rationale for splitting the function this way:**

- `δ = 24.5·M + 1` is taken **verbatim from Hansen+12 Fig 3b** — no fitting, no approximation. This is the published calibration.
- `M(γ_eff) = M_inf · (1 - exp(-γ_eff/γ_e))` is a 2-parameter sigmoid fit to Hansen+12's three (γ, M) datapoints. The functional form is the standard "saturating exponential" used in fabric-evolution models (Tommasi+2000-style, where γ_e is a characteristic strain for fabric development).
- `γ_eff(FS_AR)` makes the form work for any deformation tracked by `F`. Under simple shear, `γ_eff` is the actual γ. Under pure shear, `γ_eff` is a deterministic function of `FS_AR` that equals γ under simple shear — interpretation outside simple shear is questionable (see Risk below) but the formula remains well-defined.

**Alternatives considered:**

- *(A) Direct fit `δ = a·tanh(b·γ) + 1`* (i.e. the form I used in `hansen2012_comparison.md`'s original table): valid only for simple shear (γ has no meaning under pure shear). Rejected.
- *(B) Direct fit `M = f(FS_AR)` skipping the γ_eff intermediate*: tried `M = M_inf·(1 - FS_AR^(-c))` — doesn't fit cleanly, parameters from γ=4 vs γ=10 datapoints are inconsistent (verified during calibration). The exponential-in-γ form fits better because that's what Hansen+12's lab fabric evolution actually does in their data.
- *(C) Track M as a per-particle state variable evolved by an ODE*: real CPO modeling. Out of scope.

### D2: Calibration constants

**Choice:** Hardcoded olivine constants, calibrated by least-squares fit to the **combined Hansen-group olivine torsion dataset (38 datapoints)**:

- **Hansen+14** *EPSL* 387, 157 Table 1: 26 datapoints over γ ∈ [0, 18.7]. Hansen+14 reanalyzed previously published samples from Bystricky+00, Zhang & Karato 2005, and Hansen+12a/b/c with consistent M-index methodology (Skemer+05), so this is one unified dataset.
- **Hansen+16 Part 1** *EPSL* 445, 92 Fig. 6 + Table 1: 12 datapoints from 5 NEW samples (PT0716, PT0718, PT0751, PT0752, PT0754) not in Hansen+14.

**Mineral-specific math lives entirely in `MDLIB/FlowLaws.c`** as static helper functions, dispatched via a per-phase function pointer. This matches the existing `ReadDataPowerLaw` / `ReadDataGSE` pattern (switch on a per-phase database integer, log the citation, populate per-phase data) **but goes further** — instead of populating per-phase scalar fields and assuming a fixed functional form upstream, each case assigns a per-phase function pointer that encapsulates whatever closed form that mineral needs. `AnisotropyRoutines.c` ends up mineral-agnostic.

Two new fields on `mat_prop`:

```c
int    aniso_db[20];                              // database selector: 0=off, 1=Hansen olivine, ...
double (*aniso_delta_fn[20])(double FS_AR);       // per-phase mineral-specific δ(FS_AR) closed form
```

`ReadDataAnisotropy` currently has one entry:

```c
// File-scope static helper (constants baked into the function body):
static double anisoDelta_HansenOlivine( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.536;       // calibrated to 38-point Hansen-group dataset
    const double gamma_e   = 3.96;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return 24.5 * M + 1.0;                // Hansen+12 Fig 3b linear coefficient
}

// In ReadDataAnisotropy switch:
case 1 :
    LOG_INFO("Olivine viscous anisotropy - Hansen+12/+14/+16 (combined 38-point fit):");
    mat->aniso_delta_fn[k] = anisoDelta_HansenOlivine;
```

Future minerals (calcite, quartz, …) get their own static helpers (with whatever functional form they need — saturating exponential, slip-system-aware lookup, tanh, piecewise) and a new `case N` that assigns them. **Adding a new mineral requires zero changes to `AnisotropyRoutines.c`.**

`AnisotropyRoutines.c`'s `AnisoFactorEvolv` becomes pure dispatch:

```c
double AnisoFactorEvolv( double FS_AR, double aniso_fac_max, int ani_fstrain,
                         double (*aniso_delta_fn)(double FS_AR) ) {
  if (ani_fstrain == 2) {
    if (aniso_delta_fn == NULL) { /* exit with helpful error */ }
    return MINV(aniso_delta_fn(FS_AR), aniso_fac_max);
  }
  if (ani_fstrain == 3) return erfc_blend_delta(FS_AR, aniso_fac_max);
  return MINV(FS_AR, aniso_fac_max);
}
```

No Hansen-specific constants, no `M_inf` / `γ_e` literals, no olivine math anywhere in `AnisotropyRoutines.c`. The cap (`MINV(δ, aniso_fac_max)`) is applied uniformly by the dispatcher so mineral-specific functions need only return the uncapped δ.

Selection: per-phase input parameter `aniso_db` (default 0; auto-set to 1 when `ani_fstrain == 2` is requested without explicit selection — both at parse time in `InputOutput.c` and post-`MutateInput` in `Main_DOODZ.c`, so test harnesses that flip `ani_fstrain` after parse still get the calibration loaded).

**Asymptote — independent two-method validation**: the linear-formula extrapolation gives `24.5·0.536 + 1 = 14.13`. Hansen+16 Part 1 Eq. 3-4 reports a stress-exponent-aware direct calculation `δ = (F_s/F_w)^n = (1.39/0.73)^4.1 = 14.02`. **Agreement: 100.8%.** Two completely independent methods (statistical fit to 38 (γ, M) datapoints + Hansen+12 linear coefficient, vs ratio-of-Hill-parameters^stress-exponent from Hansen+16 Part 1) converge on `δ_∞ ≈ 14`. Capped at user-set `ani_fac_max` for safety (D5).

**Per-dataset fit comparison** (from `calibrate_hansen_olivine.py`):

| Dataset | N | M_inf | γ_e | δ_∞ (24.5·M+1) | RMS(M) |
|---|---|---|---|---|---|
| Hansen+12 Nature | 3 | 0.620 | 4.59 | 16.18 | 0.015 |
| Hansen+14 EPSL | 26 | 0.568 | 4.71 | 14.91 | 0.075 |
| Hansen+16 EPSL Part 1 | 12 | 0.500 | 3.00 | 13.24 | 0.071 |
| **COMBINED (NEW)** | **38** | **0.536** | **3.96** | **14.13** | **0.076** |

The combined fit is the most strongly constrained (largest N) and lands almost exactly on Hansen+16 Eq. 3-4's stress-aware target. Single-paper fits scatter ±5% around the combined value, reflecting modest sample-to-sample variation between independent sample populations.

**RMS of candidate δ-forms vs the 38-point combined dataset** (cap = 14.02, mapping `δ_lab = 24.5·M_obs + 1`):

| Form | RMS |
|---|---|
| `δ = min(FS_AR, 14)` (MDOODZ default) | 2.87 |
| `δ = erfc-blend(FS_AR, 14)` | 3.11 |
| Hansen+12 3-point fit | 2.07 |
| Hansen+16 12-point fit | 1.93 |
| **Combined 38-point fit (NEW)** | **1.85** |

The combined fit improves RMS by ~10% over the 12-point fit and ~36% over the existing `min` form.

**Two structural notes**:

1. **High-γ scatter is sample variability, not fit failure.** Hansen+14 Fig. 6 explicitly notes "fabric strength is roughly correlated with grain size and inversely correlated with shear stress" — explaining the ~0.2 spread in M at γ ≈ 8–10. Hansen+16 Part 2 §4.1.2 attributes this to recrystallization and grain-size feedback (the same `ε_c ≈ 1` co-evolution that drives MDOODZ's wattmeter — Hansen+16 Part 2 Eq. 10 `ḋ = ε̇·(d_ss − d)/ε_c` is the same form MDOODZ already uses for grain size). The exponential fit's smooth saturation is the right ensemble approximation.

2. **Heterogeneity within the combined dataset**: Hansen+16 Part 1's PT0718 radial profile reports `M ≈ 0.11` at γ ∈ [0.3, 1.2] because that sample retained residual fabric from a **prior extension experiment**. For fresh olivine in MDOODZ scenarios (`F = I` at t=0 so `FS_AR = 1`), `M(γ=0) = 0` is the right anchor. The exponential form forces this; Hansen+14's PT0535 (γ=0, M=0.01) confirms it directly.

### D3: γ_eff outside simple shear — caveat in code comment, not a runtime guard

**Choice:** the relation `γ_eff(FS_AR) = √FS_AR - 1/√FS_AR` is valid as a strain-magnitude proxy ONLY under simple shear. Under pure shear (`FS_AR(t) = exp(2·ε̇·t)`), the formula gives `γ_eff = exp(ε̇·t) - exp(-ε̇·t) = 2·sinh(ε̇·t)` — a deterministic function of FS_AR but not the actual shear strain (there isn't one in pure shear).

**Document but don't guard.** A runtime check on deformation type would require introspecting the velocity field, which is non-local information not available to `AnisoFactorEvolv`. Instead, add a comment in the code explaining the simple-shear assumption and note that under pure shear or general 2D flow, the form gives a defensible but interpretation-questionable saturation curve. Users running pure-shear scenarios with `ani_fstrain = 2` are signing up for "Hansen+12-shaped saturation as a function of accumulated stretch," not "Hansen+12 olivine fabric strength."

**Rationale:** the alternative — disabling `ani_fstrain = 2` outside simple shear — would force users to know which scenarios are which, breaking the per-phase opt-in cleanliness. Better to document and let users decide.

### D4: Backward compatibility — additive `ani_fstrain = 2` dispatch

**Choice:**

```c
double AnisoFactorEvolv(double FS_AR, double aniso_fac_max, int ani_fstrain) {
    if (ani_fstrain == 1) {
        return MINV(FS_AR, aniso_fac_max);    // EXISTING — byte-identical preserved
    }
    if (ani_fstrain == 2) {
        return hansen2012_delta(FS_AR, aniso_fac_max);   // NEW
    }
    return aniso_fac_max;   // ani_fstrain == 0 fallback (matches existing behaviour
                            //   when AnisoFactorEvolv isn't called)
}
```

**Migration**: the existing `AnisoFactorEvolv(FS_AR, aniso_fac_max)` signature must be extended with the `ani_fstrain` integer. All call sites in [AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) currently pass `(FS_AR, materials->ani_fac_max[p])` — they need updating to also pass `materials->ani_fstrain[p]`. Six call sites (3 in centroid loop lines 382-393, 3 in vertex loop). All under `if (materials->ani_fstrain[p] == 1)` guards, so the new form ALSO requires extending those guards to `if (materials->ani_fstrain[p] != 0)` — and the dispatch happens inside `AnisoFactorEvolv`.

**Rationale:** putting the dispatch inside `AnisoFactorEvolv` (rather than in the per-cell call sites) keeps the cell-loop logic clean: same call, behaviour controlled by the per-phase setting. This also means future phases (smoother forms, per-mineral calibrations) only edit `AnisoFactorEvolv`.

**Existing `ani_fstrain = 1` byte-identical**: yes, by direct construction — the `if (ani_fstrain == 1)` branch is the unmodified `MINV(FS_AR, aniso_fac_max)` line.

### D5: Saturation cap

**Choice:** Final returned δ is `min(hansen2012_delta(FS_AR), aniso_fac_max)`. If a user sets `ani_fac_max < 16.22` (the Hansen+12 asymptote), the cap takes precedence — preserves the meaning of `ani_fac_max` as a user-set ceiling.

**Rationale:** Without the cap, scenarios that reused our `ani_fac_max = 4` value (calcite Barnhoorn+04 calibration) and switched to `ani_fstrain = 2` would silently get olivine-asymptote saturation (~16) instead of the user's calcite asymptote (4). The cap preserves the user's intent. Documented in the code: "ani_fac_max acts as a hard ceiling regardless of dispatch form."

### D6: Erfc-blend form (`-DANISO_ERFC_BLEND`) — keep as `ani_fstrain = 3` for completeness

**Choice:** Promote the build-time toggle to a runtime dispatch arm `ani_fstrain = 3`. Keeps the work we put into prototyping it; allows comparison runs without rebuilding; removes the unusual `#ifdef` from the source.

**Rationale:** the erfc form is empirically inferior (33% worse RMS than `min`), so we wouldn't recommend it. But it's still useful for the comparison runs we already have GIFs for, and removing it would invalidate those runs' reproducibility. Cheaper to keep than delete.

**Migration:** The `#ifdef ANISO_ERFC_BLEND` block becomes `if (ani_fstrain == 3) { ... }`. The `cmake-build-erfc/` build directory becomes unnecessary — single mdoodz build supports all three dispatch values. The existing `cmake-exec/PinchSwellGSEAniso/PinchSwellGSEAniso_erfc` binary stays for backward compat with the saved overnight results, but new runs use the standard binary with `ani_fstrain = 3`.

### D7: Test strategy — analytical-reference gate at calibration points

**Choice:** New GTest case `AnisotropyBenchmark.AniFstrainHansenOlivine` reusing the [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) base fixture via `MutateInput` (sets `ani_fstrain = 2`, `ani_fac_max = 100` to disable the cap, simple-shear BCs).

The analytical reference is:
```
δ_ana(t) = 24.5 · M_inf · (1 - exp(-γ_eff(FS_AR(t)) / γ_e)) + 1
        with M_inf  = 0.536
             γ_e    = 3.96
             FS_AR(t) = 1 + γ²(t)/2 + γ(t)·√(γ²(t)/4+1)
             γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
             γ(t)   = γ_dot · t
```

Per-step assertions: `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` with the off-by-one output indexing convention from previous tests (`Output<k>` reflects F at end of step k-1).

**Rationale:** this gates the new form against ITS OWN analytical specification — the formula evaluated symbolically vs evaluated through the MDOODZ rheology pass. Catches regressions in the dispatch wiring, the calibration constants, and the `γ_eff` inversion. Doesn't gate against Hansen+16 lab data (deliberate non-goal — the lab data has 0.07 sample-scatter in M, much larger than CI threshold).

### D8: Documentation — calibration plot + updated `hansen2012_comparison.md`

**Choice:** Three artifacts in `SETS/AnisoFstrainResearch/` (the script and plot already exist as of this session):

1. **`calibrate_hansen_olivine.py`** — the calibration script. Reproducible: pure stdlib + numpy + matplotlib, no scipy (uses 2D coarse-then-fine grid search to minimise SSE on `M_inf`, `γ_e`). Loads the **combined 38-point Hansen-group dataset** (26 from Hansen+14 Table 1 + 12 from Hansen+16 Part 1 Fig. 6 + Table 1) plus the 3 Hansen+12 Nature points for comparison. Outputs the per-dataset fit summary table, the RMS comparison of all candidate δ-forms vs the combined data, and renders `hansen_olivine_calibration.png`.
2. **`hansen_olivine_calibration.png`** — two-panel output. Panel (a): M(γ) showing all 38 + 3 datapoints with the combined-fit curve in solid black plus per-dataset fits as light dashed reference curves. Panel (b): five candidate δ-forms (`min`, `erfc`, Hansen+12 3-pt, Hansen+16 12-pt, combined 38-pt) overlaid with the combined lab data and the Hansen+16 Eq. 3-4 asymptote marker at δ ≈ 14.0.
3. **Update `hansen2012_comparison.md`** — new section §7 documenting:
   - the calibration procedure (38-point combined least-squares fit),
   - the source breakdown (26 Hansen+14 Table 1 + 12 Hansen+16 Part 1 Fig. 6 + Table 1),
   - the resulting constants `M_inf = 0.536`, `γ_e = 3.96`,
   - the **two-method asymptote validation**: linear-formula extrapolation 14.13 vs Hansen+16 Eq. 3-4 stress-aware direct calc 14.02 → 100.8% agreement,
   - the embedded comparison plot,
   - the RMS comparison (`min` = 2.87, `erfc` = 3.11, Hansen+12 3-pt = 2.07, Hansen+16 12-pt = 1.93, **combined 38-pt = 1.85**),
   - the FS_AR-vs-M structural caveat,
   - the Hansen+16 Part 2 stress-orientation caveat (real δ depends on alignment between principal stress and texture's strong axis; MDOODZ's δ is a scalar magnitude / texture-aligned upper bound),
   - the Hansen+16 Part 2 ε_c ≈ 1 finding (texture and grain size co-evolve on the same strain timescale, consistent with MDOODZ's wattmeter coupling).

The `min()` vs Hansen+12 mismatch documentation in the existing markdown stays (sections 1–6) for historical context.

## Risks / Trade-offs

- **[Risk] `γ_eff(FS_AR)` outside simple shear is not a real strain.** Under pure shear, the formula gives `γ_eff = 2·sinh(ε̇·t)` — defensible but not directly the accumulated strain. → **Mitigation**: documented in code comment and in `hansen2012_comparison.md`. Users opting into `ani_fstrain = 2` are choosing "Hansen+12-shaped saturation parameterised by accumulated stretch", not "Hansen+12 olivine physics under arbitrary deformation."
- **[Risk] Calibration data has substantial sample-to-sample scatter at high γ.** Across the 38-point combined dataset, M at γ ≈ 8–10 spans roughly 0.35–0.69 (~2× spread) — Hansen+14 Fig. 6 attributes this to grain-size and shear-stress correlations, and Hansen+16 Part 2 §4.1.2 attributes it to recrystallization and grain-size feedback. → **Mitigation**: the 38-point least-squares fit smooths through ensemble scatter; the resulting form is consistent with Hansen+16 Part 2's director-method texture-evolution model curves (Fig. 7), confirming the ensemble shape is right. The asymptote agreement between linear-formula extrapolation (14.13) and Hansen+16 Eq. 3-4 stress-aware direct calc (14.02) at 100.8% is independent two-method validation that the asymptote level is correctly identified.
- **[Risk] Stress-orientation dependence of viscous anisotropy is invisible.** Hansen+16 Part 2 Fig. 10 shows real δ depends on alignment between principal stress and texture's strong axis, with a ~10× spread relative to the isotropic case. MDOODZ's scalar δ corresponds to the worst-case (texture-aligned) value. → **Mitigation**: documented in `hansen2012_comparison.md` §7 as a structural limitation. A future change could replace scalar δ with a Hill-tensor parameterisation; out of scope here.
- **[Risk] Hansen+12 is olivine.** Using `ani_fstrain = 2` for calcite/quartz scenarios would be physically wrong (different constants apply). → **Mitigation**: documented prominently in code comments and the research note. Per-mineral calibration is explicit future work.
- **[Risk] Existing scenarios that opt into `ani_fstrain = 2` post-change get olivine-asymptote δ even if they had calcite-tuned `ani_fac_max = 4`.** → **Mitigation**: D5 cap. `ani_fac_max` overrides the formula's asymptote, so `ani_fstrain = 2` with `ani_fac_max = 4` saturates at 4 (calcite-like) but with Hansen+12-shape transition. Not perfectly calibrated for calcite, but at least respecting the user's ceiling.
- **[Trade-off] Function-pointer dispatch vs scalar-field dispatch.** Function pointers add a level of indirection but make the framework genuinely mineral-agnostic — different minerals can use different functional forms without `AnisotropyRoutines.c` knowing. Scalar-field dispatch (e.g. populating `aniso_M_inf`/`aniso_gamma_e` per phase) would lock all minerals into the saturating-exponential form, which doesn't fit calcite (CPO is only 1/3 of weakening) or quartz (slip-system-dependent fabric type). Going with function pointers as the long-term architecture; constants are encapsulated inside each `static double anisoDelta_<Mineral>(double FS_AR)` helper in `FlowLaws.c`.

## Migration Plan

1. Update [`AnisoFactorEvolv`](MDLIB/AnisotropyRoutines.c#L46) signature: add `int ani_fstrain` parameter. Implement the three-arm dispatch (D4 + D6).
2. Implement `hansen2012_delta(FS_AR, aniso_fac_max)` helper using D1's formula and D2's constants.
3. Update all 6 call sites (3 centroid lines 382-393, 3 vertex lines 432-443) to pass `materials->ani_fstrain[p]`.
4. **Verify backward compat**: run full `AnisotropyBenchmarkTests` (9 tests) + `RheologyCreepTests` (13 tests) — all should be byte-identical (`ani_fstrain ∈ {0, 1}` paths unchanged).
5. Add `AnisotropyBenchmark.AniFstrainHansenOlivine` GTest case per D7. Calibrate threshold to FP precision (it's just a function evaluation — no Stokes solve required).
6. ~~Write~~ Already done: [SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py](SETS/AnisoFstrainResearch/calibrate_hansen_olivine.py) + [hansen_olivine_calibration.png](SETS/AnisoFstrainResearch/hansen_olivine_calibration.png) produced this session.
7. Update `SETS/AnisoFstrainResearch/hansen2012_comparison.md` with §7 per D8.
8. Run a quick sanity simple-shear scenario with `ani_fstrain = 2` and verify δ-evolution matches the analytical curve.
9. Manual regression smoke: `SETS/RiftingComprehensive.txt` (uses `ani_fstrain = 1`) for ≥ 2 timesteps — must be byte-identical.

Rollback: revert in reverse order. The new dispatch arm is opt-in only; reverting it has no effect on existing scenarios.

## Open Questions

- **Should `M_inf` and `γ_e` be exposed as user parameters per-phase?** Defer to demand. If only olivine scenarios use this, hardcoded is fine. If calcite users want it, they'd need their own calibration anyway and exposing parameters is the obvious extension.
- **How to render the calibration comparison plot** — same style as `hansen2012_erfc_comparison.png` (now showing 4 curves: `min`, `erfc`, `hansen2012`, ref). Format detail; defer to implementation.
- **Cap behaviour for `δ_max < user's ani_fac_max`** — currently capped at user value, but if user sets `ani_fac_max = 100` (the test config), the formula's natural asymptote of 16.22 prevails. That's fine for the test but worth documenting as the expected behaviour.
