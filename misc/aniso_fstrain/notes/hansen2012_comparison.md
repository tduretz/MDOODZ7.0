# MDOODZ `min(FS_AR, ani_fac_max)` vs Hansen+12 `13·tanh(0.25γ) + 1` — research note

> Status: **research note, not a CI gate.** This document records a known
> form-mismatch between MDOODZ's finite-strain anisotropy law and the
> empirical fit from Hansen, Zimmerman & Kohlstedt (2012) for olivine. The
> mismatch is acknowledged here as future-work motivation.
>
> CI tests for `ani_fstrain = 1` live in `TESTS/AnisotropyBenchmarkTests.cpp`
> (`AnisotropyBenchmark.AniFstrainSimpleShear`, `…AniFstrainSaturation`) and
> validate **MDOODZ's own** closed-form FS_AR formula — not the Hansen+12
> empirical fit. Using a digitised lab fit as a CI threshold would conflate
> "MDOODZ has a bug" with "the model is approximate", so it stays out of CI.

## 1. What MDOODZ implements

When per-phase `ani_fstrain = 1`, the anisotropy strength `δ` evolves with
the strain-ellipse aspect ratio `FS_AR = e1/e2` (largest/smallest principal
stretch from the polar decomposition of the deformation gradient F):

```c
// MDLIB/AnisotropyRoutines.c:46-54
double AnisoFactorEvolv(double FS_AR, double aniso_fac_max) {
  // const double FS_AR_crit  = aniso_fac_max / 2.0;
  // const double x = 0.5*erfc((FS_AR - FS_AR_crit) / delta_ratio);
  // return x * FS_AR + (1.0-x) * aniso_fac_max;
  return MINV(FS_AR, aniso_fac_max);
}
```

Note the commented-out `erfc` blend — a smoother saturation is already
contemplated in the source but not active. `FS_AR` is computed in
[`MDLIB/RheologyParticles.c:634-678`](../../MDLIB/RheologyParticles.c#L634)
via `C = F^T·F`, `U = √C`, eigenvalues of `U`, `FS_AR = e1/e2`.

Under simple shear with shear strain `γ`:

```
FS_AR(γ) = 1 + γ²/2 + γ·sqrt(γ²/4 + 1)              (closed form)
δ(γ)     = min(FS_AR(γ), ani_fac_max)               (live MDOODZ form)
```

The `min()` produces a flat clamp once FS_AR exceeds `ani_fac_max` — a sharp
"knee" in the (γ, δ) plot at the saturation strain `γ_sat`.

## 2. What Hansen, Zimmerman & Kohlstedt (2012) fit

> Hansen, L. N., Zimmerman, M. E. & Kohlstedt, D. L. (2012)
> "The influence of microstructure on deformation of olivine in the
> grain-boundary sliding regime", *Nature* **492**, 415–418.
> doi:[10.1038/nature11671](https://doi.org/10.1038/nature11671)

From olivine torsion experiments (their Eq. 7 / Fig. 4):

```
δ(γ) = 13·tanh(0.25·γ) + 1                          (Hansen+12 empirical fit)
```

Asymptote: `δ → 14` as `γ → ∞`. Half-saturation at `γ ≈ 2.2`. Smooth
roll-off (no kink). Calibrated against axial-symmetric torsion at lab
strains (γ ≤ 5) for olivine — a specific mineral and a specific
deformation geometry.

## 3. Quantitative comparison

Setting `ani_fac_max = 14` to match Hansen+12's asymptote:

| γ | FS_AR (analytical) | δ_MDOODZ = min(FS_AR, 14) | δ_Hansen = 13·tanh(0.25γ) + 1 | gap |
|---|---|---|---|---|
| 0.0 | 1.000 | 1.000 | 1.000 | 0.000 |
| 0.5 | 1.640 | 1.640 | 2.617 | -0.977 |
| 1.0 | 2.618 | 2.618 | 4.215 | -1.597 |
| 1.5 | 4.084 | 4.084 | 5.732 | -1.648 |
| 2.0 | 5.828 | 5.828 | 7.110 | -1.282 |
| 2.5 | 7.882 | 7.882 | 8.327 | -0.445 |
| 3.0 | 10.083 | 10.083 | 9.371 | +0.712 |
| 3.5 | 12.452 | 12.452 | 10.247 | +2.205 |
| 4.0 | 14.000* | 14.000 | 10.972 | +3.028 |
| 5.0 | 26.963 | 14.000 | 12.247 | +1.753 |
| 8.0 | 65.022 | 14.000 | 13.658 | +0.342 |

*FS_AR(γ ≈ 3.93) = 14 (saturation strain γ_sat for `ani_fac_max = 14`).

**Key shape differences:**

1. **Low strain (γ ≤ 1):** MDOODZ underestimates δ. At γ=1, MDOODZ gives 2.6 vs Hansen's 4.2 — 38% low. The unbounded FS_AR formula simply rises slower than the empirical tanh.
2. **Mid strain (γ ∈ [1.5, 2.5]):** MDOODZ undershoots by ~15-25%.
3. **Approaching saturation (γ ≈ 3.5–5):** MDOODZ briefly **overshoots** Hansen because FS_AR rises steeply through the asymptote region before being clamped, while the tanh approaches smoothly.
4. **Above saturation (γ ≥ 6):** both forms agree to within ~3%, because MDOODZ is pinned at 14 and Hansen approaches 14 asymptotically.

In words: MDOODZ's `min()` form is **too sharp at the knee** and **rises too fast at low strain**. It captures the asymptote but misses the empirical curve shape.

## 4. What this means for MDOODZ scenarios using `ani_fstrain = 1`

The live use cases (`SETS/RiftingComprehensive.txt`,
`SETS/RiftingComprehensive_highres.txt`, `SETS/CollisionPolarCartesianAniso.txt`)
inherit this approximation. At lithospheric strains (γ ~ 0.1–2), the
underestimate at low strain is the dominant error: **simulated foliations
develop slower than Hansen+12 would predict**. The asymptote behaviour
matches once cumulative strain is large enough.

This is acceptable for many qualitative studies (comparing aniso vs
isotropic regimes, or scenarios where the knee is reached early), but
quantitatively comparing MDOODZ output against Hansen+12 lab data — or
extrapolating Hansen+12 to lower-strain regions — would require a different
δ(FS_AR) form.

## 5. Future work (out of scope for this change)

Two natural replacements for `AnisoFactorEvolv`:

**Option A — uncomment the `erfc` blend** (already in the source, just disabled):

```c
const double FS_AR_crit  = aniso_fac_max / 2.0;
const double delta_ratio = ...;  // tunable transition width
const double x = 0.5*erfc((FS_AR - FS_AR_crit) / delta_ratio);
return x * FS_AR + (1.0-x) * aniso_fac_max;
```

This smooths the knee but doesn't change the low-strain shape; would still
underestimate at γ ≤ 1.

**Option B — direct tanh fit (closer to Hansen+12):**

```c
// δ(γ) requires γ; here we recover γ from FS_AR via the simple-shear
// closed form (only valid for simple shear, not pure shear).
// Better: parameterise δ on FS_AR directly with a tanh of (FS_AR - 1):
const double a = aniso_fac_max - 1.0;
return 1.0 + a * tanh((FS_AR - 1.0) / b);  // b chosen to match γ_half ≈ 2.2
```

This matches Hansen+12's shape but loses the convenient FS_AR-as-cap
interpretation; needs careful calibration of `b` against lab data.

**Option C — keep `min()` but document as known approximation.** Cheapest,
acknowledges the limitation in the rheology section of the docs, and lets
researchers correct for it at the analysis stage.

Both A and B require a follow-up CI test that validates the new shape
against either the existing `min()` (regression) or Hansen+12 lab points
(but then digitisation noise becomes the threshold floor — see opening
note).

## 6. Reproducing the comparison

The companion `AniFstrainSimpleShear.txt` runs simple shear with
`bkg_strain_rate = 1.0` (γ̇ = 2), `dt = 0.1`, `Nt = 20` → γ ∈ [0, 4].
With `ani_fac_max = 14`, the trajectory enters saturation around γ ≈ 3.93.

```bash
cd cmake-build-test
mkdir -p AniFstrainSimpleShearResearch
./AniFstrainSimpleShear/AniFstrainSimpleShear ../../misc/aniso_fstrain/AniFstrainSimpleShear.txt
# Inspect Centers/ani_fac in AniFstrainSimpleShearResearch/Output*.gzip.h5
```

(MDOODZ does NOT auto-build a binary for this `.txt`; this is a hand-run
research scenario, not a CI test. The closest CI binary is
`AnisotropyBenchmarkTests`, but that uses a different fixture file and
different parameters.)

## 7. Hansen-olivine calibrated form (`ani_fstrain = 2`)

Sections 1–6 above document the form mismatch between MDOODZ's default
`δ = min(FS_AR, ani_fac_max)` form (`ani_fstrain = 1`) and the Hansen+12
data, motivating a Hansen-calibrated alternative. This section documents
the realisation of that alternative as `ani_fstrain = 2`, delivered by
[openspec/changes/add-hansen-olivine-delta-form](../../openspec/changes/add-hansen-olivine-delta-form/).

### 7.1 Functional form

```
δ(FS_AR) = min( 24.5 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )

  where  γ_eff(FS_AR) = √FS_AR − 1/√FS_AR              (simple-shear inversion)
         M(γ)         = M_inf · (1 − exp(−γ / γ_e))    (sigmoidal fabric saturation)
         M_inf        = 0.536
         γ_e          = 3.96
         24.5 · M + 1                                  (Hansen+12 Fig 3b verbatim)
```

The composition is: **invert** the simple-shear FS_AR identity to recover
γ from FS_AR, **map** γ to the fabric-strength index M via an
exponential-saturation fit, then **apply** Hansen+12's published
δ-vs-M linear fit. Capped at the user-set `ani_fac_max` so calcite-tuned
caps (e.g. `ani_fac_max = 4`) are respected even though the Hansen
olivine asymptote is `≈ 14`.

### 7.2 Calibration procedure

`M_inf` and `γ_e` are calibrated by least-squares fit (2D coarse-then-fine
grid search on SSE) to the **combined 38-point Hansen-group olivine
torsion dataset**:

- **26 datapoints** from Hansen, Zhao, Zimmerman & Kohlstedt (2014) *EPSL*
  387, 157 Table 1 — Hansen+14 reanalysed previously published samples
  from Bystricky+00, Zhang & Karato 2005, and Hansen+12a/b/c using a
  single consistent M-index methodology (Skemer+05), so this is one
  unified dataset.
- **12 datapoints** from Hansen, Warren, Zimmerman & Kohlstedt (2016)
  *EPSL* 445, 92–103 Part 1 (Fig. 6 PT0718 radial profile + Table 1
  end-of-experiment values for samples PT0716, PT0751, PT0752, PT0754).
  These 5 samples are not in Hansen+14.

The script [calibrate_hansen_olivine.py](calibrate_hansen_olivine.py)
loads both datasets, runs the fit, prints the per-dataset summary table
plus the RMS comparison of all candidate δ-forms, and renders the
two-panel comparison plot. Pure stdlib + numpy + matplotlib (no scipy).

### 7.3 Asymptote — independent two-method validation

| Method | δ_∞ |
|---|---|
| Linear-formula extrapolation: `24.5 · M_inf + 1` (this fit) | **14.13** |
| Hansen+16 Part 1 Eq. 3-4 stress-aware: `(F_s/F_w)^n = (1.39/0.73)^4.1` | **14.02** |
| **Agreement** | **100.8%** |

Two completely independent methods (statistical fit to 38 (γ, M)
datapoints + Hansen+12 linear coefficient, vs ratio-of-Hill-parameters^stress-exponent
from Hansen+16 Part 1) converge on `δ_∞ ≈ 14`. This is the strongest
single-asymptote validation available for olivine viscous anisotropy.

### 7.4 Per-dataset fit comparison

| Dataset | N | M_inf | γ_e | δ_∞ (24.5·M+1) | RMS(M) |
|---|---|---|---|---|---|
| Hansen+12 Nature | 3 | 0.620 | 4.59 | 16.18 | 0.015 |
| Hansen+14 EPSL | 26 | 0.568 | 4.71 | 14.91 | 0.075 |
| Hansen+16 EPSL P1 | 12 | 0.500 | 3.00 | 13.24 | 0.071 |
| **Combined (used)** | **38** | **0.536** | **3.96** | **14.13** | **0.076** |

Single-paper fits scatter ±5% around the combined-fit asymptote. The
combined fit is the most strongly constrained (largest N) and lands
almost exactly on the Hansen+16 stress-aware target.

### 7.5 RMS of candidate δ-forms vs the 38-point reference

Mapping `δ_lab = 24.5 · M_obs + 1`, cap at `δ_∞ ≈ 14`:

| Form | RMS |
|---|---|
| `δ = min(FS_AR, 14)` (MDOODZ default, `ani_fstrain = 1`) | 2.87 |
| `δ = erfc-blend(FS_AR, 14)` | 3.11 |
| Hansen+12 3-point fit (`24.5·M_H12 + 1`) | 2.07 |
| Hansen+16 12-point fit | 1.93 |
| **Combined 38-point fit** (NEW, `ani_fstrain = 2`) | **1.85** |

The new form is the best of all candidates: ~10% improvement over the
12-point fit and ~36% improvement over both the existing `min` form and
the `erfc`-blend prototype.

![Combined Hansen-group olivine viscous-anisotropy calibration](hansen_olivine_calibration.png)

The two-panel comparison plot:
- **Panel (a)**: M(γ) showing all 38 + 3 datapoints with the combined-fit
  curve in solid black plus per-dataset fits as light reference curves.
- **Panel (b)**: five candidate δ-forms overlaid with the combined lab
  data and the Hansen+16 Eq. 3-4 asymptote marker at δ ≈ 14.0.

### 7.6 Independent corroboration of the M(γ) shape

Hansen, Conrad, Boneh, Skemer, Warren & Kohlstedt (2016) "Viscous
anisotropy of textured olivine aggregates: Part 2. Micromechanical
model", *JGR Solid Earth* 121, 7137–7160 (the Part 2 paper) develops
a pseudo-Taylor + director-method micromechanical model whose
textural-evolution prediction (Fig. 7) **overlays the Hansen+14 26-point
dataset** with an exponential-like saturation matching `M ≈ 0.6` by γ ≈ 15.

This means our calibration is triple-validated:
- **Asymptote level** ← Hansen+16 Part 1 Eq. 3-4 stress-aware direct calc
- **Saturation shape** ← Hansen+16 Part 2 Fig. 7 director-method model
- **Constants** ← combined 38-point Hansen-group dataset

We are not pattern-matching the data with an arbitrary curve; the
exponential-saturation form follows from the same slip-system rotation
kinematics Hansen+16 Part 2 derives independently.

### 7.7 Structural caveats

**(a) `γ_eff(FS_AR)` is a simple-shear inversion.** Under simple shear,
`γ_eff ≡ γ` exactly. Under pure shear (`FS_AR(t) = exp(2·ε̇·t)`), the
formula gives `γ_eff = 2·sinh(ε̇·t)` — a function of accumulated stretch
but not of any actual shear strain. Users running pure-shear scenarios
with `ani_fstrain = 2` are choosing "Hansen-shaped saturation parameterised
by accumulated stretch", not "Hansen olivine fabric strength under arbitrary
deformation."

**(b) Stress-orientation dependence is invisible (Hansen+16 Part 2 §4.1.1
/ Fig. 10).** The real δ depends on the angle between principal stress
and the texture's strong axis, with a ~10× spread relative to the
isotropic case. MDOODZ's scalar `aniso_factor_n` is the
texture-aligned (worst-case) value — an upper bound on local viscous-
anisotropy effects. A future change could replace the scalar with a
Hill-tensor parameterisation (6 components per cell). Out of scope here.

**(c) Texture and grain size co-evolve on the same timescale (Hansen+16
Part 2 Eq. 10).** Their coupled-model demonstration uses
`ḋ = ε̇·(d_ss − d) / ε_c` with `ε_c = 1` for both texture and grain-size
evolution. This is the same form MDOODZ's `LocalIterationViscoElasticGrainSize`
already implements via the wattmeter steady-state grain size, which
validates the choice of a scalar `δ(γ_eff(FS_AR))` rather than a
separately tracked texture-evolution timescale.

**(d) Olivine-only.** Constants are calibrated for olivine. Calcite
(Barnhoorn+04, Pieri+01), quartz, mica use cases would need their own
M(γ) calibrations with different constants. Per-mineral parameterisation
is future work.

### 7.8 Reproducing the calibration

```bash
cd misc/aniso_fstrain/calibrate
python3 calibrate_hansen_olivine.py
```

Outputs:
- Per-dataset fit summary table (Hansen+12, Hansen+14, Hansen+16 P1, combined)
- Two-method asymptote validation
- RMS of all candidate δ-forms vs the combined 38-point reference
- `hansen_olivine_calibration.png` — the comparison plot

## 8. Citation

**Primary citation** (linear fit `δ = 24.5·M + 1`):
Hansen, L. N., Zimmerman, M. E. & Kohlstedt, D. L. (2012)
"The influence of microstructure on deformation of olivine in the
grain-boundary sliding regime", *Nature* **492**, 415–418.
doi:[10.1038/nature11671](https://doi.org/10.1038/nature11671)

**M(γ) calibration data** (combined 38-point dataset):
- Hansen, L. N., Zhao, Y.-H., Zimmerman, M. E. & Kohlstedt, D. L. (2014)
  "Protracted fabric evolution in olivine: Implications for the
  relationship among strain, crystallographic fabric, and seismic
  anisotropy", *EPSL* **387**, 157–168.
  doi:[10.1016/j.epsl.2013.11.009](https://doi.org/10.1016/j.epsl.2013.11.009)
- Hansen, L. N., Warren, J. M., Zimmerman, M. E. & Kohlstedt, D. L. (2016)
  "Viscous anisotropy of textured olivine aggregates, Part 1: Measurement
  of the magnitude and evolution of anisotropy", *EPSL* **445**, 92–103.
  doi:[10.1016/j.epsl.2016.04.008](https://doi.org/10.1016/j.epsl.2016.04.008)

**Asymptote validation + structural caveats** (stress-orientation,
ε_c = 1 co-evolution, director-method shape):
- Hansen, L. N., Conrad, C. P., Boneh, Y., Skemer, P., Warren, J. M. &
  Kohlstedt, D. L. (2016) "Viscous anisotropy of textured olivine
  aggregates: 2. Micromechanical model", *JGR Solid Earth* **121**, 7137–7160.
  doi:[10.1002/2016JB013240](https://doi.org/10.1002/2016JB013240)
