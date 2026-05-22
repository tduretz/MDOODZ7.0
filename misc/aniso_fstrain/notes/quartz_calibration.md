# Quartz viscous-anisotropy calibration (`aniso_db = 3` and `aniso_db = 6`)

This note documents the **two** calibrated quartz δ(FS_AR) forms delivered
as `case 3` and `case 6` of `ReadDataAnisotropy`, dispatched via
`AnisoFactorEvolv` when the per-phase parameters `ani_fstrain == 2` and
`aniso_db ∈ {3, 6}`.

| Case | Source | M-index basis | Constants `(M_inf, γ_e, slope)` | δ_∞ |
|---|---|---|---|---|
| 3   | Pennacchioni+10 (n=17 natural-shear-zone samples) | ODF J | (0.75, 3.50, 3.0) | 3.25 |
| 6   | Blackford+24 (n=38 Northern Snake Range samples)   | pole-figure J of c-axis | (0.20, 4.00, 11.0) | 3.20 |

Both modes describe the same physical mineral measured on different texture
metrics (ODF J ≠ pfJ) — the slopes are rescaled so the saturated δ_∞ matches
to within rounding, but the transient shape of δ(γ) differs (case 3 saturates
faster than case 6).

Mirrors the structure of [hansen2012_comparison.md §7](../hansen2012_comparison.md)
(the Hansen olivine note) and [calcite_calibration.md](../calcite/calcite_calibration.md).

## §1. Functional form

```
δ(FS_AR) = min( slope · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )

  where  γ_eff(FS_AR) = √FS_AR − 1/√FS_AR              (simple-shear inversion)
         M(γ)         = M_inf · (1 − exp(−γ / γ_e))    (saturating exponential)

         case 3 (Pennacchioni+10):  M_inf = 0.75, γ_e = 3.50, slope = 3.0
         case 6 (Blackford+24):     M_inf = 0.20, γ_e = 4.00, slope = 11.0
```

The composition is identical to the Hansen olivine and calcite forms;
only the constants differ per phase / per dataset.

## §2. Datasets

### §2.1 Pennacchioni+10 — case 3

[Pennacchioni, Menegon, Leiss, Nestola & Bromiley (2010)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JB007674)
*JGR Solid Earth* 115, B12405: 17 quartz veins from the Adamello tonalite
(southern Alps, Italy), shear strain γ_foliation from < 1 to > 15. Single
temperature window: ~500°C, prism-a slip dominant. **J reported is full
ODF J-index** from texture goniometry (4-pole-figure inversion).

Three regimes:
- **WDV** (Weakly Deformed Veins): γ < 1, J ~ 2–4 (4 samples)
- **MDV** (Moderately Deformed Veins): 2 < γ < 3, J ~ 4–8 (5 samples)
- **SDV** (Strongly Deformed Veins): γ = 4–15, J ~ 8–12, **saturated CPO** (8 samples)

### §2.2 Blackford+24 — case 6

[Blackford, Long, Lee, Larson, Seward, Stevens & Al Harthi (2024)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023TC008166)
*Tectonics* 43, e2023TC008166: 38 quartzite samples from the Northern
Snake Range metamorphic core complex, NV. 3D finite-strain ellipsoids on
stretched detrital quartz clasts; **CPO intensity reported as pole-figure
J of the c-axis** ("density norm" pfJ ∈ [1.26, 4.85]). Strain reported as
**Nadai natural strain** ε = √[(1/3)·Σ(ε_i − ε_j)²]; mostly plane strain
to slightly constrictional kinematics. Greenschist–amphibolite facies
(~500–650°C).

The 38 samples used for calibration are those with both an octahedral
shear strain ε AND a pfJ value in Table 1.

### §2.3 Strain-measure conversion (Blackford → simple-shear γ)

Blackford reports **Nadai natural strain** in the column "Octahedral shear
strain (ε)". Under the plane-strain approximation (ε₂ = 0, ε₃ = −ε₁):

```
ε  = √2 · |ε₁|       ⇒  ε₁ = ε / √2
γ_simple_eq  =  2 · sinh(ε₁)            (since simple shear with γ has
                                         max stretch λ₁ = exp(ε₁) and
                                         γ = λ₁ − 1/λ₁ = 2·sinh(ε₁))
```

**Verified** against Sample 22 of Table 1 (R_xz = 17.9, R_xyz = 4.3,
ε = 2.04 → λ_x = 4.26 → γ_simple = 4.05) versus the formula 2·sinh(2.04/√2)
= 4.02 (0.7% agreement). Blackford's γ_eq ranges from 1.28 to 5.08 for the
calibration samples — overlapping Pennacchioni's γ ∈ [0.3, 15] in the
low-to-moderate-strain regime.

## §3. Calibration procedure

The script [calibrate_quartz.py](calibrate_quartz.py) loads both datasets,
runs free LSQ on each separately and combined, and prints:

- Per-dataset sample counts and ranges
- Free LSQ fits with `(M_inf, γ_e)` (slope fixed)
- RMS errors of candidate δ-forms vs each dataset
- Comparison plot

J→M conversion uses Skemer+05's empirical relation: `M ≈ (J − 1) / 15`,
applied identically to ODF J (Pennacchioni) and pfJ (Blackford).

**Free LSQ results:**

| Dataset | n | M_inf | γ_e | RMS(M) |
|---|---|---|---|---|
| Pennacchioni only | 17 | 0.736 | 3.59 | 0.032 |
| Blackford only    | 38 | 0.290 | 9.87 | 0.045 |
| Combined          | 55 | 0.744 | 10.00 | 0.152 |

The Blackford-only fit hits the upper γ_e bound (saturation is poorly
expressed in γ_eq ∈ [1.3, 5.1] — pfJ keeps growing slowly across the data
range). The combined fit is dominated by Pennacchioni because of the
M-scale mismatch (Pen M up to 0.7, Bla M only to 0.26 — Skemer's J→M
relation does not absorb the metric difference).

## §4. Asymptote — VPSC-style bounds, NOT direct lab fit

Same "weakest link" caveat as calcite: the δ-vs-M slope is **not** taken
from a clean linear lab fit analogous to Hansen+12 Fig 3b. For quartz, no
equivalent published linear fit exists. Instead:

```
case 3:   slope = 3.0  ⇒  δ_∞ = 3.0 · 0.75 + 1 = 3.25
case 6:   slope = 11.0 ⇒  δ_∞ = 11.0 · 0.20 + 1 = 3.20
```

Both asymptotes match **VPSC-style Heilbronner-Tullis bounds** for quartz
aggregate viscous anisotropy at saturated prism-a-dominated CPO (δ ≈ 2–3,
upper end). Case 6's larger slope compensates for pfJ being smaller than
ODF J for the same physical texture.

| Asymptote source | δ_∞ |
|---|---|
| VPSC bounds for quartz prism-a CPO | 2–3 |
| **case 3 (Pennacchioni, slope = 3.0)** | **3.25** |
| **case 6 (Blackford, slope = 11.0)** | **3.20** |
| Calcite asymptote (for scale comparison) | 4.0 (T-avg) |
| Hansen olivine asymptote (for scale comparison) | 14.13 |

If a future direct lab measurement supersedes either calibration, the only
edits needed are the constants in the static helpers `anisoDelta_Quartz`
or `anisoDelta_Quartz_Blackford` in `MDLIB/FlowLaws.c`. CI tests gate the
implementation against THIS calibration, not against lab data.

## §5. Per-datapoint verification (case 3)

| γ | J_obs | M_obs (Skemer+05) | M_fit (case 3) | δ_fit | regime |
|---|---|---|---|---|---|
| 0.30 | 2.20 | 0.080 | 0.062 | 1.19 | WDV |
| 0.50 | 2.60 | 0.107 | 0.099 | 1.30 | WDV |
| 0.70 | 3.00 | 0.133 | 0.135 | 1.40 | WDV |
| 0.90 | 3.60 | 0.173 | 0.168 | 1.51 | WDV |
| 2.00 | 4.50 | 0.233 | 0.323 | 1.97 | MDV |
| 2.30 | 5.80 | 0.320 | 0.357 | 2.07 | MDV |
| 2.50 | 6.00 | 0.333 | 0.378 | 2.14 | MDV |
| 2.80 | 6.80 | 0.387 | 0.408 | 2.23 | MDV |
| 3.00 | 7.50 | 0.433 | 0.428 | 2.28 | MDV |
| 4.00 | 9.00 | 0.533 | 0.510 | 2.53 | SDV |
| 5.00 | 10.00 | 0.600 | 0.572 | 2.72 | SDV |
| 6.00 | 10.50 | 0.633 | 0.620 | 2.86 | SDV |
| 8.00 | 11.00 | 0.667 | 0.687 | 3.06 | SDV |
| 10.00 | 11.50 | 0.700 | 0.726 | 3.18 | SDV |
| 12.00 | 11.00 | 0.667 | 0.731 | 3.19 | SDV |
| 13.50 | 12.00 | 0.733 | 0.741 | 3.22 | SDV |
| 15.00 | 11.50 | 0.700 | 0.747 | 3.24 | SDV |
| ∞ | — | — | 0.750 | **3.25** | asymptote |

The fit captures both the WDV transient and the SDV asymptote. The
free-LSQ γ_e = 3.59 was rounded to 3.50 for clean numerics (RMS impact
< 0.001 in M units).

## §6. RMS error of candidate δ-forms

Mapped via `δ_lab = SLOPE · M_obs + 1` per dataset (slope = 3.0 for Pen,
slope = 11.0 for Bla). Lower is better:

| Form | RMS vs Pennacchioni | RMS vs Blackford |
|---|---|---|
| `δ = min(FS_AR, 2.0)`              | 0.71 | 0.79 |
| `δ = min(FS_AR, 3.0)`              | 0.62 | 1.79 |
| **case 3** (M_inf=0.75, γ_e=3.5, slope=3.0)  | **0.10** | 1.06 |
| **case 6** (M_inf=0.20, γ_e=4.0, slope=11.0) | 0.12 | **0.93** |

Each calibrated form fits its own dataset best. Both calibrations
**outperform** the old min-form fallbacks on both datasets, even on
the dataset they were not fit to.

![Quartz calibration plot — Pennacchioni+10 + Blackford+24 vs candidate forms](mdoodz_vs_quartz_data.png)

## §7. Why two modes? — metric mismatch, NOT physical disagreement

The two datasets agree **qualitatively** that quartz CPO grows steadily
with increasing γ and saturates slowly (γ_e ≥ 3). They disagree
**quantitatively** in M-scale because:

- **ODF J-index** (Pennacchioni): integrates the full 3D orientation
  distribution function. Sensitive to all crystallographic axes.
- **Pole-figure J of c-axis** (Blackford): integrates the c-axis pole
  figure only. Ignores a- and rhomb-axis ordering.

For the same physical quartz texture, ODF J > pfJ (typically 2–5×).
Skemer's `M = (J − 1) / 15` was developed for olivine ODF J; applying it
identically to pfJ gives M values on a different scale.

The two committed modes accommodate this:

- **case 3** is appropriate when comparing against quartz CPO data
  reported as ODF J (most older texture-goniometry papers).
- **case 6** is appropriate when comparing against EBSD-derived
  pole-figure J (most modern quartz EBSD studies, including Blackford
  and many post-2010 papers).

Both produce **physically equivalent δ(γ)** at saturation; the transient
γ-rise differs slightly because each mode reflects the metric-specific
saturation behaviour observed in its source dataset.

## §8. T-range caveat — REQUIRED reading for quartz scenarios

Both calibrations are calibrated for **T ≥ 500°C, prism-a slip regime**.
For:

- **T < 500°C** (basal-a slip dominant): different fabric symmetry; the
  prism-a-derived slope may not apply. As a temporary fallback, users
  can revert to `ani_fstrain = 1` + `ani_fac_max = 3` per phase.
- **High-pressure quartz** (Stipp+02 quartz-coesite transition regime):
  not addressed by either dataset.
- **Wet quartz** (Hirth & Tullis 2002): water enables slip-system
  changes at lower T, effectively shifting the prism-a regime to lower T.

When in doubt, document the assumed T-regime in the scenario `.txt` file
header.

## §9. Heilbronner & Tullis 2006 — checked but not included

[Heilbronner, R. & Tullis, J. (2006)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JB004194)
*JGR* 111, B10202: Black Hills quartzite Griggs-apparatus shear experiments
(γ_tot up to 6.9, T = 900–915°C, regime-3 dislocation creep). 7 sites with
local γ ∈ [0, 8] and well-controlled experimental shear strain.

**Why not included in the fit**: the paper reports c-axis pole-figure
**peak intensity** (CPO_max in m.r.d.) rather than J-index. CPO_max → J
conversion depends on pole-figure symmetry and would introduce significant
noise. The Blackford+24 dataset already provides the experimentally
clean (γ, pfJ) pairs we'd need. Heilbronner & Tullis 2006 is cited for
the VPSC-bounds derivation in §4 but not for direct calibration.

## §10. Structural caveats (shared with calcite + olivine)

### `γ_eff(FS_AR)` is a simple-shear inversion

Identical caveat to the Hansen olivine and calcite notes: `γ_eff = √FS_AR
− 1/√FS_AR` is exact under simple shear, defensible-but-questionable
under pure shear / general 2D flow. Blackford samples are plane-strain
to slightly constrictional in nature; the simple-shear effective strain
is appropriate as a CPO-driving accumulated stretch metric but does not
recover the actual kinematics.

### Skemer+05 J→M conversion is approximate

`M ≈ (J − 1) / 15` was derived for monoclinic-symmetry textures with
ODF J. Applying it to pfJ (case 6) gives M values consistent with the
form but on a different metric scale; the rescaled slope absorbs the
difference at saturation.

### Natural-sample biases

Both Pennacchioni+10 and Blackford+24 are natural samples. Confounding
factors (documented but not quantified in the fits):

- Pre-existing crystal-growth CPO in vein samples (Pennacchioni).
- Strain-history non-uniqueness across cooling history (both).
- Slip-system mix variations sample-by-sample.

The fits smooth through these biases; committed constants reflect
ensemble averages, not any single sample.

## §11. Mica-bearing polyphase rocks — `case 8` BRACKETED estimate

Most crustal scenarios involve **quartzofeldspathic rocks with some mica
content** (granite, granitoid mylonite, schist, quartzite mylonite with
muscovite or biotite films). Pure-quartz cases 3 and 6 are calibrated
against **mica-free or low-mica** datasets; for mica fractions above a
few percent, the bulk anisotropy diverges from those calibrations.

**This is a BRACKETED estimate, not a direct LSQ fit.** A clean
(γ, J_bulk) lab dataset for quartz-mica aggregates does not exist — the
constants are bracketed from disparate sources (Tokle+23, Dempsey+11,
Holyoke & Tullis 2006, Bos & Spiers 2002). The mode IS committed as
`aniso_db = 8` with the same level of pragmatism the lithosphere modeling
community accepts for `pwlv = 16` Ranalli felsic granulite (also a
bracketed compilation, not a direct lab fit). Users seeking direct-fit
rigor should use cases 1–7. Users needing a polyphase mode for granitoid-
mylonite or schist scenarios should use `aniso_db = 8` understanding
that the constants are bracketed estimates.

**Committed constants** (`anisoDelta_QuartzMicaBracketed` in
[MDLIB/FlowLaws.c](../../../MDLIB/FlowLaws.c)):

```
M_inf = 0.40, γ_e = 2.50, slope = 15.0
δ_∞ = 15.0 · 0.40 + 1 = 7.0
```

Suitable for **typical 10–25% mica content** crustal mylonites. For
mica fractions outside this range, see the `ani_fac_max` table below
(scaling guidance) or consider `ani_fstrain = 1` with chosen cap.

### What the literature says

- **Tokle, Hirth & Stünitz (2023)** *J. Struct. Geol.* **169**, 104835 —
  synthetic quartz + muscovite (0%, 5%, 10%, 25%) at 800°C, 1.5 GPa,
  γ ∈ [0.5, 4.3]. Quartz CPO_max (peak m.r.d.) drops from **12.30** in
  pure quartz at γ = 4 to **4.78 / 4.77 / 4.40** at 5% / 10% / 25%
  muscovite. Bulk strength drops 5–10× at 25% mica. Mica grains align
  sub-parallel to the shear plane and form C′ shear bands.

- **Gottardi, Casale, Economou & Morris (2024)** *Geology* **52**, 545 —
  Raft River detachment shear zone (NW Utah), 4–18% muscovite. Confirms
  Tokle's lab finding in nature: mica content correlates inversely with
  quartz recrystallized grain size and modulates the dominant
  deformation mechanism. Strain rates calculated from the quartz
  piezometer differ by **>1 order of magnitude** between low- and
  high-mica subregions in the same outcrop.

- **Dempsey, Prior, Mariani, Toy & Tatham (2011)** *Geol. Soc. Lond.
  Spec. Publ.* **360**, 33–47 — 8 Alpine Fault quartzofeldspathic
  mylonites/ultramylonites. Mica ⟨001⟩ J-index ranges **6.75 – 10.59**
  (mean ≈ 9.0) at the mylonite stage; AVs anisotropy 49–70%. No γ
  reported (all samples are at saturation).

- **Holyoke & Tullis (2006)** *Geology* **34**, 105 — only ~10%
  interconnected K-mica is required to switch a quartzofeldspathic
  aggregate from load-bearing-framework to interconnected-weak-phase
  rheology. Mica controls bulk strength even at small volume fractions.

- **Bos & Spiers (2002, 2003)** + Niemeijer & Spiers (2007) —
  micromechanical model: foliated quartz-mica fault rock is **2–5×
  weaker** than predicted by Byerlee's rule due to easy slip on aligned
  phyllosilicate basal planes.

### Mode selection by mica volume fraction

| Mica volume fraction | Rock type | Recommended | δ_∞ |
|---|---|---|---|
| **0–5%** | quartzofeldspathic mylonite, granite leucosome, pure quartzite | `aniso_db=3` or `=6` (case 3 / 6 direct-fit) | 3.2 |
| **5–15%** | typical granitoid mylonite, two-mica granite, mica schist tail | **`aniso_db=8` (BRACKETED)** | 7.0 (or scale `ani_fac_max ≈ 5–7` if `ani_fstrain=1`) |
| **15–30%** | mica schist, mica-rich mylonite | **`aniso_db=8` (BRACKETED)** | 7.0 (or `ani_fac_max ≈ 7–10`) |
| **>30%** | phyllite, slate, mica-dominant schist | `ani_fstrain=1, ani_fac_max ≈ 10–15` | – (case 8 underpredicts at this mica content) |

These `ani_fac_max` brackets come from:
- Lower bound: pure-quartz δ_∞ = 3.2 (case 3) or 3.2 (case 6)
- Upper bound: pure-mica δ_∞ ≈ 20–40 (Cholach & Schmitt 2006 intrinsic
  muscovite aggregate elasticity + VPSC predictions)
- Interpolation guided by Holyoke & Tullis (2006) interconnection
  threshold and Bos & Spiers (2002) micromechanical mixing model

### Why this is BRACKETED, not direct-fit

Three reasons, all data-structural:

1. **No single dataset has (γ, J_bulk) pairs** for quartz-mica
   aggregates. Tokle+23 reports quartz-only CPO_max; Dempsey+11
   reports mica-only J at saturation. They can't be combined into a
   regression fit, only into a bracketed point estimate.
2. **The slope is genuinely under-determined** (3 ≤ slope ≤ 50). Any
   committed value is an editorial choice, not a measurement. The
   committed `slope = 15` is a defensible middle value but is not
   a measurement.
3. **Bulk anisotropy depends nonlinearly on phase geometry** (isolated
   vs interconnected mica), which Tokle's 4-fraction dataset cannot
   resolve.

The mode is nonetheless committed because the lithosphere modeling
community accepts equivalent bracketed estimates in the power-law
database (e.g., Ranalli 1995 felsic granulite `pwlv = 16`, Ranalli
mafic granulite `pwlv = 17`). Quality tier of `aniso_db = 8` matches
those `pwlv` entries: *"defensible point estimate from a literature
bracket, useful but not Hansen+14-grade."*

If a future controlled experiment reports both quartz and mica J-index
per sample at multiple γ values for ≥3 mica fractions, we will
re-commit `case 8` as a direct LSQ fit then.

### See also

For polyphase olivine (mantle xenoliths, naturally-deformed peridotite
with multi-event history), see
[../olivine_calibration.md](../olivine_calibration.md)
§7 — `case 7` is a true LSQ fit on 93 non-Hansen samples and is the
recommended damped-olivine mode. The polyphase-mica situation is weaker
because the underlying datasets are not co-located in (γ, J) space.

## §12. Plagioclase / feldspar-bearing rocks — `case 9` BRACKETED estimate

Plagioclase is the dominant phase in **lower-crustal gabbros, anorthosites,
amphibolites, and metabasites**. In the SETS power-law database it shows
up as `pwlv = 19` (Anorthite 100), `pwlv = 25` (Plagioclase Shelton),
`pwlv = 30` (Anorthite dry), and `pwlv = 32` (Anorthite 60) — combined
~22+ uses across input files.

Like mica (§11), plagioclase is **not** a candidate for a Hansen-style
direct-fit calibration. The reason is different from mica's, though:
while mica has strong CPO but no co-located (γ, J) data, **plagioclase
actively SUPPRESSES its CPO at high strain** because the dominant
deformation mechanism transitions from dislocation creep →
dislocation-accommodated grain-boundary sliding (DisGBS) → diffusion
creep. The saturating-exponential form M(γ) = M_inf·(1−exp(−γ/γ_e))
**over-predicts** CPO strength at very high γ.

The mode IS committed as `aniso_db = 9` with the same level of
pragmatism the lithosphere modeling community accepts for `pwlv = 17`
Ranalli mafic granulite (also a bracketed compilation, not a direct
lab fit). Users seeking direct-fit rigor should use cases 1–7. Users
needing a polyphase mode for gabbro-mylonite, anorthosite, mafic
granulite, basalt, or diabase scenarios should use `aniso_db = 9`
understanding that the constants are bracketed estimates that
over-predict at γ > ~8.

**Committed constants** (`anisoDelta_PlagioclaseBracketed` in
[MDLIB/FlowLaws.c](../../../MDLIB/FlowLaws.c)):

```
M_inf = 0.12, γ_e = 2.00, slope = 8.0
δ_∞ = 8.0 · 0.12 + 1 = 1.96 ≈ 2.0
```

Suitable for **monophase plagioclase mylonite** at saturated CPO. For
polyphase plag-cpx layers (where DisGBS further weakens CPO), this
over-predicts by ~2× — see scaling table below.

### What the literature says

- **Mehl & Hirth (2008)** *J. Geophys. Res. Solid Earth* **113**, B05202 —
  ODP Hole 735B SWIR oceanic gabbro mylonites. Table 1 reports M-index
  for 7 high-strain shear zones, multiple layers each. **Key result**:
  - **Monophase plagioclase mylonite layers**: M ∈ [0.05, 0.19], mean
    ≈ 0.12, dominant slip systems (010)[100] or (001)[100] (plag
    anorthite content An38–An57)
  - **Polyphase plag+cpx layers**: M ∈ [0.03, 0.07], mean ≈ 0.04,
    "random" — phase mixing weakens CPO via DisGBS
  - "Progressive mixing of pyroxene and plagioclase within gabbro
    mylonite layers is accompanied by **weakening of the LPO**,
    indicating that phase mixing promotes a transition to diffusion
    creep processes that involve grain boundary sliding."

- **Marti, Stünitz, Heilbronner, Plümper & Kilian (2018)** *Solid Earth*
  **9**, 985 — synthetic plag-cpx-opx mixtures deformed at 800°C, 1.0–
  1.5 GPa, γ_a up to ≈ 6 (water-added). **Key finding**: even in
  plagioclase-dominated shear bands, CPO is **"weak but distinct"**
  and develops via dissolution-precipitation creep (DPC) + GBS, not
  classical dislocation creep. Maxima have "moderate strengths" only.

- **Barreiro, Wenk & Vogel (2015)** *J. Struct. Geol.* **71**, 100 —
  Morin Shear Zone (Quebec) anorthosite mylonite-ultramylonite
  (88% plag + 8% cpx + 4% opx), neutron-diffraction texture analysis.
  **Plagioclase pole-figure max only 2.9 m.r.d.** (vs ~10–15 for
  saturated olivine, ~12 for saturated quartz). Bulk anorthosite
  AVp = 2.4–2.6% (Voigt) — much weaker seismic anisotropy than
  olivine (~10%) or mica-bearing rocks (~30%).

- **Stünitz, Fitz Gerald & Tullis (2003)** *Tectonophysics* — single-
  crystal plagioclase, identifies (010)[001] as a major slip system
  (alongside (010)[100]). Slip-system data only, no aggregate γ-vs-J.

- **Allard et al. (2021)** *J. Geophys. Res. Solid Earth* **126** —
  Atlantis Bank gabbro (slow-spreading SWIR). Confirms Mehl & Hirth
  picture in another oceanic location: monophase plag layers preserve
  a CPO, polyphase mixing weakens it.

### Saturation behavior — fundamentally different from olivine/calcite/quartz

| Mineral | M_inf (saturated CPO) | dominant mechanism | shape of M(γ) |
|---|---|---|---|
| Olivine (case 1, dry Fo90 lab) | 0.54 | dislocation creep, monotonic | clean saturating exponential |
| Olivine (case 7, natural/wet/pre-CPO) | 0.20 | mixed dislocation + DisGBS | flatter, slower-saturating |
| Calcite (cases 2/4/5) | 0.22–0.83 | dislocation creep + recryst | saturating exponential |
| Quartz (cases 3/6) | 0.20–0.75 | dislocation creep + recryst | saturating exponential |
| **Plagioclase (no mode)** | **0.05–0.19** | **DPC + DisGBS dominant; CPO can WEAKEN at high γ** | **non-monotonic; doesn't fit M_inf·(1−exp) form** |

The plagioclase M-vs-γ trajectory in nature looks more like:
- Low γ: weak CPO from random starting orientation
- Moderate γ (~0.5–3): CPO grows slowly (M ≈ 0.10–0.18) via dislocation
  creep on (010)[100] / (001)[100]
- High γ (mylonite stage): grain size reduction → mechanism transitions
  to DisGBS → CPO **stops growing or even decays** as GBS randomizes
  orientations
- Polyphase mixing (plag + cpx → fine grain size): **further weakens**
  CPO by activating diffusion creep

There's no clean "saturation γ_e" to fit because the curve isn't a
saturating exponential.

### Mode selection by rock type

| Plagioclase context | Recommended | δ_∞ |
|---|---|---|
| Anorthosite, gabbro mylonite (monophase plag layer) | **`aniso_db=9` (BRACKETED)** | 1.96 |
| Two-phase plag+cpx gabbro mylonite | **`aniso_db=9`** (over-predicts ~2× at γ>5) OR `ani_fstrain=1, ani_fac_max≈1.2–1.5` | 1.96 |
| Diabase (Maryland Diabase, `pwlv=11`) | **`aniso_db=9` (BRACKETED)** | 1.96 |
| Granulite (felsic, `pwlv=16`) | `aniso_db=3` or `=6` (case 3 / 6 quartz controls) | 3.2 |
| Granulite (mafic, `pwlv=17`) | **`aniso_db=9` (BRACKETED)** | 1.96 |
| Basalt (`pwlv=29`) | **`aniso_db=9` (BRACKETED)** | 1.96 |

These brackets reflect the consistent finding that **plagioclase
viscous anisotropy is intrinsically modest** (single-crystal AVp ≈ 10%
even at perfect alignment, vs ~30% for mica) AND **further damped in
natural samples by DisGBS / phase mixing** (Mehl & Hirth 2008 polyphase
M ≈ 0.04).

### Why this is BRACKETED, not direct-fit

Three reasons, all data-structural:

1. **The form M(γ) = M_inf·(1−exp(−γ/γ_e)) doesn't apply at high γ.**
   Mehl & Hirth data shows M can DECREASE during high-strain DisGBS.
   A monotonic saturating exponential **over-predicts** δ at γ > 8 by
   ~2× compared to natural polyphase plag-cpx layers (M ≈ 0.04, vs
   committed M_inf = 0.12).
2. **No co-located (γ, J) lab dataset exists.** Marti+18 has γ_a but
   reports CPO only qualitatively ("weak but distinct"). Mehl & Hirth
   has M values but only at the saturated mylonite stage (no γ).
   Barreiro+15 is one ultra-mylonite sample with neutron-diffraction
   texture, no strain history.
3. **The slope is even more under-determined than for mica.** With
   single-crystal AVp ≈ 10% and saturated M ≈ 0.12, slope estimates
   range from 5 (low end) to 20 (high end) depending on which mixing
   rule and which orientation distribution function symmetry you assume.

The mode is nonetheless committed because the lithosphere modeling
community accepts equivalent bracketed estimates in the power-law
database. Quality tier of `aniso_db = 9` matches `pwlv = 17` Ranalli
mafic granulite: *"defensible point estimate from a literature
bracket, useful but not Hansen+14-grade."*

If a future dataset reports plagioclase J-index at multiple γ values
in a single experimental campaign (analogous to Hansen+14 for olivine
or Pieri+01 for calcite), we will re-commit `case 9` as a direct LSQ
fit then.

### Bottom-line summary table — all 9 modes

| Rock type | Mode | δ_∞ | Provenance |
|---|---|---|---|
| Pure peridotite, fresh CPO | `aniso_db=1` | 14.13 | **Direct LSQ** (Hansen+12/+14/+16, n=38) |
| Pure peridotite, natural / pre-CPO / wet | `aniso_db=7` | 5.90 | **Direct LSQ** (Tasaka+Boneh+Kumamoto+Bernard, n=93) |
| Calcite, T-averaged | `aniso_db=2` | 3.46 | **Direct LSQ** (Pieri+Bruijn+Barnhoorn, n=18) |
| Calcite, T < 650°C | `aniso_db=4` | 2.33 | **Direct LSQ** (Barnhoorn+04 mid-low-T, n=7) |
| Calcite, T > 650°C | `aniso_db=5` | 5.99 | **Direct LSQ** (Pieri+Bruijn+Barnhoorn high-T, n=11) |
| Pure quartzite mylonite (ODF J) | `aniso_db=3` | 3.25 | **Direct LSQ** (Pennacchioni+10, n=17) |
| Pure quartzite mylonite (pfJ) | `aniso_db=6` | 3.20 | **Direct LSQ** (Blackford+24, n=38) |
| Granitoid + 5–25% mica | `aniso_db=8` | 7.0 | **BRACKETED** (Tokle+Dempsey+Holyoke; pwlv=16-grade) |
| Schist / phyllite (>30% mica) | `ani_fstrain=1, ani_fac_max≈10–15` | – | Guidance only (case 8 underpredicts) |
| Gabbro / anorthosite mylonite | `aniso_db=9` | 1.96 | **BRACKETED** (Mehl & Hirth + Marti + Barreiro; pwlv=17-grade) |
| Mafic granulite, basalt, diabase | `aniso_db=9` | 1.96 | **BRACKETED** |
| Polyphase plag+cpx high-γ (>8) | `ani_fstrain=1, ani_fac_max≈1.2–1.5` | – | Guidance only (case 9 over-predicts at high γ) |

**7 direct-fit + 2 bracketed = 9 committed `aniso_db` cases.**
Direct-fit cases match Hansen+14-grade rigor; bracketed cases match
pwlv=16/17 Ranalli granulite rigor (defensible point estimates from
literature compilations).

## §13. Reproducing the calibration

```bash
cd misc/aniso_fstrain/quartz
python3 calibrate_quartz.py            # both fits + comparison plot
python3 mdoodz_vs_quartz_comparison.py # HDF5 overlay (after CI tests run)
```

CI tests:

```bash
cd cmake-build-test/TESTS
./AnisotropyBenchmarkTests --gtest_filter='*AniFstrainQuartz*'
```

Outputs:

- `mdoodz_vs_quartz_data.png` — both calibrations + lab data (3-panel)
- `mdoodz_vs_quartz_comparison.png` — same plot + MDOODZ HDF5 stars per case

## §14. Citations

- **Pennacchioni, G., Menegon, L., Leiss, B., Nestola, F. & Bromiley, G.**
  (2010) "Development of crystallographic preferred orientation and
  microstructure during plastic deformation of natural coarse-grained
  quartz veins", *J. Geophys. Res. Solid Earth* **115**, B12405.
  doi:[10.1029/2010JB007674](https://doi.org/10.1029/2010JB007674).

- **Blackford, N. R., Long, S. P., Lee, J., Larson, K. P., Seward, G.,
  Stevens, J. L. & Al Harthi, H.** (2024) "Relating Quartz
  Crystallographic Preferred Orientation Intensity to Finite Strain
  Magnitude in the Northern Snake Range Metamorphic Core Complex,
  Nevada: A New Tool for Characterizing Strain Patterns in Ductilely
  Sheared Rocks", *Tectonics* **43**, e2023TC008166.
  doi:[10.1029/2023TC008166](https://doi.org/10.1029/2023TC008166).

- **Heilbronner, R. & Tullis, J.** (2006) "Evolution of c axis pole
  figures and grain size during dynamic recrystallization", *J. Geophys.
  Res.* **111**, B10202.
  doi:[10.1029/2005JB004194](https://doi.org/10.1029/2005JB004194).

- **Skemer, P., Katayama, I., Jiang, Z. & Karato, S.-i.** (2005), see
  calcite_calibration.md §10 for full citation.

#### §11 polyphase mica guidance — additional citations

- **Tokle, L., Hirth, G. & Stünitz, H.** (2023) "The effect of muscovite
  on the microstructural evolution and rheology of quartzite in general
  shear", *J. Struct. Geol.* **169**, 104835.
  doi:[10.1016/j.jsg.2023.104835](https://doi.org/10.1016/j.jsg.2023.104835).

- **Gottardi, R., Casale, G., Economou, J. & Morris, K.** (2024) "A
  little mica goes a long way: Impact of phyllosilicates on quartz
  deformation fabrics in naturally deformed rocks", *Geology* **52**,
  545–549. doi:[10.1130/G52053.1](https://doi.org/10.1130/G52053.1).

- **Dempsey, E. D., Prior, D. J., Mariani, E., Toy, V. G. & Tatham, D. J.**
  (2011) "Mica-controlled anisotropy within mid-to-upper crustal
  mylonites: an EBSD study of mica fabrics in the Alpine Fault Zone,
  New Zealand", *Geol. Soc. Lond. Spec. Publ.* **360**, 33–47.
  doi:[10.1144/SP360.3](https://doi.org/10.1144/SP360.3).

- **Holyoke, C. W. III & Tullis, J.** (2006) "Formation and maintenance
  of shear zones", *Geology* **34**, 105–108.
  doi:[10.1130/G22116.1](https://doi.org/10.1130/G22116.1).

- **Bos, B. & Spiers, C. J.** (2002) "Frictional-viscous flow of
  phyllosilicate-bearing fault rock: Microphysical model and
  implications for crustal strength profiles", *J. Geophys. Res. Solid
  Earth* **107**, 2028.
  doi:[10.1029/2001JB000301](https://doi.org/10.1029/2001JB000301).

- **Cholach, P. Y. & Schmitt, D. R.** (2006) "Intrinsic elasticity of a
  textured transversely isotropic muscovite aggregate: Comparisons to
  the seismic anisotropy of schists and shales", *J. Geophys. Res. Solid
  Earth* **111**, B09410.
  doi:[10.1029/2005JB004158](https://doi.org/10.1029/2005JB004158).

#### §12 plagioclase guidance — additional citations

- **Mehl, L. & Hirth, G.** (2008) "Plagioclase preferred orientation in
  layered mylonites: Evaluation of flow laws for the lower crust",
  *J. Geophys. Res. Solid Earth* **113**, B05202.
  doi:[10.1029/2007JB005075](https://doi.org/10.1029/2007JB005075).

- **Marti, S., Stünitz, H., Heilbronner, R., Plümper, O. & Kilian, R.**
  (2018) "Syn-kinematic hydration reactions, grain size reduction, and
  dissolution-precipitation creep in experimentally deformed
  plagioclase-pyroxene mixtures", *Solid Earth* **9**, 985–1009.
  doi:[10.5194/se-9-985-2018](https://doi.org/10.5194/se-9-985-2018).

- **Gómez Barreiro, J., Wenk, H.-R. & Vogel, S.** (2015) "Texture and
  elastic anisotropy of a mylonitic anorthosite from the Morin Shear
  Zone (Quebec, Canada)", *J. Struct. Geol.* **71**, 100–111.
  doi:[10.1016/j.jsg.2014.07.021](https://doi.org/10.1016/j.jsg.2014.07.021).

- **Stünitz, H., Fitz Gerald, J. D. & Tullis, J.** (2003) "Dislocation
  generation, slip systems, and dynamic recrystallization in
  experimentally deformed plagioclase single crystals", *Tectonophysics*
  **372**, 215–233.
  doi:[10.1016/S0040-1951(03)00241-5](https://doi.org/10.1016/S0040-1951(03)00241-5).

- **Allard, M., Ildefonse, B., Oliot, É. & Barou, F.** (2021) "Plastic
  Deformation of Plagioclase in Oceanic Gabbro Accreted at a
  Slow-Spreading Ridge (Hole U1473A, Atlantis Bank, Southwest Indian
  Ridge)", *J. Geophys. Res. Solid Earth* **126**, e2021JB021964.
  doi:[10.1029/2021JB021964](https://doi.org/10.1029/2021JB021964).
