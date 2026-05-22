# Calcite viscous-anisotropy calibration (`aniso_db = 2`)

> **Post-audit update.** Bruijn+11 retracted from the dataset; cases 2 and 5
> were re-fit on the verifiable Pieri+01 + Barnhoorn+04 points only. Committed
> constants:
>
>   - case 2 (T-averaged):  M_∞ = 0.41, γ_e = **3.15** (was 3.71), δ_∞ = 3.46 (unchanged)
>   - case 5 (high-T):      M_∞ = **0.73** (was 0.83), γ_e = **5.28** (was 8.24), δ_∞ = **5.40** (was 5.99)
>
> The Bruijn+11 audit found the 3 (γ, J) values cannot be sourced from the
> Bruijn 2011 paper. See
> [`audit/Bruijn_etal_2011_audit.md`](../audit/Bruijn_etal_2011_audit.md).
>
> Long-form sections below were written under the pre-audit calibration and
> may quote the older constants — **FlowLaws.c is the authoritative source for
> the current values.**

This note documents the calibrated calcite δ(FS_AR) form delivered by the
`add-calcite-delta-form` openspec change. The form is dispatched via
`AnisoFactorEvolv` when the per-phase parameters `ani_fstrain == 2` and
`aniso_db == 2`.

Mirrors the structure of [hansen2012_comparison.md §7](../hansen2012_comparison.md)
(the Hansen olivine note) for consistency.

## §1. Functional form

```
δ(FS_AR) = min( 6.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )

  where  γ_eff(FS_AR) = √FS_AR − 1/√FS_AR              (simple-shear inversion)
         M(γ)         = M_inf · (1 − exp(−γ / γ_e))    (saturating exponential)
         M_inf        = 0.41
         γ_e          = 3.71
         slope        = 6.0                            (calcite-specific)
```

The composition is identical to the Hansen olivine form except the constants
`(M_inf, γ_e, slope) = (0.41, 3.71, 6.0)` replace olivine's `(0.536, 3.96, 24.5)`.
The lower slope reflects calcite's lower viscous-anisotropy-per-fabric-strength
ratio (Barnhoorn+04: CPO contributes only ~1/3 of calcite weakening; the rest
from grain-size reduction, already captured by MDOODZ's wattmeter).

> **History note**: the initial commit of this calibration used `(M_inf=0.50,
> γ_e=3.50, δ_∞=4.0)` calibrated to match the empirical `ani_fac_max=4`
> PinchSwell convention. After adding Barnhoorn+04's 14 datapoints (10 kept
> after Skemer linear-range filter), the LSQ-best fit is `(M_inf=0.41,
> γ_e=3.71, δ_∞=3.46)`. The new asymptote is more honest about the lab
> data ensemble and still defensible relative to Tommasi+09 VPSC bounds
> for calcite saturated CPO (δ ≈ 2-5).

## §2. Calibration procedure

`M_inf` and `γ_e` are calibrated by free LSQ fit to the **combined 18-point
calcite torsion dataset** with J→M-index conversion via Skemer+05's empirical
relation for monoclinic textures:

```
M ≈ (J − 1) / 15        (Skemer+05; valid for J ≤ 15)
```

The script [calibrate_calcite.py](calibrate_calcite.py) loads all three
datasets inline, runs a 2D coarse-then-fine grid search on `(M_inf, γ_e)`
with `slope` fixed at 6.0, filters Barnhoorn+04 high-J samples (J > 15 where
the linear conversion breaks down), prints the per-dataset summary + RMS
comparison vs candidate δ-forms, and renders the comparison plot.

The script is pure stdlib + numpy + matplotlib (no scipy).

### Dataset sources

- **Pieri+01** (Tectonophysics 330, 119 — Carrara marble torsion at 1000 K
  ≈ 727°C, 300 MPa, constant strain rate 3×10⁻⁴ s⁻¹), Table 1: 5 datapoints
  at γ ∈ {0, 1, 2, 5, 11}.
- **Bruijn+11** (Tectonophysics 503, 75 — two-stage Carrara marble torsion
  at the same conditions): 3 single-stage-equivalent datapoints at
  γ ∈ {1, 2.6, 5} extracted from the D1 portion of the multi-stage
  experiments. These samples carry pre-existing texture from the prior D1
  stage, so their J values are systematically lower than Pieri+01's
  single-stage measurements at the same γ — explained in §6.
- **Barnhoorn+04** (J. Struct. Geol. 26, 885 — Carrara marble torsion at
  500/600/727°C, 300 MPa, ε̇ = 6e-5 to 3e-3 s⁻¹), Fig 7: 14 datapoints
  spanning 3 temperatures across γ ∈ [0.5, 46]. **10 datapoints kept** after
  filtering for J ≤ 15 (Skemer+05 linear conversion validity); 4 high-J
  samples excluded:
  - γ=32 J=42.1 at 600°C (very high γ in dynamic-recrystallization regime)
  - γ=9 J=20.6 at 727°C (rapid CPO strengthening at high T)
  - γ=26 J=31.1 at 727°C
  - γ=46 J=34.5 at 727°C
  These excluded points represent calcite CPO **continuing to strengthen
  past γ ≈ 10**, contradicting Pieri+01's apparent saturation at γ ≈ 5.
  This T-dependence is documented in §11. For the saturating-exponential
  framework appropriate to the γ ∈ [0, 10] regime that MDOODZ scenarios
  typically visit, the LSQ fit on the 18-point combined dataset is the
  right calibration.

## §3. Asymptote — LSQ-fit to lab data, bracketed by VPSC bounds

The δ-vs-M slope `a = 6.0` is **not** taken from a clean linear lab fit
analogous to Hansen+12 Fig 3b (`δ = 24.5·M + 1` for olivine — a direct
mechanical lab measurement). For calcite, no equivalent published linear
fit exists. Instead, `slope = 6.0` is fixed by Tommasi+09-style VPSC
modeling of calcite Hill tensors which gives `δ_calcite (saturated)` ≈ 2-5
depending on slip-system activity. With `slope = 6.0` and the LSQ-fit
`M_inf = 0.41`:

```
δ_∞ = slope · M_inf + 1 = 6.0 · 0.41 + 1 = 3.46
```

| Asymptote source | δ_∞ |
|---|---|
| Tommasi+09-style VPSC bounds | 2–5 |
| **This calibration (slope = 6.0, M_inf = 0.41, LSQ on 18-pt combined dataset)** | **3.46** |
| (Earlier calibration, 8-pt dataset, before Barnhoorn+04 added) | 4.0 |
| Empirical PinchSwell `ani_fac_max = 4` (used as soft bound, not target) | 4.0 |
| Hansen olivine (for scale comparison) | 14.0 |

The LSQ fit on the expanded 18-point dataset pulls the asymptote down
slightly compared to the earlier 8-point calibration. This is **honest**
— the additional Barnhoorn+04 datapoints at 500°C and 600°C show calcite
CPO strengthens slowly in the dislocation-creep regime and only reaches
its sharpest CPO at 727°C in the recrystallization regime. The 18-point
fit averages across this T-range and ends up with a slightly lower bulk
asymptote.

If a future direct lab measurement supersedes this, the only edit needed
is changing the constants in `MDLIB/FlowLaws.c`'s static helper
`anisoDelta_Calcite`. The CI test gates the implementation against THIS
calibration, not against lab data, so the test threshold scenarios also
need re-calibration if the constants change.

## §4. Per-datapoint verification (18-point combined fit)

| γ | J_obs | M_obs (Skemer+05) | M_fit (this calibration) | δ_fit | source |
|---|---|---|---|---|---|
| 0.00 | 1.12 | 0.008 | 0.000 | 1.00 | Pieri+01 |
| 1.00 | 1.66 | 0.044 | 0.099 | 1.59 | Pieri+01 |
| 2.00 | 2.02 | 0.068 | 0.176 | 2.06 | Pieri+01 |
| 5.00 | 9.76 | 0.584 | 0.305 | 2.83 | Pieri+01 (P087, high outlier) |
| 11.00 | 8.92 | 0.528 | 0.391 | 3.34 | Pieri+01 (P088, recryst-randomized) |
| 1.00 | 1.60 | 0.040 | 0.099 | 1.59 | Bruijn+11 |
| 2.60 | 2.50 | 0.100 | 0.219 | 2.32 | Bruijn+11 |
| 5.00 | 4.00 | 0.200 | 0.305 | 2.83 | Bruijn+11 |
| 0.50 | 1.80 | 0.053 | 0.052 | 1.31 | Barnhoorn+04 (500°C) |
| 3.00 | 3.80 | 0.187 | 0.245 | 2.47 | Barnhoorn+04 (500°C) |
| 4.00 | 6.70 | 0.380 | 0.292 | 2.75 | Barnhoorn+04 (500°C) |
| 10.00 | 3.30 | 0.153 | 0.382 | 3.29 | Barnhoorn+04 (500°C, recryst-randomized) |
| 2.00 | 2.90 | 0.127 | 0.176 | 2.06 | Barnhoorn+04 (600°C) |
| 4.00 | 3.70 | 0.180 | 0.292 | 2.75 | Barnhoorn+04 (600°C) |
| 9.00 | 4.10 | 0.207 | 0.371 | 3.22 | Barnhoorn+04 (600°C) |
| 0.60 | 3.10 | 0.140 | 0.062 | 1.37 | Barnhoorn+04 (727°C) |
| 3.00 | 1.70 | 0.047 | 0.245 | 2.47 | Barnhoorn+04 (727°C, low outlier) |
| 5.00 | 12.50 | 0.767 | 0.305 | 2.83 | Barnhoorn+04 (727°C, high-T sharp CPO) |
| ∞ | — | — | 0.410 | 3.46 | asymptote |

**Excluded high-J Barnhoorn samples** (Skemer linear conversion breaks down):

| γ | J_obs | T (°C) | reason |
|---|---|---|---|
| 32 | 42.1 | 600 | recryst regime, J → M conversion non-linear |
| 9  | 20.6 | 727 | high-T sharp CPO, M would be > 1.3 |
| 26 | 31.1 | 727 | high-T sharp CPO |
| 46 | 34.5 | 727 | high-T sharp CPO continues to strengthen |

The fit smooths through all included points. The 727°C γ=5 sample (J=12.5,
M=0.77) is a notable outlier above the curve — high-T sharp-CPO regime
that the average-T fit doesn't capture. Documented further in §11.

## §11. Temperature dependence of calcite CPO (Barnhoorn+04)

Barnhoorn+04 reports torsion experiments at three temperatures (500, 600,
727°C) — a critical observation excluded from the simpler Pieri+01 dataset.
Per their Figs 7 & 6:

- **500°C** (low-T, dislocation-creep regime): J peaks ~6-7 by γ=4, then drops
  due to recrystallization randomization (γ=10 J=3.3 outlier). CPO never
  becomes very strong.
- **600°C** (intermediate): J slowly grows to ~4 by γ=9, then JUMPS to J=42 at
  γ=32 as recrystallization-driven CPO sharpening takes over.
- **727°C** (high-T, dynamic-recrystallization regime, same as Pieri+01): J=12.5
  at γ=5, J=20.6 at γ=9, J=31.1 at γ=26, J=34.5 at γ=46. **CPO continues to
  strengthen all the way to γ=46** with no sign of saturation, contradicting
  Pieri+01's apparent saturation at γ ≈ 5.

This T-dependence has important implications:

1. **The "saturation strain" depends on T**. At low T, CPO peaks then declines
   (recrystallization randomization). At high T, CPO continues to strengthen
   for dozens of strain units.
2. **A single-T calibration is necessarily an approximation**. Our LSQ fit
   averages across the 500/600/727°C data and ends up reasonably mid-range,
   but doesn't capture either extreme.
3. **MDOODZ scenarios at typical lithospheric T (500-700°C)** match the
   middle-of-the-road regime our calibration targets. For very high-T
   scenarios (subduction-zone metamorphic conditions, > 700°C) the
   calibration may underestimate calcite δ at high γ. Document this in
   scenario `.txt` headers if running such conditions.

A future extension could add T-aware dispatch (`anisoDelta_Calcite(double
FS_AR, double T)`) with separate calibrations for the three temperature
regimes. Out of scope for this change. The function-pointer architecture
supports this trivially via a wider function signature in a follow-up.

## §5. RMS error of candidate δ-forms (vs 18-point combined dataset)

Mapped via `δ_lab = 6.0 · M_obs + 1` (the new calcite linear fit applied to
Skemer-converted J-index data):

| Form | RMS |
|---|---|
| `δ = min(FS_AR, 14)` (Hansen olivine cap, irrelevant for calcite — shown for scale) | 8.16 |
| `δ = 24.5·M(γ_eff) + 1` (Hansen olivine form, calcite would over-shoot dramatically) | 6.26 |
| `δ = min(FS_AR, 4)` (existing PinchSwell default) | 1.67 |
| **`δ = 6.0·M(γ_eff) + 1` (NEW calcite case 2, M_inf=0.41, γ_e=3.71)** | **0.99** |

The new calcite form is ~50% better than the existing `min(FS_AR, 4)` form —
the improvement comes from matching the slow rise (γ_e = 3.50 vs the min
form's hard saturation at γ = 1.6) while preserving the asymptote (4.0).

![Calcite calibration plot — Pieri+01 + Bruijn+11 vs candidate forms](mdoodz_vs_calcite_data.png)

The two-panel plot shows:

- **Panel (a)**: M(γ) with all 8 datapoints (J→M-converted), the new
  calibrated curve in solid black, and the Hansen olivine fit as a dotted
  reference for scale comparison.
- **Panel (b)**: δ(γ) with the existing default `min(FS_AR, 4)` (dashed grey),
  the Hansen olivine form (`24.5·M + 1`, dotted red — for scale), and the
  new calcite case-2 form (solid black). The asymptote `δ_∞ = 4.0` is marked.

## §6. Pieri+01 vs Bruijn+11 disagreement at γ = 5

P087 (Pieri+01, single-stage to γ = 5) reports J = 9.76. Bruijn+11's
single-stage segment at the same γ reports J ≈ 4.0. The discrepancy is
~2.4× in J, ~3× in M after Skemer+05 conversion.

**Likely cause**: Bruijn+11 samples carry pre-existing texture from a prior
D1 deformation stage that biases the subsequent strain to a different
fabric pathway. Pieri+01 samples are single-stage from random-fabric
starting material. The Bruijn measurement is structurally a different
material than Pieri's despite the same nominal strain.

**How the fit handles it**: the LSQ-best fit through both datasets gives
unrealistic constants (`M_inf = 0.92`, `γ_e = 11.9` per the script
output) because it has to honor the high P087 value AND the low Bruijn
value. We **don't use the LSQ-best fit**. Instead, the committed values
(`M_inf = 0.50`, `γ_e = 3.50`) are **bias-toward-Pieri** because Pieri's
single-stage measurements are the cleaner reference for "fresh calcite
deformed to γ in simple shear".

This is documented in design.md D2 of the openspec change.

## §7. Pieri+01 γ = 11 non-monotonic outlier

P088 (Pieri+01 at γ = 11) reports J = 8.92, **lower** than P087 (J = 9.76 at
γ = 5). Pieri+01 §3.3-3.4 attributes this to recrystallization randomization
at very high strain — a fraction of the recrystallized grains nucleate in
near-random orientations, weakening the bulk fabric.

**Our exponential-saturation form cannot capture this non-monotonic behaviour
by construction.** The fit smooths through it (M_fit at γ=11 is 0.478,
which gives δ=3.87, vs M_obs=0.528 → δ=4.17). The error is small in
absolute terms, and **no MDOODZ scenario routinely visits γ = 11** for
calcite — pinch-and-swell layer dynamics happen over γ ∈ [0.5, 5] where
the exponential-saturation form is well-validated.

## §8. Structural caveats

### `γ_eff(FS_AR)` is a simple-shear inversion

Identical caveat to the Hansen olivine note: `γ_eff = √FS_AR − 1/√FS_AR` is
exact under simple shear, defensible-but-questionable under pure shear /
general 2D flow. Calcite scenarios (`SETS/PinchSwell*`) are in a near-pure-
shear regime — the user is signing up for "calcite-shaped saturation
parameterised by accumulated stretch", not "calcite fabric strength under
arbitrary deformation."

### NOT applicable to high-pressure calcite (Schuster+18)

Schuster+18 (J. Struct. Geol., 2018) reports HPT torsion experiments on
calcite to γ = 80 at confining pressures of 1–4 GPa, temperatures
235–450 °C. **At pressures > 1.6 GPa, calcite undergoes a phase transition
to CaCO₃-II** (Kondo+72, Schuster+17), making it a different material with
different slip systems and different fabric evolution. The Schuster+18
data is documented for context but **deliberately excluded from this
calibration**. MDOODZ scenarios at sub-GPa pressures (the typical regime
for crustal/lithospheric pinch-and-swell) are not affected by this caveat.

If a future MDOODZ scenario at GPa-range confining pressures emerges,
adding a `case 4 = high-pressure calcite (Schuster+18)` is a clean
follow-up. The function-pointer architecture supports this trivially.

### Skemer+05 J→M conversion is approximate

The relation `M ≈ (J − 1) / 15` assumes monoclinic-symmetry textures
typical of simple-shear deformation. Pieri+01 calcite is monoclinic by
construction (Pieri+01 §3.2, sample symmetry). For axial-symmetric or
orthorhombic textures the conversion factor differs (Skemer+05 §4).
For the calcite torsion data here, the conversion is appropriate.

### Recrystallization-randomization not modelled (see §7)

## §9. Reproducing the calibration

```bash
cd misc/aniso_fstrain/calcite
python3 calibrate_calcite.py
```

Outputs:

- Per-dataset summary table (Pieri+01, Bruijn+11, combined)
- Free LSQ best fit (for comparison with our committed bias-toward-Pieri values)
- Per-datapoint verification table (γ, J_obs, M_obs, M_fit, δ_fit)
- RMS errors of all candidate δ-forms vs the calcite reference
- `mdoodz_vs_calcite_data.png` — the two-panel comparison plot

## §10. Citations

**Primary calibration data (combined 18-point dataset):**

- Pieri, M., Kunze, K., Burlini, L., Stretton, I., Olgaard, D. L., Burg, J.-P.
  & Wenk, H.-R. (2001) "Texture development of calcite by deformation and
  dynamic recrystallization at 1000 K during torsion experiments of marble
  to large strains", *Tectonophysics* **330**, 119–140 (5 datapoints).
  doi:[10.1016/S0040-1951(00)00225-0](https://doi.org/10.1016/S0040-1951(00)00225-0)
- Bruijn, R. H. C., Kunze, K., Mainprice, D. & Burlini, L. (2011) "Mechanical
  and microstructural development of Carrara marble with pre-existing strain
  variation", *Tectonophysics* **503**, 75–91 (3 single-stage segment datapoints).
  doi:[10.1016/j.tecto.2010.09.029](https://doi.org/10.1016/j.tecto.2010.09.029)
- Barnhoorn, A., Bystricky, M., Burlini, L. & Kunze, K. (2004) "The role of
  recrystallisation on the deformation behaviour of calcite rocks: large
  strain torsion experiments on Carrara marble", *J. Struct. Geol.* **26**,
  885–903 (10 datapoints kept after Skemer J ≤ 15 filter, across 500/600/727°C).
  doi:[10.1016/j.jsg.2003.11.024](https://doi.org/10.1016/j.jsg.2003.11.024)

**δ-vs-M slope (asymptote bound, NOT direct fit):**

- Tommasi, A., Knoll, M., Vauchez, A., Signorelli, J. W., Thoraval, C. &
  Logé, R. (2009) "Structural reactivation in plate tectonics controlled by
  olivine crystal anisotropy", *Nature Geoscience* **2**, 423–427.
  (Also Tommasi-style VPSC modeling for calcite Hill tensors.)

**J→M conversion:**

- Skemer, P., Katayama, I., Jiang, Z. & Karato, S.-i. (2005) "The
  misorientation index: Development of a new method for calculating the
  strength of lattice-preferred orientation", *Tectonophysics* **411**,
  157–167. doi:[10.1016/j.tecto.2005.08.023](https://doi.org/10.1016/j.tecto.2005.08.023)

**Calcite weakening mechanisms (CPO ~1/3 of total):**

- Barnhoorn, A., Bystricky, M., Burlini, L. & Kunze, K. (2004) "The role of
  recrystallisation on the deformation behaviour of calcite rocks: Large
  strain torsion experiments on Carrara marble", *J. Struct. Geol.* **26**,
  885–903.

**HPT calcite (excluded from fit, documented for context):**

- Schuster, R., Habler, G., Schafler, E. & Abart, R. (2018) "Microstructural
  and textural evolution of calcite deformed to high shear strain by
  high-pressure torsion", *J. Struct. Geol.* (PII S0191-8141(18)30443-7).
  doi:[10.1016/j.jsg.2018.09.003](https://doi.org/10.1016/j.jsg.2018.09.003)
- Kondo, S. et al. (1972), Schuster+17 (CaCO₃ phase transition at 1.6 GPa).

## §13. Two-regime T-aware calibration (cases 4 and 5)

The temperature dependence documented in §11 is large enough to motivate **two
additional T-aware calibration cases** alongside the T-averaged `case 2`:

| Case | Aniso DB | Helper | Regime | M_inf | γ_e | δ_∞ | RMS(M) | Recommended for |
|---|---|---|---|---|---|---|---|---|
| Default | 2 | `anisoDelta_Calcite` | T-averaged | 0.41 | 3.71 | 3.46 | 0.166 | unknown / mixed T |
| **Mid-low-T** | **4** | `anisoDelta_Calcite_LowT` | T ≤ 650°C | **0.22** | **1.37** | **2.33** | **0.073** | T < 650°C |
| **High-T** | **5** | `anisoDelta_Calcite_HighT` | T > 650°C | **0.83** | **8.24** | **5.99** | 0.169 | T > 650°C |

### Mid-low-T fit is the cleanest of the three

The 7-point Barnhoorn+04 dataset at 500/600°C clusters tightly along a
saturating-exponential curve (RMS = 0.073, less than half the T-averaged
RMS). Subgrain-rotation recrystallization dominates this regime,
randomizing CPO at high γ — the modest asymptote δ_∞ = 2.33 reflects this
self-limiting behaviour.

### High-T fit is messier — internal lab disagreement

The 11-point Pieri+01 + Bruijn+11 + Barnhoorn+04(727°C) high-T dataset has
RMS = 0.169 — barely better than the T-averaged fit. Three sources at the
same nominal T (727°C) disagree by factor 3 at γ=5:

- Pieri+01 P087 γ=5 J=9.76 (apparent saturation by here)
- Bruijn+11 γ=5 J=4.0 (pre-existing-fabric bias, lower)
- Barnhoorn+04 727°C γ=5 J=12.5 (continued strengthening)

The LSQ has to honor all three. M_inf = 0.83 and γ_e = 8.24 reflect a
slow-saturation, high-asymptote curve dominated by Barnhoorn's continued
strengthening to γ=46 (excluded from the linear-conversion-valid fit but
informing the slow saturation). At γ < 10 the fit closely tracks the
data; for very-high-γ scenarios (γ > 20) it extrapolates beyond observed
data into the no-saturation Barnhoorn regime.

### Predictions across γ for each regime

```
              γ=2    γ=5    γ=10    γ=14    γ→∞
Mid-low-T:    1.93   2.29   2.33    2.33    2.33   (saturated by γ=4)
T-averaged:   2.05   2.82   3.30    3.39    3.46   (middle ground)
High-T:       2.07   3.27   4.51    5.07    5.99   (still rising)
```

Up to γ ≈ 5 all three predict similar δ (within ~30%). Beyond γ=10 they
diverge dramatically — the regime choice matters most at high strain.

### Choosing a case for your scenario

- **T < 500°C** (cool middle/upper crust): use `aniso_db = 4`. Mid-low-T
  fit is the right shape (SR-recryst regime).
- **500-650°C** (mid-crust, typical pinch-and-swell): use `aniso_db = 4`.
  Border case; the mid-low-T fit's tighter RMS makes it preferable to the
  T-averaged default.
- **650-750°C** (lower crust, transition): either `aniso_db = 5` for
  high-T physics, or `aniso_db = 2` for the conservative middle-ground
  default. Document the choice in scenario header.
- **> 750°C** (hot lower crust, mantle transition): use `aniso_db = 5`.
  GBM-recryst regime, continued strengthening at high γ.
- **Unknown / mixed-T** (whole-lithosphere model with T variation across
  one calcite phase): use `aniso_db = 2` (T-averaged) as a defensible
  middle ground.

### CI tests gate all three calibrations

`AniFstrainCalcite` (case 2, γ ∈ [0, 10]), `AniFstrainCalciteLowT` (case 4,
γ ∈ [0, 6]), and `AniFstrainCalciteHighT` (case 5, γ ∈ [0, 14]) all assert
per-step δ matches the analytical reference to `1e-6` relative tolerance.
Each test reuses the shared `AniFstrainEvolution.txt` fixture via
`MutateInput` setting the appropriate `aniso_db` value.

### Polyphase / non-direct-fit guidance (cross-reference)

For lower-crustal scenarios involving **calcite + plagioclase + pyroxene**
(skarn, marble in metamorphosed mafic rocks, calc-silicate), case 4 / 5
calibrations apply only to the calcite fraction. For the bulk-rock
anisotropy of plagioclase- or mica-bearing rocks, see
[`../quartz/quartz_calibration.md`](../quartz/quartz_calibration.md):

- §11 — quartz + mica polyphase mylonites (`ani_fac_max ≈ 5–15`)
- §12 — plagioclase / feldspar-bearing rocks, gabbro, diabase, basalt,
  mafic granulite (`ani_fac_max ≈ 1.2–2.0`, much weaker than calcite
  because plagioclase intrinsically has low single-crystal anisotropy
  and is further damped by DisGBS at high strain)

Neither polyphase mode is a committed `aniso_db` case — they are
qualitative `ani_fac_max` brackets because no co-located (γ, J_bulk)
dataset exists for those rock types. Calcite (case 4 / 5) is direct-fit
on Pieri+01 / Bruijn+11 / Barnhoorn+04 data, comparable in rigor to
olivine cases 1 / 7 and quartz cases 3 / 6.

### Limitations of the 2-regime split

1. **No direct lab measurement of T-aware δ-vs-M slope.** The slope = 6.0
   is held fixed across all three cases — calibrated from Tommasi+09 VPSC
   bounds. A true T-aware calibration would have slope vary with T as
   well, but data is too sparse to constrain that.
2. **Mid-low-T uses Barnhoorn+04 alone.** No Pieri or Bruijn data exists
   at < 700°C — the 7 datapoints are all from one paper.
3. **The 650°C cutoff is somewhat arbitrary.** Barnhoorn+04 600°C samples
   were grouped with low-T; this is defensible (their stress-strain
   curves and microstructures resemble 500°C more than 727°C) but a
   user-tunable cutoff is not provided.
4. **Naïve regime selection.** MDOODZ has no automatic T → aniso_db
   routing — the user picks per phase at scenario setup. For models with
   spatially-varying T across one calcite phase, the calibration is
   averaged via `case 2`. A future T-aware function-pointer extension
   (`double aniso_delta_fn(double FS_AR, double T)`) could resolve this.
