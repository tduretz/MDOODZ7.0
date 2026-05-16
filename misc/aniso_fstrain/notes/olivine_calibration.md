# Olivine viscous-anisotropy calibration (`aniso_db = 1` and `aniso_db = 7`)

This note documents the **two** calibrated olivine δ(FS_AR) forms delivered
as `case 1` and `case 7` of `ReadDataAnisotropy`, dispatched via
`AnisoFactorEvolv` when the per-phase parameters `ani_fstrain == 2` and
`aniso_db ∈ {1, 7}`.

| Case | Source | Conditions | `(M_inf, γ_e, slope)` | δ_∞ |
|---|---|---|---|---|
| 1 | Hansen+12/+14/+16 (n=38) | dry Fo90, lab torsion, fresh CPO | (0.536, 3.96, 24.5) | **14.13** |
| 7 | Tasaka+Boneh+Kumamoto+Bernard (n=93) | wet / pre-CPO / natural | (0.20, 1.30, 24.5) | **5.90** |

Both modes use the same δ-vs-M slope (24.5) because Hansen+12 Fig 3b's
linear M↔δ relation is universal across olivine fabric strength — set by
single-crystal anisotropy of olivine, not by M_inf magnitude. The two
modes differ only in the (M_inf, γ_e) saturation parameters.

## §1. Why two modes?

The olivine literature contains **two distinct populations** of γ-vs-M data:

**Population A — fresh-CPO lab experiments** (clean saturating exponential, M_inf ≈ 0.5):
- Hansen, Zhao, Zimmerman & Kohlstedt (2014) EPSL 387:157 — 26 Paterson torsion samples, dry Fo90, γ ∈ [0, 18.7]
- Hansen, Warren, Zimmerman & Kohlstedt (2016) Part 1, EPSL 445:92 — 12 additional samples
- Bystricky, Kunze, Burlini & Burg (2000) Science 290:1564 — 4 lab torsion samples (also dry Fo90, consistent with Hansen)

**Population B — damped CPO** (M_inf < 0.30, no clean saturating shape):
- Tasaka, Zimmerman & Kohlstedt (2016) JGR-SE 121:92 — 12 wet Fo50 lab torsion samples
- Boneh & Skemer (2014) EPSL 406:213 — 16 Aheim dunite samples with **pre-existing CPO**
- Kumamoto, Warren & Hansen (2019) JGR-SE 124:12763 — 30 natural Josephine Peridotite shear-zone samples
- Bernard, Behr, Becker & Young (2019) G3 20:3469 — 35 natural xenoliths with X/Z aspect ratios

The two populations share the same physics (olivine slip systems, single-
crystal anisotropy) but differ in how cleanly CPO grows with strain:

| Mechanism damping CPO growth (Population B) | Source |
|---|---|
| Hydrous + Fe-rich olivine: water + Fe shifts active slip systems | Tasaka 2016 |
| Pre-existing CPO disruption: starting fabric in unfavorable orientation can DECREASE M before new CPO grows | Boneh & Skemer 2014 |
| Complex natural strain history: multiple deformation events, polyphase strain, refoliation | Kumamoto+19, Bernard+19 |

Three independent mechanisms, but all converge on the same observational
signature: **M caps at 0.20–0.30 across all γ**, even at γ = 20+ in
Kumamoto's natural samples.

![Olivine M-vs-γ survey across all 6 datasets](olivine_all_datasets.png)

The plot shows clear separation: Hansen + Bystricky (Population A) trace
the solid line up to M ≈ 0.65 at γ = 10+; the four Population B datasets
stay below M = 0.30 across the full strain range.

## §2. Functional form

Both cases share the same parametric form:

```
δ(FS_AR) = min( 24.5 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )

  where  γ_eff(FS_AR) = √FS_AR − 1/√FS_AR              (simple-shear inversion)
         M(γ)         = M_inf · (1 − exp(−γ / γ_e))    (saturating exponential)

         case 1 (fresh):    M_inf = 0.536, γ_e = 3.96
         case 7 (damped):   M_inf = 0.20,  γ_e = 1.30
```

The Hansen+12 Fig 3b slope of 24.5 in `δ = 24.5·M + 1` is shared between
both cases — it reflects the intrinsic single-crystal anisotropy of
olivine and is independent of CPO strength.

## §3. Strain-measure conversions for Population B

Since the four Population B datasets report different strain measures,
they are converted to a "simple-shear-equivalent" γ for the combined
LSQ fit:

| Dataset | Reported strain | Conversion to γ_eq |
|---|---|---|
| Tasaka 2016 | γ (lab simple shear) | **identity** (already simple-shear γ) |
| Boneh & Skemer 2014 | ε (axial natural strain, Griggs) | γ_eq = 2·sinh(ε/2)  (≈ ε for small ε) |
| Kumamoto+19 | γ (foliation deflection in shear zone) | **identity** |
| Bernard+19 | X/Z aspect ratio (3D ellipsoid) | γ_eq = √(X/Z) − 1/√(X/Z) (plane-strain inversion) |

These conversions are documented in
[`survey_olivine_datasets.py`](survey_olivine_datasets.py).

## §4. Calibration procedure

The script [`survey_olivine_datasets.py`](survey_olivine_datasets.py)
loads all six datasets, applies the strain conversions, runs a free LSQ
fit on each population separately, and renders the comparison plot.

**Free LSQ results:**

| Population | n | M_inf | γ_e | RMS(M) | δ_∞ |
|---|---|---|---|---|---|
| A (Hansen + Bystricky) | 42 | 0.533 | 3.72 | 0.078 | 14.06 |
| B (Tasaka + Boneh + Kumamoto + Bernard) | 93 | 0.200 | 1.27 | 0.071 | 5.90 |

The Population A free LSQ is consistent with the existing `case 1` commit
(0.536, 3.96, 14.13) — adding Bystricky's 4 points doesn't shift the
fit meaningfully because they fall on the existing trend.

The Population B free LSQ is committed to `case 7` (rounded to 0.20, 1.30
for clean numerics; RMS(M) ≈ 0.073 either way).

## §5. Asymptote — physical interpretation

Both asymptotes are derived from the SAME slope (24.5):

```
case 1:   δ_∞ = 24.5 · 0.536 + 1 = 14.13
case 7:   δ_∞ = 24.5 · 0.20  + 1 =  5.90
```

| Asymptote source | δ_∞ |
|---|---|
| Single-crystal olivine bound (Tommasi+09 VPSC) | ~20 |
| **case 1 fresh CPO (lab Fo90)** | **14.13** |
| Hansen+16 Part 1 stress-aware direct calc | 14.02 |
| **case 7 damped CPO (natural / wet)** | **5.90** |
| Calcite high-T (case 5) | 5.99 — coincidentally close |
| Quartz cases 3 / 6 | ~3.2 |

The damped olivine's δ_∞ ≈ 6 is roughly **half** of fresh-CPO olivine's
14.13, reflecting the ~50% lower M_inf observed in natural / pre-CPO /
wet samples. It is also coincidentally close to calcite high-T's 5.99 —
suggesting that for natural lithospheric scenarios, olivine and calcite
have comparable viscous anisotropy magnitudes.

## §6. Choosing between case 1 and case 7

| Scenario | Recommended case | Why |
|---|---|---|
| Asthenosphere, MOR, post-melting fresh peridotite | **case 1** | Fresh CPO grows monotonically, lab calibration applies |
| Lithospheric mantle xenoliths, mantle-shear-zone mylonites | **case 7** | Pre-existing CPO + multiple deformation events damp fabric |
| Subduction zones (wet, water-bearing slab/wedge) | **case 7** | Wet conditions activate alternate slip systems, damping CPO |
| Lab simple-shear comparison (validating against torsion) | **case 1** | Direct match to Hansen lab data |

**Default behavior**: `ani_fstrain = 2 && aniso_db == 0` auto-routes to
`aniso_db = 1` (Hansen fresh CPO) — this is the **olivine-specific
default** preserved from the original implementation. Users wanting the
damped mode MUST set `aniso_db = 7` explicitly per phase.

## §7. Per-population dataset details

### §7.1 Hansen+12/+14/+16 (case 1, n=38)

All samples are dry Fo90, Paterson torsion at 1200°C, 300 MPa. Hansen+14
Table 1 reanalyzes 26 previously-published samples (Bystricky 2000, Zhang
& Karato 2005, Hansen+12a/b/c) using consistent Skemer-2005 M-index
methodology. Hansen+16 Part 1 adds 12 datapoints from 5 new samples
(PT0716, PT0718, PT0751, PT0752, PT0754).

### §7.2 Bystricky 2000 (case 1, n=4)

Original Paterson torsion experiments, San Carlos olivine, 1200°C, 300 MPa,
dry, γ ∈ [0.5, 5]. ODF J-index reported per sample (Fig 4). Converted to
M via Skemer's linear relation `M ≈ (J−1)/15`. Already incorporated into
Hansen+14 Table 1, so listed separately here only for completeness.

### §7.3 Tasaka 2016 (case 7, n=12)

Wet Fo50 (iron-rich) torsion, 1200°C, 300 MPa, γ ∈ [0.6, 5]. Both M and J
reported per sample (Fig 5). Damping mechanism: water + Fe content shifts
active slip systems and slows CPO development relative to dry Fo90.

### §7.4 Boneh & Skemer 2014 (case 7, n=16)

Aheim dunite (with strong starting CPO M_initial = 0.13) deformed in
triaxial Griggs apparatus at 1 GPa, 1200°C. Three sample geometries
(perpendicular, oblique, parallel to pre-existing foliation) → different
evolution. Strains ε ∈ [0, 0.72] (small, lab-limited). Damping mechanism:
pre-existing CPO is disrupted before new CPO grows.

### §7.5 Kumamoto+19 (case 7, n=30)

Three natural Josephine Peridotite shear zones (SZA, SZG, SZP). Strain
γ from foliation deflection (Ramsay & Graham 1970). γ ∈ [0, 65]; M ∈
[0.03, 0.34]. Paper §4.1: *"In all three shear zones, the texture
strength does not increase systematically with increasing strain."*

### §7.6 Bernard+19 (case 7, n=35)

Compilation of natural peridotite xenoliths and continental massifs. We
use the 35 Table 2 entries with X/Z aspect ratios reported (samples
without aspect ratios are excluded since strain magnitude can't be
estimated). γ_eq derived from X/Z under plane-strain assumption. Paper
§5.2: *"There did not appear to be any strong correlation between olivine
CPO type and fabric strength."*

## §8. Datasets we checked but did not include

- **Zhang & Karato (1995) Nature 375:774** — 13 simple-shear olivine
  experiments at 1473–1573 K, but the paper reports only **angle of [100]
  peak vs flow direction** (kinematic), not J or M (intensity). Not
  usable for our γ-vs-J calibration.
- **Ohuchi+17 (Contrib Min Petrol 172:65)** — high-P hydrous olivine,
  focus on rheology / flow behavior, no per-sample (γ, J, M) data.
- **Skemer & Hansen 2016 (Tectonophys 668:1)** — review paper, cites the
  primary datasets above.
- **Bernard+19 supporting info Table S5** — additional 445 literature
  samples without uniform strain reporting; aspect ratios available for a
  subset but not consistently across the compilation.

### Polyphase crustal rocks (cross-reference)

Olivine cases 1 and 7 both assume an essentially pure-olivine aggregate
(>95% olivine, as in lab dunites or natural peridotites). For crustal-
scale models that involve **other rock types**, see
[`quartz_calibration.md`](../quartz/quartz_calibration.md):

- §11 — quartzofeldspathic + mica mylonites (granitoid mylonite,
  schist, phyllite). `ani_fac_max ≈ 5–15` depending on mica fraction.
- §12 — plagioclase / feldspar-bearing rocks (gabbro, anorthosite,
  diabase, basalt, mafic granulite). `ani_fac_max ≈ 1.2–2.0` because
  plagioclase has intrinsically modest single-crystal anisotropy AND
  is further damped in nature by DisGBS / phase mixing.

We deliberately did **not** commit polyphase modes for either case
because no single dataset reports (γ, J_bulk) pairs in the form a
saturating exponential M(γ) = M_inf·(1−exp(−γ/γ_e)) requires. For
plagioclase specifically, the form doesn't apply at all because the
deformation mechanism transitions from dislocation creep to DisGBS at
high strain, and CPO can WEAKEN rather than saturate. The olivine
case-7 calibration, by contrast, IS a true LSQ fit on 93 co-located
(γ, M) samples and is fully committed.

## §9. Reproducing the calibration

```bash
cd misc/aniso_fstrain/olivine_extended
python3 survey_olivine_datasets.py            # both fits + comparison plot

cd cmake-build-test/TESTS
./AnisotropyBenchmarkTests --gtest_filter='*Olivine*'
```

## §10. Citations

- **Hansen, L. N., Zhao, Y.-H., Zimmerman, M. E. & Kohlstedt, D. L.** (2014)
  "Protracted fabric evolution in olivine: Implications for the
  relationship among strain, crystallographic fabric, and seismic
  anisotropy", *Earth Planet. Sci. Lett.* **387**, 157–168.
  doi:[10.1016/j.epsl.2013.11.009](https://doi.org/10.1016/j.epsl.2013.11.009).

- **Hansen, L. N., Warren, J. M., Zimmerman, M. E. & Kohlstedt, D. L.** (2016)
  "Viscous anisotropy of textured olivine aggregates, Part 1: Measurement
  of the magnitude and evolution of anisotropy", *Earth Planet. Sci. Lett.*
  **445**, 92–103. doi:[10.1016/j.epsl.2016.04.008](https://doi.org/10.1016/j.epsl.2016.04.008).

- **Bystricky, M., Kunze, K., Burlini, L. & Burg, J.-P.** (2000) "High shear
  strain of olivine aggregates: Rheological and seismic consequences",
  *Science* **290**, 1564–1567.
  doi:[10.1126/science.290.5496.1564](https://doi.org/10.1126/science.290.5496.1564).

- **Tasaka, M., Zimmerman, M. E. & Kohlstedt, D. L.** (2016) "Evolution of
  the rheological and microstructural properties of olivine aggregates
  during dislocation creep under hydrous conditions", *J. Geophys. Res.
  Solid Earth* **121**, 92–113.
  doi:[10.1002/2015JB012134](https://doi.org/10.1002/2015JB012134).

- **Boneh, Y. & Skemer, P.** (2014) "The effect of deformation history on
  the evolution of olivine CPO", *Earth Planet. Sci. Lett.* **406**,
  213–222. doi:[10.1016/j.epsl.2014.09.018](https://doi.org/10.1016/j.epsl.2014.09.018).

- **Kumamoto, K. M., Warren, J. M. & Hansen, L. N.** (2019) "Evolution of
  the Josephine Peridotite Shear Zones: 2. Influences on Olivine CPO
  Evolution", *J. Geophys. Res. Solid Earth* **124**, 12763–12781.
  doi:[10.1029/2019JB017968](https://doi.org/10.1029/2019JB017968).

- **Bernard, R. E., Behr, W. M., Becker, T. W. & Young, D. J.** (2019)
  "Relationships Between Olivine CPO and Deformation Parameters in
  Naturally Deformed Rocks and Implications for Mantle Seismic
  Anisotropy", *Geochem. Geophys. Geosyst.* **20**, 3469–3494.
  doi:[10.1029/2019GC008289](https://doi.org/10.1029/2019GC008289).

- **Skemer, P., Katayama, I., Jiang, Z. & Karato, S.-i.** (2005) "The
  misorientation index: development of a new method for calibrating the
  strength of lattice-preferred orientation", *Tectonophysics* **411**,
  157–167. doi:[10.1016/j.tecto.2005.08.023](https://doi.org/10.1016/j.tecto.2005.08.023).
