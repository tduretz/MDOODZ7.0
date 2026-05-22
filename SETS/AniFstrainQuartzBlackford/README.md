# AniFstrainQuartzBlackford — Blackford+24 quartz calibration set

Homogeneous simple-shear MDOODZ run that traces δ(γ) under the Blackford+24
(Tectonics 43, e2023TC008166) Northern Snake Range quartzite calibration —
`aniso_db = 6` (case 6, `anisoDelta_Quartz_Blackford`, T = 500-650°C,
M_∞ = 0.20, γ_e = 4.00, slope = 11.0, asymptote δ_∞ = 3.20).

The 38-sample (γ_eq, pfJ → M) corpus is overlaid on the MDOODZ HDF5 output by
[mdoodz_vs_blackford.py](mdoodz_vs_blackford.py), producing the
datapoint-vs-MDOODZ plot in `img/`.

## Why slope = 11.0 (vs 3.0 for case 3)

Blackford+24 reports pole-figure J (pfJ ∈ [1.3, 4.85]), which is on a
different scale to Pennacchioni+10's ODF J (case 3): for the same physical
texture, pfJ < ODF J, so M = (J−1)/15 gives ~5× smaller M values.  The
δ-vs-M slope is rescaled to 11.0 (vs 3.0 for case 3) so the saturation
asymptote δ_∞ = 11.0·0.20 + 1 = 3.20 matches case 3's δ_∞ = 3.25 to within
rounding — both modes describe the same physical mineral, just measured on
different texture metrics.

## Why γ_e = 4.0 is a "centered defensible point"

Free LSQ on Blackford+24 alone produces a flat trough where M_∞ ≈ 0.17–0.29
across γ_e ∈ [3, 10] — saturation is poorly constrained because Blackford's
γ_eq spans only [1.28, 5.08], well below the saturation γ-range.  Committed
(M_∞ = 0.20, γ_e = 4.0) is centered in the LSQ envelope; the calibration is
honest about this caveat (see FlowLaws.c case-6 docstring); γ_e ≥ 4 is
generic for natural-quartz saturation regimes.

Unlike case 3 (Pennacchioni), case 6 IS data-driven on real per-sample MTEX —
σ/M_∞ = 0.191 is well within the real-data floor σ/M_∞ ∈ [0.13, 0.20].
The Blackford CSV is real-MTEX-grade.

## How to run

```bash
# 1. build (from MDOODZ7.0 repo root, after cmake configure)
cmake --build cmake-exec --target AniFstrainQuartzBlackford

# 2. run the simulation
cd cmake-exec/AniFstrainQuartzBlackford
mkdir -p output                          # writer_subfolder must exist
./AniFstrainQuartzBlackford              # writes 12 HDF5 outputs into output/

# 3. produce the comparison plot
cd <repo>/SETS/AniFstrainQuartzBlackford
python3 mdoodz_vs_blackford.py           # writes img/mdoodz_vs_blackford.png
```

Or use the convenience wrapper:

```bash
cd <repo>/SETS/AniFstrainQuartzBlackford
./run.sh                                 # build + run + plot in one shot
```

## What the plot shows

- 38 Blackford+24 datapoints (γ_eq ∈ [1.28, 5.08], pfJ → M via Skemer+05),
  coloured by Northern Snake Range strain domain (1-5)
- Calibrated case-6 saturating-exp curve δ(γ) = 11.0·0.20·(1 − exp(−γ/4.00)) + 1
- 11 MDOODZ HDF5 stars (per-step δ_mean) sampled at γ ∈ {1, 2, ..., 11}

The MDOODZ stars trace the case-6 saturating-exp curve essentially exactly
(rel err ~ 1e-6), confirming `aniso_db = 6` evaluates correctly. The
Blackford datapoints scatter naturally around the curve at σ/M_∞ ≈ 0.19 —
the legitimate real-MTEX scatter signature. Datapoints concentrated in
γ_eq ∈ [1.5, 5] cover the rising part of the curve; the saturation tail
(γ > 6) is shown by the curve + MDOODZ stars beyond Blackford's coverage.

## References

- **Blackford+24** Tectonics 43, e2023TC008166 — Northern Snake Range
  quartzite, 38 samples with 3D strain ellipsoids on detrital quartz
  porphyroclasts. PDF in
  [`~/bib/anisotropy_calibration_2026/quartz/Blackford_etal_2024_NSR_quartzite_3D_strain_Tectonics.pdf`](../../../bib/anisotropy_calibration_2026/quartz/Blackford_etal_2024_NSR_quartzite_3D_strain_Tectonics.pdf).
- Derived 38-point CSV at
  [`misc/aniso_fstrain/data/quartz/derived_Blackford24_38pts.csv`](../../misc/aniso_fstrain/data/quartz/derived_Blackford24_38pts.csv).
- Calibration script that fitted case 6:
  [`misc/aniso_fstrain/calibrate/calibrate_quartz.py`](../../misc/aniso_fstrain/calibrate/calibrate_quartz.py).
