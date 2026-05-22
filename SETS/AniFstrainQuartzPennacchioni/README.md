# AniFstrainQuartzPennacchioni — Pennacchioni+10 quartz calibration set

Homogeneous simple-shear MDOODZ run that traces δ(γ) under the Pennacchioni+10
(JGR 115, B12405) natural quartz shear-zone calibration —
`aniso_db = 3` (case 3, `anisoDelta_Quartz`, T ~ 500°C, M_∞ = 0.75,
γ_e = 3.50, slope = 3.0, asymptote δ_∞ = 3.25).

The 17-point (γ, J) corpus is overlaid on the MDOODZ HDF5 output by
[mdoodz_vs_pennacchioni.py](mdoodz_vs_pennacchioni.py), producing the
datapoint-vs-MDOODZ plot in `img/`.

## Audit caveat — case 3 is physical-prior, not data-driven

The 17 (γ, J) values used to fit case 3 are **NOT real per-sample MTEX data**.
Pennacchioni+10 reports no numerical J-index in tables, figures, or supplementary
material; the values were originally attributed to "(estimated from Fig 7)" but
Fig 7 is porphyroclast COIs, not a J-vs-γ plot.

Residual smoothness on the 17 points is σ/M_∞ = 0.045, **3–4× below the
real-data floor** seen in comparable corpora (Hansen lab combined 0.142,
Blackford+24 natural quartz 0.191, case-7 multi-source 0.430). Real per-sample
MTEX produces σ/M_∞ ∈ [0.13, 0.20]. The values are most likely a hand-fitted
smooth curve drawn to match Pennacchioni+10's qualitative WDV/MDV/SDV regime
descriptions. See
[Pennacchioni_etal_2010_audit.md](../../misc/aniso_fstrain/audit/Pennacchioni_etal_2010_audit.md)
for the full audit.

The case-3 calibration is still useful — M_∞ = 0.75, γ_e = 3.50 are consistent
with natural-quartz-mylonite expectations — but its defensibility is "calibrated
to match qualitative regime descriptions" rather than "fit to per-sample J-data".
For data-driven natural quartz, use case 6 (Blackford+24).

## How to run

```bash
# 1. build (from MDOODZ7.0 repo root, after cmake configure)
cmake --build cmake-exec --target AniFstrainQuartzPennacchioni

# 2. run the simulation
cd cmake-exec/AniFstrainQuartzPennacchioni
mkdir -p output                            # writer_subfolder must exist
./AniFstrainQuartzPennacchioni             # writes 17 HDF5 outputs into output/

# 3. produce the comparison plot
cd <repo>/SETS/AniFstrainQuartzPennacchioni
python3 mdoodz_vs_pennacchioni.py          # writes img/mdoodz_vs_pennacchioni.png
```

Or use the convenience wrapper:

```bash
cd <repo>/SETS/AniFstrainQuartzPennacchioni
./run.sh                                   # build + run + plot in one shot
```

## What the plot shows

- 17 Pennacchioni+10 datapoints (γ ∈ [0.3, 15] at T ~ 500°C, ODF J → M via
  Skemer+05) coloured by deformation regime (WDV / MDV / SDV)
- Calibrated case-3 saturating-exp curve δ(γ) = 3.0·0.75·(1 − exp(−γ/3.50)) + 1
- 16 MDOODZ HDF5 stars (per-step δ_mean) sampled at γ ∈ {1, 2, ..., 16}

If the MDOODZ stars track the calibrated curve to within ~1e-6 relative
error, `ani_fstrain = 2` with `aniso_db = 3` is correctly evaluating
`anisoDelta_Quartz(FS_AR)`. The unusual smoothness of the Pennacchioni
datapoints around the curve (σ/M_∞ = 0.045, vs Blackford+24 real-data
scatter at σ/M_∞ = 0.191) is the audit signature noted above.

## References

- **Pennacchioni+10** JGR 115, B12405 — natural quartz shear zones, ODF J
  values audit-suspect (see caveat above). PDF in
  [`~/bib/anisotropy_calibration_2026/quartz/Pennacchioni_etal_2010_quartz_CPO_natural_shear_zones_JGR.pdf`](../../../bib/anisotropy_calibration_2026/quartz/Pennacchioni_etal_2010_quartz_CPO_natural_shear_zones_JGR.pdf).
- Derived 17-point CSV at
  [`misc/aniso_fstrain/data/quartz/derived_Pennacchioni10_17pts.csv`](../../misc/aniso_fstrain/data/quartz/derived_Pennacchioni10_17pts.csv).
- Calibration script:
  [`misc/aniso_fstrain/calibrate/calibrate_quartz.py`](../../misc/aniso_fstrain/calibrate/calibrate_quartz.py).
- Audit memo:
  [`misc/aniso_fstrain/audit/Pennacchioni_etal_2010_audit.md`](../../misc/aniso_fstrain/audit/Pennacchioni_etal_2010_audit.md).
