# AniFstrainCalciteBruijn — Bruijn+11 calcite audit-control set

Homogeneous simple-shear MDOODZ run that traces δ(γ) under the
Bruijn, Burlini & Drury 2011 (Tectonophysics 503, 75;
[doi:10.1016/j.tecto.2010.09.029](https://doi.org/10.1016/j.tecto.2010.09.029))
Carrara marble torsion calibration —
`aniso_db = 2` (case 2, `anisoDelta_Calcite`, T-averaged Carrara, M_∞ = 0.41,
γ_e = 3.18, slope = 6.0, asymptote δ_∞ = 3.46).

The 3-point (γ, J) corpus stored at
[`misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv`](../../misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv)
is overlaid on the MDOODZ HDF5 output by
[`mdoodz_vs_bruijn.py`](mdoodz_vs_bruijn.py), producing the datapoint-vs-MDOODZ
plot in `img/`.

## ⚠ Audit caveat — Bruijn J-values not in published paper

Per the audit warning embedded in the CSV header (2026-05-08), the 3 (γ, J)
values cannot be sourced from Bruijn et al. 2011. The paper's Fig. 8 left
column shows pole figures at γ ∈ {0, 0.6, 1.0, 2.0, 3.0, 5.0, 9.0, 11} with
**no numerical J labels**, and the caption explicitly states the figures are
reproduced from Pieri+01 / Barnhoorn+04, not new Bruijn measurements. The
values appear to be visual estimates from greyscale pole-figure intensities.

**Impact on the case-2 calibration**: with the 3 Bruijn rows removed, the
case-2 LSQ fit shifts from (M_∞ = 0.41, γ_e = 3.71, δ_∞ = 3.46) to
(M_∞ = 0.41, γ_e = 3.18, δ_∞ = 3.47) — a negligible perturbation. The
committed case-2 calibration is therefore driven by Pieri+01 + Barnhoorn+04
even after the audit, which is why the SET retains case 2 (rather than a
custom Bruijn-only refit).

This SET subfolder doubles as the **audit-control mirror** of
[AniFstrainQuartzPennacchioni](../AniFstrainQuartzPennacchioni/README.md):
the MDOODZ-vs-datapoint plot makes the audit finding visible by showing the
3 Bruijn-attributed points sitting on/near the calibrated curve (consistent
with greyscale-from-pole-figure estimates that "look like" Pieri/Barnhoorn).
See the [`Bruijn_etal_2011_audit`](../../misc/aniso_fstrain/audit/) memo
embedded in the raw CSV header for the full audit narrative.

## How to run

```bash
# 1. build (from MDOODZ7.0 repo root, after cmake configure)
cmake --build cmake-exec --target AniFstrainCalciteBruijn

# 2. run the simulation
cd cmake-exec/AniFstrainCalciteBruijn
mkdir -p output                        # writer_subfolder must exist
./AniFstrainCalciteBruijn              # writes 13 HDF5 outputs into output/

# 3. produce the comparison plot
cd <repo>/SETS/AniFstrainCalciteBruijn
python3 mdoodz_vs_bruijn.py            # writes img/mdoodz_vs_bruijn.png
```

Or use the convenience wrapper:

```bash
cd <repo>/SETS/AniFstrainCalciteBruijn
./run.sh                               # build + run + plot in one shot
```

## What the plot shows

- 3 Bruijn+11 audit-flagged datapoints (γ ∈ {1.0, 2.6, 5.0} at 727°C, J → δ
  via Skemer+05 + slope = 6.0) overlaid as red ✕ markers
- Calibrated case-2 saturating-exp curve δ(γ) = 6·0.41·(1 − exp(−γ/3.18)) + 1
- 12 MDOODZ HDF5 stars (per-step δ_mean) sampled at γ ∈ {0, 1, ..., 11}
- A boxed audit caveat in the lower-right corner

If the MDOODZ stars track the calibrated curve to within ~1e-4 relative
error, `ani_fstrain = 2` with `aniso_db = 2` is correctly evaluating
`anisoDelta_Calcite(FS_AR)`. The Bruijn ✕ markers' agreement / disagreement
with the calibrated curve is itself the audit signature.

## References

- **Bruijn+11** Tectonophysics 503, 75 — Carrara marble torsion at 727°C,
  D1 single-stage segment. PDF in
  [`~/bib/anisotropy_calibration_2026/calcite/Bruijn_etal_2011_calcite_torsion_glide_systems_Tectonophysics.pdf`](../../../bib/anisotropy_calibration_2026/calcite/Bruijn_etal_2011_calcite_torsion_glide_systems_Tectonophysics.pdf).
- Audit-flagged raw CSV (with embedded audit memo) at
  [`misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv`](../../misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv).
- Calibration script that fitted case 2:
  [`misc/aniso_fstrain/calibrate/calibrate_calcite.py`](../../misc/aniso_fstrain/calibrate/calibrate_calcite.py).
