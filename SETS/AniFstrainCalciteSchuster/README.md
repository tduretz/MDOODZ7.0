# AniFstrainCalciteSchuster — Schuster+18 calcite γ-coverage annotation set

Homogeneous simple-shear MDOODZ run that traces δ(γ) under the
Schuster, Habler, Schafler & Abart 2018 (J. Struct. Geol.;
[doi:10.1016/j.jsg.2018.09.003](https://doi.org/10.1016/j.jsg.2018.09.003))
Carrara marble high-pressure-torsion (HPT) experiments at 450°C, 1-4 GPa —
`aniso_db = 4` (case 4, `anisoDelta_Calcite_LowT`, M_∞ = 0.22, γ_e = 1.37,
slope = 6.0, asymptote δ_∞ = 2.32).

The 16-row γ-coverage CSV at
[`misc/aniso_fstrain/data/calcite/raw_Schuster_etal_2018_JSG.csv`](../../misc/aniso_fstrain/data/calcite/raw_Schuster_etal_2018_JSG.csv)
is overlaid on the MDOODZ HDF5 output by
[`mdoodz_vs_schuster.py`](mdoodz_vs_schuster.py) as a horizontal rug at the
bottom of the δ(γ) panel, color-coded by confining pressure (1.2 / 2 / 3 /
4 GPa). Output plot in `img/`.

## ⚠ Annotation-only — Schuster+18 reports no J/M values

Per the audit memo
[`misc/aniso_fstrain/audit/Schuster_etal_2018_annotation.md`](../../misc/aniso_fstrain/audit/Schuster_etal_2018_annotation.md):

- Schuster+18 reports γ ∈ [3, 63] in Table 2 (16 EBSD maps from 450°C
  samples) but **does NOT tabulate J-index, M-index, or peak-m.r.d. for any
  sample**. Texture is shown qualitatively as pole figures (Figs 9, 13).
- All 16 rows in the raw CSV have `fit_status = annotation-only`; they are
  excluded from the LSQ fits in `calibrate_calcite.py`.
- The forward paths to make Schuster's γ-coverage quantitatively usable are
  (a) MTEX re-extraction from digitised pole figures (~2-3 days), (b) author
  contact for raw EBSD `.ctf` files, (c) the realised path 3 — annotation
  only.

### Regime caveat (must accompany every mention)

Schuster's HPT environment is **off-spec for MDOODZ lithospheric calibration**:

| Parameter | Schuster+18 | MDOODZ-relevant (Pieri+01 / Barnhoorn+04) |
|---|---|---|
| Temperature | 235 / 450 °C | 500-727 °C |
| Confining pressure | 1-4 GPa | 300 MPa |
| Strain rate | 1.7e-3 to 8.8e-2 s⁻¹ | 3e-4 s⁻¹ |
| Polymorph | calcite + CaCO₃-II (above 1.6 GPa) | calcite only |

Schuster's §5.2 explicitly notes that HPT pressures produce "low-temperature
characteristics" — narrower e-twins, smaller recrystallized grains — that
mimic deformation at sub-GPa pressures and lower T. So **the same γ under
HPT does not yield the same fabric** as Pieri's Paterson rig at 727°C. Even
if J/M values are eventually extracted, they would not be a straightforward
extension of case 4 — they would either need a new high-pressure case or a
documented regime caveat.

## How to run

```bash
# 1. build (from MDOODZ7.0 repo root, after cmake configure)
cmake --build cmake-exec --target AniFstrainCalciteSchuster

# 2. run the simulation
cd cmake-exec/AniFstrainCalciteSchuster
mkdir -p output                        # writer_subfolder must exist
./AniFstrainCalciteSchuster            # writes 17 HDF5 outputs into output/

# 3. produce the comparison plot
cd <repo>/SETS/AniFstrainCalciteSchuster
python3 mdoodz_vs_schuster.py          # writes img/mdoodz_vs_schuster.png
```

Or use the convenience wrapper:

```bash
cd <repo>/SETS/AniFstrainCalciteSchuster
./run.sh                               # build + run + plot in one shot
```

## What the plot shows

- Calibrated case-4 saturating-exp curve δ(γ) = 6·0.22·(1 − exp(−γ/1.37)) + 1
  reaching δ_∞ = 2.32 by γ ≈ 5
- 16 MDOODZ HDF5 stars (per-step δ_mean) sampled at γ ∈ {0, 1, ..., 16}
- 16 Schuster+18 γ-rug ticks at the bottom of the panel, color-coded by
  confining pressure, spanning γ ∈ [3, 63]
- A boxed annotation in the lower-right corner

If the MDOODZ stars track the calibrated curve to within ~1e-4 relative
error, `ani_fstrain = 2` with `aniso_db = 4` is correctly evaluating
`anisoDelta_Calcite_LowT(FS_AR)`. The Schuster γ-rug is purely
illustrative — it shows where Schuster's high-strain HPT data live in
γ-space without claiming numerical agreement.

## References

- **Schuster+18** JSG — Carrara marble HPT 450°C, 1-4 GPa. PDF in
  [`~/bib/anisotropy_calibration_2026/calcite/Schuster_etal_2018_calcite_HPT_high_shear_strain_JSG.pdf`](../../../bib/anisotropy_calibration_2026/calcite/Schuster_etal_2018_calcite_HPT_high_shear_strain_JSG.pdf).
- Annotation-only raw CSV at
  [`misc/aniso_fstrain/data/calcite/raw_Schuster_etal_2018_JSG.csv`](../../misc/aniso_fstrain/data/calcite/raw_Schuster_etal_2018_JSG.csv).
- Annotation audit memo:
  [`misc/aniso_fstrain/audit/Schuster_etal_2018_annotation.md`](../../misc/aniso_fstrain/audit/Schuster_etal_2018_annotation.md).
- Calibration script that fitted case 4 (Schuster excluded):
  [`misc/aniso_fstrain/calibrate/calibrate_calcite.py`](../../misc/aniso_fstrain/calibrate/calibrate_calcite.py).
- Research question driving this annotation:
  
