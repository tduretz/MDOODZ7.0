# AniFstrainCalcitePieri — Pieri+01 calcite calibration set

Homogeneous simple-shear MDOODZ run that traces δ(γ) under the Pieri+01
(Tectonophysics 330, 119) Carrara marble torsion calibration —
`aniso_db = 5` (case 5, `anisoDelta_Calcite_HighT`, T = 727°C, M_∞ = 0.73,
γ_e = 5.28, slope = 6.0, asymptote δ_∞ = 5.40).

The lab corpus (γ ∈ {0, 1, 2, 5, 11}, J converted to M via Skemer+05) is
overlaid on the MDOODZ HDF5 output by [mdoodz_vs_pieri.py](mdoodz_vs_pieri.py),
producing the datapoint-vs-MDOODZ plot in `img/`.

## How to run

```bash
# 1. build (from MDOODZ7.0 repo root, after cmake configure)
cmake --build cmake-exec --target AniFstrainCalcitePieri

# 2. run the simulation
cd cmake-exec/AniFstrainCalcitePieri
mkdir -p output                        # writer_subfolder must exist
./AniFstrainCalcitePieri               # writes 12 HDF5 outputs into output/

# 3. produce the comparison plot
cd <repo>/SETS/AniFstrainCalcitePieri
python3 mdoodz_vs_pieri.py             # writes img/mdoodz_vs_pieri.png
```

Or use the convenience wrapper:

```bash
cd <repo>/SETS/AniFstrainCalcitePieri
./run.sh                               # build + run + plot in one shot
```

## What the plot shows

- 5 Pieri+01 datapoints (γ ∈ {0, 1, 2, 5, 11} at 727°C) as δ = 6·(J−1)/15 + 1
- Calibrated case-5 saturating-exp curve δ(γ) = 6·0.73·(1 − exp(−γ/5.28)) + 1
- 12 MDOODZ HDF5 stars (per-step δ_mean) sampled at γ ∈ {0, 1, ..., 11}

If the MDOODZ stars track the calibrated curve to within ~1e-4 relative
error, `ani_fstrain = 2` with `aniso_db = 5` is correctly evaluating
`anisoDelta_Calcite_HighT(FS_AR)`. Disagreement between lab points and
MDOODZ stars is the calibration residual (Pieri-vs-Barnhoorn 727°C
disagreement at γ = 5 contributes ~17% scatter — see the docstring of
[`MDLIB/FlowLaws.c::anisoDelta_Calcite_HighT`](../../MDLIB/FlowLaws.c)).

## References

- **Pieri+01** Tectonophysics 330, 119 — Carrara marble torsion at 727°C,
  300 MPa, γ ∈ {0, 1, 2, 5, 11}. PDF in
  `~/bib/anisotropy_calibration_2026/calcite/Pieri_etal_2001_Carrara_marble_torsion_727C_Tectonophysics.pdf`.
- Combined raw data CSV at
  [`misc/aniso_fstrain/data/calcite/raw_Pieri_etal_2001_Tectonophysics.csv`](../../misc/aniso_fstrain/data/calcite/raw_Pieri_etal_2001_Tectonophysics.csv).
- Calibration script that fitted case 5:
  [`misc/aniso_fstrain/calibrate/calibrate_calcite.py`](../../misc/aniso_fstrain/calibrate/calibrate_calcite.py).
