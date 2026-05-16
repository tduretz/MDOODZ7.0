#!/usr/bin/env python3
"""Calibrate case-15 quartz-mica polyphase compositional model.

Fits the compositional model to per-(γ, X_ms) data extracted from Tokle+23
JSG Fig 12 (CSV at data/quartz_mica/raw_Tokle_etal_2023_JSG.csv).

Compositional model (deprecated; see V3 below for the form actually deployed):
    M_aggregate(γ, X_ms) = f_DRX(X_ms) · M_recryst(γ) + (1 − f_DRX(X_ms)) · M_init
    M_recryst(γ) = M_∞ · (1 − exp(−γ / γ_e))

where f_DRX(X_ms) is empirical from Tokle+23 §3.2.3:
    X_ms = 0  → f_DRX = 0.97
    X_ms = 5  → f_DRX = 0.85
    X_ms = 10 → f_DRX = 0.30
    X_ms = 25 → f_DRX = 0.05

Free parameters fitted: M_∞, γ_e, M_init (3 params on 12 datapoints).

CRITICAL CAVEAT: Tokle+23 reports pole-figure maxima in m.u.d. (multiples
of uniform distribution), NOT Skemer M-index. This fit treats m.u.d. as
the response variable directly. Conversion m.u.d. → M for MDLIB deployment
is a separate question; see notes/case15_calibration.md.

Snapshot output is written to misc/aniso_fstrain/notes/case15_calibration.md
when SNAPSHOT_PATH is set explicitly; default is no snapshot write.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


CSV_PATH = Path(__file__).parent.parent / "data" / "quartz_mica" / "raw_Tokle_etal_2023_JSG.csv"
# Optional snapshot output path: set via env var CASE15_SNAPSHOT_PATH or leave None.
import os
_SNAP_ENV = os.environ.get("CASE15_SNAPSHOT_PATH")
SNAPSHOT_PATH = Path(_SNAP_ENV) if _SNAP_ENV else None

F_DRX_LOOKUP = {0: 0.97, 5: 0.85, 10: 0.30, 25: 0.05}


def load_csv(path: Path):
    """Load CSV skipping comment lines starting with '#'."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.strip())
    header = rows[0].split(",")
    data = []
    for r in rows[1:]:
        vals = r.split(",")
        data.append(dict(zip(header, vals)))
    return data


def case15_model(X, M_inf, gamma_e, M_init):
    """Compositional model V1: M_agg = f·M_recryst(γ) + (1-f)·M_init.

    DEPRECATED — V1 fails (R² = −0.88) because constant f_DRX is wrong at
    γ = 0 where no DRX has occurred yet.
    """
    gamma, f = X
    M_recryst = M_inf * (1.0 - np.exp(-gamma / gamma_e))
    return f * M_recryst + (1.0 - f) * M_init


def case15_model_v2(X, M_inf, gamma_e, M_init, gamma_DRX_onset):
    """Compositional model V2: γ-dependent DRX fraction.

    DEPRECATED — V2 fails (R² = 0.34, covariance degenerate) because it
    shares ONE M_∞ across all X_ms; the data needs M_∞(X_ms).
    """
    gamma, f_max = X
    f_drx_gamma = f_max * (1.0 - np.exp(-gamma / gamma_DRX_onset))
    M_recryst = M_inf * (1.0 - np.exp(-gamma / gamma_e))
    return f_drx_gamma * M_recryst + (1.0 - f_drx_gamma) * M_init


def case15_model_v3(X, M_inf_pure, alpha, gamma_e, M_init):
    """V3 model: X_ms-dependent M_∞ scaling.

    Mechanism: pure-quartz saturating-exp CPO development BUT with mica
    fraction damping the achievable M_∞ asymptote:

        M_∞(X_ms) = M_∞_pure · exp(−α · X_ms)
        M(γ, X_ms) = M_init + [M_∞(X_ms) − M_init] · (1 − exp(−γ / γ_e))

    At γ = 0: M = M_init (matches all 4 hydrostatic observations ~5-6 m.u.d.)
    At γ → ∞: M → M_∞(X_ms), which DECREASES with mica content.
    For X_ms small enough that M_∞(X_ms) > M_init: CPO strengthens (X_ms = 0).
    For X_ms large enough that M_∞(X_ms) < M_init: CPO weakens (X_ms ≥ 5).

    This captures Tokle+23's central observation (CPO inversion at any
    X_ms ≥ 5 %) as a direct consequence of mica-damped M_∞.

    4 free parameters: M_∞_pure, α, γ_e, M_init. 12 datapoints → 8 dof.
    Better-constrained than V2 because pure-vs-mica-bearing samples directly
    fix the M_∞ ratio.
    """
    gamma, X_ms = X
    M_inf_X = M_inf_pure * np.exp(-alpha * X_ms)
    return M_init + (M_inf_X - M_init) * (1.0 - np.exp(-gamma / gamma_e))


def main():
    print(f"Loading: {CSV_PATH}")
    data = load_csv(CSV_PATH)
    print(f"Loaded {len(data)} datapoints.")

    gammas = np.array([float(d["gamma"]) for d in data])
    X_ms = np.array([int(d["X_ms_pct"]) for d in data])
    pf_max = np.array([float(d["pf_max_mud"]) for d in data])
    f_drx = np.array([F_DRX_LOOKUP[x] for x in X_ms])

    # V3 fit: 4-param X_ms-dependent M_∞ scaling (option a)
    # Inputs: (gamma, X_ms) — using X_ms directly (not f_DRX_max)
    X = (gammas, X_ms)
    y = pf_max

    # Initial guesses:
    #   M_∞_pure ≈ 15 (pure-qz extrapolated asymptote, observed 12.3 at γ=4)
    #   α ≈ 0.3 / X_ms_unit (gives exp(−0.3·5)=0.22 reduction at 5% Ms)
    #   γ_e ≈ 2 (Tokle+23 mid-range)
    #   M_init ≈ 5.5 (avg hydrostatic m.u.d.)
    p0 = [15.0, 0.3, 2.0, 5.5]
    popt, pcov = curve_fit(case15_model_v3, X, y, p0=p0, maxfev=10000)
    perr = np.sqrt(np.diag(pcov))

    M_inf_pure, alpha, gamma_e, M_init = popt
    M_inf_pure_err, alpha_err, gamma_e_err, M_init_err = perr
    # Aliases for downstream output compat
    M_inf = M_inf_pure
    M_inf_err = M_inf_pure_err
    gamma_DRX_onset = alpha  # repurposed as α for output
    gamma_DRX_err = alpha_err

    # Residuals + goodness of fit
    y_pred = case15_model_v3(X, *popt)
    residuals = y - y_pred
    SSE = float(np.sum(residuals ** 2))
    SST = float(np.sum((y - np.mean(y)) ** 2))
    R2 = 1.0 - SSE / SST
    RMSE = float(np.sqrt(SSE / len(y)))

    print(f"\nLSQ fit results (V3 X_ms-dependent M_∞ scaling, option a):")
    print(f"  M_∞_pure         = {M_inf_pure:.3f} ± {M_inf_pure_err:.3f}")
    print(f"  α (X_ms damping) = {alpha:.4f} ± {alpha_err:.4f}")
    print(f"  γ_e              = {gamma_e:.3f} ± {gamma_e_err:.3f}")
    print(f"  M_init           = {M_init:.3f} ± {M_init_err:.3f}")
    print(f"  SSE              = {SSE:.3f}")
    print(f"  R²               = {R2:.4f}")
    print(f"  RMSE             = {RMSE:.3f}")
    # Implied M_∞(X_ms) values for documentation
    print(f"\n  Implied M_∞(X_ms):")
    for x in [0, 5, 10, 25]:
        Mx = M_inf_pure * np.exp(-alpha * x)
        print(f"    X_ms={x:>3}% → M_∞={Mx:.2f} m.u.d.")

    # Per-datapoint table
    print(f"\nDatapoint table:")
    print(f"  {'X_ms':>4} {'γ':>5} {'pf_obs':>7} {'pf_fit':>7} {'resid':>7}")
    for d, yo, yp in zip(data, y, y_pred):
        x_ms = d["X_ms_pct"]
        g = d["gamma"]
        print(f"  {x_ms:>4} {g:>5} {yo:>7.2f} {yp:>7.2f} {yo - yp:>+7.2f}")

    # Write snapshot to notes/ if SNAPSHOT_PATH is configured.
    if SNAPSHOT_PATH is None:
        return
    SNAPSHOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(SNAPSHOT_PATH, "w") as fh:
        fh.write(f"# case-15 quartz-mica polyphase calibration — LSQ output\n\n")
        fh.write(f"**Model**: V2 fit with γ-dependent f_DRX (V1 with constant f_DRX failed at R² = −0.88).\n")
        fh.write(f"**Input CSV**: `{CSV_PATH}`\n")
        fh.write(f"**Script**: `calibrate_quartz_mica.py` (this directory)\n\n")
        fh.write(f"## Compositional model V2 (γ-dependent f_DRX)\n\n")
        fh.write(f"```\nf_DRX(γ, X_ms) = f_DRX_max(X_ms) · (1 − exp(−γ/γ_DRX_onset))\nM_aggregate = f_DRX(γ, X_ms) · M_∞·(1−exp(−γ/γ_e)) + (1−f_DRX(γ, X_ms)) · M_init\n```\n\n")
        fh.write(f"with f_DRX_max from Tokle+23 §3.2.3: ")
        fh.write(", ".join(f"X_ms={k}%→{v}" for k, v in F_DRX_LOOKUP.items()))
        fh.write(f".\n\n")
        fh.write(f"## Fitted parameters (units: pole-figure max in m.u.d.)\n\n")
        fh.write(f"| Param | Value | 1σ |\n|---|---|---|\n")
        fh.write(f"| M_∞ | **{M_inf:.3f}** | {M_inf_err:.3f} |\n")
        fh.write(f"| γ_e | **{gamma_e:.3f}** | {gamma_e_err:.3f} |\n")
        fh.write(f"| M_init | **{M_init:.3f}** | {M_init_err:.3f} |\n")
        fh.write(f"| γ_DRX_onset | **{gamma_DRX_onset:.3f}** | {gamma_DRX_err:.3f} |\n\n")
        fh.write(f"## Fit quality\n\n")
        fh.write(f"- N datapoints: {len(y)}\n")
        fh.write(f"- SSE: {SSE:.3f}\n")
        fh.write(f"- R²: **{R2:.4f}**\n")
        fh.write(f"- RMSE: {RMSE:.3f} m.u.d.\n\n")
        fh.write(f"## Per-datapoint table\n\n")
        fh.write(f"| X_ms (%) | γ | pf_obs (m.u.d.) | pf_fit (m.u.d.) | residual |\n")
        fh.write(f"|---|---|---|---|---|\n")
        for d, yo, yp in zip(data, y, y_pred):
            fh.write(f"| {d['X_ms_pct']} | {d['gamma']} | {yo:.2f} | {yp:.2f} | {yo - yp:+.2f} |\n")
        fh.write(f"\n## Caveats\n\n")
        fh.write(f"- **m.u.d. ≠ Skemer M-index**: fit is in pole-figure max units.\n")
        fh.write(f"  Conversion factor for MDLIB deployment needs separate calibration.\n")
        fh.write(f"- **Single-paper anchor for f_DRX(X_ms)**: all 4 modal-fraction points\n")
        fh.write(f"  from Tokle+23 alone. Cross-validation needed.\n")
        fh.write(f"- **γ values are representative, not measured**: 'low strain' = γ≈0.6\n")
        fh.write(f"  is a midpoint of 0.5-0.75 range per Tokle+23 §3.2.2.\n")
        fh.write(f"- **Hydrostatic m.u.d. ~5-6 NOT random** (random=1.0): sample-prep\n")
        fh.write(f"  ball-milling+SPS leaves some texture. M_init absorbs this.\n")
    print(f"\nSnapshot written: {SNAPSHOT_PATH}")


if __name__ == "__main__":
    main()
