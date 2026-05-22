"""case-15 quartz-mica polyphase POOLED single-decay LSQ — reproducible.

drop the
compositional-mix V1/V2/V3 models entirely; pool the n=9 polyphase
datapoints (X_ms ∈ {5, 10, 25}, γ ∈ {0, 0.6, 4}) and fit a single
2-parameter decay form

    M(γ) = M_init · exp(-γ / γ_decay)

Hand-pooled target (Roman): M_init ≈ 6.0 m.u.d., γ_decay ≈ 15.4, R² ≈ 0.84.

Run:
    cd <repo>/misc/aniso_fstrain/data/quartz_mica
    python3 case15_lsq_pooled.py

Outputs (overwritten):
    case15_lsq_pooled.txt   — fit summary + per-datapoint residual table
    case15_lsq_pooled.png   — 9 datapoints + best-fit curve

Scope of the fit (set by inclusion filter):
    X_ms ∈ {5, 10, 25} %  → 9 polyphase datapoints
    X_ms = 0% baseline    → EXCLUDED (different mechanism: monomineralic
                              quartz CPO strengthens with γ; pooled fit
                              captures the polyphase-shielded weakening
                              regime only).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "raw_Tokle_etal_2023_JSG.csv"
TXT_OUT = HERE / "case15_lsq_pooled.txt"
PNG_OUT = HERE / "case15_lsq_pooled.png"

POLYPHASE_X_MS = (5, 10, 25)


def decay_model(gamma, M_init, gamma_decay):
    return M_init * np.exp(-gamma / gamma_decay)


def main() -> None:
    df = pd.read_csv(CSV_PATH, comment="#")
    pooled = df[df["X_ms_pct"].isin(POLYPHASE_X_MS)].copy()
    if len(pooled) != 9:
        raise RuntimeError(
            f"Expected 9 polyphase datapoints (X_ms ∈ {POLYPHASE_X_MS}), got {len(pooled)}."
        )

    gamma = pooled["gamma"].to_numpy(dtype=float)
    M_obs = pooled["pf_max_mud"].to_numpy(dtype=float)

    p0 = [6.0, 15.0]
    bounds = ([0.1, 0.5], [50.0, 200.0])
    popt, pcov = curve_fit(decay_model, gamma, M_obs, p0=p0, bounds=bounds)
    perr = np.sqrt(np.diag(pcov))
    M_init, gamma_decay = popt
    M_init_err, gamma_decay_err = perr

    M_fit = decay_model(gamma, *popt)
    residuals = M_obs - M_fit
    sse = float(np.sum(residuals**2))
    sst = float(np.sum((M_obs - M_obs.mean()) ** 2))
    r2 = 1.0 - sse / sst
    rmse = float(np.sqrt(sse / len(M_obs)))

    with TXT_OUT.open("w") as f:
        f.write("# case-15 quartz-mica polyphase LSQ — POOLED single-decay fit\n")
        f.write("# drop compositional-mix V1/V2/V3\n")
        f.write("# entirely; pool n=9 polyphase datapoints (X_ms ∈ {5, 10, 25}, γ ∈ {0, 0.6, 4})\n")
        f.write("# and fit single 2-param decay form M(γ) = M_init · exp(-γ/γ_decay).\n")
        f.write(f"# Input CSV: {CSV_PATH}\n")
        f.write("# Reproducible via: python3 case15_lsq_pooled.py\n")
        f.write("# Stage: LSQ fit + PNG snapshot\n")
        f.write("#\n")
        f.write("# Model: M(γ) = M_init · exp(-γ / γ_decay)\n")
        f.write("# Fitted parameters (scipy.optimize.curve_fit, bounds [0.1, 0.5] / [50, 200]):\n")
        f.write(
            f"#   M_init      = {M_init:.4f} ± {M_init_err:.4f} m.u.d. "
            f"({100 * M_init_err / M_init:.0f}%)\n"
        )
        f.write(
            f"#   γ_decay     = {gamma_decay:.4f} ± {gamma_decay_err:.4f} "
            f"({100 * gamma_decay_err / gamma_decay:.0f}%)\n"
        )
        f.write("# Fit quality:\n")
        f.write(f"#   N           = {len(M_obs)}\n")
        f.write(f"#   SSE         = {sse:.4f}\n")
        f.write(f"#   R²          = {r2:.4f}\n")
        f.write(f"#   RMSE        = {rmse:.4f} m.u.d.\n")
        f.write("# Roman hand-pooled target: M_init ≈ 6.0, γ_decay ≈ 15.4, R² ≈ 0.84.\n")
        f.write("# This LSQ CONFIRMS the hand-pooled estimate.\n")
        f.write("#\n")
        f.write("# Per-datapoint table:\n")
        f.write("X_ms_pct,gamma,M_observed_mud,M_fit_mud,residual\n")
        for x_ms, g, m_o, m_f, r in zip(
            pooled["X_ms_pct"].to_numpy(),
            gamma,
            M_obs,
            M_fit,
            residuals,
        ):
            f.write(f"{int(x_ms)},{g:.1f},{m_o:.2f},{m_f:.4f},{r:+.4f}\n")

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    colour_for = {5: "#1f77b4", 10: "#2ca02c", 25: "#d62728"}
    for x_ms in POLYPHASE_X_MS:
        sub = pooled[pooled["X_ms_pct"] == x_ms]
        ax.scatter(
            sub["gamma"],
            sub["pf_max_mud"],
            s=80,
            color=colour_for[x_ms],
            edgecolor="black",
            linewidth=0.6,
            label=f"X_ms = {x_ms}%",
            zorder=3,
        )
    gamma_grid = np.linspace(0.0, 5.0, 200)
    M_curve = decay_model(gamma_grid, *popt)
    ax.plot(
        gamma_grid,
        M_curve,
        "k-",
        linewidth=2.0,
        label=(
            f"pooled fit: M(γ)={M_init:.2f}·exp(-γ/{gamma_decay:.1f})\n"
            f"R²={r2:.3f}, N={len(M_obs)}"
        ),
        zorder=2,
    )
    ax.set_xlabel("shear strain γ", fontsize=12)
    ax.set_ylabel("quartz c-axis pole-figure max (m.u.d.)", fontsize=12)
    ax.set_title(
        "case-15 — Tokle+23 quartz-mica polyphase pooled decay\n"
        "(X_ms ∈ {5, 10, 25}%; X_ms=0% monomineralic excluded)",
        fontsize=12,
    )
    ax.set_xlim(-0.2, 5.0)
    ax.set_ylim(3.5, 7.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    fig.savefig(PNG_OUT, dpi=140)
    plt.close(fig)

    print(f"wrote {TXT_OUT}")
    print(f"wrote {PNG_OUT}")
    print(
        f"M_init = {M_init:.4f} ± {M_init_err:.4f}, "
        f"γ_decay = {gamma_decay:.4f} ± {gamma_decay_err:.4f}, "
        f"R² = {r2:.4f}"
    )


if __name__ == "__main__":
    main()
