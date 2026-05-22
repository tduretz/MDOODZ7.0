#!/usr/bin/env python3
"""Calibrate quartz M(γ) → δ for cases 3 and 6.

Reproduces MDLIB/FlowLaws.c:
    case 3  anisoDelta_Quartz            (M_inf = 0.75, γ_e = 3.50, slope = 3.0)
    case 6  anisoDelta_Quartz_Blackford  (M_inf = 0.20, γ_e = 4.00, slope = 11.0)

Reads from misc/aniso_fstrain/data/quartz/:
    derived_Pennacchioni10_17pts.csv
    derived_Blackford24_38pts.csv
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aniso_data import load_raw

OUT_DIR   = Path(__file__).resolve().parent.parent / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_PATH = OUT_DIR / "quartz_calibration.png"

SLOPE_C3  = 3.0
SLOPE_C6  = 11.0


def m_of_gamma(g, M_inf, gamma_e):
    return M_inf * (1.0 - np.exp(-g / gamma_e))


def delta_quartz(g, M_inf, gamma_e, slope):
    return slope * m_of_gamma(g, M_inf, gamma_e) + 1.0


def grid_search(g, m, mi_r, ge_r, n):
    mi_grid = np.linspace(*mi_r, n)
    ge_grid = np.linspace(*ge_r, n)
    best = (None, None, float("inf"))
    for mi in mi_grid:
        for ge in ge_grid:
            sse = float(np.sum((m_of_gamma(g, mi, ge) - m) ** 2))
            if sse < best[2]:
                best = (float(mi), float(ge), sse)
    return best


def fit(g, m, mi_lo=0.05, mi_hi=0.95, ge_lo=0.5, ge_hi=8.0):
    coarse = grid_search(g, m, (mi_lo, mi_hi), (ge_lo, ge_hi), n=60)
    cmi, cge, _ = coarse
    fine = grid_search(g, m,
                       (max(0.01, cmi - 0.15), min(0.99, cmi + 0.15)),
                       (max(0.10, cge - 1.5), cge + 1.5),
                       n=80)
    fmi, fge, _ = fine
    return grid_search(g, m,
                       (max(0.01, fmi - 0.05), min(0.99, fmi + 0.05)),
                       (max(0.10, fge - 0.5), fge + 0.5),
                       n=200)


def main() -> None:
    pen = load_raw("quartz/derived_Pennacchioni10_17pts.csv")
    g_p, M_p = pen["gamma_simple"], pen["M_skemer"]

    bla = load_raw("quartz/derived_Blackford24_38pts.csv")
    g_b, M_b = bla["gamma_simple_eq"], bla["M_skemer_pfJ"]

    mi_p, ge_p, sse_p = fit(g_p, M_p)
    # case 6: γ_e is poorly constrained because Blackford samples span only
    # γ ∈ [1.3, 5.1] — saturation is not strongly expressed.  Free LSQ
    # produces a flat trough where M_inf ≈ 0.17–0.29 over γ_e ∈ [3, 10].
    # We commit M_inf = 0.20, γ_e = 4.0 as a centered, defensible point in
    # that envelope (matches FlowLaws.c case 6 comment).
    mi_b, ge_b, sse_b = fit(g_b, M_b, mi_lo=0.05, mi_hi=0.50, ge_lo=2.0, ge_hi=8.0)
    print(f"  (Note: case-6 γ_e is poorly constrained by Blackford γ-range; "
          f"committed value γ_e = 4.0 is a centered choice within the LSQ envelope.)")

    print("=== Quartz calibration summary ===")
    print(f"  {'case':<22s} {'N':>3s}  {'M_inf':>7s}  {'γ_e':>7s}  "
          f"{'slope':>5s}  {'δ_∞':>7s}  {'RMS(M)':>7s}")
    for label, n, mi, ge, sse, slope in [
        ("case 3 Pennacchioni+10",  len(g_p), mi_p, ge_p, sse_p, SLOPE_C3),
        ("case 6 Blackford+24",     len(g_b), mi_b, ge_b, sse_b, SLOPE_C6),
    ]:
        print(f"  {label:<22s} {n:>3d}  {mi:>7.4f}  {ge:>7.4f}  "
              f"{slope:>5.1f}  {slope * mi + 1.0:>7.3f}  {math.sqrt(sse / n):>7.4f}")

    # ----- 2-panel plot -----
    grid = np.linspace(0.0, 20.0, 400)
    fig, (ax_p, ax_b) = plt.subplots(1, 2, figsize=(14.0, 5.4),
                                      constrained_layout=True)

    # Pennacchioni — color by regime
    regimes = {0: ("WDV", "C0"), 1: ("MDV", "C1"), 2: ("SDV", "C3")}
    for r_code, (lbl, c) in regimes.items():
        mask = pen["regime"] == r_code
        if mask.any():
            ax_p.scatter(g_p[mask], M_p[mask], s=70,
                          facecolors="none", edgecolors=c, linewidths=1.3,
                          label=f"Pennacchioni+10 {lbl} (n={mask.sum()})")
    ax_p.plot(grid, m_of_gamma(grid, mi_p, ge_p), color="black", lw=2.0,
               label=f"case 3 fit  M_inf={mi_p:.3f}, γ_e={ge_p:.2f}")
    ax_p.plot(grid, m_of_gamma(grid, 0.75, 3.50), color="black", lw=1.0, ls="--",
               label="committed (rounded) 0.75 / 3.50")
    ax_p.set_xlabel("shear strain γ")
    ax_p.set_ylabel("M-index (Skemer of ODF J)")
    ax_p.set_xlim(-0.5, 17)
    ax_p.set_ylim(0, 0.85)
    ax_p.set_title(f"case 3 — Pennacchioni+10 ODF-J fit (n={len(g_p)})")
    ax_p.grid(alpha=0.3)
    ax_p.legend(loc="lower right", fontsize=8)

    # Blackford — color by strain domain
    domains = {1: "C0", 2: "C1", 3: "C2", 4: "C3", 5: "C4"}
    for d_code, c in domains.items():
        mask = bla["strain_domain"] == d_code
        if mask.any():
            ax_b.scatter(g_b[mask], M_b[mask], s=55,
                          facecolors="none", edgecolors=c, linewidths=1.2,
                          label=f"domain {d_code} (n={mask.sum()})")
    ax_b.plot(grid, m_of_gamma(grid, mi_b, ge_b), color="black", lw=2.0,
               label=f"case 6 fit  M_inf={mi_b:.3f}, γ_e={ge_b:.2f}")
    ax_b.plot(grid, m_of_gamma(grid, 0.20, 4.00), color="black", lw=1.0, ls="--",
               label="committed (centered) 0.20 / 4.00")
    ax_b.set_xlabel("shear strain γ_simple_eq")
    ax_b.set_ylabel("M-index (Skemer of pfJ)")
    ax_b.set_xlim(0, 6)
    ax_b.set_ylim(0, 0.35)
    ax_b.set_title(f"case 6 — Blackford+24 pfJ fit (n={len(g_b)})")
    ax_b.grid(alpha=0.3)
    ax_b.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"Quartz calibrations — case 3 (ODF J): δ_∞ = {SLOPE_C3*mi_p+1:.2f}; "
        f"case 6 (pfJ): δ_∞ = {SLOPE_C6*mi_b+1:.2f}",
        fontsize=11,
    )
    fig.savefig(PLOT_PATH, dpi=140)
    print(f"\nWrote {PLOT_PATH} ({PLOT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
