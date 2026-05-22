#!/usr/bin/env python3
"""Calibrate calcite M(γ) → δ for cases 2, 4, 5.

Each case loads the relevant derived CSV and refits M(γ) = M_inf · (1 −
exp(−γ/γ_e)) with slope = 6.0 (calcite-specific Tommasi+09 VPSC bound).

Reproduces MDLIB/FlowLaws.c:
    case 2  anisoDelta_Calcite          (M_inf = 0.41,  γ_e = 3.18)
    case 4  anisoDelta_Calcite_LowT     (M_inf = 0.22,  γ_e = 1.37)
    case 5  anisoDelta_Calcite_HighT    (M_inf = 0.73,  γ_e = 5.25)

Reads from misc/aniso_fstrain/data/calcite/:
    derived_combined_15pts_filtered.csv     (Pieri+01 + Barnhoorn+04, no Bruijn)
    derived_LowT_Barnhoorn04_7pts.csv       (Barnhoorn 500 + 600°C, J ≤ 15)
    derived_HighT_8pts.csv                  (Pieri 727°C + Barnhoorn 727°C J ≤ 15)

History: cases 2 and 5 were re-fit on 2026-05-08 after the Bruijn+11
audit (audit/Bruijn_etal_2011_audit.md) found the 3 (γ, J) values
attributed to Bruijn 2011 are not in the paper.  Pre-audit committed
constants were  (case 2: 0.41/3.71; case 5: 0.83/8.24).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aniso_data import load_raw, fs_ar_from_gamma

OUT_DIR   = Path(__file__).resolve().parent.parent / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_PATH = OUT_DIR / "calcite_calibration.png"

SLOPE = 6.0   # calcite-specific δ-vs-M slope (Tommasi+09 VPSC bracket)


def m_of_gamma(g, M_inf, gamma_e):
    return M_inf * (1.0 - np.exp(-g / gamma_e))


def delta_calcite(g, M_inf, gamma_e):
    return SLOPE * m_of_gamma(g, M_inf, gamma_e) + 1.0


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


def fit(g, m):
    coarse = grid_search(g, m, (0.05, 0.95), (0.3, 12.0), n=80)
    cmi, cge, _ = coarse
    fine = grid_search(g, m,
                       (max(0.02, cmi - 0.15), cmi + 0.15),
                       (max(0.10, cge - 2.0), cge + 2.0),
                       n=120)
    fmi, fge, _ = fine
    return grid_search(g, m,
                       (max(0.02, fmi - 0.04), fmi + 0.04),
                       (max(0.10, fge - 0.5), fge + 0.5),
                       n=200)


def load_kept(rel_path: str) -> tuple[np.ndarray, np.ndarray]:
    """Load a derived calcite CSV and return (γ_kept, M_kept)."""
    d = load_raw(f"calcite/{rel_path}")
    if "fit_status" in d:
        keep = d["fit_status"] == "kept"
        return d["gamma_simple"][keep], d["M_skemer"][keep]
    return d["gamma_simple"], d["M_skemer"]


def load_schuster_gamma_only() -> np.ndarray:
    """Schuster+18 HPT calcite: γ values for 450 °C samples only.
    Annotation-only — Schuster reports no J/M, so these never enter the LSQ.
    See audit/Schuster_etal_2018_annotation.md for regime caveats.
    """
    d = load_raw("calcite/raw_Schuster_etal_2018_JSG.csv")
    return np.asarray(d["gamma"], dtype=float)


def main() -> None:
    # ----- case 2 — combined T-averaged -----
    g_c2, m_c2 = load_kept("derived_combined_15pts_filtered.csv")
    mi2, ge2, sse2 = fit(g_c2, m_c2)

    # ----- case 4 — mid-low-T (T ≤ 650°C) -----
    g_c4, m_c4 = load_kept("derived_LowT_Barnhoorn04_7pts.csv")
    mi4, ge4, sse4 = fit(g_c4, m_c4)

    # ----- case 5 — high-T (T > 650°C) -----
    g_c5, m_c5 = load_kept("derived_HighT_8pts.csv")
    mi5, ge5, sse5 = fit(g_c5, m_c5)

    print("=== Calcite calibration summary (M(γ) = M_inf · (1 − exp(−γ/γ_e))) ===")
    print(f"  {'case':<22s} {'N':>3s}  {'M_inf':>7s}  {'γ_e':>7s}  "
          f"{'δ_∞':>7s}  {'RMS(M)':>7s}")
    for label, n, mi, ge, sse in [
        ("case 4 mid-low-T (≤650°C)", len(g_c4), mi4, ge4, sse4),
        ("case 2 T-averaged",         len(g_c2), mi2, ge2, sse2),
        ("case 5 high-T (>650°C)",    len(g_c5), mi5, ge5, sse5),
    ]:
        print(f"  {label:<22s} {n:>3d}  {mi:>7.4f}  {ge:>7.4f}  "
              f"{SLOPE * mi + 1.0:>7.3f}  {math.sqrt(sse / n):>7.4f}")

    g_schuster = load_schuster_gamma_only()

    # ----- plot -----
    grid = np.linspace(0.0, 16.0, 400)
    fig, (ax_lo, ax_hi) = plt.subplots(1, 2, figsize=(14.0, 5.4),
                                        constrained_layout=True)

    def plot_panel(ax, gamma_max):
        # Lab data (case 2 = full combined kept = case 4 + case 5 + Bruijn)
        ax.scatter(g_c4, SLOPE * m_c4 + 1.0, s=70, marker="^",
                    facecolors="none", edgecolors="darkorange", linewidths=1.3,
                    zorder=3, label=f"Barnhoorn ≤650°C (n={len(g_c4)})")
        ax.scatter(g_c5, SLOPE * m_c5 + 1.0, s=70, marker="s",
                    facecolors="none", edgecolors="firebrick", linewidths=1.3,
                    zorder=3, label=f"727°C (Pieri+Bruijn+Barnhoorn, n={len(g_c5)})")

        # Schuster+18 annotation: γ-only rug at top of panel (no M/J)
        in_view = g_schuster[g_schuster <= gamma_max]
        if in_view.size:
            y_rug = 6.6
            ax.scatter(in_view, np.full_like(in_view, y_rug), marker="|",
                        s=180, color="purple", linewidths=1.4, zorder=4,
                        label=f"Schuster+18 HPT 450°C γ-coverage (n={len(in_view)}/{len(g_schuster)}, J not tabulated)")

        # Three calibrated curves
        ax.plot(grid, delta_calcite(grid, mi4, ge4), color="darkorange", lw=2.0,
                 label=f"case 4 mid-low-T (M_inf={mi4:.3f}, γ_e={ge4:.2f}, δ_∞={SLOPE*mi4+1:.2f})")
        ax.plot(grid, delta_calcite(grid, mi2, ge2), color="black", lw=2.0,
                 label=f"case 2 T-avg (M_inf={mi2:.3f}, γ_e={ge2:.2f}, δ_∞={SLOPE*mi2+1:.2f})")
        ax.plot(grid, delta_calcite(grid, mi5, ge5), color="firebrick", lw=2.0,
                 label=f"case 5 high-T (M_inf={mi5:.3f}, γ_e={ge5:.2f}, δ_∞={SLOPE*mi5+1:.2f})")

        ax.set_xlabel("shear strain γ")
        ax.set_ylabel("δ (= 6·M + 1)")
        ax.set_xlim(-0.3, gamma_max)
        ax.set_ylim(0, 7)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8)

    plot_panel(ax_lo, 5.0)
    ax_lo.set_title("Low-strain detail γ ∈ [0, 5]")
    plot_panel(ax_hi, 16.0)
    ax_hi.set_title("Full saturation range γ ∈ [0, 16]")

    fig.suptitle(
        "Calcite δ(γ) — three regime calibrations (cases 4, 2, 5)  "
        f"[mid-low: δ_∞={SLOPE*mi4+1:.2f}, T-avg: {SLOPE*mi2+1:.2f}, "
        f"high-T: {SLOPE*mi5+1:.2f}]",
        fontsize=11,
    )

    fig.savefig(PLOT_PATH, dpi=140)
    print(f"\nWrote {PLOT_PATH} ({PLOT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
