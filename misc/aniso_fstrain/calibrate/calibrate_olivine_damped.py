#!/usr/bin/env python3
"""Calibrate the damped-olivine M(γ) → δ mapping (case 7).

Free LSQ on the combined non-Hansen olivine dataset (Tasaka 2016 wet Fo50 +
Boneh & Skemer 2014 pre-CPO + Kumamoto+19 Josephine peridotite + Bernard+19
natural xenoliths = 93 pts).  Constants must reproduce
MDLIB/FlowLaws.c anisoDelta_HansenOlivine_Damped (M_inf = 0.20, γ_e = 1.30).

Plus a panel-style survey overlay showing all 6 olivine datasets with the
case-1 and case-7 curves.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aniso_data import (
    load_raw,
    j_to_m_skemer,
    gamma_from_axial_eps,
    gamma_from_xz,
)

OUT_DIR   = Path(__file__).resolve().parent.parent / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_PATH = OUT_DIR / "olivine_damped_calibration.png"
SURVEY    = OUT_DIR / "olivine_all_datasets.png"

SLOPE = 24.5  # Hansen+12 Fig 3b — universal across olivine fabric strength

# Reference: case 1 committed values for overlay
CASE_1 = dict(M_inf=0.536, gamma_e=3.96)
# Target (will be recomputed): case 7 committed values
CASE_7 = dict(M_inf=0.20,  gamma_e=1.30)


def m_of_gamma(g, M_inf, gamma_e):
    return M_inf * (1.0 - np.exp(-g / gamma_e))


def grid_search(g, m, mi_range, ge_range, n):
    mi_grid = np.linspace(*mi_range, n)
    ge_grid = np.linspace(*ge_range, n)
    best = (None, None, float("inf"))
    for mi in mi_grid:
        for ge in ge_grid:
            sse = float(np.sum((m_of_gamma(g, mi, ge) - m) ** 2))
            if sse < best[2]:
                best = (float(mi), float(ge), sse)
    return best


def fit(g, m):
    coarse = grid_search(g, m, (0.05, 0.60), (0.2, 8.0), n=80)
    cmi, cge, _ = coarse
    fine = grid_search(g, m,
                       (max(0.02, cmi - 0.10), cmi + 0.10),
                       (max(0.05, cge - 1.5), cge + 1.5),
                       n=120)
    fmi, fge, _ = fine
    return grid_search(g, m,
                       (max(0.02, fmi - 0.03), fmi + 0.03),
                       (max(0.05, fge - 0.5), fge + 0.5),
                       n=200)


def main() -> None:
    # ----- per-dataset (γ, M) tuples -----
    h14 = load_raw("olivine/raw_Hansen_etal_2014_EPSL.csv")
    h16 = load_raw("olivine/raw_Hansen_etal_2016_EPSL.csv")
    g_hansen = np.concatenate([h14["gamma"], h16["gamma"]])
    M_hansen = np.concatenate([h14["M"],     h16["M"]])

    b00 = load_raw("olivine/raw_Bystricky_etal_2000_Science.csv")
    g_byst, M_byst = b00["gamma"], j_to_m_skemer(b00["J"])

    t16 = load_raw("olivine/raw_Tasaka_etal_2016_JGR.csv")
    g_tas, M_tas = t16["gamma"], t16["M"]

    bn = load_raw("olivine/raw_Boneh_Skemer_2014_EPSL.csv")
    g_bon = gamma_from_axial_eps(bn["epsilon"])
    M_bon = bn["M"]
    bon_geom = bn["geometry"].astype(int)

    k19 = load_raw("olivine/raw_Kumamoto_etal_2019_JGR.csv")
    g_kum, M_kum = k19["gamma"], k19["M"]

    bj = load_raw("olivine/raw_Bernard_etal_2019_G3.csv")
    g_brn = gamma_from_xz(bj["xz"])
    M_brn = bj["M"]

    # ----- damped combined fit (case 7) -----
    g_damped = np.concatenate([g_tas, g_bon, g_kum, g_brn])
    M_damped = np.concatenate([M_tas, M_bon, M_kum, M_brn])
    print(f"Damped (non-Hansen) combined dataset: "
          f"Tasaka {len(g_tas)} + Boneh {len(g_bon)} + "
          f"Kumamoto {len(g_kum)} + Bernard {len(g_brn)} = "
          f"{len(g_damped)} pts")

    mi_d, ge_d, sse_d = fit(g_damped, M_damped)
    rms_d = math.sqrt(sse_d / len(g_damped))
    print(f"  free LSQ:  M_inf = {mi_d:.4f}  γ_e = {ge_d:.4f}  "
          f"RMS(M) = {rms_d:.4f}  → δ_∞ = {SLOPE * mi_d + 1.0:.3f}")
    print(f"  committed: M_inf = {CASE_7['M_inf']}  γ_e = {CASE_7['gamma_e']}  "
          f"(rounded for clean numerics)")

    # ----- summary stats -----
    print("\n=== Per-dataset survey ===")
    fmt = "  {:<32s} n={:>3d}  γ ∈ [{:>5.2f}, {:>5.2f}]  M ∈ [{:>5.3f}, {:>5.3f}]"
    print(fmt.format("Hansen+12/+14/+16 (lab dry Fo90)", len(g_hansen),
                      g_hansen.min(), g_hansen.max(), M_hansen.min(), M_hansen.max()))
    print(fmt.format("Bystricky 2000   (lab dry Fo90)",  len(g_byst),
                      g_byst.min(), g_byst.max(), M_byst.min(), M_byst.max()))
    print(fmt.format("Tasaka 2016      (lab WET Fo50)",  len(g_tas),
                      g_tas.min(), g_tas.max(), M_tas.min(), M_tas.max()))
    print(fmt.format("Boneh & Skemer14 (lab pre-CPO)",   len(g_bon),
                      g_bon.min(), g_bon.max(), M_bon.min(), M_bon.max()))
    print(fmt.format("Kumamoto+19   (natural Josephine)", len(g_kum),
                      g_kum.min(), g_kum.max(), M_kum.min(), M_kum.max()))
    print(fmt.format("Bernard+19   (natural xenoliths)", len(g_brn),
                      g_brn.min(), g_brn.max(), M_brn.min(), M_brn.max()))

    # ----- 2-panel SURVEY plot (full range + low-γ zoom) -----
    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(18.0, 7.0),
                                            constrained_layout=True)
    grid = np.linspace(0.0, 22.0, 400)

    def plot_all(ax):
        ax.plot(grid, m_of_gamma(grid, CASE_1["M_inf"], CASE_1["gamma_e"]),
                color="black", lw=2.0,
                label=f"case 1 fresh: M_inf={CASE_1['M_inf']}, γ_e={CASE_1['gamma_e']}, "
                      f"δ_∞={SLOPE * CASE_1['M_inf'] + 1.0:.2f}")
        ax.plot(grid, m_of_gamma(grid, CASE_7["M_inf"], CASE_7["gamma_e"]),
                color="purple", lw=2.0,
                label=f"case 7 damped: M_inf={CASE_7['M_inf']}, γ_e={CASE_7['gamma_e']}, "
                      f"δ_∞={SLOPE * CASE_7['M_inf'] + 1.0:.2f}")
        ax.scatter(g_hansen, M_hansen, s=55, marker="o",
                    facecolors="lightgray", edgecolors="black",
                    linewidths=0.7, alpha=0.9, zorder=3,
                    label=f"Hansen+12/+14/+16 (n={len(g_hansen)}, lab dry Fo90)")
        ax.scatter(g_byst, M_byst, s=110, marker="*",
                    facecolors="C0", edgecolors="black",
                    linewidths=0.8, zorder=5,
                    label=f"Bystricky 2000 (n={len(g_byst)}, lab dry Fo90)")
        ax.scatter(g_tas, M_tas, s=70, marker="s",
                    facecolors="none", edgecolors="C1",
                    linewidths=1.5, zorder=4,
                    label=f"Tasaka 2016 (n={len(g_tas)}, lab WET Fo50)")
        boneh_colors = {0: "k", 1: "C2", 2: "C3", 3: "C4"}
        boneh_labels = {0: "undef", 1: "perp", 2: "oblique", 3: "parallel"}
        for code in [0, 1, 2, 3]:
            mask = bon_geom == code
            if mask.any():
                ax.scatter(g_bon[mask], M_bon[mask], s=55, marker="^",
                            facecolors="none", edgecolors=boneh_colors[code],
                            linewidths=1.2, zorder=4,
                            label=f"Boneh&Skemer14 {boneh_labels[code]} (n={mask.sum()})")
        ax.scatter(g_kum, M_kum, s=60, marker="D",
                    facecolors="none", edgecolors="C5",
                    linewidths=1.3, zorder=4,
                    label=f"Kumamoto+19 (n={len(g_kum)}, natural Josephine)")
        ax.scatter(g_brn, M_brn, s=45, marker="v",
                    facecolors="none", edgecolors="C6",
                    linewidths=1.0, zorder=4,
                    label=f"Bernard+19 (n={len(g_brn)}, natural xenoliths)")
        ax.grid(alpha=0.3)
        ax.set_xlabel("strain γ (or simple-shear-equivalent)")
        ax.set_ylabel("M-index")
        ax.legend(loc="upper left", fontsize=8)

    plot_all(ax_full)
    ax_full.set_title("Full range γ ∈ [0, 22]")
    ax_full.set_xlim(-0.5, 22)
    ax_full.set_ylim(-0.02, 0.75)
    plot_all(ax_zoom)
    ax_zoom.set_title("Low-strain zoom γ ∈ [0, 6]")
    ax_zoom.set_xlim(-0.2, 6.0)
    ax_zoom.set_ylim(-0.02, 0.55)
    fig.suptitle(
        "Olivine M(γ) — full survey (Hansen-group + non-Hansen).  "
        f"Case 7 LSQ-best:  M_inf = {mi_d:.3f}, γ_e = {ge_d:.2f}",
        fontsize=11,
    )
    fig.savefig(SURVEY, dpi=140)
    print(f"\nWrote {SURVEY} ({SURVEY.stat().st_size / 1024:.0f} KB)")
    plt.close(fig)

    # ----- single-panel CASE 7 calibration plot -----
    fig2, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.scatter(g_tas, M_tas, s=70, marker="s",
                facecolors="none", edgecolors="C1", linewidths=1.5,
                label=f"Tasaka 2016 (n={len(g_tas)}, wet Fo50)")
    ax.scatter(g_bon, M_bon, s=55, marker="^",
                facecolors="none", edgecolors="C2", linewidths=1.2,
                label=f"Boneh & Skemer 2014 (n={len(g_bon)}, pre-CPO)")
    ax.scatter(g_kum, M_kum, s=60, marker="D",
                facecolors="none", edgecolors="C5", linewidths=1.3,
                label=f"Kumamoto+19 (n={len(g_kum)}, natural)")
    ax.scatter(g_brn, M_brn, s=45, marker="v",
                facecolors="none", edgecolors="C6", linewidths=1.0,
                label=f"Bernard+19 (n={len(g_brn)}, natural)")
    ax.plot(grid, m_of_gamma(grid, mi_d, ge_d), color="purple", lw=2.4,
             label=f"case 7 LSQ-best  M_inf={mi_d:.3f}, γ_e={ge_d:.2f}, "
                   f"δ_∞={SLOPE*mi_d+1:.2f}")
    ax.plot(grid, m_of_gamma(grid, CASE_7["M_inf"], CASE_7["gamma_e"]),
             color="purple", lw=1.4, ls="--",
             label=f"committed (rounded)  M_inf={CASE_7['M_inf']}, γ_e={CASE_7['gamma_e']}")
    ax.set_xlabel("shear strain γ (or simple-shear equiv.)")
    ax.set_ylabel("M-index")
    ax.set_xlim(-0.5, 22)
    ax.set_ylim(0, 0.45)
    ax.set_title(f"Damped olivine LSQ fit (case 7) — n = {len(g_damped)}, RMS(M) = {rms_d:.3f}")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig2.savefig(PLOT_PATH, dpi=140)
    print(f"Wrote {PLOT_PATH} ({PLOT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
