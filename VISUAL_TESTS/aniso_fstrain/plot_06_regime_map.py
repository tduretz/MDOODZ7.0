#!/usr/bin/env python3
"""Test 6 plot — 2D T × ε̇ regime map for ani_fstrain=3 final δ.

Reads regime_map_runs/results.json, plots heatmap with regime annotations.
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/06_regime_map"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/06_regime_map.png"

TS = [500, 700, 900, 1100, 1300, 1500]
EDOTS = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10]

# Boneh DiSRX kinetics constants (olivine)
Q = 133e3; M0 = 2e-11; mu = 50e9; b = 0.6e-9; R = 8.314
L_relax = 2e-3
DRHO_MIN = 1e10
DRHO_MAX = 1e13


def tau_min(T):
    """Lower bound — fully strained marker (drho_max)."""
    V = M0 * math.exp(-Q / (R * T)) * mu * b * b * DRHO_MAX
    return L_relax / V


def main():
    res = json.loads((RUNS / "results.json").read_text())
    # Build matrix
    M = np.full((len(TS), len(EDOTS)), np.nan)
    for r in res:
        i = TS.index(int(r["T"]))
        j = EDOTS.index(float(r["edot"]))
        M[i, j] = r["delta"]

    fig = plt.figure(figsize=(11, 7.5))
    ax = fig.add_axes([0.08, 0.10, 0.65, 0.78])

    log_e = np.log10(np.array(EDOTS))
    extent = [log_e[0] - 0.5, log_e[-1] + 0.5,
              TS[0] - 100, TS[-1] + 100]
    # T axis in C
    T_C = [T - 273 for T in TS]
    im = ax.imshow(M, origin="lower", aspect="auto",
                   extent=[log_e[0] - 0.5, log_e[-1] + 0.5, T_C[0] - 100, T_C[-1] + 100],
                   cmap="viridis", vmin=1.0, vmax=13.13)

    # Annotate each cell with δ value
    for i, Tc in enumerate(T_C):
        for j, le in enumerate(log_e):
            d = M[i, j]
            color = "white" if d < 7 else "black"
            ax.text(le, Tc, f"{d:.1f}", ha="center", va="center",
                    fontsize=8.5, color=color)

    # Regime corner annotations
    corners = [
        (-11, T_C[0] - 50, "FROZEN / SATURATED\n(δ → Hansen ceiling)", "white"),
        (-15.5, T_C[-1] + 30, "ERASED\n(hot + quiescent)", "black"),
        (-13, T_C[-1] + 30, "COMPETITION\n(reset-strain)", "black"),
    ]
    for lx, ly, txt, col in corners:
        ax.annotate(txt, xy=(lx, ly), ha="center", fontsize=9, color=col,
                    fontweight="bold")

    # Analytic contour τ_min(T)·ε̇ = 1
    # Equivalent: log10(ε̇) = -log10(τ_min(T))
    T_dense = np.linspace(500, 1500, 200)
    log_edot_eq1 = np.array([-math.log10(tau_min(T)) for T in T_dense])
    ax.plot(log_edot_eq1, T_dense - 273, "k--", lw=2,
            label=r"$\tau_{\rm relax}(T)\cdot \dot\epsilon = 1$ (drho_max)")

    # Gate threshold ε̇ = 1e-13
    ax.axvline(-13, color="red", ls=":", lw=2, label=r"$\dot\epsilon_{\rm max} = 10^{-13}$ s$^{-1}$ (gate threshold; OFF here)")

    # Rifting reference point
    ax.plot(np.log10(3e-15), 1100 - 273, "*", ms=22, color="red",
            markeredgecolor="white", markeredgewidth=1.5, zorder=10,
            label="rifting reference (T=1100K, ε̇=3e-15)")

    ax.set_xlabel(r"$\log_{10}\dot\epsilon$ (s$^{-1}$)")
    ax.set_ylabel(r"T (°C)")
    ax.set_title("Test 6 — 2D regime map: final δ vs T and ε̇\n"
                 "Olivine simple shear to γ=5, ani_fstrain=3, ani_relax_eps_max=−1, aniso_factor=4 (init)",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.92)
    ax.set_xticks(log_e)
    ax.set_xticklabels([f"{int(x)}" for x in log_e])
    ax.set_yticks(T_C)

    # Colorbar
    cax = fig.add_axes([0.78, 0.10, 0.025, 0.78])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(r"final $\delta$ at $\gamma=5$", fontsize=10)

    # Inset — 1D cut at T=1100K
    axi = fig.add_axes([0.10, 0.62, 0.18, 0.18])
    j_1100 = TS.index(1100)
    axi.semilogx(EDOTS, M[j_1100, :], "o-", color="darkred", lw=2)
    axi.set_xlabel(r"$\dot\epsilon$ (s$^{-1}$)", fontsize=8)
    axi.set_ylabel(r"$\delta_{\rm final}$", fontsize=8)
    axi.set_title("cut at T = 1100 K (1827°C)", fontsize=8)
    axi.tick_params(labelsize=7)
    axi.grid(alpha=0.3, which="both")
    axi.set_ylim(0, 14)

    cap = ("FROZEN: cold T → τ_relax much larger than t_total. SATURATED: high ε̇ → FS_AR drives δ to Hansen ceiling faster than relaxation.\n"
           "ERASED: hot + slow → relaxation outraces FS_AR growth, δ drops. COMPETITION: both rates comparable — the rifting-relevant regime.\n"
           "Dashed black: kinetic boundary τ_relax·ε̇ = 1 (with drho_max). Red dotted: gate threshold (informational; gate is OFF in this map).")
    fig.text(0.5, 0.01, cap, ha="center", fontsize=9, style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
