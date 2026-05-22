#!/usr/bin/env python3
"""Test 3 plot — δ round-trip across aniso_factor.

Two stacked panels:
  Top: extracted aniso_delta vs input aniso_factor with y=x identity line.
  Bottom: Hansen analytic δ(γ_eff) curve with 9 red dots at (γ_eff_pred, δ_extracted).
"""
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/03_finit_factor"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/03_finit_factor_roundtrip.png"

FACTORS = [1, 2, 3, 4, 6, 8, 10, 12, 14]

# Hansen olivine constants (from FlowLaws.c::ReadDataAnisotropy)
C2 = 13.132
GAMMA_E = 3.96
DELTA_MAX = 14.13


def hansen_delta(g):
    return C2 * (1.0 - np.exp(-g / GAMMA_E)) + 1.0


def hansen_inv(delta):
    # γ_eff = −γ_e · ln(1 − (δ−1)/c2)
    arg = 1.0 - (delta - 1.0) / C2
    arg = np.clip(arg, 1e-30, None)
    return -GAMMA_E * np.log(arg)


def fs_ar_from_gamma(g):
    """FS_AR for simple shear at γ_eff: F = [[1,γ],[0,1]] → λ1/λ2."""
    return (g + np.sqrt(g * g + 4.0)) / 2.0 / ((-g + np.sqrt(g * g + 4.0)) / 2.0)


def main():
    extracted = []
    fs_ars = []
    for f in FACTORS:
        fn = RUNS / f"f{f:02d}/output/Output00000.gzip.h5"
        with h5py.File(fn, "r") as fh:
            d = fh["Centers/aniso_delta"][:].mean()
            Fxx = fh["Centers/Fxx"][:].mean()
            Fxz = fh["Centers/Fxz"][:].mean()
            Fzx = fh["Centers/Fzx"][:].mean()
            Fzz = fh["Centers/Fzz"][:].mean()
        F = np.array([[Fxx, Fxz], [Fzx, Fzz]])
        s = np.linalg.svd(F, compute_uv=False)
        extracted.append(d)
        fs_ars.append(s[0] / s[1])

    extracted = np.array(extracted)
    fs_ars = np.array(fs_ars)
    g_eff_pred = hansen_inv(np.array(FACTORS, dtype=float))

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 9))

    # Top — y=x round-trip
    ax_top.plot([0, 15], [0, 15], "k--", lw=1, alpha=0.6, label="y = x identity")
    ax_top.axhline(DELTA_MAX, color="gray", ls=":", lw=1)
    ax_top.text(0.5, DELTA_MAX + 0.2, f"Hansen ceiling δ_max = {DELTA_MAX}", fontsize=8, color="gray")
    ax_top.scatter(FACTORS, extracted, c="red", s=60, zorder=5, label="extracted aniso_delta")
    ax_top.set_xlabel("input aniso_factor")
    ax_top.set_ylabel("extracted Centers/aniso_delta at step 0")
    ax_top.set_title("Round-trip: extracted δ vs input aniso_factor")
    ax_top.set_xlim(0, 15); ax_top.set_ylim(0, 15)
    ax_top.grid(alpha=0.3)
    ax_top.legend(loc="lower right", fontsize=9)
    # Annotate each point with (g_eff, FS_AR)
    for k, f in enumerate(FACTORS):
        ax_top.annotate(f"f={f}\nγ_eff={g_eff_pred[k]:.2f}",
                        (f, extracted[k]),
                        xytext=(4, -10), textcoords="offset points", fontsize=7,
                        color="darkred")

    # Bottom — Hansen analytic + dots
    g_grid = np.linspace(0, 25, 400)
    ax_bot.plot(g_grid, hansen_delta(g_grid), "b-", lw=1.5,
                label=fr"$\delta(\gamma_{{\rm eff}}) = {C2}\,(1-e^{{-\gamma/{GAMMA_E}}})+1$")
    ax_bot.axhline(DELTA_MAX, color="gray", ls=":", lw=1)
    ax_bot.scatter(g_eff_pred, extracted, c="red", s=60, zorder=5,
                   label="9 SET runs at (γ_eff_pred, δ_extracted)")
    for k, f in enumerate(FACTORS):
        ax_bot.annotate(f"f={f}",
                        (g_eff_pred[k], extracted[k]),
                        xytext=(6, -4), textcoords="offset points", fontsize=8,
                        color="darkred")
    ax_bot.set_xlabel(r"$\gamma_{\rm eff}$ (= aniso_delta_inv_hansen(aniso_factor))")
    ax_bot.set_ylabel(r"$\delta$ at step 0")
    ax_bot.set_title("Hansen analytic curve with extracted δ values overlaid")
    ax_bot.set_xlim(0, 25); ax_bot.set_ylim(0, 15.5)
    ax_bot.grid(alpha=0.3)
    ax_bot.legend(loc="lower right", fontsize=9)

    # Twin x: FS_AR for simple shear
    ax_bot_twin = ax_bot.twiny()
    g_ticks = np.array([0, 1, 2, 3, 4, 6, 10, 20])
    fs_ticks = fs_ar_from_gamma(g_ticks)
    ax_bot_twin.set_xticks(g_ticks)
    ax_bot_twin.set_xticklabels([f"{x:.1f}" for x in fs_ticks])
    ax_bot_twin.set_xlim(ax_bot.get_xlim())
    ax_bot_twin.set_xlabel(r"corresponding simple-shear FS_AR = $\lambda_1/\lambda_2$")

    fig.suptitle("Test 3 — δ round-trip across aniso_factor (Hansen olivine inverse)",
                 fontsize=12, fontweight="bold")
    cap = ("F_init is constructed from γ_eff = aniso_delta_inv_hansen(aniso_factor) at θ = aniso_angle−90°.\n"
           "Re-extracting δ from F via the Hansen forward law closes the loop exactly. "
           "Saturation lurks at δ_max=14.13.")
    fig.text(0.5, 0.005, cap, ha="center", fontsize=9, style="italic")

    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
