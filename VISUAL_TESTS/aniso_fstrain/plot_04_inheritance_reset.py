#!/usr/bin/env python3
"""Test 4 plot — Inheritance reset-strain sweep.

δ(γ) curves for aniso_factor ∈ {1, 2, 4, 8} overlaid, with the Hansen analytic
curve (the f=1 case must overlay it). Annotated grey vertical band marks the
reset-strain window from Warren+08 / Boneh+15.
"""
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/04_inheritance"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/04_inheritance_reset.png"

FACTORS = [1, 2, 4, 8]
NT = 21
GAMMA_DOT = 2.0
DT = 0.5

C2 = 13.132
GAMMA_E = 3.96
DELTA_MAX = 14.13


def hansen(g):
    return C2 * (1.0 - np.exp(-g / GAMMA_E)) + 1.0


def collect(setname):
    outdir = RUNS / setname / "output"
    gs, ds = [], []
    for step in range(NT):
        fn = outdir / f"Output{step:05d}.gzip.h5"
        if not fn.exists():
            continue
        with h5py.File(fn, "r") as f:
            af = f["Centers/ani_fac"][:].mean()
        g = step * GAMMA_DOT * DT
        gs.append(g); ds.append(af)
    return np.array(gs), np.array(ds)


def main():
    fig, ax = plt.subplots(figsize=(10, 6.5))
    cmap = cm.viridis
    series = {}
    for k, f in enumerate(FACTORS):
        g, d = collect(f"f{f}")
        series[f] = (g, d)
        c = cmap(k / max(1, len(FACTORS) - 1))
        ax.plot(g, d, marker="o", lw=2.0, ms=6, color=c,
                label=fr"aniso_factor = {f} ($\delta_{{\rm init}}={f}$)")

    # Hansen analytic
    g_grid = np.linspace(0, max(series[1][0].max(), 20), 400)
    ax.plot(g_grid, hansen(g_grid), "k--", lw=1.5, label="Hansen analytic δ(γ)")
    ax.axhline(DELTA_MAX, color="black", ls=":", lw=1, alpha=0.6)
    ax.text(0.5, DELTA_MAX + 0.15, f"Hansen ceiling δ_max = {DELTA_MAX}",
            fontsize=8, color="black")

    # Reset-strain band
    ax.axvspan(1, 5, alpha=0.12, color="gray")
    ax.text(3, 1.5, "reset strain γ ∈ [1, 5]\n(Warren+08, Boneh+15)",
            ha="center", fontsize=9, color="dimgray", style="italic")

    ax.set_xlabel(r"$\gamma$ (cumulative shear strain)")
    ax.set_ylabel(r"$\delta$ = ani_fac (mesh)")
    ax.set_title("Test 4 — Inheritance reset-strain sweep (gate OFF: ani_relax_eps_max = −1)\n"
                 "Olivine simple shear, T=300K, aniso_angle=80°, ani_fstrain=3, aniso_db=1",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xlim(0, max(series[1][0])); ax.set_ylim(0, DELTA_MAX + 1.5)

    # Inset — log-scale residual |δ_initN − δ_init1|
    axins = ax.inset_axes([0.55, 0.32, 0.42, 0.35])
    g_ref, d_ref = series[1]
    for k, f in enumerate(FACTORS[1:], start=1):
        g, d = series[f]
        residual = np.abs(d - d_ref)
        # avoid log(0)
        c = cmap(k / max(1, len(FACTORS) - 1))
        axins.semilogy(g, np.maximum(residual, 1e-6), "o-", ms=4, color=c,
                       label=f"f={f}")
    axins.set_xlabel(r"$\gamma$", fontsize=8)
    axins.set_ylabel(r"$|\delta_{{\rm initN}} - \delta_{{\rm init1}}|$", fontsize=8)
    axins.set_title("inheritance residual (log)", fontsize=8)
    axins.tick_params(labelsize=7)
    axins.grid(alpha=0.3, which="both")
    axins.legend(loc="upper right", fontsize=6)

    cap = ("With gate OFF, all four cases ride the Hansen analytic curve at large γ. "
           "The δ_init=1 (f=1) trajectory tracks Hansen exactly from γ=0.\n"
           "Higher δ_init starts above the curve and converges as the strain history is "
           "rewritten — the 'inheritance bonus' decays with reset strain.")
    fig.text(0.5, -0.02, cap, ha="center", fontsize=9, style="italic")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
