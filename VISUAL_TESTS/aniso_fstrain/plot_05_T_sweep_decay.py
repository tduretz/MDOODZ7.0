#!/usr/bin/env python3
"""Test 5 plot — T-sweep pure decay (Arrhenius staircase).

mechanical=1 with bkg_strain_rate≈1e-20 (effectively frozen F).
δ_init = 4 at step 0; evolves via the Boneh DiSRX exponential relaxation only.
Six curves for T ∈ {500, 700, 900, 1100, 1300, 1500} K.
"""
import math
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/05_T_sweep"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/05_T_sweep_decay.png"

Ts = [500, 700, 900, 1100, 1300, 1500]
NT = 400
DT_SI = 1e11  # s
TSCALE = DT_SI * 1  # 1e11 s per step → 400 steps × 1e11 s = 4e13 s ≈ 1.27 Myr
YR = 3.15576e7

# Olivine relaxation constants
Q = 133e3; M0 = 2e-11; mu = 50e9; b = 0.6e-9; R = 8.314
L_relax = 2e-3  # gs_ref default
DRHO_MIN = 1e10


def tau_relax_analytic(T_K):
    V = M0 * math.exp(-Q / (R * T_K)) * mu * b * b * DRHO_MIN
    return L_relax / V


def collect(T):
    outdir = RUNS / f"T{T:04d}" / "output"
    ts, ds = [], []
    for step in range(NT + 1):
        fn = outdir / f"Output{step:05d}.gzip.h5"
        if not fn.exists():
            continue
        with h5py.File(fn, "r") as f:
            ds.append(f["Centers/aniso_delta"][:].mean())
        ts.append(step * DT_SI)
    return np.array(ts), np.array(ds)


def main():
    fig, ax = plt.subplots(figsize=(11, 6.5))
    cmap = cm.coolwarm
    for k, T in enumerate(Ts):
        t, d = collect(T)
        c = cmap(k / max(1, len(Ts) - 1))
        tau = tau_relax_analytic(T)
        label = fr"$T = {T}$ K, $\tau_{{\rm relax}}={tau:.1e}$ s ({tau/YR:.1e} yr)"
        ax.semilogx(np.maximum(t / YR, 1e-6), d, "-", color=c, lw=2, label=label)

    ax.axhline(1.0, color="black", ls=":", lw=1, alpha=0.5)
    ax.axhline(4.0, color="black", ls=":", lw=1, alpha=0.5)
    ax.text(1e-6, 4.05, r"$\delta_{\rm init} = 4$", fontsize=9)
    ax.text(1e-6, 1.05, r"$\delta_{\rm eq} = 1$", fontsize=9)

    ax.set_xlabel("Time (years)")
    ax.set_ylabel(r"$\delta$ = ani_fac (mesh)")
    ax.set_ylim(0.9, 4.2)
    ax.set_xlim(1e-3, 1.5e7)
    ax.set_title("Test 5 — T-sweep pure δ-decay (Arrhenius staircase)\n"
                 "mechanical=1 with bkg_strain_rate≈0 (frozen F), ani_relax_eps_max=−1, "
                 "L_relax = gs_ref = 2 mm, Δρ = drho_min = 1e10 m⁻²",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=8.5, framealpha=0.92)
    ax.grid(alpha=0.3, which="both")

    cap = ("Caveat (per design): with mechanical=1 + tiny BSR, strain_pwl≈0 → Δρ=drho_min=1e10 m⁻² (NOT drho_max=1e13).\n"
           "This shifts the Arrhenius staircase UP in T by ~3 decades. The half-decay point sits at T≈1000 K (not 850 K). "
           "At T=500 K δ is frozen; at T=1300 K relaxation is faster than dt → near-instantaneous.")
    fig.text(0.5, -0.02, cap, ha="center", fontsize=9, style="italic")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
