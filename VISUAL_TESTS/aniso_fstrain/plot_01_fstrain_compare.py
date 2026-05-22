#!/usr/bin/env python3
"""Test 1 plot — ani_fstrain 1 vs 2 vs 3 comparison.

Reads HDF5 outputs from
    VISUAL_TESTS/aniso_fstrain/runs/01_fstrain_compare/f{1,2,3}/output/
and produces a two-panel PNG showing δ(γ) for each ani_fstrain value
and FS_AR(γ) (shared kinematics).
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/01_fstrain_compare"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/01_fstrain_compare.png"

CASES = {
    "F1 (ani_fstrain=1)\nDefault MDOODZ — raw FS_AR capped": "f1",
    "F2 (ani_fstrain=2)\nHansen-strict — δ = aniso_delta_fn(FS_AR)": "f2",
    "F3 (ani_fstrain=3)\nHansen + Boneh δ-relaxation":          "f3",
}

HANSEN_CEILING = 14.13
GAMMA_DOT = 2.0  # 2·bkg_strain_rate
DT = 0.5
NT = 17


def collect(setname):
    outdir = RUNS / setname / "output"
    gammas, deltas, fsars = [], [], []
    for step in range(NT):
        fn = outdir / f"Output{step:05d}.gzip.h5"
        if not fn.exists():
            continue
        with h5py.File(fn, "r") as f:
            ani_fac = f["Centers/ani_fac"][:].mean()
            Fxx = f["Centers/Fxx"][:].mean()
            Fxz = f["Centers/Fxz"][:].mean()
            Fzx = f["Centers/Fzx"][:].mean()
            Fzz = f["Centers/Fzz"][:].mean()
        F = np.array([[Fxx, Fxz], [Fzx, Fzz]])
        s = np.linalg.svd(F, compute_uv=False)
        FS_AR = s[0] / max(s[1], 1e-30)
        gamma = step * GAMMA_DOT * DT  # γ = γ̇·t = γ̇·step·dt
        gammas.append(gamma)
        deltas.append(ani_fac)
        fsars.append(FS_AR)
    return np.array(gammas), np.array(deltas), np.array(fsars)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax_d, ax_fs = axes
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    markers = ["o", "s", "^"]
    for (label, setname), c, m in zip(CASES.items(), colors, markers):
        g, d, fs = collect(setname)
        ax_d.plot(g, d, marker=m, color=c, label=label, lw=1.5, ms=6, alpha=0.9)
        if "F1" in label:
            ax_fs.plot(g, fs, marker=m, color="black", lw=1.5, ms=6,
                       label="FS_AR(γ) shared kinematics")

    ax_d.axhline(HANSEN_CEILING, ls="--", color="gray", lw=1)
    ax_d.text(0.05, HANSEN_CEILING - 0.4, f"Hansen ceiling δ_max ≈ {HANSEN_CEILING}",
              transform=ax_d.get_yaxis_transform(),
              ha="left", va="top", fontsize=8, color="gray")
    ax_d.set_xlabel(r"$\gamma$ (shear strain)")
    ax_d.set_ylabel(r"$\delta$ = ani_fac (mesh)")
    ax_d.set_title(r"Aniso factor $\delta(\gamma)$ — three ani_fstrain modes")
    ax_d.set_yscale("symlog", linthresh=1)
    ax_d.legend(loc="lower right", fontsize=8, frameon=False)
    ax_d.grid(alpha=0.3)

    ax_fs.set_yscale("log")
    ax_fs.set_xlabel(r"$\gamma$")
    ax_fs.set_ylabel(r"FS_AR = $\lambda_1/\lambda_2$ from F")
    ax_fs.set_title("Finite-strain aspect ratio (identical across F1/F2/F3)")
    ax_fs.grid(alpha=0.3, which="both")
    ax_fs.legend(loc="lower right", fontsize=8, frameon=False)

    fig.suptitle("Test 1 — ani_fstrain 1 vs 2 vs 3 (Olivine, cold T=300K, simple shear γ̇=2)",
                 fontsize=12, fontweight="bold")
    cap = ("F1 grows unbounded (capped at ani_fac_max=100). F2 saturates at the Hansen ceiling.\n"
           "F3 overlays F2 exactly under cold T + strain-rate gate (DiSRX relaxation inactive — the "
           "operator split collapses to F2). All three share the same FS_AR(γ).")
    fig.text(0.5, -0.02, cap, ha="center", va="top", fontsize=9, style="italic")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
