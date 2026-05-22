#!/usr/bin/env python3
"""Test 2 plot — F-init geometry across aniso_angle.

2x3 panel grid showing, per aniso_angle ∈ {0, 30, 60, 80, 90, 120}:
  - imshow of Centers/aniso_delta (should be uniform ≈ 4)
  - principal-axis arrows of F at every 3rd cell
  - red dashed line: director (aniso_angle)
  - blue solid line: foliation (aniso_angle − 90°) → F-major should align
"""
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/this_script.py → ../..
ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "VISUAL_TESTS/aniso_fstrain/runs/02_finit_angle"
OUT = ROOT / "VISUAL_TESTS/img/aniso_fstrain/02_finit_angle_panel.png"

ANGLES = [0, 30, 60, 80, 90, 120]


def load_step0(ang):
    fn = RUNS / f"a{ang:03d}" / "output/Output00000.gzip.h5"
    with h5py.File(fn, "r") as f:
        d = f["Centers/aniso_delta"][:]
        Fxx = f["Centers/Fxx"][:]
        Fxz = f["Centers/Fxz"][:]
        Fzx = f["Centers/Fzx"][:]
        Fzz = f["Centers/Fzz"][:]
        # X is in Centers — gives the center coordinates flat
        try:
            X = f["Centers/X"][:]
        except KeyError:
            X = None
    return d, Fxx, Fxz, Fzx, Fzz


def principal_axes(Fxx, Fxz, Fzx, Fzz):
    """Return (major_vec, minor_vec, lambda1, lambda2) for F at each cell.
    Diagonalize FᵀF; eigenvectors are right-singular vectors → axes of strain ellipse in reference.
    Actually we want the deformed ellipse: eigenvectors of FFᵀ (left-singular vectors)
    point along the major/minor axes of the strain ellipse in deformed coords.
    """
    n = Fxx.shape[0]
    major = np.zeros((n, 2))
    minor = np.zeros((n, 2))
    lam1 = np.zeros(n)
    lam2 = np.zeros(n)
    for i in range(n):
        F = np.array([[Fxx[i], Fxz[i]], [Fzx[i], Fzz[i]]])
        B = F @ F.T
        w, V = np.linalg.eigh(B)
        # w is sorted ascending; major = V[:,1] (largest eigenvalue)
        major[i] = V[:, 1]
        minor[i] = V[:, 0]
        lam1[i] = np.sqrt(max(w[1], 0))
        lam2[i] = np.sqrt(max(w[0], 0))
    return major, minor, lam1, lam2


def main():
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9.0))
    Nx = 20  # Centers grid Nx-1 = 20 (Nx=21 in .txt)
    Nz = 20

    for ax, ang in zip(axes.flat, ANGLES):
        d, Fxx, Fxz, Fzx, Fzz = load_step0(ang)
        d2 = d.reshape((Nz, Nx))
        # imshow with vmin=3, vmax=5 — should be 4 everywhere
        im = ax.imshow(d2, origin="lower", extent=[-0.5, 0.5, -0.5, 0.5],
                       vmin=3.0, vmax=5.0, cmap="viridis", aspect="equal")

        # Principal axes — sample every 3rd cell
        xs = np.linspace(-0.5 + 0.5 / Nx, 0.5 - 0.5 / Nx, Nx)
        zs = np.linspace(-0.5 + 0.5 / Nz, 0.5 - 0.5 / Nz, Nz)
        XS, ZS = np.meshgrid(xs, zs)
        XS = XS.flatten(); ZS = ZS.flatten()
        maj, mn, l1, l2 = principal_axes(Fxx, Fxz, Fzx, Fzz)
        # arrow length scale
        arrL = 0.03
        step = 3
        idxs = np.arange(0, len(XS))
        # filter to every (step, step) grid pos
        keep = []
        for r in range(0, Nz, step):
            for c in range(0, Nx, step):
                keep.append(r * Nx + c)
        keep = np.array(keep)
        ax.quiver(XS[keep], ZS[keep], maj[keep, 0], maj[keep, 1],
                  color="black", scale=1/arrL, width=0.005, pivot="mid",
                  headlength=0, headaxislength=0, headwidth=0)
        ax.quiver(XS[keep], ZS[keep], mn[keep, 0]*0.5, mn[keep, 1]*0.5,
                  color="gray", scale=1/arrL, width=0.003, pivot="mid",
                  headlength=0, headaxislength=0, headwidth=0, alpha=0.5)

        # Director line — at aniso_angle from horizontal
        th_dir = np.radians(ang)
        ax.plot([-0.45 * np.cos(th_dir), 0.45 * np.cos(th_dir)],
                [-0.45 * np.sin(th_dir), 0.45 * np.sin(th_dir)],
                "r--", lw=1.5, label=f"director (aniso_angle = {ang}°)")
        # Foliation = aniso_angle − 90°
        th_fol = np.radians(ang - 90)
        ax.plot([-0.45 * np.cos(th_fol), 0.45 * np.cos(th_fol)],
                [-0.45 * np.sin(th_fol), 0.45 * np.sin(th_fol)],
                "b-", lw=2.0, label=f"foliation (aniso_angle − 90° = {ang - 90}°)")

        # F-major mean angle
        meanFxx = Fxx.mean(); meanFxz = Fxz.mean()
        meanFzx = Fzx.mean(); meanFzz = Fzz.mean()
        B = np.array([[meanFxx*meanFxx + meanFxz*meanFxz, meanFxx*meanFzx + meanFxz*meanFzz],
                      [meanFxx*meanFzx + meanFxz*meanFzz, meanFzx*meanFzx + meanFzz*meanFzz]])
        wB, VB = np.linalg.eigh(B)
        theta_F = np.degrees(np.arctan2(VB[1, 1], VB[0, 1]))

        ax.set_title(f"aniso_angle = {ang}°  |  θ_F = {theta_F:.1f}°  |  FS_AR ≈ {np.sqrt(wB[1]/wB[0]):.2f}",
                     fontsize=10)
        ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(loc="lower right", fontsize=7, framealpha=0.85)

    fig.suptitle(("Test 2 — F-init geometry across aniso_angle\n"
                  "(Hansen olivine, aniso_factor=4 → γ_eff≈1.03, FS_AR≈2.68).  "
                  "Black bars = F-major principal axis; gray = F-minor.\n"
                  "Visual acceptance: black bars lie on the BLUE foliation line, not the red director."),
                 fontsize=11, fontweight="bold")

    # Shared colorbar
    cbar = fig.colorbar(axes[0, 0].images[0], ax=axes, shrink=0.6,
                        fraction=0.025, pad=0.02)
    cbar.set_label(r"$\delta$ = Centers/aniso_delta")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
