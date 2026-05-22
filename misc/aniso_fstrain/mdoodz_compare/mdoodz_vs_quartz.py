#!/usr/bin/env python3
"""Overlay MDOODZ HDF5 output (from AniFstrainQuartz + AniFstrainQuartzBlackford
GTests) on the quartz calibration plot. CSV-driven port.

Reads HDF5 files written by:
    AnisotropyBenchmark.AniFstrainQuartz           (aniso_db = 3, Pennacchioni+10)
    AnisotropyBenchmark.AniFstrainQuartzBlackford  (aniso_db = 6, Blackford+24)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "calibrate"))

from aniso_data import (  # noqa: E402
    load_raw, j_to_m_skemer, fs_ar_from_gamma, gamma_oct_to_gamma_simple,
)

# Mirror MDLIB constants
M_INF       = 0.75
GAMMA_E     = 3.50
SLOPE       = 3.0
M_INF_BLA   = 0.20
GAMMA_E_BLA = 4.00
SLOPE_BLA   = 11.0

OUT = HERE.parent / "img"
OUT.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT / "mdoodz_vs_quartz.png"


def m_of_gamma(g, mi, ge):
    return mi * (1.0 - np.exp(-g / ge))


def delta_quartz(g, mi, ge, slope):
    return slope * m_of_gamma(g, mi, ge) + 1.0


def delta_min(fs_ar, ani_fac_max):
    return np.minimum(fs_ar, ani_fac_max)


def read_mdoodz_outputs(d: Path):
    outs = sorted(d.glob("Output*.gzip.h5"))
    if not outs:
        return []
    rows = []
    for p in outs:
        step = int(p.name.replace("Output", "").replace(".gzip.h5", ""))
        if step == 0:
            continue
        gamma = float(step - 1)
        with h5py.File(p, "r") as f:
            ani = np.array(f["Centers"]["ani_fac"]).ravel()
            ani = ani[ani > 0.0]
            rows.append((gamma, float(ani.mean()),
                         float(ani.min()), float(ani.max())))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    base = HERE.parent.parent.parent / "cmake-build-test" / "TESTS"
    ap.add_argument("--pen-build-dir", type=Path,
                    default=base / "AniFstrainQuartz")
    ap.add_argument("--bla-build-dir", type=Path,
                    default=base / "AniFstrainQuartzBlackford")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows_pen = read_mdoodz_outputs(args.pen_build_dir)
    rows_bla = read_mdoodz_outputs(args.bla_build_dir)

    print("MDOODZ outputs:")
    print(f"  case 3 Pennacchioni: {len(rows_pen):>2d} pts")
    print(f"  case 6 Blackford:    {len(rows_bla):>2d} pts")

    for label, rows, mi, ge, sl in [
        ("case 3 Pennacchioni", rows_pen, M_INF,     GAMMA_E,     SLOPE),
        ("case 6 Blackford",    rows_bla, M_INF_BLA, GAMMA_E_BLA, SLOPE_BLA),
    ]:
        if not rows:
            continue
        g = np.array([r[0] for r in rows]); d = np.array([r[1] for r in rows])
        d_form = delta_quartz(g, mi, ge, sl)
        rel_err = (d - d_form) / np.maximum(d_form, 1e-12)
        print(f"  {label} max rel error: {np.max(np.abs(rel_err)):.2e}")

    # ----- Lab data via CSVs -----
    pen = load_raw("quartz/raw_Pennacchioni_etal_2010_JGR.csv")
    g_pen = pen["gamma"]; J_pen = pen["J"]
    regimes = pen["regime"].astype(int)
    delta_pen = SLOPE * j_to_m_skemer(J_pen) + 1.0

    bla = load_raw("quartz/raw_Blackford_etal_2024_Tectonics.csv")
    g_bla_eq = gamma_oct_to_gamma_simple(bla["gamma_oct"])
    delta_bla = SLOPE_BLA * j_to_m_skemer(bla["pfJ"]) + 1.0

    grid = np.linspace(0.0, 16.0, 400)
    fs_ar_grid = fs_ar_from_gamma(grid)
    regime_labels = {0: "WDV (γ<1)", 1: "MDV (2<γ<3)", 2: "SDV (γ=4-15)"}
    regime_markers = {0: "o", 1: "s", 2: "^"}

    fig, (ax_lo, ax_hi) = plt.subplots(1, 2, figsize=(14.0, 5.4),
                                        constrained_layout=True)

    def plot_panel(ax, gamma_max):
        for r in [0, 1, 2]:
            mask = regimes == r
            if mask.any():
                ax.scatter(g_pen[mask], delta_pen[mask], s=60,
                            marker=regime_markers[r],
                            facecolors="none", edgecolors="C3",
                            linewidths=1.3, zorder=3,
                            label=f"Pennacchioni+10 {regime_labels[r]} (n={mask.sum()})")
        ax.scatter(g_bla_eq, delta_bla, s=45, marker="D",
                    facecolors="none", edgecolors="C0",
                    linewidths=1.0, zorder=3,
                    label=f"Blackford+24 (n={len(g_bla_eq)}, γ_eq from γ_oct)")

        ax.plot(grid, delta_min(fs_ar_grid, 3.0), color="gray", lw=1.0,
                 linestyle=":", label="Old: δ = min(FS_AR, 3)")
        ax.plot(grid, delta_quartz(grid, M_INF, GAMMA_E, SLOPE),
                 color="C3", lw=2.0,
                 label=f"case 3 Pennacchioni (M_inf={M_INF}, γ_e={GAMMA_E}, "
                       f"slope={SLOPE}, δ_∞={SLOPE*M_INF+1:.2f})")
        ax.plot(grid, delta_quartz(grid, M_INF_BLA, GAMMA_E_BLA, SLOPE_BLA),
                 color="C0", lw=2.0,
                 label=f"case 6 Blackford (M_inf={M_INF_BLA}, γ_e={GAMMA_E_BLA}, "
                       f"slope={SLOPE_BLA}, δ_∞={SLOPE_BLA*M_INF_BLA+1:.2f})")

        if rows_pen:
            g = np.array([r[0] for r in rows_pen]); d = np.array([r[1] for r in rows_pen])
            ax.scatter(g, d, s=140, marker="*", color="C3",
                        edgecolors="black", linewidths=0.6, zorder=5,
                        label=f"MDOODZ case 3 (n={len(rows_pen)})")
        if rows_bla:
            g = np.array([r[0] for r in rows_bla]); d = np.array([r[1] for r in rows_bla])
            ax.scatter(g, d, s=140, marker="*", color="C0",
                        edgecolors="black", linewidths=0.6, zorder=5,
                        label=f"MDOODZ case 6 (n={len(rows_bla)})")

        ax.set_xlabel("shear strain γ")
        ax.set_ylabel("δ (viscous-anisotropy magnitude)")
        ax.set_xlim(-0.3, gamma_max)
        ax.set_ylim(0, 4.5)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=7)

    plot_panel(ax_lo, 6.0)
    ax_lo.set_title("Low-strain detail γ ∈ [0, 6]")
    plot_panel(ax_hi, 16.0)
    ax_hi.set_title(f"Full saturation γ ∈ [0, 16]  "
                     f"(δ_∞: case 3 = {SLOPE*M_INF+1:.2f}, case 6 = {SLOPE_BLA*M_INF_BLA+1:.2f})")
    fig.suptitle(
        "MDOODZ ani_fstrain=2 vs quartz lab data — case 3 (Pennacchioni ODF J) + case 6 (Blackford pfJ)",
        fontsize=11,
    )
    fig.savefig(args.out, dpi=140)
    print(f"\nWrote {args.out} ({args.out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
