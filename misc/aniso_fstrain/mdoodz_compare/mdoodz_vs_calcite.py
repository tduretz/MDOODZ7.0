#!/usr/bin/env python3
"""Overlay MDOODZ HDF5 output from the 3 calcite CI tests on the calibration
plot. CSV-driven port of the old calcite/mdoodz_vs_calcite_comparison.py.

Reads HDF5 files written by:
    AnisotropyBenchmark.AniFstrainCalcite       (aniso_db = 2, T-averaged)
    AnisotropyBenchmark.AniFstrainCalciteLowT   (aniso_db = 4, T ≤ 650°C)
    AnisotropyBenchmark.AniFstrainCalciteHighT  (aniso_db = 5, T > 650°C)

Plots each per-step δ_mean as colored stars on its respective calibrated
curve, alongside Pieri+01 + Barnhoorn+04 lab data (Bruijn+11 retracted
post-2026-05-08 audit) and the existing min(FS_AR, 4) form for comparison.
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
    load_raw, j_to_m_skemer, fs_ar_from_gamma,
)

# Mirror MDLIB constants (post-2026-05-08 audit; Bruijn+11 retracted)
M_INF        = 0.41
GAMMA_E      = 3.15
M_INF_LOW    = 0.22
GAMMA_E_LOW  = 1.37
M_INF_HIGH   = 0.73
GAMMA_E_HIGH = 5.28
SLOPE        = 6.0
J_FILTER     = 15.0

OUT = HERE.parent / "img"
OUT.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT / "mdoodz_vs_calcite.png"


def m_of_gamma(g, mi, ge):
    return mi * (1.0 - np.exp(-g / ge))


def delta_calcite(g, mi, ge):
    return SLOPE * m_of_gamma(g, mi, ge) + 1.0


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
    ap.add_argument("--avg-build-dir",  type=Path, default=base / "AniFstrainCalcite")
    ap.add_argument("--low-build-dir",  type=Path, default=base / "AniFstrainCalciteLowT")
    ap.add_argument("--high-build-dir", type=Path, default=base / "AniFstrainCalciteHighT")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows_avg  = read_mdoodz_outputs(args.avg_build_dir)
    rows_low  = read_mdoodz_outputs(args.low_build_dir)
    rows_high = read_mdoodz_outputs(args.high_build_dir)

    print("MDOODZ outputs:")
    print(f"  case 2 (T-avg):     {len(rows_avg):>2d} pts")
    print(f"  case 4 (mid-low-T): {len(rows_low):>2d} pts")
    print(f"  case 5 (high-T):    {len(rows_high):>2d} pts")

    for label, rows, mi, ge in [
        ("case 2 T-avg",      rows_avg,  M_INF,      GAMMA_E),
        ("case 4 mid-low-T",  rows_low,  M_INF_LOW,  GAMMA_E_LOW),
        ("case 5 high-T",     rows_high, M_INF_HIGH, GAMMA_E_HIGH),
    ]:
        if not rows:
            continue
        g = np.array([r[0] for r in rows]); d = np.array([r[1] for r in rows])
        rel_err = (d - delta_calcite(g, mi, ge)) / np.maximum(delta_calcite(g, mi, ge), 1e-12)
        print(f"  {label} max rel error: {np.max(np.abs(rel_err)):.2e}")

    # ----- Lab data via CSVs (Bruijn+11 retracted post-2026-05-08 audit) -----
    pieri = load_raw("calcite/raw_Pieri_etal_2001_Tectonophysics.csv")
    barn  = load_raw("calcite/raw_Barnhoorn_etal_2004_JSG.csv")
    # Schuster+18 — γ-only annotation; no J/M tabulated (see audit/Schuster_etal_2018_annotation.md)
    schuster = load_raw("calcite/raw_Schuster_etal_2018_JSG.csv")
    g_schuster = np.asarray(schuster["gamma"], dtype=float)

    delta_p = SLOPE * j_to_m_skemer(pieri["J"]) + 1.0
    barn_kept = barn["J"] <= J_FILTER
    delta_n = SLOPE * j_to_m_skemer(barn["J"][barn_kept]) + 1.0

    fig, (ax_lo, ax_hi) = plt.subplots(1, 2, figsize=(14.0, 5.4),
                                        constrained_layout=True)
    grid = np.linspace(0.0, 16.0, 400)
    fs_ar_grid = fs_ar_from_gamma(grid)

    def plot_panel(ax, gamma_max):
        ax.scatter(pieri["gamma"], delta_p, s=70, marker="s",
                    facecolors="none", edgecolors="C0", linewidths=1.3, zorder=3,
                    label=f"Pieri+01 (n={len(pieri['gamma'])}, 727°C)")
        in_view = g_schuster[g_schuster <= gamma_max]
        if in_view.size:
            ax.scatter(in_view, np.full_like(in_view, 6.6), marker="|",
                        s=180, color="purple", linewidths=1.4, zorder=4,
                        label=f"Schuster+18 HPT 450°C γ-coverage (n={len(in_view)}/{len(g_schuster)}, J not tabulated)")
        for T_label, col in [(500, "C2"), (600, "C4"), (727, "C5")]:
            mask = (barn["T_C"][barn_kept] == T_label)
            if mask.any():
                gT = barn["gamma"][barn_kept][mask]
                dT = delta_n[mask]
                ax.scatter(gT, dT, s=60, marker="^",
                            facecolors="none", edgecolors=col, linewidths=1.2, zorder=3,
                            label=f"Barnhoorn+04 {T_label}°C (n={mask.sum()})")

        ax.plot(grid, delta_min(fs_ar_grid, 4.0), color="gray", lw=1.0,
                 linestyle=":", label="Old: δ = min(FS_AR, 4)")
        ax.plot(grid, delta_calcite(grid, M_INF_LOW,  GAMMA_E_LOW),
                 color="darkorange", lw=2.0,
                 label=f"case 4 mid-low-T (M_inf={M_INF_LOW}, γ_e={GAMMA_E_LOW}, δ_∞={SLOPE*M_INF_LOW+1:.2f})")
        ax.plot(grid, delta_calcite(grid, M_INF, GAMMA_E),
                 color="black", lw=2.0,
                 label=f"case 2 T-avg (M_inf={M_INF}, γ_e={GAMMA_E}, δ_∞={SLOPE*M_INF+1:.2f})")
        ax.plot(grid, delta_calcite(grid, M_INF_HIGH, GAMMA_E_HIGH),
                 color="firebrick", lw=2.0,
                 label=f"case 5 high-T (M_inf={M_INF_HIGH}, γ_e={GAMMA_E_HIGH}, δ_∞={SLOPE*M_INF_HIGH+1:.2f})")

        if rows_low:
            g = np.array([r[0] for r in rows_low]); d = np.array([r[1] for r in rows_low])
            ax.scatter(g, d, s=120, marker="*", color="darkorange",
                        edgecolors="black", linewidths=0.6, zorder=5,
                        label=f"MDOODZ case 4 (n={len(rows_low)})")
        if rows_avg:
            g = np.array([r[0] for r in rows_avg]); d = np.array([r[1] for r in rows_avg])
            ax.scatter(g, d, s=120, marker="*", color="black",
                        edgecolors="white", linewidths=0.6, zorder=5,
                        label=f"MDOODZ case 2 (n={len(rows_avg)})")
        if rows_high:
            g = np.array([r[0] for r in rows_high]); d = np.array([r[1] for r in rows_high])
            ax.scatter(g, d, s=120, marker="*", color="firebrick",
                        edgecolors="black", linewidths=0.6, zorder=5,
                        label=f"MDOODZ case 5 (n={len(rows_high)})")

        ax.set_xlabel("shear strain γ")
        ax.set_ylabel("δ (= 6·M + 1)")
        ax.set_xlim(-0.3, gamma_max)
        ax.set_ylim(0, 7)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=7.5)

    plot_panel(ax_lo, 5.0)
    ax_lo.set_title("Low-strain detail γ ∈ [0, 5]")
    plot_panel(ax_hi, 16.0)
    ax_hi.set_title("Full saturation range γ ∈ [0, 14]  (δ_∞: 2.33 / 3.46 / 5.40)")

    fig.suptitle(
        "MDOODZ ani_fstrain=2 vs calcite lab data — three regime calibrations",
        fontsize=11,
    )
    fig.savefig(args.out, dpi=140)
    print(f"\nWrote {args.out} ({args.out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
