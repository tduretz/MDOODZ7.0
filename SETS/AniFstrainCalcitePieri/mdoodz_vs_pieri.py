#!/usr/bin/env python3
"""Overlay MDOODZ AniFstrainCalcitePieri HDF5 output on Pieri+01 (Tecto 330, 119)
Carrara marble torsion datapoints at 727°C.

Reads HDF5 files from ../../cmake-exec/AniFstrainCalcitePieri/AniFstrainCalcitePieri/
and writes img/mdoodz_vs_pieri.png. Mirrors the structure of
misc/aniso_fstrain/mdoodz_compare/mdoodz_vs_olivine.py.
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
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "misc" / "aniso_fstrain" / "calibrate"))

from aniso_data import load_raw, j_to_m_skemer, fs_ar_from_gamma  # noqa: E402

# Mirror MDLIB/FlowLaws.c::anisoDelta_Calcite_HighT — must match aniso_db = 5.
M_INF   = 0.73
GAMMA_E = 5.28
SLOPE   = 6.0

OUT_DIR = HERE / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT_DIR / "mdoodz_vs_pieri.png"


def m_of_gamma(g, mi, ge):
    return mi * (1.0 - np.exp(-g / ge))


def delta_of_gamma(g, mi, ge):
    return SLOPE * m_of_gamma(g, mi, ge) + 1.0


def read_mdoodz_outputs(d: Path):
    outs = sorted(d.glob("Output*.gzip.h5"))
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
    default_build = (ROOT / "cmake-exec" / "AniFstrainCalcitePieri"
                     / "output")
    ap.add_argument("--build-dir", type=Path, default=default_build,
                    help="Directory containing Output*.gzip.h5")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows = read_mdoodz_outputs(args.build_dir)
    print(f"MDOODZ outputs: {len(rows)} from {args.build_dir}")
    if rows:
        g = np.array([r[0] for r in rows])
        d = np.array([r[1] for r in rows])
        d_form = delta_of_gamma(g, M_INF, GAMMA_E)
        rel_err = (d - d_form) / np.maximum(d_form, 1e-12)
        print(f"  max relative error vs case-5 formula: "
              f"{np.max(np.abs(rel_err)):.2e}")

    pieri = load_raw("calcite/raw_Pieri_etal_2001_Tectonophysics.csv")
    g_pieri = np.asarray(pieri["gamma"], dtype=float)
    delta_pieri = SLOPE * j_to_m_skemer(pieri["J"]) + 1.0

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.4), constrained_layout=True)
    grid = np.linspace(0.0, 12.0, 400)
    fs_ar_grid = fs_ar_from_gamma(grid)
    delta_inf = SLOPE * M_INF + 1.0

    ax.scatter(g_pieri, delta_pieri, s=110, marker="s",
               facecolors="none", edgecolors="C0", linewidths=1.6, zorder=3,
               label=f"Pieri+01 (n={len(g_pieri)}, Carrara 727°C)")
    ax.plot(grid, delta_of_gamma(grid, M_INF, GAMMA_E),
            color="firebrick", lw=2.0,
            label=f"case 5 high-T  (M_∞={M_INF}, γ_e={GAMMA_E}, "
                  f"δ_∞={delta_inf:.2f})")
    ax.axhline(delta_inf, color="firebrick", lw=0.5, ls=":", alpha=0.6)

    if rows:
        g_md = np.array([r[0] for r in rows])
        d_md = np.array([r[1] for r in rows])
        ax.scatter(g_md, d_md, s=170, marker="*", color="firebrick",
                   edgecolors="black", linewidths=0.7, zorder=6,
                   label=f"MDOODZ AniFstrainCalcitePieri (n={len(rows)})")

    ax.set_xlabel("shear strain γ")
    ax.set_ylabel("δ (= 6·M + 1)")
    ax.set_xlim(-0.5, 12.0)
    ax.set_ylim(0, 6)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title(
        "MDOODZ ani_fstrain=2 + aniso_db=5 vs Pieri+01 Carrara marble (727°C)",
        fontsize=11,
    )

    fig.savefig(args.out, dpi=140)
    size_kb = args.out.stat().st_size / 1024
    print(f"\nWrote {args.out} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
