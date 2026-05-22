#!/usr/bin/env python3
"""Overlay MDOODZ AniFstrainQuartzBlackford HDF5 output on Blackford+24
(Tectonics 43, e2023TC008166) Northern Snake Range quartzite datapoints.

Reads HDF5 files from ../../cmake-exec/AniFstrainQuartzBlackford/output/
and writes img/mdoodz_vs_blackford.png. Mirrors mdoodz_vs_pieri.py.
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

from aniso_data import load_raw, fs_ar_from_gamma  # noqa: E402

# Mirror MDLIB/FlowLaws.c::anisoDelta_Quartz_Blackford — must match aniso_db = 6.
M_INF   = 0.20
GAMMA_E = 4.00
SLOPE   = 11.0

OUT_DIR = HERE / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT_DIR / "mdoodz_vs_blackford.png"


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
    default_build = (ROOT / "cmake-exec" / "AniFstrainQuartzBlackford"
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
        print(f"  max relative error vs case-6 formula: "
              f"{np.max(np.abs(rel_err)):.2e}")

    bla = load_raw("quartz/derived_Blackford24_38pts.csv")
    g_bla = np.asarray(bla["gamma_simple_eq"], dtype=float)
    M_bla = np.asarray(bla["M_skemer_pfJ"], dtype=float)
    delta_bla = SLOPE * M_bla + 1.0
    domain = np.asarray(bla["strain_domain"], dtype=int)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.4), constrained_layout=True)
    grid = np.linspace(0.0, 12.0, 400)
    delta_inf = SLOPE * M_INF + 1.0

    domain_colors = {1: "C0", 2: "C1", 3: "C2", 4: "C3", 5: "C4"}
    for d_code, c in domain_colors.items():
        mask = domain == d_code
        if mask.any():
            ax.scatter(g_bla[mask], delta_bla[mask], s=70, marker="s",
                       facecolors="none", edgecolors=c, linewidths=1.4, zorder=3,
                       label=f"NSR domain {d_code} (n={mask.sum()})")

    ax.plot(grid, delta_of_gamma(grid, M_INF, GAMMA_E),
            color="darkgreen", lw=2.0,
            label=f"case 6 quartz  (M_∞={M_INF}, γ_e={GAMMA_E}, "
                  f"δ_∞={delta_inf:.2f})")
    ax.axhline(delta_inf, color="darkgreen", lw=0.5, ls=":", alpha=0.6)

    if rows:
        g_md = np.array([r[0] for r in rows])
        d_md = np.array([r[1] for r in rows])
        ax.scatter(g_md, d_md, s=170, marker="*", color="darkgreen",
                   edgecolors="black", linewidths=0.7, zorder=6,
                   label=f"MDOODZ AniFstrainQuartzBlackford (n={len(rows)})")

    ax.set_xlabel("shear strain γ_simple_eq")
    ax.set_ylabel("δ (= 11·M + 1)")
    ax.set_xlim(-0.5, 12.0)
    ax.set_ylim(0, 3.6)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(
        "MDOODZ ani_fstrain=2 + aniso_db=6 vs Blackford+24 NSR quartzite",
        fontsize=11,
    )

    fig.savefig(args.out, dpi=140)
    size_kb = args.out.stat().st_size / 1024
    print(f"\nWrote {args.out} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
