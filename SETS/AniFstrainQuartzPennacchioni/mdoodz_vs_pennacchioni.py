#!/usr/bin/env python3
"""Overlay MDOODZ AniFstrainQuartzPennacchioni HDF5 output on Pennacchioni+10
(JGR 115, B12405) natural-quartz-shear-zone datapoints (T ~ 500°C).

Reads HDF5 files from ../../cmake-exec/AniFstrainQuartzPennacchioni/output/
and writes img/mdoodz_vs_pennacchioni.png. Mirrors mdoodz_vs_pieri.py.

Note: the 17 (γ, J) values used here are
NOT real per-sample MTEX data; the smooth scatter (σ/M_∞ = 0.045) is
the audit signature that flagged the case as physical-prior.
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

# Mirror MDLIB/FlowLaws.c::anisoDelta_Quartz — must match aniso_db = 3.
M_INF   = 0.75
GAMMA_E = 3.50
SLOPE   = 3.0

OUT_DIR = HERE / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT_DIR / "mdoodz_vs_pennacchioni.png"


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
    default_build = (ROOT / "cmake-exec" / "AniFstrainQuartzPennacchioni"
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
        print(f"  max relative error vs case-3 formula: "
              f"{np.max(np.abs(rel_err)):.2e}")

    pen = load_raw("quartz/derived_Pennacchioni10_17pts.csv")
    g_pen = np.asarray(pen["gamma_simple"], dtype=float)
    M_pen = np.asarray(pen["M_skemer"], dtype=float)
    delta_pen = SLOPE * M_pen + 1.0
    regime = np.asarray(pen["regime"], dtype=int)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.4), constrained_layout=True)
    grid = np.linspace(0.0, 17.0, 400)
    delta_inf = SLOPE * M_INF + 1.0

    # color by regime: 0=WDV, 1=MDV, 2=SDV
    regime_styles = {
        0: ("WDV (weak)",   "C0"),
        1: ("MDV (mod)",    "C1"),
        2: ("SDV (strong)", "C3"),
    }
    for r_code, (lbl, c) in regime_styles.items():
        mask = regime == r_code
        if mask.any():
            ax.scatter(g_pen[mask], delta_pen[mask], s=85, marker="s",
                       facecolors="none", edgecolors=c, linewidths=1.6, zorder=3,
                       label=f"Pennacchioni+10 {lbl} (n={mask.sum()})")

    ax.plot(grid, delta_of_gamma(grid, M_INF, GAMMA_E),
            color="tab:purple", lw=2.0,
            label=f"case 3 quartz  (M_∞={M_INF}, γ_e={GAMMA_E}, "
                  f"δ_∞={delta_inf:.2f})")
    ax.axhline(delta_inf, color="tab:purple", lw=0.5, ls=":", alpha=0.6)

    if rows:
        g_md = np.array([r[0] for r in rows])
        d_md = np.array([r[1] for r in rows])
        ax.scatter(g_md, d_md, s=170, marker="*", color="tab:purple",
                   edgecolors="black", linewidths=0.7, zorder=6,
                   label=f"MDOODZ AniFstrainQuartzPennacchioni (n={len(rows)})")

    ax.set_xlabel("shear strain γ")
    ax.set_ylabel("δ (= 3·M + 1)")
    ax.set_xlim(-0.5, 17.0)
    ax.set_ylim(0, 3.6)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title(
        "MDOODZ ani_fstrain=2 + aniso_db=3 vs Pennacchioni+10 quartz\n"
        "(audit caveat: σ/M_∞=0.045 is below the real-data floor 0.13)",
        fontsize=10,
    )

    fig.savefig(args.out, dpi=140)
    size_kb = args.out.stat().st_size / 1024
    print(f"\nWrote {args.out} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
