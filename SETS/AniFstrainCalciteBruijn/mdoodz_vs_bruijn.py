#!/usr/bin/env python3
"""Overlay MDOODZ AniFstrainCalciteBruijn HDF5 output on Bruijn+11
(Tectonophysics 503, 75; doi:10.1016/j.tecto.2010.09.029) Carrara marble
torsion datapoints at 727°C.

Reads HDF5 files from ../../cmake-exec/AniFstrainCalciteBruijn/output/ and
writes img/mdoodz_vs_bruijn.png. Mirrors mdoodz_vs_pieri.py.

Audit caveat (per audit warning embedded in
misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv):
the 3 (γ, J) points are NOT sourced in Bruijn's published paper. Bruijn's
Fig. 8 reproduces Pieri+01 / Barnhoorn+04 pole figures with no numerical J
labels. The values appear to be greyscale visual estimates and do NOT enter
the case-2 LSQ fit. This SET visualises the audit finding by overlaying the
3 audit-flagged points on the case-2 calibrated curve.
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

from aniso_data import load_raw, j_to_m_skemer  # noqa: E402

# Mirror MDLIB/FlowLaws.c::anisoDelta_Calcite — must match aniso_db = 2.
M_INF   = 0.41
GAMMA_E = 3.18
SLOPE   = 6.0

OUT_DIR = HERE / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT_DIR / "mdoodz_vs_bruijn.png"


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
    default_build = (ROOT / "cmake-exec" / "AniFstrainCalciteBruijn"
                     / "output")
    ap.add_argument("--build-dir", type=Path, default=default_build)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows = read_mdoodz_outputs(args.build_dir)
    print(f"MDOODZ outputs: {len(rows)} from {args.build_dir}")
    if rows:
        g = np.array([r[0] for r in rows])
        d = np.array([r[1] for r in rows])
        d_form = delta_of_gamma(g, M_INF, GAMMA_E)
        rel_err = (d - d_form) / np.maximum(d_form, 1e-12)
        print(f"  max relative error vs case-2 formula: "
              f"{np.max(np.abs(rel_err)):.2e}")

    bru = load_raw("calcite/raw_Bruijn_etal_2011_Tectonophysics.csv")
    g_bru = np.asarray(bru["gamma"], dtype=float)
    delta_bru = SLOPE * j_to_m_skemer(bru["J"]) + 1.0

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.4), constrained_layout=True)
    grid = np.linspace(0.0, 13.0, 400)
    delta_inf = SLOPE * M_INF + 1.0

    ax.scatter(g_bru, delta_bru, s=140, marker="x", color="darkred",
               linewidths=2.0, zorder=4,
               label=f"Bruijn+11 (n={len(g_bru)}, AUDIT-FLAGGED — see caveat)")
    ax.plot(grid, delta_of_gamma(grid, M_INF, GAMMA_E),
            color="black", lw=2.0,
            label=f"case 2 T-avg  (M_∞={M_INF}, γ_e={GAMMA_E}, "
                  f"δ_∞={delta_inf:.2f})")
    ax.axhline(delta_inf, color="black", lw=0.5, ls=":", alpha=0.6)

    if rows:
        g_md = np.array([r[0] for r in rows])
        d_md = np.array([r[1] for r in rows])
        ax.scatter(g_md, d_md, s=170, marker="*", color="black",
                   edgecolors="white", linewidths=0.7, zorder=6,
                   label=f"MDOODZ AniFstrainCalciteBruijn (n={len(rows)})")

    ax.text(0.98, 0.04,
            "⚠ Audit: Bruijn+11 J-values not sourced in published paper.\n"
            "   See raw CSV header + audit/Bruijn_etal_2011_audit (memo).",
            transform=ax.transAxes, fontsize=8, color="darkred",
            ha="right", va="bottom",
            bbox=dict(facecolor="white", edgecolor="darkred", alpha=0.85))

    ax.set_xlabel("shear strain γ")
    ax.set_ylabel("δ (= 6·M + 1)")
    ax.set_xlim(-0.5, 13.0)
    ax.set_ylim(0, 5.5)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(
        "MDOODZ ani_fstrain=2 + aniso_db=2 vs Bruijn+11 Carrara marble (727°C)\n"
        "(audit-control SET — Bruijn J-values not in paper, σ-fit driven by Pieri/Barnhoorn)",
        fontsize=10,
    )

    fig.savefig(args.out, dpi=140)
    size_kb = args.out.stat().st_size / 1024
    print(f"\nWrote {args.out} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
