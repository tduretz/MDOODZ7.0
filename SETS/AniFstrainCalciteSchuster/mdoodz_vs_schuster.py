#!/usr/bin/env python3
"""Overlay MDOODZ AniFstrainCalciteSchuster HDF5 output on Schuster+18
(JSG; doi:10.1016/j.jsg.2018.09.003) Carrara marble HPT γ-coverage at 450°C.

Reads HDF5 files from ../../cmake-exec/AniFstrainCalciteSchuster/output/ and
writes img/mdoodz_vs_schuster.png. Mirrors mdoodz_vs_pieri.py.

Annotation-only data — Schuster+18 reports γ ∈ [3, 63] in Table 2 but does
NOT tabulate J/M for any sample. The plot therefore shows:
  - case-4 saturating-exp curve (M_∞ = 0.22, γ_e = 1.37, δ_∞ = 2.32)
  - 16 MDOODZ HDF5 stars from the simple-shear run (γ ≤ 16)
  - Schuster's γ-coverage as a horizontal rug at the bottom of the panel
    (so the high-strain saturation regime Schuster sampled is visible even
    though no J/M is plottable)
See audit/Schuster_etal_2018_annotation.md for the full annotation memo +
regime caveat (Schuster's HPT 1-4 GPa is off-spec for MDOODZ lithospheric
calibration; same γ does not produce same fabric as Pieri's Paterson rig).
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

from aniso_data import load_raw  # noqa: E402

# Mirror MDLIB/FlowLaws.c::anisoDelta_Calcite_LowT — must match aniso_db = 4.
M_INF   = 0.22
GAMMA_E = 1.37
SLOPE   = 6.0

OUT_DIR = HERE / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT_DIR / "mdoodz_vs_schuster.png"


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
    default_build = (ROOT / "cmake-exec" / "AniFstrainCalciteSchuster"
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
        print(f"  max relative error vs case-4 formula: "
              f"{np.max(np.abs(rel_err)):.2e}")

    sch = load_raw("calcite/raw_Schuster_etal_2018_JSG.csv")
    g_sch = np.asarray(sch["gamma"], dtype=float)
    P_sch = np.asarray(sch["P_GPa"], dtype=float)

    fig, ax = plt.subplots(1, 1, figsize=(9.0, 5.4), constrained_layout=True)
    grid = np.linspace(0.0, 65.0, 600)
    delta_inf = SLOPE * M_INF + 1.0

    ax.plot(grid, delta_of_gamma(grid, M_INF, GAMMA_E),
            color="darkorange", lw=2.0,
            label=f"case 4 LowT  (M_∞={M_INF}, γ_e={GAMMA_E}, "
                  f"δ_∞={delta_inf:.2f})")
    ax.axhline(delta_inf, color="darkorange", lw=0.5, ls=":", alpha=0.6)

    if rows:
        g_md = np.array([r[0] for r in rows])
        d_md = np.array([r[1] for r in rows])
        ax.scatter(g_md, d_md, s=170, marker="*", color="darkorange",
                   edgecolors="black", linewidths=0.7, zorder=6,
                   label=f"MDOODZ AniFstrainCalciteSchuster (n={len(rows)})")

    # γ-rug for Schuster: one tick per (sample, scan), color by P_GPa
    rug_y = 0.08
    P_unique = sorted(set(P_sch.tolist()))
    P_color = {1.2: "C0", 2.0: "C1", 3.0: "C2", 4.0: "C3"}
    for P in P_unique:
        mask = P_sch == P
        ax.scatter(g_sch[mask], np.full(mask.sum(), rug_y),
                   marker="|", s=200,
                   color=P_color.get(P, "gray"), linewidths=1.5, zorder=4,
                   label=f"Schuster+18 P={P:.1f} GPa (n={mask.sum()}, J not tabulated)")

    ax.text(0.98, 0.04,
            "γ-coverage rug: Schuster+18 reports γ ∈ [3, 63] in Table 2\n"
            "   but no J/M values. HPT 1-4 GPa is off-spec for lithospheric\n"
            "   calibration (see annotation memo §Regime caveat).",
            transform=ax.transAxes, fontsize=8, color="black",
            ha="right", va="bottom",
            bbox=dict(facecolor="white", edgecolor="gray", alpha=0.85))

    ax.set_xlabel("shear strain γ")
    ax.set_ylabel("δ (= 6·M + 1)")
    ax.set_xlim(-1.0, 65.0)
    ax.set_ylim(0, 3.5)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(
        "MDOODZ ani_fstrain=2 + aniso_db=4 vs Schuster+18 Carrara marble HPT (450°C)\n"
        "(annotation-only SET — Schuster reports γ-coverage; J/M not tabulated)",
        fontsize=10,
    )

    fig.savefig(args.out, dpi=140)
    size_kb = args.out.stat().st_size / 1024
    print(f"\nWrote {args.out} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
