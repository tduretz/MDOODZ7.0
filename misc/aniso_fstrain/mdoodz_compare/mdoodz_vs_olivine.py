#!/usr/bin/env python3
"""Overlay MDOODZ HDF5 output from the olivine CI tests on the combined
6-dataset olivine survey.  CSV-driven (raw data under ../data/olivine/).

Reads HDF5 files written by:
    AnisotropyBenchmark.AniFstrainHansenOlivine        (aniso_db = 1, fresh CPO)
    AnisotropyBenchmark.AniFstrainHansenOlivineDamped  (aniso_db = 7, damped)

Plots each per-step δ_mean as colored stars on its respective calibrated
curve, alongside all 6 olivine datasets converted to δ via δ = 24.5·M + 1.

Usage:
    python3 mdoodz_vs_olivine.py
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

# Hook the calibrate package on sys.path so we can import aniso_data.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "calibrate"))

from aniso_data import (  # noqa: E402
    load_raw,
    j_to_m_skemer,
    fs_ar_from_gamma,
    gamma_from_axial_eps,
    gamma_from_xz,
)

# Mirror MDLIB constants — must match MDLIB/FlowLaws.c.
M_INF_FRESH    = 0.536
GAMMA_E_FRESH  = 3.96
M_INF_DAMPED   = 0.16
GAMMA_E_DAMPED = 0.76
SLOPE          = 24.5

OUT = HERE.parent / "img"
OUT.mkdir(parents=True, exist_ok=True)
DEFAULT_OUT = OUT / "mdoodz_vs_olivine.png"


def m_of_gamma(g, mi, ge):
    return mi * (1.0 - np.exp(-g / ge))


def delta_of_gamma(g, mi, ge):
    return SLOPE * m_of_gamma(g, mi, ge) + 1.0


def delta_min(fs_ar, ani_fac_max):
    return np.minimum(fs_ar, ani_fac_max)


def read_mdoodz_outputs(build_dir: Path):
    outs = sorted(build_dir.glob("Output*.gzip.h5"))
    if not outs:
        return []
    rows = []
    for p in outs:
        m = p.name.replace("Output", "").replace(".gzip.h5", "")
        step = int(m)
        if step == 0:
            continue
        gamma = float(step - 1)  # γ̇ = 2, dt = 0.5 → γ at step k = (k-1)
        with h5py.File(p, "r") as f:
            ani = np.array(f["Centers"]["ani_fac"]).ravel()
            ani = ani[ani > 0.0]
            rows.append((gamma, float(ani.mean()),
                         float(ani.min()), float(ani.max())))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    base = HERE.parent.parent.parent / "cmake-build-test" / "TESTS"
    ap.add_argument("--fresh-build-dir",  type=Path,
                    default=base / "AniFstrainHansenOlivine")
    ap.add_argument("--damped-build-dir", type=Path,
                    default=base / "AniFstrainHansenOlivineDamped")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows_fresh  = read_mdoodz_outputs(args.fresh_build_dir)
    rows_damped = read_mdoodz_outputs(args.damped_build_dir)

    print("MDOODZ outputs:")
    print(f"  case 1 fresh : {len(rows_fresh):>2d} pts from {args.fresh_build_dir.name}")
    print(f"  case 7 damped: {len(rows_damped):>2d} pts from {args.damped_build_dir.name}")

    # MDOODZ vs analytical formula per case
    for label, rows, mi, ge in [
        ("case 1 fresh", rows_fresh, M_INF_FRESH, GAMMA_E_FRESH),
        ("case 7 damped", rows_damped, M_INF_DAMPED, GAMMA_E_DAMPED),
    ]:
        if not rows:
            continue
        g = np.array([r[0] for r in rows]); d = np.array([r[1] for r in rows])
        d_form = delta_of_gamma(g, mi, ge)
        rel_err = (d - d_form) / np.maximum(d_form, 1e-12)
        print(f"  {label} max rel error vs formula: {np.max(np.abs(rel_err)):.2e}")

    # Lab data via CSVs
    h14 = load_raw("olivine/raw_Hansen_etal_2014_EPSL.csv")
    h16 = load_raw("olivine/raw_Hansen_etal_2016_EPSL.csv")
    g_hansen = np.concatenate([h14["gamma"], h16["gamma"]])
    M_hansen = np.concatenate([h14["M"],     h16["M"]])
    d_hansen = SLOPE * M_hansen + 1.0

    b00 = load_raw("olivine/raw_Bystricky_etal_2000_Science.csv")
    g_byst = b00["gamma"]; d_byst = SLOPE * j_to_m_skemer(b00["J"]) + 1.0

    t16 = load_raw("olivine/raw_Tasaka_etal_2016_JGR.csv")
    g_tas = t16["gamma"]; d_tas = SLOPE * t16["M"] + 1.0

    bn = load_raw("olivine/raw_Boneh_Skemer_2014_EPSL.csv")
    g_bon = gamma_from_axial_eps(bn["epsilon"]); d_bon = SLOPE * bn["M"] + 1.0

    k19 = load_raw("olivine/raw_Kumamoto_etal_2019_JGR.csv")
    g_kum = k19["gamma"]; d_kum = SLOPE * k19["M"] + 1.0

    bj = load_raw("olivine/raw_Bernard_etal_2019_G3.csv")
    g_brn = gamma_from_xz(bj["xz"]); d_brn = SLOPE * bj["M"] + 1.0

    # ----- 2-panel plot -----
    fig, (ax_lo, ax_hi) = plt.subplots(1, 2, figsize=(15.5, 6.0),
                                        constrained_layout=True)
    grid = np.linspace(0.0, 22.0, 500)
    fs_ar_grid = fs_ar_from_gamma(grid)
    delta_inf_fresh  = SLOPE * M_INF_FRESH  + 1.0
    delta_inf_damped = SLOPE * M_INF_DAMPED + 1.0

    def plot_panel(ax, gamma_max):
        ax.scatter(g_hansen, d_hansen, s=40, marker="o",
                    facecolors="lightgray", edgecolors="black",
                    linewidths=0.7, alpha=0.85, zorder=3,
                    label=f"Hansen+12/+14/+16 (n={len(g_hansen)}, lab dry Fo90)")
        ax.scatter(g_byst, d_byst, s=110, marker="*",
                    facecolors="C0", edgecolors="black",
                    linewidths=0.8, zorder=3,
                    label=f"Bystricky 2000 (n={len(g_byst)}, lab dry Fo90)")
        ax.scatter(g_tas, d_tas, s=55, marker="s",
                    facecolors="none", edgecolors="C1",
                    linewidths=1.3, zorder=3,
                    label=f"Tasaka 2016 (n={len(g_tas)}, wet Fo50)")
        ax.scatter(g_bon, d_bon, s=45, marker="^",
                    facecolors="none", edgecolors="C2",
                    linewidths=1.2, zorder=3,
                    label=f"Boneh&Skemer14 (n={len(g_bon)}, pre-CPO)")
        ax.scatter(g_kum, d_kum, s=50, marker="D",
                    facecolors="none", edgecolors="C5",
                    linewidths=1.3, zorder=3,
                    label=f"Kumamoto+19 (n={len(g_kum)}, natural Josephine)")
        ax.scatter(g_brn, d_brn, s=40, marker="v",
                    facecolors="none", edgecolors="C6",
                    linewidths=1.0, zorder=3,
                    label=f"Bernard+19 (n={len(g_brn)}, natural xenoliths)")

        ax.plot(grid, delta_min(fs_ar_grid, 14.02), color="C7", lw=1.0,
                 linestyle=":", label="OLD: δ = min(FS_AR, 14.02)")
        ax.plot(grid, delta_of_gamma(grid, M_INF_FRESH, GAMMA_E_FRESH),
                 color="black", lw=2.0,
                 label=f"case 1 fresh  (M_inf={M_INF_FRESH}, γ_e={GAMMA_E_FRESH}, "
                       f"δ_∞={delta_inf_fresh:.2f})")
        ax.plot(grid, delta_of_gamma(grid, M_INF_DAMPED, GAMMA_E_DAMPED),
                 color="purple", lw=2.0,
                 label=f"case 7 damped (M_inf={M_INF_DAMPED}, γ_e={GAMMA_E_DAMPED}, "
                       f"δ_∞={delta_inf_damped:.2f})")

        if rows_fresh:
            g_md = np.array([r[0] for r in rows_fresh])
            d_md = np.array([r[1] for r in rows_fresh])
            ax.scatter(g_md, d_md, s=160, marker="*", color="black",
                        edgecolors="white", linewidths=0.7, zorder=6,
                        label=f"MDOODZ case 1 HDF5 (n={len(rows_fresh)})")
        if rows_damped:
            g_md = np.array([r[0] for r in rows_damped])
            d_md = np.array([r[1] for r in rows_damped])
            ax.scatter(g_md, d_md, s=160, marker="*", color="purple",
                        edgecolors="white", linewidths=0.7, zorder=6,
                        label=f"MDOODZ case 7 HDF5 (n={len(rows_damped)})")

        ax.axhline(delta_inf_fresh,  color="black",  lw=0.4, ls=":", alpha=0.5)
        ax.axhline(delta_inf_damped, color="purple", lw=0.4, ls=":", alpha=0.5)

        ax.set_xlabel("shear strain γ  (Boneh: γ_eq from ε; Bernard: γ_eq from X/Z)")
        ax.set_ylabel("δ (viscous-anisotropy magnitude)")
        ax.set_xlim(-0.5, gamma_max)
        ax.set_ylim(0, 18)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=7.5)

    plot_panel(ax_lo, 6.0)
    ax_lo.set_title("Low-strain detail γ ∈ [0, 6]  (MDOODZ HDF5 stars)")
    plot_panel(ax_hi, 22.0)
    ax_hi.set_title(f"Full saturation range γ ∈ [0, 22]  "
                     f"(δ_∞: case 1 = {delta_inf_fresh:.2f}, "
                     f"case 7 = {delta_inf_damped:.2f})")

    fig.suptitle(
        "MDOODZ ani_fstrain=2 vs olivine literature — "
        "case 1 (Hansen lab) + case 7 (Tasaka+Boneh+Kumamoto+Bernard)",
        fontsize=11,
    )
    fig.savefig(args.out, dpi=140)
    print(f"\nWrote {args.out} ({args.out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
