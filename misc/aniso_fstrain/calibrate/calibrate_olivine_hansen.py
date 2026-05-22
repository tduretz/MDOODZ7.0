#!/usr/bin/env python3
"""Calibrate the Hansen-group olivine M(γ) → δ mapping (case 1).

Fits M(γ) = M_inf · (1 − exp(−γ/γ_e)) against the combined Hansen+14 +
Hansen+16 dataset (38 pts).  Constants must reproduce MDLIB/FlowLaws.c
anisoDelta_HansenOlivine (M_inf = 0.536, γ_e = 3.96, δ_∞ = 14.13).

Reads from the CSV-driven data tree under misc/aniso_fstrain/data/olivine/:
    raw_Hansen_etal_2014_EPSL.csv
    raw_Hansen_etal_2016_EPSL.csv
    raw_Hansen_etal_2012_Nature.csv  (historic comparison only)

The Hansen+16 Part 1 stress-aware reference δ_∞ = (1.39/0.73)^4.1 ≈ 14.0
serves as an independent check on the Hansen+12 Fig 3b linear δ-vs-M slope
(24.5).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aniso_data import load_raw, fs_ar_from_gamma

OUT_DIR   = Path(__file__).resolve().parent.parent / "img"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_PATH = OUT_DIR / "olivine_hansen_calibration.png"

SLOPE = 24.5   # Hansen+12 Fig 3b linear δ-vs-M coefficient


# ---------------------------------------------------------------------------
# Functional forms (mirror MDLIB/FlowLaws.c case 1 exactly)
# ---------------------------------------------------------------------------
def m_of_gamma(gamma: np.ndarray, m_inf: float, gamma_e: float) -> np.ndarray:
    return m_inf * (1.0 - np.exp(-gamma / gamma_e))


def delta_hansen(gamma: np.ndarray, m_inf: float, gamma_e: float) -> np.ndarray:
    return SLOPE * m_of_gamma(gamma, m_inf, gamma_e) + 1.0


def delta_min(fs_ar: np.ndarray, ani_fac_max: float) -> np.ndarray:
    return np.minimum(fs_ar, ani_fac_max)


def delta_erfc(fs_ar: np.ndarray, ani_fac_max: float) -> np.ndarray:
    """erfc-blend prototype (gated by -DANISO_ERFC_BLEND in legacy builds)."""
    delta_ratio = ani_fac_max / 4.0
    fs_crit = ani_fac_max / 2.0
    erfc = np.vectorize(math.erfc)
    x = 0.5 * erfc((fs_ar - fs_crit) / delta_ratio)
    return x * fs_ar + (1.0 - x) * ani_fac_max


# ---------------------------------------------------------------------------
# Fit: 2D coarse-then-fine grid search on (M_inf, γ_e)
# ---------------------------------------------------------------------------
def grid_search(gammas, m_obs, m_inf_range, gamma_e_range, n):
    mi_grid = np.linspace(*m_inf_range, n)
    ge_grid = np.linspace(*gamma_e_range, n)
    best = (None, None, float("inf"))
    for mi in mi_grid:
        for ge in ge_grid:
            sse = float(np.sum((m_of_gamma(gammas, mi, ge) - m_obs) ** 2))
            if sse < best[2]:
                best = (float(mi), float(ge), sse)
    return best


def fit(gammas, m_obs):
    coarse = grid_search(gammas, m_obs, (0.30, 0.95), (0.5, 15.0), n=50)
    cmi, cge, _ = coarse
    fine = grid_search(
        gammas, m_obs,
        (max(0.20, cmi - 0.20), min(0.99, cmi + 0.20)),
        (max(0.10, cge - 4.0), cge + 4.0),
        n=80,
    )
    fmi, fge, _ = fine
    return grid_search(
        gammas, m_obs,
        (max(0.20, fmi - 0.05), min(0.99, fmi + 0.05)),
        (max(0.10, fge - 1.0), fge + 1.0),
        n=200,
    )


def rms(pred, obs):
    return float(np.sqrt(np.mean((pred - obs) ** 2)))


# ---------------------------------------------------------------------------
def main() -> None:
    h14 = load_raw("olivine/raw_Hansen_etal_2014_EPSL.csv")
    h16 = load_raw("olivine/raw_Hansen_etal_2016_EPSL.csv")
    h12 = load_raw("olivine/raw_Hansen_etal_2012_Nature.csv")

    g14, m14 = h14["gamma"], h14["M"]
    g16, m16 = h16["gamma"], h16["M"]
    g12, m12 = h12["gamma"], h12["M"]

    g_all = np.concatenate([g14, g16])
    m_all = np.concatenate([m14, m16])
    print(f"Combined Hansen-group dataset: {len(g14)}+{len(g16)} = {len(g_all)} pts")

    mi_all, ge_all, sse_all = fit(g_all, m_all)
    mi_h14, ge_h14, sse_h14 = fit(g14, m14)
    mi_h16, ge_h16, sse_h16 = fit(g16, m16)
    mi_h12, ge_h12, sse_h12 = fit(g12, m12)

    # Hansen+16 Part 1 Eq. 3-4 stress-aware δ_∞
    F_w, F_s, n_exp = 0.73, 1.39, 4.1
    delta_inf_eq34 = (F_s / F_w) ** n_exp

    print("\n=== Calibration summary (M(γ) = M_inf · (1 − exp(−γ/γ_e))) ===")
    print(f"  {'dataset':<28s} {'N':>3s}  {'M_inf':>7s}  {'γ_e':>7s}  "
          f"{'δ_∞':>9s}  {'RMS(M)':>7s}")
    for name, n, mi, ge, sse in [
        ("Hansen+12 Nature",  len(g12), mi_h12, ge_h12, sse_h12),
        ("Hansen+14 EPSL",    len(g14), mi_h14, ge_h14, sse_h14),
        ("Hansen+16 P1 EPSL", len(g16), mi_h16, ge_h16, sse_h16),
        ("COMBINED 38 pts",   len(g_all), mi_all, ge_all, sse_all),
    ]:
        print(f"  {name:<28s} {n:>3d}  {mi:>7.4f}  {ge:>7.4f}  "
              f"{SLOPE * mi + 1.0:>9.3f}  {math.sqrt(sse / n):>7.4f}")
    print(f"\n  Independent: Hansen+16 Eq. 3-4 stress-aware δ_∞ = "
          f"(1.39/0.73)^4.1 = {delta_inf_eq34:.3f}")

    # δ-form RMS comparison
    delta_lab = SLOPE * m_all + 1.0
    fs_ar_all = fs_ar_from_gamma(g_all)
    cap = delta_inf_eq34
    print("\n=== RMS error of candidate δ-forms vs the COMBINED 38-pt lab data ===")
    print(f"  (cap = {cap:.2f})")
    print(f"  δ = min(FS_AR, cap)               RMS = {rms(delta_min(fs_ar_all, cap), delta_lab):.3f}")
    print(f"  δ = erfc-blend(FS_AR, cap)        RMS = {rms(delta_erfc(fs_ar_all, cap), delta_lab):.3f}")
    print(f"  δ = 24.5·M_H12(γ) + 1   (3-pt)    RMS = {rms(delta_hansen(g_all, mi_h12, ge_h12), delta_lab):.3f}")
    print(f"  δ = 24.5·M_H16(γ) + 1   (12-pt)   RMS = {rms(delta_hansen(g_all, mi_h16, ge_h16), delta_lab):.3f}")
    print(f"  δ = 24.5·M_ALL(γ) + 1   (38-pt)   RMS = {rms(delta_hansen(g_all, mi_all, ge_all), delta_lab):.3f}")

    # ----- 2-panel calibration plot -----
    grid = np.linspace(0.0, 20.0, 500)
    fs_ar_grid = fs_ar_from_gamma(grid)

    fig, (ax_m, ax_d) = plt.subplots(1, 2, figsize=(14.0, 5.4),
                                      constrained_layout=True)

    ax_m.scatter(g14, m14, s=42, facecolors="none", edgecolors="C1",
                  linewidths=1.0, label=f"Hansen+14 EPSL (n={len(g14)})")
    ax_m.scatter(g16, m16, s=55, marker="s", facecolors="none",
                  edgecolors="C0", linewidths=1.2,
                  label=f"Hansen+16 EPSL P1 (n={len(g16)})")
    ax_m.scatter(g12, m12, s=130, marker="x", color="C3", linewidths=2.5,
                  label=f"Hansen+12 Nature (n={len(g12)})")
    ax_m.plot(grid, m_of_gamma(grid, mi_all, ge_all), color="black", lw=2.4,
               label=f"COMBINED fit  M = {mi_all:.3f}·(1−exp(−γ/{ge_all:.2f}))")
    ax_m.plot(grid, m_of_gamma(grid, mi_h16, ge_h16), color="C0", lw=1.4,
               linestyle="--", alpha=0.8,
               label=f"Hansen+16-only (M_inf={mi_h16:.3f}, γ_e={ge_h16:.2f})")
    ax_m.plot(grid, m_of_gamma(grid, mi_h12, ge_h12), color="C3", lw=1.2,
               linestyle=":", alpha=0.8,
               label=f"Hansen+12-only (M_inf={mi_h12:.3f}, γ_e={ge_h12:.2f})")
    ax_m.set_xlabel("shear strain γ")
    ax_m.set_ylabel("M-index")
    ax_m.set_title(f"Olivine M(γ) — Hansen-group dataset (n={len(g_all)})")
    ax_m.set_xlim(-0.5, 20)
    ax_m.set_ylim(0, 0.78)
    ax_m.legend(loc="lower right", fontsize=8)
    ax_m.grid(alpha=0.3)

    ax_d.scatter(g_all, delta_lab, s=40, facecolors="none", edgecolors="C7",
                  zorder=4,
                  label=f"Combined lab data (n={len(g_all)}) → δ = 24.5·M + 1")
    ax_d.plot(grid, delta_min(fs_ar_grid, cap), color="C7", lw=1.4,
               linestyle="--", label=f"δ = min(FS_AR, {cap:.1f})  (MDOODZ default)")
    ax_d.plot(grid, delta_erfc(fs_ar_grid, cap), color="C8", lw=1.4,
               linestyle="-.", label="δ = erfc-blend (prototype)")
    ax_d.plot(grid, delta_hansen(grid, mi_h12, ge_h12), color="C3", lw=1.4,
               linestyle=":",
               label="δ = 24.5·M_H12(γ_eff) + 1  (3-pt fit)")
    ax_d.plot(grid, delta_hansen(grid, mi_h16, ge_h16), color="C0", lw=1.4,
               linestyle="--",
               label="δ = 24.5·M_H16(γ_eff) + 1  (12-pt fit)")
    ax_d.plot(grid, delta_hansen(grid, mi_all, ge_all), color="black", lw=2.4,
               label="δ = 24.5·M_ALL(γ_eff) + 1  (38-pt fit, NEW)")
    ax_d.axhline(delta_inf_eq34, color="black", lw=0.6, linestyle=":")
    ax_d.text(19.5, delta_inf_eq34 + 0.3,
               f"Hansen+16 Eq. 3-4 δ_∞ ≈ {delta_inf_eq34:.1f}",
               fontsize=8, ha="right", va="bottom", color="black")
    ax_d.set_xlabel("shear strain γ")
    ax_d.set_ylabel("δ (viscous-anisotropy magnitude)")
    ax_d.set_title("Candidate δ-forms vs Hansen-group reference")
    ax_d.set_xlim(-0.5, 20)
    ax_d.set_ylim(0, 20)
    ax_d.legend(loc="lower right", fontsize=8)
    ax_d.grid(alpha=0.3)

    fig.suptitle(
        "Hansen-group olivine viscous-anisotropy calibration — case 1.  "
        f"Combined 38-pt fit:  M_inf = {mi_all:.3f}, γ_e = {ge_all:.2f}, "
        f"δ_∞ = {SLOPE * mi_all + 1.0:.2f}  (Hansen+16 Eq. 3-4: {delta_inf_eq34:.2f})",
        fontsize=11,
    )
    fig.savefig(PLOT_PATH, dpi=140)
    print(f"\nWrote {PLOT_PATH} ({PLOT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
