#!/usr/bin/env python3
"""Calibrate amphibolite/hornblende saturating-exp form.

Fits the saturating-exp form to per-sample data extracted from Ko-Jung 2015
Table 2 (CSV at data/amphibole/raw_Ko_Jung_2015_NatComms_table2.csv); avoids
prior literature-estimated values (M_∞ = 0.45, γ_e = 7 from Tatham+08 with
n = 3 plus Vp-anisotropy proxy).

Substrate: 14 deformed amphibolite samples, Griggs simple shear at P=1 GPa,
T ∈ [480, 700] °C, γ ∈ [0.9, 5.7]. Three fabric types (I/II/III) developed
depending on (T, σ). Response: AVp_max (%), the max P-wave seismic anisotropy
computed via Voigt-Reuss-Hill of single-crystal hornblende Cij × EBSD ODF.

Model: AVp(γ) = AVp_∞ · (1 − exp(−γ / γ_e))

Free parameters: AVp_∞, γ_e (2 params).

CRITICAL CAVEATS:
1. AVp(%) is a CPO-strength proxy for fixed mineralogy, NOT a direct
   Skemer M-index. Conversion AVp → M for MDLIB deployment requires
   separate calibration. For now this fit operates in AVp space.
2. Fabric types I/II/III are categorically distinct CPO geometries.
   A unified saturating-exp may smear that distinction; stratified
   fits (per-fabric-type) are run as secondary diagnostic.
3. Bounded range AVp ∈ [9.0, 14.6] % suggests near-saturation across
   most of the dataset — γ_e may be poorly constrained because no γ = 0
   anchor exists in Ko-Jung's deformed sample set. The "v2_anchored"
   variant below augments with a synthetic γ = 0 anchor to test this.

Optional snapshot path: set env var CASE16_SNAPSHOT_PATH to write a markdown
summary alongside the fit.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


CSV_PATH = Path(__file__).parent.parent / "data" / "amphibole" / "raw_Ko_Jung_2015_NatComms_table2.csv"
# Optional snapshot output path: set via env var CASE16_SNAPSHOT_PATH or leave None.
import os
_SNAP_ENV = os.environ.get("CASE16_SNAPSHOT_PATH")
SNAPSHOT_PATH = Path(_SNAP_ENV) if _SNAP_ENV else None
RUN_MODE = "compare"   # "v1": just the 2-param baseline fit;
                       # "compare": also report v2_anchored + v2_3param variants.


def load_csv(path: Path):
    """Load CSV skipping comment lines starting with '#'."""
    rows = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            rows.append(stripped)
    if not rows:
        raise RuntimeError(f"No data rows in {path}")
    header = rows[0].split(",")
    data = []
    for line in rows[1:]:
        fields = line.split(",")
        record = dict(zip(header, fields))
        data.append(record)
    return data


def to_float(records, key):
    return np.array([float(r[key]) for r in records])


def saturating_exp(gamma, avp_inf, gamma_e):
    """AVp(γ) = AVp_∞ · (1 − exp(−γ / γ_e))."""
    return avp_inf * (1.0 - np.exp(-gamma / gamma_e))


def saturating_exp_init(gamma, avp_inf, gamma_e, avp_init):
    """V2 form with free initial CPO: AVp(γ) = AVp_init + (AVp_∞ − AVp_init)·(1−exp(−γ/γ_e)).

    Tests whether adding a γ = 0 hydrostatic anchor (Ko-Jung Fig 4 starting
    material max m.u.d. ≈ 2, AVp ≈ 2 %) constrains γ_e — i.e. whether the
    "near-saturation" character of the substrate is responsible for the
    γ_e uncertainty.
    """
    return avp_init + (avp_inf - avp_init) * (1.0 - np.exp(-gamma / gamma_e))


def fit_subset(records, label, model="v1"):
    """Run LSQ on a subset of records; return (popt, perr, sse, r2, rmse).

    model="v1": saturating_exp 2-param baseline.
    model="v2_anchored": augment data with synthetic γ=0, AVp=2 % hydrostatic
        anchor (Ko-Jung Fig 4 starting material) then fit 2-param V1 form.
    model="v2_3param": fit 3-param saturating_exp_init with AVp_init free
        (no synthetic anchor); tests whether 3-param LSQ recovers the V2
        intuition without forcing a value.
    """
    gamma = to_float(records, "gamma")
    avp = to_float(records, "AVp_max_pct")
    if len(gamma) < 3:
        return None

    if model == "v2_anchored":
        # Augment with synthetic γ=0 anchor; AVp=2% from Ko-Jung Fig 4
        # starting-material pole-fig max m.u.d. ≈ 2 (page 4 Fig 4 contour scale).
        gamma_aug = np.concatenate([[0.0], gamma])
        avp_aug = np.concatenate([[2.0], avp])
        fit_gamma, fit_avp = gamma_aug, avp_aug
    else:
        fit_gamma, fit_avp = gamma, avp

    if model == "v2_3param":
        # 3 free params: AVp_∞, γ_e, AVp_init.
        p0 = [float(np.max(avp)), float(np.median(gamma)), 2.0]
        try:
            popt, pcov = curve_fit(
                saturating_exp_init, fit_gamma, fit_avp, p0=p0,
                bounds=([0.0, 0.01, 0.0], [50.0, 100.0, 15.0]),
                maxfev=10000,
            )
        except RuntimeError as exc:
            print(f"  {label}: 3-param fit failed — {exc}", file=sys.stderr)
            return None
        avp_fit = saturating_exp_init(gamma, *popt)
    else:
        p0 = [float(np.max(avp)), float(np.median(gamma))]
        try:
            popt, pcov = curve_fit(
                saturating_exp, fit_gamma, fit_avp, p0=p0,
                bounds=([0.0, 0.01], [50.0, 100.0]),
                maxfev=10000,
            )
        except RuntimeError as exc:
            print(f"  {label}: fit failed — {exc}", file=sys.stderr)
            return None
        avp_fit = saturating_exp(gamma, *popt)

    perr = np.sqrt(np.diag(pcov))
    sse = float(np.sum((avp - avp_fit) ** 2))
    ss_tot = float(np.sum((avp - np.mean(avp)) ** 2))
    r2 = 1.0 - sse / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(sse / len(avp)))
    return {
        "label": label,
        "n": len(gamma),
        "model": model,
        "popt": popt,
        "perr": perr,
        "sse": sse,
        "r2": r2,
        "rmse": rmse,
        "gamma": gamma,
        "avp_obs": avp,
        "avp_fit": avp_fit,
    }


def format_result(r):
    if r is None:
        return "(fit failed)\n"
    model = r.get("model", "v1")
    lines = [f"### {r['label']} (n={r['n']}, model={model})", ""]
    if model == "v2_3param":
        avp_inf, gamma_e, avp_init = r["popt"]
        avp_inf_err, gamma_e_err, avp_init_err = r["perr"]
        lines += [
            f"- AVp_∞     = **{avp_inf:.3f}** ± {avp_inf_err:.3f} % ({100.0*avp_inf_err/max(abs(avp_inf),1e-9):.0f}%)",
            f"- γ_e       = **{gamma_e:.3f}** ± {gamma_e_err:.3f} ({100.0*gamma_e_err/max(abs(gamma_e),1e-9):.0f}%)",
            f"- AVp_init  = **{avp_init:.3f}** ± {avp_init_err:.3f} % ({100.0*avp_init_err/max(abs(avp_init),1e-9):.0f}%)",
        ]
    else:
        avp_inf, gamma_e = r["popt"]
        avp_inf_err, gamma_e_err = r["perr"]
        lines += [
            f"- AVp_∞ = **{avp_inf:.3f}** ± {avp_inf_err:.3f} % ({100.0*avp_inf_err/max(abs(avp_inf),1e-9):.0f}%)",
            f"- γ_e   = **{gamma_e:.3f}** ± {gamma_e_err:.3f} ({100.0*gamma_e_err/max(abs(gamma_e),1e-9):.0f}%)",
        ]
    lines += [
        f"- SSE   = {r['sse']:.3f}",
        f"- R²    = **{r['r2']:.4f}**",
        f"- RMSE  = {r['rmse']:.3f} %",
        "",
    ]
    return "\n".join(lines)


def main():
    records = load_csv(CSV_PATH)
    print(f"Loaded {len(records)} amphibole records from {CSV_PATH}")

    # Stratifications to test:
    all_records = records
    type_i = [r for r in records if r["fabric"] == "Type-I"]
    type_ii = [r for r in records if r["fabric"] == "Type-II"]
    type_iii = [r for r in records if r["fabric"] == "Type-III"]
    low_T = [r for r in records if float(r["T_C"]) < 600.0]
    high_T = [r for r in records if float(r["T_C"]) >= 600.0]
    dry = [r for r in records if r["water"] == "dry"]
    wet = [r for r in records if r["water"] == "wet"]

    results = []
    if RUN_MODE == "compare":
        # Compare three model variants on all-samples + Type-I + Type-II,
        # testing whether a γ = 0 anchor (forced or fitted) constrains γ_e
        # versus the 2-param baseline. R² is computed on the original n = 14
        # observed AVp values (NOT augmented data) so cross-model comparison
        # is on the same observations.
        for label, subset in [
            ("All samples", all_records),
            ("Type-I (low-T / [001]-fabric)", type_i),
            ("Type-II (high-σ / (010)-fabric)", type_ii),
        ]:
            for model in ("v1", "v2_anchored", "v2_3param"):
                r = fit_subset(subset, f"{label} [{model}]", model=model)
                results.append(r)
    else:
        for label, subset in [
            ("All samples", all_records),
            ("Type-I (low-T / [001]-fabric)", type_i),
            ("Type-II (high-σ / (010)-fabric)", type_ii),
            ("Type-III (high-T / girdle)", type_iii),
            ("T < 600 °C", low_T),
            ("T ≥ 600 °C", high_T),
            ("Dry samples", dry),
            ("Wet samples", wet),
        ]:
            r = fit_subset(subset, label)
            results.append(r)

    # Build per-datapoint table for the all-samples fit.
    all_r = results[0]
    per_dp_lines = ["| sample | T (°C) | γ | fabric | AVp_obs (%) | AVp_fit (%) | residual |",
                    "|---|---|---|---|---|---|---|"]
    for i, rec in enumerate(all_records):
        obs = all_r["avp_obs"][i]
        fit = all_r["avp_fit"][i]
        per_dp_lines.append(
            f"| {rec['sample']} | {rec['T_C']} | {rec['gamma']} | {rec['fabric']} | "
            f"{obs:.2f} | {fit:.2f} | {obs - fit:+.2f} |"
        )

    if RUN_MODE == "compare":
        md_parts = [
            "# case-16 amphibole — model-variant comparison (γ = 0 anchor counter-check)",
            "",
            f"**Input CSV**: `{CSV_PATH}`",
            f"**Script**: `calibrate_amphibole.py` (RUN_MODE = compare)",
            "",
            "## Model variants compared",
            "",
            "- **v1**: `AVp(γ) = AVp_∞ · (1 − exp(−γ / γ_e))` — 2-param baseline.",
            "- **v2_anchored**: same form, but data augmented with synthetic γ = 0, AVp = 2 %",
            "  anchor (Ko-Jung Fig 4 starting-material pole-fig max m.u.d. ≈ 2 → AVp ≈ 2 %).",
            "- **v2_3param**: 3-param `AVp(γ) = AVp_init + (AVp_∞ − AVp_init) · (1 − exp(−γ / γ_e))`",
            "  with AVp_init free (no augmentation).",
            "",
            "**Interpretation rule**: if v2_anchored or v2_3param γ_e converges with bounded",
            "uncertainty AND R² improves materially over v1, the γ-anchor hypothesis is",
            "supported. If γ_e still collapses or R² stays ≈ 0, the substrate is too",
            "saturated / regime-categorical for a unified saturating-exp.",
            "R² computed on n = 14 observed values (NOT augmented data) for cross-model parity.",
            "",
            "## Fit results across stratifications",
            "",
        ]
    else:
        md_parts = [
            "# case-16 amphibolite calibration — LSQ output",
            "",
            f"**Input CSV**: `{CSV_PATH}`",
            f"**Script**: `calibrate_amphibole.py`",
            "",
            "## Model",
            "",
            "`AVp(γ) = AVp_∞ · (1 − exp(−γ / γ_e))`",
            "",
            "Response variable: AVp_max (%) from Ko-Jung Table 2.",
            "AVp is a CPO-strength proxy via VRH of single-crystal hornblende Cij × ODF.",
            "AVp → Skemer M conversion left to future calibration.",
            "",
            "## Fit results across stratifications",
            "",
        ]
    for r in results:
        md_parts.append(format_result(r))

    md_parts += [
        "## Per-datapoint table (All samples fit)",
        "",
        *per_dp_lines,
        "",
        "## Interpretation hooks",
        "",
        "- Compare R² across stratifications. If per-fabric or per-T-segment fits",
        "  outperform the unified fit by Δ R² > 0.1, case-16 likely needs",
        "  fabric- or T-stratified MDLIB deployment.",
        "- Check γ_e uncertainty. With γ_min = 0.9 (no hydrostatic anchor),",
        "  expect γ_e poorly constrained if AVp_∞ is in the saturating regime",
        "  across most data.",
        "- If unified fit has acceptable R² (> 0.5), report as the primary",
        "  case-16 calibration with explicit γ ∈ [0.9, 5.7] scope and an",
        "  AVp → M conversion caveat.",
    ]

    if SNAPSHOT_PATH is None:
        print("(SNAPSHOT_PATH unset; skipping markdown write — set CASE16_SNAPSHOT_PATH env var to enable.)")
    else:
        SNAPSHOT_PATH.parent.mkdir(parents=True, exist_ok=True)
        SNAPSHOT_PATH.write_text("\n".join(md_parts))
        print(f"Snapshot written to {SNAPSHOT_PATH}")

    # Also dump to stdout for quick check.
    for r in results:
        if r is None:
            continue
        model = r.get("model", "v1")
        avp_inf = r["popt"][0]
        gamma_e = r["popt"][1]
        extra = f", AVp_init={r['popt'][2]:.2f}" if model == "v2_3param" else ""
        print(f"  {r['label']} (n={r['n']}): AVp_∞={avp_inf:.2f}, "
              f"γ_e={gamma_e:.2f}{extra}, R²={r['r2']:.3f}, RMSE={r['rmse']:.3f}%")


if __name__ == "__main__":
    main()
