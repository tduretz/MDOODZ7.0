#!/usr/bin/env python3
"""Generate compact derived CSVs from the raw_*.csv files.

For each mineral, the derived CSV reduces every dataset to the (γ_simple, M)
pair that the calibration scripts actually consume — applying Skemer J→M
conversions, axial-ε / X-Z / γ_oct → γ_simple inversions, and any documented
filtering (e.g., Barnhoorn+04 J ≤ 15 cutoff).

Re-run this script whenever a raw_*.csv is edited.  The committed calibration
constants in MDLIB/FlowLaws.c should be reproduced byte-identically by the
calibrate_*.py scripts after `make_derived.py`.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from aniso_data import (
    DATA_ROOT,
    load_raw,
    j_to_m_skemer,
    j_to_m_natural,
    gamma_from_axial_eps,
    gamma_from_xz,
    gamma_oct_to_gamma_simple,
)

OUT = DATA_ROOT


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_csv(path: Path, header: list[str], rows: list[list], comments: list[str]) -> None:
    """Write a CSV with `## ` comment header lines, then header, then rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        for c in comments:
            fh.write(f"## {c}\n")
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path.relative_to(OUT.parent)}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Olivine
# ---------------------------------------------------------------------------
def derive_olivine() -> None:
    print("[olivine]")

    # ----- HANSEN GROUP (Hansen+14 + Hansen+16, 38 pts) ---------------------
    h14 = load_raw("olivine/raw_Hansen_etal_2014_EPSL.csv")
    h16 = load_raw("olivine/raw_Hansen_etal_2016_EPSL.csv")
    rows = []
    for g, m, src in zip(h14["gamma"], h14["M"], h14["sample"]):
        rows.append([f"{g:.4g}", f"{m:.4g}", f"H14:{src}"])
    for g, m, src in zip(h16["gamma"], h16["M"], h16["sample"]):
        rows.append([f"{g:.4g}", f"{m:.4g}", f"H16:{src}"])
    write_csv(
        OUT / "olivine" / "derived_HansenGroup_38pts.csv",
        ["gamma_simple", "M", "source"],
        rows,
        [
            "Combined Hansen+14 + Hansen+16 olivine torsion dataset (38 pts).",
            "Used for case 1 (anisoDelta_HansenOlivine) calibration.",
            "gamma_simple = engineering shear strain γ as reported (lab torsion).",
            "M = Skemer M-index, no further conversion.",
            "Generated from raw_Hansen_etal_2014_EPSL.csv + raw_Hansen_etal_2016_EPSL.csv.",
        ],
    )

    # ----- BYSTRICKY (J → M Skemer) ----------------------------------------
    b00 = load_raw("olivine/raw_Bystricky_etal_2000_Science.csv")
    rows = []
    for g, J in zip(b00["gamma"], b00["J"]):
        rows.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}"])
    write_csv(
        OUT / "olivine" / "derived_Bystricky00_Mskemer.csv",
        ["gamma_simple", "J", "M_skemer"],
        rows,
        [
            "Bystricky+00 with M obtained via Skemer+05 linear J→M = max((J-1)/15, 0).",
            "Used in olivine survey (overlay panel) — NOT a calibration target.",
            "Generated from raw_Bystricky_etal_2000_Science.csv.",
        ],
    )

    # ----- TASAKA, KUMAMOTO -------------------------------------------------
    t16 = load_raw("olivine/raw_Tasaka_etal_2016_JGR.csv")
    rows = [[f"{g:.4g}", f"{m:.4g}", str(s)] for g, m, s in
            zip(t16["gamma"], t16["M"], t16["sample"])]
    write_csv(
        OUT / "olivine" / "derived_Tasaka16_gamma_M.csv",
        ["gamma_simple", "M", "sample"],
        rows,
        ["Tasaka+16 reduced to (γ, M); Fo50 wet — case 7 calibration input.",
         "Generated from raw_Tasaka_etal_2016_JGR.csv."],
    )

    k19 = load_raw("olivine/raw_Kumamoto_etal_2019_JGR.csv")
    rows = [[f"{g:.4g}", f"{m:.4g}", str(s), str(z)] for g, m, s, z in
            zip(k19["gamma"], k19["M"], k19["sample"], k19["shear_zone"])]
    write_csv(
        OUT / "olivine" / "derived_Kumamoto19_gamma_M.csv",
        ["gamma_simple", "M", "sample", "shear_zone"],
        rows,
        ["Kumamoto+19 reduced to (γ, M); natural Josephine peridotite — case 7 input.",
         "Generated from raw_Kumamoto_etal_2019_JGR.csv."],
    )

    # ----- BONEH (axial ε → γ_eq) ------------------------------------------
    b14 = load_raw("olivine/raw_Boneh_Skemer_2014_EPSL.csv")
    rows = []
    geom_label = {0: "undeformed", 1: "perpendicular", 2: "oblique", 3: "parallel"}
    for s, e, m, g in zip(b14["sample"], b14["epsilon"], b14["M"], b14["geometry"]):
        gamma_eq = float(gamma_from_axial_eps(e))
        rows.append([f"{gamma_eq:.4g}", f"{m:.4g}", str(s), geom_label.get(int(g), str(int(g)))])
    write_csv(
        OUT / "olivine" / "derived_Boneh14_gammaEq_M.csv",
        ["gamma_simple_eq", "M", "sample", "geometry"],
        rows,
        ["Boneh & Skemer 2014 with γ_simple_eq = 2·sinh(ε/2) from axial natural strain.",
         "Used in case 7 (damped) calibration as one of 4 dataset inputs.",
         "Generated from raw_Boneh_Skemer_2014_EPSL.csv."],
    )

    # ----- BERNARD (X/Z → γ_eq) --------------------------------------------
    bn19 = load_raw("olivine/raw_Bernard_etal_2019_G3.csv")
    rows = []
    for m, J, xz in zip(bn19["M"], bn19["J"], bn19["xz"]):
        gamma_eq = float(gamma_from_xz(xz))
        rows.append([f"{gamma_eq:.4g}", f"{m:.4g}", f"{J:.4g}", f"{xz:.4g}"])
    write_csv(
        OUT / "olivine" / "derived_Bernard19_gammaEq_M.csv",
        ["gamma_simple_eq", "M", "J", "xz"],
        rows,
        ["Bernard+19 with γ_simple_eq = √(X/Z) − 1/√(X/Z).",
         "Used in case 7 (damped) calibration as one of 4 dataset inputs.",
         "Generated from raw_Bernard_etal_2019_G3.csv."],
    )

    # ----- DAMPED COMBINED (case 7) — 12+16+30+35 = 93 pts -----------------
    # Bernard+19 M re-derived from J via j_to_m_natural() = (J−1)/24 per the Bernard+19 natural-sample regression
    # (research/aniso_fstrain) — natural-sample regression on Bernard's own
    # n=65 M/J pairs gives R²=0.825 for the (J−1)/24 form vs R²=−0.69 for
    # Skemer's lab (J−1)/15 on the same data. Tasaka16, Boneh14, Kumamoto19
    # keep their as-reported M because (a) Tasaka16 / Boneh14 are lab-deformed
    # synthetics for which Skemer's slope is presumed valid, and (b)
    # Kumamoto19 reports M only (no J in the raw CSV), so we cannot re-derive.
    rows = []
    for g, m, s in zip(t16["gamma"], t16["M"], t16["sample"]):
        rows.append([f"{g:.4g}", f"{m:.4g}", f"Tasaka16:{s}"])
    for s, e, m, g in zip(b14["sample"], b14["epsilon"], b14["M"], b14["geometry"]):
        gamma_eq = float(gamma_from_axial_eps(e))
        rows.append([f"{gamma_eq:.4g}", f"{m:.4g}", f"Boneh14:{s}"])
    for s, g, m in zip(k19["sample"], k19["gamma"], k19["M"]):
        rows.append([f"{g:.4g}", f"{m:.4g}", f"Kumamoto19:{s}"])
    for m_rep, J, xz in zip(bn19["M"], bn19["J"], bn19["xz"]):
        gamma_eq = float(gamma_from_xz(xz))
        m_corr = float(j_to_m_natural(J))
        rows.append([f"{gamma_eq:.4g}", f"{m_corr:.4g}", "Bernard19"])
    write_csv(
        OUT / "olivine" / "derived_DampedOlivine_93pts.csv",
        ["gamma_simple", "M", "source"],
        rows,
        [f"Combined damped olivine dataset for case 7 calibration ({len(rows)} pts).",
         "Tasaka16 wet Fo50 + Boneh14 pre-CPO + Kumamoto19 Josephine + Bernard19 natural xenoliths.",
         "γ_simple is either as-reported (Tasaka, Kumamoto) or converted from ε_axial (Boneh) / X-Z (Bernard).",
         "M for Tasaka16/Boneh14/Kumamoto19 is as-reported (lab or no-J).",
         "M for Bernard19 is re-derived from J via j_to_m_natural() = max((J−1)/24, 0),",
         "per the Bernard+19 natural-sample regression, R²=0.825) — supersedes Bernard's",
         "as-reported M which inherits Skemer's lab (J−1)/15 bias on natural samples."],
    )


# ---------------------------------------------------------------------------
# Calcite
# ---------------------------------------------------------------------------
def derive_calcite() -> None:
    print("[calcite]")

    p01 = load_raw("calcite/raw_Pieri_etal_2001_Tectonophysics.csv")
    b04 = load_raw("calcite/raw_Barnhoorn_etal_2004_JSG.csv")

    # Bruijn+11: removed from the calibration after the 2026-05-08 audit
    # (audit/Bruijn_etal_2011_audit.md) found the (γ, J) values cannot be
    # sourced from the Bruijn 2011 paper — its Fig. 8 has no numerical J
    # labels and the caption says the pole figures are reproduced from
    # Pieri+01 / Barnhoorn+04.  CSV header in
    # data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv is flagged
    # but kept for historical traceability; the data is no longer used.

    J_FILTER = 15.0   # Skemer+05 linear J→M valid range
    rows_full = []
    for s, g, J, T in zip(p01["sample"], p01["gamma"], p01["J"], p01["T_C"]):
        rows_full.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}",
                          f"{T:.0f}", f"Pieri01:{s}", "kept"])
    for g, J, T, n in zip(b04["gamma"], b04["J"], b04["T_C"], b04["note"]):
        kept = "kept" if J <= J_FILTER else f"filtered (J>{J_FILTER:.0f})"
        rows_full.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}",
                          f"{T:.0f}", f"Barnhoorn04 ({n})" if n else "Barnhoorn04",
                          kept])
    n_kept = sum(1 for r in rows_full if r[-1] == "kept")
    write_csv(
        OUT / "calcite" / "derived_combined_15pts_filtered.csv",
        ["gamma_simple", "J", "M_skemer", "T_C", "source", "fit_status"],
        rows_full,
        ["Combined Pieri+01 / Barnhoorn+04 calcite torsion dataset.",
         "Bruijn+11 removed after the 2026-05-08 audit (audit/Bruijn_etal_2011_audit.md).",
         f"Skemer+05 J→M conversion applied. J ≤ {J_FILTER:g} kept for case-2 LSQ fit ({n_kept} pts).",
         "Used for case 2 (anisoDelta_Calcite combined T-averaged fit).",
         "Generated from raw_Pieri_etal_2001 + raw_Barnhoorn_etal_2004."],
    )

    # ----- LOW-T (case 4) — Barnhoorn 500 + 600 °C, J ≤ 15 only ------------
    rows_lo = []
    for g, J, T, n in zip(b04["gamma"], b04["J"], b04["T_C"], b04["note"]):
        if T <= 650.0 and J <= J_FILTER:
            rows_lo.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}",
                            f"{T:.0f}", f"Barnhoorn04 ({n})" if n else "Barnhoorn04"])
    write_csv(
        OUT / "calcite" / "derived_LowT_Barnhoorn04_7pts.csv",
        ["gamma_simple", "J", "M_skemer", "T_C", "source"],
        rows_lo,
        ["Calcite mid-low-T regime (T ≤ 650°C, subgrain-rotation recryst).",
         "Used for case 4 (anisoDelta_Calcite_LowT) calibration.",
         f"Barnhoorn+04 500°C and 600°C samples filtered to J ≤ {J_FILTER:g} ({len(rows_lo)} pts)."],
    )

    # ----- HIGH-T (case 5) — Pieri + Barnhoorn 727 (J≤15) ------------------
    # Bruijn+11 removed after the 2026-05-08 audit.
    rows_hi = []
    for s, g, J, T in zip(p01["sample"], p01["gamma"], p01["J"], p01["T_C"]):
        rows_hi.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}",
                        f"{T:.0f}", f"Pieri01:{s}"])
    for g, J, T, n in zip(b04["gamma"], b04["J"], b04["T_C"], b04["note"]):
        if T == 727.0 and J <= J_FILTER:
            rows_hi.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}",
                            f"{T:.0f}", f"Barnhoorn04 727°C ({n})" if n else "Barnhoorn04 727°C"])
    write_csv(
        OUT / "calcite" / "derived_HighT_8pts.csv",
        ["gamma_simple", "J", "M_skemer", "T_C", "source"],
        rows_hi,
        ["Calcite high-T regime (T > 650°C, grain-boundary-migration recryst).",
         "Pieri+01 (727°C) + Barnhoorn+04 727°C (J ≤ 15) — 8 pts.",
         "Bruijn+11 removed after the 2026-05-08 audit (audit/Bruijn_etal_2011_audit.md).",
         "Used for case 5 (anisoDelta_Calcite_HighT) calibration."],
    )


# ---------------------------------------------------------------------------
# Quartz
# ---------------------------------------------------------------------------
def derive_quartz() -> None:
    print("[quartz]")

    p10 = load_raw("quartz/raw_Pennacchioni_etal_2010_JGR.csv")
    b24 = load_raw("quartz/raw_Blackford_etal_2024_Tectonics.csv")

    rows = []
    for g, J, r in zip(p10["gamma"], p10["J"], p10["regime"]):
        rows.append([f"{g:.4g}", f"{J:.4g}", f"{j_to_m_skemer(J):.4g}", str(int(r))])
    write_csv(
        OUT / "quartz" / "derived_Pennacchioni10_17pts.csv",
        ["gamma_simple", "J_ODF", "M_skemer", "regime"],
        rows,
        ["Pennacchioni+10 ODF J → Skemer M.  17 pts.",
         "Used for case 3 (anisoDelta_Quartz, ODF-J-based) calibration."],
    )

    rows = []
    for s, g_oct, pfJ, dom in zip(b24["sample"], b24["gamma_oct"], b24["pfJ"],
                                    b24["strain_domain"]):
        gamma_eq = float(gamma_oct_to_gamma_simple(g_oct))
        rows.append([f"{gamma_eq:.4g}", f"{g_oct:.4g}", f"{pfJ:.4g}",
                     f"{j_to_m_skemer(pfJ):.4g}", str(int(dom)), str(s)])
    write_csv(
        OUT / "quartz" / "derived_Blackford24_38pts.csv",
        ["gamma_simple_eq", "gamma_oct", "pfJ", "M_skemer_pfJ", "strain_domain", "sample"],
        rows,
        ["Blackford+24 with γ_simple_eq from γ_oct via plane-strain inversion.",
         "Skemer+05 J→M conversion applied to pfJ (note: pfJ < ODF J, so M values",
         "are ~5× smaller than Pennacchioni+10's; the case-6 slope=11 compensates).",
         "Used for case 6 (anisoDelta_Quartz_Blackford) calibration."],
    )


# ---------------------------------------------------------------------------
def main() -> None:
    derive_olivine()
    derive_calcite()
    derive_quartz()
    print("\nAll derived CSVs regenerated under misc/aniso_fstrain/data/.")


if __name__ == "__main__":
    main()
