"""Tiny CSV loader + strain-measure conversions for the aniso_fstrain calibrations.

All raw datasets in `../data/{olivine,calcite,quartz}/raw_*.csv` use the same
convention: lines starting with `##` are comment headers (paper citation,
field documentation), then a single header row, then the data.  This module
exposes:

  load_raw(path)                — DataFrame-like read (returns dict-of-arrays)
  j_to_m_skemer(J)              — folkloric M = (J-1)/15 conversion; NOT
                                   actually published in Skemer+05. See
                                   function docstring for derivation.
  fs_ar_from_gamma(gamma)       — strain-ellipse aspect ratio under simple shear
  gamma_from_fs_ar(fs_ar)       — inverse of fs_ar_from_gamma
  gamma_from_axial_eps(eps)     — axial natural strain ε → γ_simple_eq
                                   (uses 2·sinh(ε/2) — see Boneh+14 / Bernard+19 docs)
  gamma_from_xz(xz)             — X/Z aspect ratio → γ_simple_eq (plane strain)
  gamma_oct_to_gamma_simple(oct)— Nadai octahedral shear strain ε_oct → γ_simple_eq
                                   (used for Blackford+24)

DATA_ROOT points at `misc/aniso_fstrain/data/`.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE      = Path(__file__).resolve().parent
DATA_ROOT = (HERE.parent / "data").resolve()


def load_raw(rel_path: str | Path) -> dict[str, np.ndarray]:
    """Read a raw_*.csv file.  `rel_path` is relative to misc/aniso_fstrain/data
    (e.g., 'olivine/raw_Hansen_etal_2014_EPSL.csv') OR an absolute path.

    Returns a dict keyed by column name; each value is a numpy array.
    Numeric columns are float; non-numeric columns are str (object array)."""
    path = Path(rel_path)
    if not path.is_absolute():
        path = DATA_ROOT / path
    if not path.exists():
        raise FileNotFoundError(f"Raw dataset not found: {path}")

    cols: dict[str, list[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.reader(_skip_comments(fh))
        header = next(reader)
        for h in header:
            cols[h] = []
        for row in reader:
            for h, v in zip(header, row):
                cols[h].append(v)

    out: dict[str, np.ndarray] = {}
    for h, vals in cols.items():
        try:
            out[h] = np.array([float(v) for v in vals])
        except ValueError:
            out[h] = np.array(vals, dtype=object)
    return out


def _skip_comments(fh):
    for line in fh:
        if line.lstrip().startswith("##"):
            continue
        if not line.strip():
            continue
        yield line


# ---------------------------------------------------------------------------
# J → M-index conversion (linear, J ≲ 15) — folkloric approximation,
# not Skemer+05 despite common attribution.
# ---------------------------------------------------------------------------
def j_to_m_skemer(J: np.ndarray) -> np.ndarray:
    """Linear J → M conversion: M = max((J-1)/15, 0).

    Citation provenance: this formula is NOT in Skemer, Katayama, Jiang &
    Karato (2005) Tectonophysics 411, despite being widely attributed to
    them. The closest text in Skemer+05 reads "M ≈ (J-1)/X" with X
    unspecified; X = 15 appears in downstream practice as a folkloric
    approximation. (Confirmed by direct read of Skemer+05, by inspection
    of the Skemer & Hansen 2016 review, and by literature web sweep.)

    On natural-olivine xenoliths the formula FAILS: Bernard+19 (G3, n=65
    Table 2) reports M and J side-by-side and the free linear fit is
    M = 0.028·J + 0.021 (R² = 0.71). The Skemer-(J-1)/15 fit on the same
    data has R² = -0.69 — i.e. anti-correlation, worse than the mean.
    Bias of (J-1)/15 vs the Bernard fit is roughly +0.05 in M.

    On lab calcite / quartz / calcite-mica corpora the relationship is
    UNTESTED — no published paper reports M and J side-by-side for those
    minerals.

    Use this function with the awareness that any δ_∞ value derived from
    its output inherits a regime-dependent bias; for natural-olivine
    workflows prefer Bernard+19's fit. Kept here for compatibility with
    the existing case-2/4/5 (calcite) and case-3 (quartz Pennacchioni)
    calibrations whose source CSVs report J only.
    """
    return np.maximum((np.asarray(J) - 1.0) / 15.0, 0.0)


def j_to_m_natural(J: np.ndarray) -> np.ndarray:
    """Natural-sample linear J → M conversion: M = max((J-1)/24, 0).

    Empirical fit derived from Bernard+19 (G3, n=65 Table 2) where M and J
    are reported side-by-side for naturally-deformed olivine xenoliths. A
    free linear fit gives M = 0.028·J + 0.021 (R²≈0.71); the physically-
    anchored form M=(J−1)/24 (zero-pinned at J=1, slope 1/24 ≈ 0.042) gives
    R² = 0.825 on the same n=35 calibration subset. The natural-sample 1/24
    slope is 1.6× lower than Skemer's lab (J−1)/15 slope, reflecting the
    well-documented mismatch between lab-deformed olivine and natural
    xenolith CPO strengths at the same J-index.

    Use this for natural-sample J → M conversion (case 7 damped olivine,
    Bernard+19, any future natural-CPO calibration). Use j_to_m_skemer for
    lab-deformed samples (Hansen group, calcite Pieri/Barnhoorn, quartz
    Pennacchioni). The cut-over is regime-dependent.
    """
    return np.maximum((np.asarray(J) - 1.0) / 24.0, 0.0)


# ---------------------------------------------------------------------------
# Strain-measure conversions
# ---------------------------------------------------------------------------
def fs_ar_from_gamma(gamma: np.ndarray) -> np.ndarray:
    """Forward simple-shear strain-ellipse aspect ratio:

        FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4+1)

    Standard structural-geology result: derive from F=[[1,γ],[0,1]] and
    Cauchy-Green C=F^T F; the larger eigenvalue of C is FS_AR (= λ_max²
    when normalised). See Means (1976) *Stress and Strain* eq. 4.6 or
    Passchier & Trouw (2005) *Microtectonics* §3.3 (textbooks; not in bib).
    """
    g = np.asarray(gamma)
    return 1.0 + 0.5 * g * g + g * np.sqrt(0.25 * g * g + 1.0)


def gamma_from_fs_ar(fs_ar: np.ndarray) -> np.ndarray:
    """Inverse of fs_ar_from_gamma: γ = √FS_AR − 1/√FS_AR.

    Algebraic inverse; same source citations as fs_ar_from_gamma.
    """
    s = np.sqrt(np.asarray(fs_ar))
    return s - 1.0 / s


def gamma_from_axial_eps(eps: np.ndarray) -> np.ndarray:
    """Axial natural (true) strain ε → simple-shear γ that produces the
    same maximum stretch λ_max = exp(ε/2) under plane-strain volume
    conservation:

        γ_eq = 2 · sinh(ε / 2)

    Inline derivation, NOT from Boneh & Skemer 2014 directly: the source
    paper reports ε (their fig. 4 caption labels strain as ε, line 464
    of the PDF text extract) but does not provide a γ_eq conversion. The
    formula above is the standard Hencky-strain ↔ simple-shear-equivalent
    relation under plane-strain volume conservation; see Passchier & Trouw
    (2005) *Microtectonics* §3.3 (textbook; not in bib).

    Used for Boneh & Skemer 2014 Aheim-dunite Griggs-apparatus samples
    (case-7 calibration input).
    """
    return 2.0 * np.sinh(0.5 * np.asarray(eps))


def gamma_from_xz(xz: np.ndarray) -> np.ndarray:
    """Finite-strain ellipsoid X/Z aspect ratio → simple-shear-equivalent
    γ under plane-strain volume conservation:

        γ_eq = √(X/Z) − 1/√(X/Z)

    Inline derivation, NOT from Bernard et al. 2019 directly: Bernard+19
    reports X/Z directly (Table 2 column "X/Z", their text §3.3 line 416
    "X/Y ratio calculated from (X/Z)/(Y/Z)") but does not provide a γ_eq
    conversion. The formula above is the algebraic inverse of
    fs_ar_from_gamma — equivalent to assuming the X/Z ratio is the
    plane-strain principal-stretch ratio.

    Used for Bernard et al. 2019 natural xenolith samples (case-7
    calibration input).
    """
    s = np.sqrt(np.asarray(xz))
    return s - 1.0 / s


def gamma_oct_to_gamma_simple(eps_oct: np.ndarray) -> np.ndarray:
    """Convert Nadai octahedral shear strain ε_oct = √[(1/3) · Σ(ε_i − ε_j)²]
    to simple-shear γ_simple_eq under plane strain (ε₂ = 0, ε₃ = −ε₁):

        ε_oct = √2 · |ε₁|        →     ε₁ = ε_oct / √2
        γ_eq  = 2 · sinh(ε₁)

    Inline derivation. Octahedral-shear-strain definition is canonical
    (Nadai 1950 *Theory of Flow and Fracture of Solids*, §15; textbook
    not in bib). The principal-stretch → simple-shear-equivalent step
    follows the same plane-strain Hencky-strain inversion as
    gamma_from_axial_eps.
    """
    eps1 = np.asarray(eps_oct) / np.sqrt(2.0)
    return 2.0 * np.sinh(eps1)
