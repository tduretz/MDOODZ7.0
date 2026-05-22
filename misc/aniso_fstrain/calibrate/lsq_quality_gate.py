"""LSQ pre-flight quality-gate helper — paired cardinal-arity + bilateral-straddle.

Two orthogonal pre-LSQ checks for the saturating-exponential form
`M(γ) = M_∞ · (1 − exp(−γ / γ_e))`:

**bilateral-straddle**: the form fits cleanly only when both ratios satisfy
  γ_min / γ_e ≪ 1   (rising part observable)
  γ_max / γ_e > 1   (saturation regime sampled)

**cardinal-arity**: an LSQ with k free parameters requires
  n_unique_γ ≥ k + 1   (necessary — LSQ feasibility)
  n_unique_γ ≥ k + 3   (sufficient — residual dof for fit quality)
  n_unique_γ ≥ k + 5   (strong — F-test for lack-of-fit reliable)

Usage in a calibrate_*.py preamble:

    from lsq_quality_gate import paired_pre_lsq_check
    paired_pre_lsq_check(gammas, k=2, gamma_e_guess=3.0, label="case-3")

prints a multi-line verdict; returns a dict for programmatic decisions.
"""
from __future__ import annotations

from typing import Iterable


def cardinal_arity_check(
    gammas: Iterable[float],
    k: int,
    label: str = "",
    eps: float = 0.05,
) -> dict:
    """Cardinal-arity LSQ precondition check.

    Counts unique γ values (modulo ε-tolerance) and compares to k+1 (necessary)
    and k+3 (sufficient) thresholds.
    """
    g_sorted = sorted(set(round(float(g) / eps) for g in gammas))
    n_unique = len(g_sorted)
    n_min_necessary = k + 1
    n_min_sufficient = k + 3
    n_min_strong = k + 5

    if n_unique < n_min_necessary:
        tier = "FAIL"
        verdict = (f"LSQ underdetermined — n_unique_γ={n_unique} < k+1={n_min_necessary}; "
                   f"expect substrate-gap or degenerate fit.")
    elif n_unique < n_min_sufficient:
        tier = "MARGINAL"
        verdict = (f"LSQ feasible but residual dof tight — n_unique_γ={n_unique} "
                   f"in [k+1={n_min_necessary}, k+3={n_min_sufficient}); expect bounded params "
                   f"but high uncertainty on at least one param.")
    elif n_unique < n_min_strong:
        tier = "SUFFICIENT"
        verdict = (f"LSQ well-constrained — n_unique_γ={n_unique} ≥ k+3={n_min_sufficient}; "
                   f"residual dof = {n_unique - k}.")
    else:
        tier = "STRONG"
        verdict = (f"LSQ strong — n_unique_γ={n_unique} ≥ k+5={n_min_strong}; "
                   f"F-test for lack-of-fit reliable.")

    return {
        "rule": "cardinal_arity",
        "n_unique_gamma": n_unique,
        "k": k,
        "k_plus_1": n_min_necessary,
        "k_plus_3": n_min_sufficient,
        "k_plus_5": n_min_strong,
        "tier": tier,
        "verdict": verdict,
        "label": label,
    }


def bilateral_straddle_check(
    gammas: Iterable[float],
    gamma_e_guess: float,
    label: str = "",
) -> dict:
    """Bilateral-straddle LSQ precondition check.

    Computes γ_min/γ_e and γ_max/γ_e and classifies whether the substrate
    spans both the rising and saturating parts of the curve.
    """
    gs = [float(g) for g in gammas]
    g_min, g_max = min(gs), max(gs)
    r_min = g_min / gamma_e_guess
    r_max = g_max / gamma_e_guess

    lower_ok = r_min < 0.1
    upper_ok = r_max > 1.0

    if lower_ok and upper_ok:
        tier = "PASS"
        verdict = (f"Bilateral straddle ok — γ_min/γ_e={r_min:.3f} < 0.1 AND "
                   f"γ_max/γ_e={r_max:.2f} > 1; both rising and saturating regimes sampled.")
    elif lower_ok and not upper_ok:
        tier = "UPPER_FAIL"
        verdict = (f"γ_max too low — γ_max/γ_e={r_max:.2f} ≤ 1; "
                   f"saturation regime not reached; expect γ_e poorly constrained from above.")
    elif not lower_ok and upper_ok:
        tier = "LOWER_FAIL"
        verdict = (f"γ_min too high — γ_min/γ_e={r_min:.3f} ≥ 0.1; "
                   f"rising regime not sampled; expect γ_e poorly constrained from below "
                   f"(or M_∞ shifted upward).")
    else:
        tier = "FAIL"
        verdict = (f"Both ends fail — γ_min/γ_e={r_min:.3f} (need < 0.1) AND "
                   f"γ_max/γ_e={r_max:.2f} (need > 1); substrate covers neither rising nor "
                   f"saturating regime; LSQ likely degenerate.")

    return {
        "rule": "bilateral_straddle",
        "gamma_min": g_min,
        "gamma_max": g_max,
        "gamma_e_guess": gamma_e_guess,
        "ratio_min": r_min,
        "ratio_max": r_max,
        "tier": tier,
        "verdict": verdict,
        "label": label,
    }


def paired_pre_lsq_check(
    gammas: Iterable[float],
    k: int,
    gamma_e_guess: float,
    label: str = "",
    silent: bool = False,
) -> dict:
    """Run both cardinal-arity AND bilateral-straddle checks.

    Cardinal-arity is the "necessary" filter (must satisfy LSQ-feasibility before
    γ-range analysis matters). Bilateral-straddle is the "sufficient" filter (given
    feasible LSQ, does the γ-range straddle γ_e well enough?).

    Returns dict with both sub-results + overall verdict (PASS / MARGINAL / FAIL).
    Prints multi-line summary unless silent=True.
    """
    arity = cardinal_arity_check(gammas, k=k, label=label)
    straddle = bilateral_straddle_check(gammas, gamma_e_guess=gamma_e_guess, label=label)

    if arity["tier"] == "FAIL":
        overall = "FAIL"
        reason = "cardinal-arity FAIL — LSQ underdetermined regardless of γ-range."
    elif straddle["tier"] == "FAIL":
        overall = "FAIL"
        reason = "bilateral-straddle FAIL — γ-range covers neither rising nor saturating regime."
    elif arity["tier"] == "MARGINAL" or straddle["tier"] in ("UPPER_FAIL", "LOWER_FAIL"):
        overall = "MARGINAL"
        reason = (f"Arity={arity['tier']}, straddle={straddle['tier']} — "
                  f"LSQ will run but expect at least one poorly-constrained param.")
    else:
        overall = "PASS"
        reason = "Both rules pass — expect clean fit with bounded params."

    result = {
        "label": label,
        "arity": arity,
        "straddle": straddle,
        "overall": overall,
        "reason": reason,
    }

    if not silent:
        prefix = f"[{label}] " if label else ""
        print(f"{prefix}=== LSQ pre-flight (cardinal-arity + bilateral-straddle) ===")
        print(f"{prefix}n_total = {len(list(gammas))}, k = {k}, γ_e_guess = {gamma_e_guess}")
        print(f"{prefix}cardinal-arity ({arity['tier']}): {arity['verdict']}")
        print(f"{prefix}bilateral-straddle ({straddle['tier']}): {straddle['verdict']}")
        print(f"{prefix}OVERALL: {overall} — {reason}")

    return result


if __name__ == "__main__":
    # Demo: run on Tokle+23 polyphase (case-15 pooled n=9 substrate)
    # γ ∈ {0, 0.6, 4} × X_ms ∈ {5, 10, 25} = 9 datapoints, 3 unique γ values.
    polyphase_gammas = [0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 4.0, 4.0, 4.0]
    print("=== Demo 1: case-15 pooled decay (k=2, γ_e_guess=15) ===")
    paired_pre_lsq_check(polyphase_gammas, k=2, gamma_e_guess=15.0, label="case-15 pooled")
    print()

    # Demo: case-15 V3 attempted on same substrate with k=4
    print("=== Demo 2: case-15 V3 X_ms-dependent M_∞ on same data (k=4) ===")
    paired_pre_lsq_check(
        [0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.6, 4.0, 4.0, 4.0, 4.0],  # full 12 with X_ms=0
        k=4, gamma_e_guess=20.0, label="case-15 V3",
    )
    print()

    # Demo: case-10 Marti+18 (substrate-gap pilot)
    print("=== Demo 3: case-10 Marti+18 (5 dp all at γ=4.0) ===")
    paired_pre_lsq_check(
        [4.0, 4.0, 4.0, 4.0, 4.0], k=2, gamma_e_guess=3.0, label="case-10 Marti+18",
    )
    print()

    # Demo: case-3 Pennacchioni+10 (clean-fit baseline)
    print("=== Demo 4: case-3 Pennacchioni+10 (n=17, γ ∈ [0.3, 15]) ===")
    pen_gammas = [0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 13.0, 15.0]
    paired_pre_lsq_check(pen_gammas, k=2, gamma_e_guess=3.59, label="case-3 Pennacchioni+10")
