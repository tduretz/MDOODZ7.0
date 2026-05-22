# case-15 quartz-mica polyphase calibration — quartz CPO decay perspective

**Status**: current. Earlier V1/V2/V3 compositional-mix calibrations produced
R² = 0.92 but with 52–95 % uncertainty on three of four free parameters; the
committed form here is a 2-parameter pooled decay fit on the polyphase
substrate only.

**MDLIB function**: `anisoDelta_QuartzMicaPolyphase` in
`MDLIB/FlowLaws.c` (aniso_db = 15).

**Companion case**: case-8 (`anisoDelta_QuartzMicaBracketed`) is the
*mica* perspective on the same rock; case-15 is the *quartz* perspective.

## Substrate

Tokle, Hirth, Stünitz (2023) *JSG* 169, 104835 Fig 12: quartz c-axis
pole-figure maxima in m.u.d. (multiples-of-uniform-distribution) units,
in general-shear Griggs experiments at T=800°C, P=1.5 GPa, 0.1 wt% H₂O.

Raw CSV: `data/quartz_mica/raw_Tokle_etal_2023_JSG.csv`.

Polyphase substrate used for case-15 (X_ms=0% monomineralic excluded):

| X_ms (%) | γ ≈ 0 | γ ≈ 0.6 | γ ≈ 4 |
|---|---|---|---|
| 5  | 6.24 | 6.35 | 4.78 |
| 10 | 6.00 | 5.54 | 4.77 |
| 25 | 5.94 | 5.35 | 4.40 |

n = 9 datapoints, 3 unique γ values × 3 X_ms values. Across all X_ms ∈
{5, 10, 25}% the polyphase qz pole-figure max **decreases** with γ —
mica-mediated shielding (Tokle+23 §3.2.3) suppresses quartz CPO
development. The pooled fit collapses X_ms ∈ {5, 10, 25}% into a single
γ-dependent decay.

## Model — pooled single-decay

```
M_mud(γ) = M_init · exp(−γ / γ_decay)
```

Fit via `scipy.optimize.curve_fit`, bounds [0.1, 0.5] / [50, 200], p0 =
[6.0, 15.0]. Reproducible by `python3 data/quartz_mica/case15_lsq_pooled.py`.

| Param | Value | 1σ | Rel. err |
|---|---|---|---|
| M_init | **6.0229** | 0.1442 | 2% |
| γ_decay | **15.3533** | 2.8565 | 19% |

Fit quality:

| Metric | Value |
|---|---|
| N | 9 |
| SSE | 0.7188 |
| **R²** | **0.8202** |
| RMSE | 0.2826 m.u.d. |

Hand-pooled target (hand-pooled estimate): M_init ≈ 6.0,
γ_decay ≈ 15.4, R² ≈ 0.84. LSQ confirms hand-pooled estimate at all
three metrics.

Snapshot: `data/quartz_mica/case15_lsq_pooled.{py,txt,png}`.

## m.u.d. → Skemer M-index conversion (gate-b documentation)

Tokle+23 reports pole-figure maxima in m.u.d.; MDLIB's δ-form
(`δ = slope · M + 1`) takes Skemer M-index. A conversion is required.

Three paths considered (spec):

1. **Find another paper that reports both metrics on quartz c-axis
   fabrics (cleanest).** Searched bib + web sweep —
   **NEGATIVE**. No paper reports both Skemer M-index AND m.u.d. pole-
   figure max on the same quartz c-axis dataset. (Bunge J-index ↔ Skemer
   M is the only documented relation, via Skemer 2005 lab calibration.)

2. **Use Bunge's J-index → M-index Skemer formula M = (J−1)/15 if
   Tokle+23 also reports J alongside m.u.d.** Inspected raw CSV — only
   `pf_max_mud` column. **NOT AVAILABLE.**

3. **Defensible proxy `(mud − 1) / 14`.** iter
   #674 explicit fallback. Rationale: pole-figure max for a fully
   developed quartz c-axis CPO is ~15 m.u.d. (Skemer 2005 / Hansen+14
   lab data on quartzite), and for a random fabric is 1.0 m.u.d.; the
   linear normalisation `M ∈ [0, 1]` follows. **CHOSEN.**

   - M_init_Skemer = (6.0229 − 1) / 14 ≈ **0.3588 ≈ 0.359**
   - Watson-MC theoretical cross-check (research notes
     `mud_to_M_conversion_watson.md`): for the same κ that
     produces m.u.d. ≈ 6.0, Watson Monte-Carlo gives Skemer M ≈ 0.26.
   - The 0.36 / 0.26 spread (factor 1.4×) is the conversion uncertainty
     class; **(mud − 1)/14 is the documented choice** per the model-selection spec.

This conversion is **approximate** and should be revisited if a paper
reporting both metrics surfaces for quartz c-axis fabrics. Until then,
the case-15 calibration carries the conversion uncertainty as scope
caveat.

## Implied MDLIB parameters

```
M_init      = 0.359   (Skemer M-index, via (mud−1)/14)
γ_decay     = 15.35
slope       = 15.0    (mirrors case-8 BRACKETED; defensible mid-range
                       for quartz CPO; Hansen+14 olivine 24.5 over-
                       predicts for quartz)
δ(γ = 0)    = 15.0 · 0.359 + 1 = 6.39
δ(γ → ∞)    = 1.0     (full CPO erasure; unphysical floor — flag as
                       scope caveat for γ > 5 extrapolation)
```

## Scope of validity

- **X_ms ∈ [5%, 25%]**: case-15 BRACKETED valid; use the decay form.
- **X_ms < 5%**: NOT case-15. Use case 3 (Pennacchioni+10 audit-control)
  or case 6 (Blackford+24 lab) for near-monophase quartz, saturating-exp.
- **X_ms > 25%**: NOT case-15. Fall back to case 8 BRACKETED (mica-
  dominated regime; mica CPO is the relevant rheological signal).
- **γ ∈ [0, ~5]**: calibrated range. Tokle+23 data span γ ∈ [0, 4];
  per the bilateral-straddle criterion, γ_max/γ_decay = 4/15.35 = 0.26 ≪ 1
  → saturation regime not sampled. **γ_decay is constrained by the
  early-γ slope, not by an observed asymptote.** Extrapolation beyond
  γ ≈ 5 is uncertain and the δ(γ → ∞) = 1.0 asymptote is unphysical
  (should saturate at a non-zero floor).

Scope enforcement: documentation contract only. MDLIB has no per-marker
X_ms field today; user picks `aniso_db = 15` based on their geological
setup.

## LSQ-discipline pre-flight (cardinal-arity + bilateral-straddle)

Applied via `data/calibrate/lsq_quality_gate.py paired_pre_lsq_check`:

- **cardinal-arity check**: n_unique_γ = 3, k = 2 → necessary (k+1=3) ✓
  but not sufficient (need k+3=5). **MARGINAL** — residual dof tight,
  expect γ_decay to be poorly constrained from above.
- **bilateral-straddle check**: γ_min/γ_decay = 0/15.35 = 0 < 0.1 ✓ AND
  γ_max/γ_decay = 4/15.35 = 0.26 < 1.0 ✗. **UPPER_FAIL** — saturation
  regime not reached; γ_decay constrained from below only.
- **Overall**: MARGINAL. LSQ runs but γ_decay carries 19% uncertainty
  as predicted. M_init is well-constrained (2%) because the substrate
  spans the rising part.

This MARGINAL outcome is *expected and acceptable* at the BRACKETED tier
(case-15 is not a Tier-1 LSQ-tier case like cases 1-7). The 19%
uncertainty on γ_decay is documented and the scope caveat above flags
the high-γ extrapolation risk.

## Caveats and known limitations

1. **m.u.d. → M conversion is approximate** (Path 3). Conversion
   uncertainty class ≈ 1.4× (0.26 vs 0.36 across the two routes). Will
   tighten if a dual-metric paper surfaces.
2. **High-γ extrapolation unphysical** (decay-to-zero asymptote). Use
   case-15 only within γ ∈ [0, ~5]. For γ > 5 in qz-mica polyphase
   without case-8 fallback, treat as out-of-scope.
3. **Pooled fit averages across X_ms ∈ {5, 10, 25}%**. Per-X_ms γ_decay
   would refine further (cost: ~1h LSQ per fraction). Deferred until
   case-15 has independent validation against a non-Tokle dataset.
4. **Single-paper anchor**. Only Tokle+23 substrate. No cross-validation
   from another qz-mica polyphase shear experiment with comparable
   metric reporting (Holyoke & Tullis 2006 has 13% Bt natural samples
   but reports neither m.u.d. nor M-index in a directly extractable
   form).
5. **NOT a Tier-1 LSQ-tier case** (those are cases 1-7). case-15 is
   BRACKETED — comparable rigor class to case-8 / case-9 / pwlv=16-17
   Ranalli granulite/mafic-granulite.

## Methodology negative-result note (§S1.5.6 candidate)

The compositional-mix V1 → V2 → V3 mechanism iterations are a real *negative result*:

> The compositional-mix mechanism `M_aggregate(γ, X_ms) = f_DRX(X_ms) ·
> M_recryst(γ) + (1 − f_DRX(X_ms)) · M_init` was attempted in three
> variants (V1 constant f_DRX → R² = −0.88; V2 γ-dependent f_DRX → R²
> = 0.34 with covariance degenerate; V3 X_ms-dependent M_∞ scaling →
> R² = 0.92 but 3/4 parameters with 52-95% uncertainty). The
> pooled-LSQ analysis revealed the misattribution: the
> compositional-mix framework conflates two separate phase-perspectives
> (mica CPO building up with γ vs quartz CPO breaking down with γ).
> The simpler 2-parameter pooled-decay form for the quartz perspective
> (current case-15; M_init = 6.02 m.u.d., γ_decay = 15.35, R² = 0.82,
> only 2 free params) is the correct framework. Methodology lesson: a
> high-R² compositional-mix fit can be a *false positive* if the
> mechanism conflates phase-perspectives; LSQ-discipline on pooled data
> discriminates.

The obsolete V3 calibration is preserved at
(internal research notes) as a
historical record. **This file (`MDOODZ7.0/misc/aniso_fstrain/notes/
case15_calibration.md`) is the canonical case-15 calibration doc**
that MDLIB's `anisoDelta_QuartzMicaPolyphase` docstring points at.

## References

- Tokle, Hirth, Stünitz (2023) *JSG* 169, 104835 —
  `~/bib/anisotropy_calibration_2026/mica_polyphase/Tokle_etal_2023_quartz_muscovite_polyphase_general_shear_JSG.pdf`
- Holyoke, Tullis (2006) — interconnected mica threshold
- Tokle+23 §3.2.3 shielding mechanism
- Skemer 2005 lab M-index calibration; Bunge J-index references in
  research notes `notes/mud_to_M_conversion_watson.md`
