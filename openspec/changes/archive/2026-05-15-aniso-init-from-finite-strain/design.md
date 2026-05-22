## Context

The `ani_fstrain == 3` viscous-anisotropy delta-relaxation (capability `aniso-delta-relaxation`, archived change `ani-fstrain3-delta-relaxation`, archive `openspec/changes/archive/2026-05-14-ani-fstrain3-delta-relaxation/`) was designed around an operator split

```
δ_FS_n        = aniso_delta_fn(FS_AR_n)
δ_prod_n      = aniso_delta_{n-1} + (δ_FS_n − aniso_delta_fs_prev_{n-1})
δ_new         = 1 + (δ_prod_n − 1)·exp(−Δt / τ_eff)
aniso_delta_n     = clamp(δ_new, [1, ani_fac_max])
aniso_delta_fs_prev_n = δ_FS_n
```

whose intended cold-limit (τ_eff → ∞, equivalently `Δt/τ_eff → 0`) behaviour is documented (`MDLIB/RheologyParticles.c:755-757`, archived `design.md:49`) as

> "Cold limit (tau → ∞): exp → 1 ⇒ delta_new = delta_prod ⇒ the increments telescope to aniso_delta = aniso_delta_fn(FS_AR), exactly the calibrated ani_fstrain==2 form."

The full algebra (verified independently in the 10-agent code audit, Agent 7) shows that the cold-limit telescope actually evaluates to

```
aniso_delta_N → aniso_delta_fn(FS_AR_N) + (aniso_delta_0 − aniso_delta_fs_prev_0)
```

With the current initialisation in `InitialiseAnisoDeltaParticles` (`RheologyParticles.c:55-69`) which sets `aniso_delta_0 = aniso_factor[phase]`, combined with `PartInit` setting `aniso_delta_fs_prev_0 = 1.0` (`ParticleRoutines.c:908-909`), and the marker `F`-tensor at identity, the cold-limit offset is `aniso_factor − 1`. For the rifting "inherited cratonic fabric" use case (`aniso_factor = 4`), the offset is `+3`, and δ can exceed the Hansen 38-point fit ceiling `24.5·0.536 + 1 = 14.13`. Empirically the v2 rifting run reached `δ_max = 16.04` at 5 Ma.

Three competing fixes were considered (Agent 7, Roman discussion):
1. Documentation-only: rewrite the comment to acknowledge the offset and treat it as "best-case additive inherited fabric." Behaviour-preserving but Hansen-violating.
2. `ani_fac_max = 14.13` cap: clip the overshoot but does not remove the operator-split inconsistency; cold static regions still drift toward Hansen + (`aniso_factor`−1).
3. Init-from-finite-strain: make the three boundary conditions consistent — set `F_init` such that `FS_AR_init` produces `δ_FS_init = aniso_factor` via the existing analytic `aniso_delta_fn`, and set `aniso_delta_fs_prev_init = aniso_factor`. The cold-limit telescope then collapses to `aniso_delta → aniso_delta_fn(FS_AR)` exactly, the Hansen ceiling holds, and the model interprets "inherited fabric of magnitude δ" as "this rock carries `γ_eq(δ)` of past simple-shear-equivalent strain in the `aniso_angle` direction."

This change pursues option 3.

## Goals / Non-Goals

**Goals:**

- Restore the cold-limit invariant `aniso_delta → aniso_delta_fn(FS_AR)` for all `aniso_factor ≥ 1` and `aniso_angle ∈ [0°, 180°)` at every marker, exactly (i.e. independent of the initialisation, no +offset).
- Make the Hansen 38-point fit ceiling `δ ≤ aniso_delta_fn(FS_AR_CAP) ≤ 14.13` for `aniso_db = 1` (olivine), and the analogous mineral-specific ceilings for `aniso_db ∈ {2, 6, 7, 9}`, hold *by construction* under the analytic operator split — without relying on `ani_fac_max` clamping.
- Preserve the "inherited cratonic fabric persists in cold quiescent regions" feature: cold + static markers continue to read `δ = aniso_factor` indefinitely (via the existing strain-rate-gate × Boneh-DiSRX mechanism).
- Keep `aniso_factor = 1` initialisation behaviourally identical to today (`F_init = I`, `δ = 1` everywhere at step 0). No regression on `init1`-style runs.
- Provide a documented physical interpretation: "inherited fabric of magnitude `aniso_factor` and orientation `aniso_angle` ↔ a marker that has experienced `γ_eq` of past simple-shear-equivalent strain in the `aniso_angle` direction."

**Non-Goals:**

- Change downstream consumers of `F` (advection, stress-tensor rotation, post-processing). The change is a per-marker init only; the advection / Cauchy-Green / SVD code paths are unchanged.
- Track lattice orientation as a separate field. The inheritance is encoded in the `F`-tensor's principal directions; the rheology continues to consume only the scalar `aniso_factor_n/_s` and the `nx/nz` director, which are computed downstream from `F` and `aniso_angle` exactly as before.
- Generalise to `ani_fstrain ∈ {1, 4, …}` (other Roman/Annelore variants) beyond what the existing `aniso_delta_fn` dispatch already does. The init mechanism uses whatever `aniso_delta_fn` is configured for the phase.
- Add a new dislocation-density or stored-energy field. Δρ proxy is unchanged.
- Change `aniso_angle` semantics or the rotation convention between marker-frame `F` and lab-frame stress. The change reuses the existing `aniso_angle` definition.

## Decisions

### D1. Invert `aniso_delta_fn` to obtain γ_eq from `aniso_factor`

For `aniso_db = 1` (Hansen 38-point olivine, FlowLaws.c:1006-1013):

```
δ = 24.5 · M + 1
M = 0.536 · (1 − exp(−γ_eff / 3.96))
γ_eff = √FS_AR − 1/√FS_AR
```

The analytic inverse for the magnitude-to-strain map is

```
M(δ) = (δ − 1) / 24.5
γ_eff(δ) = −3.96 · ln(1 − M(δ) / 0.536)
```

defined for `1 ≤ δ ≤ 1 + 24.5·0.536 = 14.13`. For `δ = 4` this gives `γ_eff ≈ 1.027`. The `FS_AR(γ_eff)` inverse is

```
FS_AR(γ_eff) = ((γ_eff + √(γ_eff² + 4)) / 2)²
```

For `γ_eff = 1.027` we get `FS_AR ≈ 2.69`.

For other `aniso_db` cases (calcite Pieri/Bruijn/Schuster, quartz Pennacchioni/Blackford) the same shape `δ = c2·(1 − exp(−γ_eff/γ_e)) + 1` holds with mineral-specific `c2, γ_e` from `FlowLaws.c::aniso_delta_fn_table`. The analytic inverse exists for all currently-active `aniso_db` cases.

Out-of-range handling: if a user configures `aniso_factor > aniso_delta_fn(FS_AR_CAP)` (over-saturation; e.g. `aniso_factor = 20` for olivine where `δ_max = 14.13`), the helper SHALL clamp to the achievable ceiling and emit an `LOG_WARN` rather than producing a non-finite γ_eq.

### D2. Build `F_canonical(γ_eq)` as pure-shear in the principal frame

The simplest `F` producing a target `FS_AR = α²` is pure shear

```
F_canonical(γ_eq) = diag(exp(γ_eq/2), exp(−γ_eq/2))
```

with principal stretches `(α, 1/α)` where `α = exp(γ_eq/2)`. Substituting into the existing closed-form SVD at `RheologyParticles.c:712-732`:

```
E  = (Fxx + Fzz)/2 = cosh(γ_eq/2)
Fd = (Fxx − Fzz)/2 = sinh(γ_eq/2)
G = H = 0
Q = cosh(γ_eq/2),  R = sinh(γ_eq/2)
σ_max = exp(+γ_eq/2),  σ_min = exp(−γ_eq/2)
FS_AR = exp(γ_eq)
```

Substituting in `γ_eff = √FS_AR − 1/√FS_AR = exp(γ_eq/2) − exp(−γ_eq/2) = 2·sinh(γ_eq/2)`. So the `γ_eq` from D1 and the `γ_eq` used to build `F_canonical` are related by `γ_eff = 2·sinh(γ_eq/2)`, i.e. for small γ_eq, `γ_eff ≈ γ_eq`; for larger γ_eq the two diverge.

To avoid this confusion, the helper SHALL solve `γ_eq` directly from the target `FS_AR`:

```
FS_AR = exp(γ_eq)   ⇒   γ_eq = ln(FS_AR)
```

So the full pipeline is `aniso_factor → δ → γ_eff (analytic) → FS_AR (analytic) → γ_eq = ln(FS_AR) → F_canonical`. All steps are closed-form for `aniso_db = 1`.

Alternative (simple shear) `F_canonical(γ_s) = [[1, γ_s], [0, 1]]` is rejected because it has a non-trivial rotation component (vorticity), which would make the "principal frame" rotation in D3 ambiguous. Pure-shear `F` has zero spin and clean principal directions along the eigenvectors; the easy-slip plane direction is then rotated explicitly in D3.

### D3. Rotate by `aniso_angle − π/2` to align with the configured easy-slip direction

**Convention** (verified against `MDLIB/AnisotropyRoutines.c:1209-1210` and `:924,1081`): `aniso_angle` (read from `.txt`, `materials.aniso_angle[phase]`, in degrees, converted to radians at `InputOutput.c:1456`) is the **director** angle — the unit vector `(nx, nz) = (cos(angle), sin(angle))` is **normal to the easy-slip plane**. The easy-slip plane itself lies at angle `aniso_angle − π/2` (lab-frame). The downstream rheology computes `lay_ang = angle − π/2` (`AnisotropyRoutines.c:924,1081`; `RheologyDensity.c:263`) and consumes `(lx, lz) = (cos(lay_ang), sin(lay_ang))` as the easy-slip plane direction.

For the inherited-fabric `F`, the principal extension direction of the past simple-shear-equivalent strain SHALL align with the **easy-slip plane direction** (foliation), which is consistent with A-type olivine fabric where the [100] axis lies in the foliation plane (Hansen 12/14/16; Karato 2008 §15). So the F principal axis must be at `aniso_angle − π/2` from +x, not at `aniso_angle`. The rotation is

```
θ_F = aniso_angle − π/2                              (in radians)
R(θ_F) = [[cos(θ_F), −sin(θ_F)], [sin(θ_F), cos(θ_F)]]
F_init = R(θ_F) · F_canonical(γ_eq) · R(θ_F)^T
```

Equivalently, since `R(θ − π/2) · diag(a, b) · R(θ − π/2)^T = R(θ) · diag(b, a) · R(θ)^T`, the same result is reached by *swapping* the principal stretches in `F_canonical` and using `R(aniso_angle)`. Pick whichever is clearer in implementation; the verification audit (Agent 1) reproduces step-0 δ to machine precision under either form.

**Closed form expansion** (with `θ_F = aniso_angle − π/2`):

```
c = cos(θ_F),   s = sin(θ_F),    a = exp(γ_eq/2),    b = exp(−γ_eq/2)
Fxx =  c²·a + s²·b
Fzz =  s²·a + c²·b
Fxz =  c·s·(a − b)        (= 2·c·s·sinh(γ_eq/2))
Fzx =  c·s·(a − b)        (symmetric; pure-shear F is symmetric)
```

**Step-0 δ is rotation-invariant** under R·F·R^T conjugation (Agent 1 + Agent 4 audits), so the +offset fix (D4) does not depend on the choice of `θ_F`. The 90° correction matters only for the post-step evolution: with `θ_F = aniso_angle` (legacy), the F principal extension would lie *along the director* (perpendicular to the slip plane); subsequent flow vorticity would rotate the inherited fabric the wrong handedness relative to natural CPO. With the corrected `θ_F = aniso_angle − π/2`, the inherited fabric is foliation-parallel from t=0, as observed in cratonic xenoliths (Boneh+17, Tommasi-Vauchez 15).

Note `F_init = F_init^T` for pure-shear (no vorticity). The existing SVD at `RheologyParticles.c:712-732` consumes general F (no symmetry assumed), so no downstream change.

### D4. Set `aniso_delta_fs_prev_init = aniso_factor`

After `F_init` is set, the FS_AR computed by the SVD at step 0 yields exactly `aniso_delta_fn(FS_AR_init) = aniso_factor` (by construction of D1-D3). Setting `aniso_delta_fs_prev[k] = aniso_factor` then makes the step-1 increment

```
(δ_FS_1 − aniso_delta_fs_prev_0) = (δ_FS_1 − aniso_factor)
```

— purely the change in Hansen δ from step 0 to step 1, no permanent +offset. The cold-limit telescope yields `aniso_delta_N → aniso_delta_fn(FS_AR_N)` exactly. The hot-limit pure relaxation continues to drive `aniso_delta → 1 + (aniso_factor − 1)·exp(−Σ·Δt/τ)` toward 1, because both `aniso_delta_0` and `aniso_delta_fs_prev_0` start at `aniso_factor` — the "excess over isotropy" that decays is `aniso_factor − 1`, exactly as the existing hot-limit physics requires.

### D5. Where the init helper runs in the timestep loop

`InitialiseAnisoDeltaParticles` runs once, after `PartInit`, before the time loop begins (`Main_DOODZ.c:374` per the archived change summary). The new behaviour is also invoked once at the same call site. Reseeding / inflow / free-surface marker insertion is NOT modified by this change — new mid-run markers continue to inherit local cell-mean `F` via the existing G2P at the marker insertion site (`ParticleRoutines.c:710-711, 784-785`, `FreeSurface.c:155`). That cell-mean `F` is already the strain-accumulated `F` at insertion time, so the new markers carry the correct inherited fabric automatically.

### D6. `aniso_db` dispatch for the inverse

The forward `aniso_delta_fn[phase]` is set in `FlowLaws.c::ReadDataAnisotropy` per `aniso_db` value. The actual MDLIB dispatch covers **11 currently-active cases**: `{1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15}` (verified by Agent 3 via direct read of `FlowLaws.c:1685-1805`).

| `aniso_db` | Forward function | Form | `c2` | `γ_e` | `δ_max` | Inverse |
|---|---|---|---|---|---|---|
| 1 | `anisoDelta_HansenOlivine` | saturating | 13.132 | 3.96 | **14.132** | analytic |
| 2 | `anisoDelta_Calcite_Combined` (post-2026-05-08 substrate: Pieri+01 ∪ Barnhoorn+04 T-avg; Bruijn 2011 retracted) | saturating | 2.460 | 3.15 | 3.460 | analytic |
| 3 | `anisoDelta_Quartz` (Pennacchioni+10) | saturating | 2.250 | 3.50 | 3.250 | analytic |
| 4 | `anisoDelta_Calcite_LowT` (Barnhoorn ≤650°C) | saturating | 1.320 | 1.37 | 2.320 | analytic |
| 5 | `anisoDelta_Calcite_HighT` (Pieri+Barnhoorn 727°C) | saturating | 4.380 | 5.28 | 5.380 | analytic |
| 6 | `anisoDelta_Quartz_Blackford` | saturating | 2.200 | 4.00 | 3.200 | analytic |
| 7 | `anisoDelta_HansenOlivine_Damped` (Tasaka+Boneh+Kumamoto+Bernard 93-pt) | saturating | 3.920 | 0.76 | 4.920 | analytic |
| 8 | `anisoDelta_QuartzMicaBracketed` (mica-perspective) | saturating | 6.000 | 2.50 | 7.000 | analytic |
| 9 | `anisoDelta_PlagioclaseBracketed` | saturating | 0.960 | 2.00 | 1.960 | analytic |
| 13 | `anisoDelta_Calcite_Superplastic` (GBS BRACKETED) | saturating | 0.120 | 1.00 | 1.120 | analytic |
| 15 | `anisoDelta_QuartzMicaPolyphase` (quartz-perspective) | **decay** | 5.385 | 15.35 | **6.385 at γ=0** | analytic, **reversed monotonicity** |

For saturating cases the analytic inverse is

```
γ_eff(δ) = − γ_e · ln(1 − (δ − 1)/c2)              valid for δ ∈ [1, δ_max)
```

For **case 15 (decay form)** the forward is `δ(γ) = c2 · exp(−γ/γ_e) + 1`, so

```
γ_eff(δ) = − γ_e · ln((δ − 1)/c2)                  valid for δ ∈ (1, c2+1]
```

Note the inverted semantics for case 15: a HIGH `aniso_factor` corresponds to LOW past strain (the strain has not yet erased the inherited CPO), opposite to cases 1–13. A `LOG_WARN` SHALL surface this at init time when `ani_fstrain == 3 && aniso_db == 15 && aniso_factor > 1`.

A symmetric per-phase inverse `aniso_delta_fn_inv[phase]` SHALL be populated alongside the existing forward `aniso_delta_fn[phase]`, with the same dispatch:

```c
typedef double (*aniso_delta_inv_t)(double delta);
double aniso_delta_inv_hansen(double delta);           // case 1
double aniso_delta_inv_calcite_combined(double delta); // case 2
double aniso_delta_inv_quartz(double delta);           // case 3
/* ... cases 4, 5, 6, 7, 8, 9, 13 ... */
double aniso_delta_inv_quartzmica_polyphase(double delta); // case 15 (decay form, reversed)
```

The case-9 plagioclase forward has a downstream d-coupling modifier (`aniso_d_threshold`, `aniso_d_decay`, applied in `AnisotropyRoutines.c::AnisoFactorEvolv:117-121`); the init inverse maps the per-mineral δ-law only — the d-coupling is a forward-only modifier and is NOT deconvolved at init (Agent 3 finding).

For `aniso_db = 0` (no `ani_fstrain`-driven anisotropy) the inverse is NULL and the init helper SHALL skip. The helper `InverseAnisoDeltaFn(δ, phase)` returns the equivalent `γ_eff`, and the caller proceeds with the `γ_eq = ln(FS_AR(γ_eff))` step from D2.

Out-of-range δ (above the analytic saturation, e.g. δ = 20 for olivine when δ_max = 14.13): clamp to `δ_max − ε` (ε = 1e-9), warn, continue with the saturated `γ_eff`. This preserves robustness; misconfigured `aniso_factor` does not crash the model.

### D7. Backwards compatibility — `aniso_factor = 1` is a no-op

When `aniso_factor[phase] = 1`, D1 gives `δ = 1 → γ_eff = 0 → FS_AR = 1 → F_canonical = I`, and D4 gives `aniso_delta_fs_prev = 1`. So the existing `init1` behaviour (`F = I`, `aniso_delta_fs_prev = 1`, `aniso_delta = 1`) is exactly recovered — no per-phase opt-in needed; the helper is unconditionally safe. When `ani_fstrain = 0` on a phase, the helper SHALL skip that phase entirely (no F or aniso_delta_fs_prev write), preserving the legacy isotropic init.

## Risks / Trade-offs

- **R1 — Existing inherited-fabric runs change behaviour.** Any prior run with `aniso_factor > 1` and `ani_fstrain ≥ 1` was producing the +offset behaviour; runs after this change will be Hansen-strict. Documented as **BREAKING** in proposal. Mitigated by: this is a bug fix; the prior behaviour was empirically observable as δ > 14.13 (not what the user wants); the strict version is the documented intent.
- **R2 — Inverse of `aniso_delta_fn` requires per-`aniso_db` cases.** Each existing `aniso_db ∈ {1, 2, 3, 4, 5, 6, 7, 9}` needs its analytic inverse. The mineral-specific calibrations all share the same functional shape `δ = c2·(1 − exp(−γ_eff/γ_e)) + 1`, so the inverse is one closed-form per case. Bisection fallback is available but should not be needed.
- **R3 — Pure-shear F has zero vorticity; if a downstream routine assumes `F = I + ε·G` for small ε, the non-identity init may surprise it.** Mitigation: audit the F-consumer call sites. The known consumers are (a) `FiniteStrainAspectRatio` (SVD, robust by construction), (b) the F advection (`Main_DOODZ.c` velocity-gradient update; robust for any F), (c) HDF5 output (numeric copy; robust). No issues expected, but the audit is on the tasks list.
- **R4 — Symmetric F (no vorticity) means the marker has zero "initial finite rotation."** Consequence: if the simulation begins with `aniso_angle ≠ 0` and the kinematic flow has non-zero vorticity, the F principal direction will be rotated by the subsequent vorticity. This is exactly the desired physics: the inherited fabric is rotated by ongoing deformation, just as in nature. No issue.
- **R5 — `aniso_delta_fs_prev` overwriting may regress `ani_fstrain == 1` (`Beghein`/older) phases.** Mitigated by the conditional `if (materials->ani_fstrain[phase] == 3) { ... }` guard preserved from the existing `InitialiseAnisoDeltaParticles`. The new init applies only on `ani_fstrain == 3` phases. (Future scope to extend to `ani_fstrain == 1, 2, 4` if desired.)
- **R6 — Initial F-tensor at non-identity may make the `FS_AR_CAP` (`= 1e12`) clamp in `FiniteStrainAspectRatio` reachable sooner.** The clamp protects against ill-conditioned F; `FS_AR_init = 2.69` for `δ = 4` is nowhere near the cap. For `aniso_factor = 14` (effectively at saturation), `FS_AR_init = exp(γ_eff(δ=14)) = exp(−3.96·ln(1 − 13/13.13)) ≈ exp(18.3) ≈ 8.9e7` — within the cap but high. Acceptable; users targeting `aniso_factor` near the ceiling are intentionally encoding near-saturated inherited fabric.
- **R7 — Restart from a pre-fix breakpoint will read the legacy `F = I, aniso_delta = aniso_factor, aniso_delta_fs_prev = 1` state and run forward with the post-fix telescope, which still has the +offset (because the state is already inconsistent on disk).** Mitigated by: documenting the change as forward-only (a fresh run, not a restart from a legacy breakpoint, will get the strict-Hansen behaviour). Restart-of-a-post-fix-breakpoint is clean.
- **R8 — Output HDF5 at step 0 changes (F ≠ I, FS_AR > 1).** Plot tooling that asserts F = I at step 0 (none known) would need updating.

## Open Questions

- **OQ1**: Should the helper be invoked unconditionally for every `ani_fstrain ≥ 1` phase (not just `== 3`)? Pro: cleaner reset on user setting `aniso_factor > 1` for any of the legacy `ani_fstrain` modes. Con: changes init behaviour for existing `ani_fstrain ∈ {1, 2}` SETs that intentionally relied on the legacy `F = I, aniso_delta_fs_prev = 1` (e.g. the `AniFstrainQuartz*` calibration SETs). Default: gate strictly on `ani_fstrain == 3`; widen later if asked.
- **OQ2**: For `aniso_db = 0` (no analytic forward) with `aniso_factor > 1`: what should the init do? Default: skip (legacy `F = I`). The legacy +offset would only manifest if someone runs `ani_fstrain = 3` with no `aniso_db`, which is already caught by the `aniso_delta_fn != NULL` guard at `RheologyParticles.c:771-773`.
- **OQ3**: HDF5 output of `aniso_delta_fs_prev` for post-mortem audit — currently this field is per-marker but not written. Adding it as a centroid P2G would make the init-from-FS bug class testable from offline data. Scope: optional follow-up; not in this change.
