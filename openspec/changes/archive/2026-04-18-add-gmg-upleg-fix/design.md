## Context

The GMG path (`lin_solver = 3`) shipped in `add-gmg-stokes-solver` and was measured in `add-gmg-stokes-defence`. On production-sized grids it stalls: the per-iteration FGMRES residual drops 687× on the first restart (the smoother is working) and then only 1.77× across the next 540 iterations (the V-cycle coarse correction has stopped doing work). Defence §8.4 item 1 localised this to the V-cycle upleg — the recursive return path from coarsest → finest — and specifically identified a **20–40× per-level residual amplification** compounding to ~20× per V-cycle, which saturates FGMRES restart=30. Four pieces of evidence pin this to the upleg rather than any of the other V-cycle components:

1. The first restart drops residual by 687× — so the smoother and the fine-level operator are both correct.
2. `gmg_levels = 3` converges at 201² today — receipt at `benchmarks/ec2/results/vcycle_probe_201_lvl3/mdoodz.log` — so the bug appears as the hierarchy deepens past 3 levels, not on any single operator.
3. The `gmg_dump_vcycle = 1` HDF5 snapshots show the residual field growing between the restrict step and the post-smooth step on levels 3 and 4 — not on the downleg, not on the coarse solve.
4. The identical hierarchy at 41² SolCx produces the spec-promised ρ ≈ 0.1 per V-cycle — so the bug is grid-scale-dependent, consistent with a missing `dx²` factor or an unscaled boundary-row contribution somewhere on the way up.

Everything else about the GMG path is healthy: the stencil bridge matches `StokesAssemblyDecoupled.c` to 1e-12 (golden test), the Vanka smoother reduces residual on every test, the coarse UMFPACK solve is bit-exact to independent invocation, the 41² convergence factor is grid-independent. This change is a focused fix to the upleg plus two companion tidy-ups (the hardcoded `max_restarts = 20` and the misaligned convergence-logging predicate), not a redesign.

## Goals / Non-Goals

**Goals**

- Localise the offending upleg operator with enough precision to explain the 20–40× per-level amplification factor — a single operator, a single file, a specific code path.
- Fix it with ≤ 50 LOC of core-code change, per the `add-gmg-stokes-solver` D8 budget.
- Add a regression assertion that upleg amplification stays bounded on constant-viscosity fixtures, so a future refactor cannot silently reintroduce the bug.
- Add a production-grid convergence fixture (201² at default `gmg_levels`) that currently would fail and passes after the fix.
- Promote the hardcoded FGMRES `max_restarts = 20` to `gmg_fgmres_max_restarts` in the `.txt` parameter family.
- Align the FGMRES convergence predicate and its log message so they report the same quantity relative to `gmg_fgmres_tol`.
- Ship an interim-workaround path (`gmg_levels ≤ 3` clamp for grids ≥ 161²) if the root-cause fix exceeds the 50-LOC budget, so production users have a usable GMG path either way.

**Non-Goals**

- Fine-level Vanka smoother tuning (ω sweep, damping per-level). That is defence §8.4 #2 and is deferred to `add-gmg-smoother-tuning` so its measurements aren't polluted by upleg behaviour.
- Drucker-Prager Newton-mode NaN. That is STATUS.md §1 / defence §8.4 #3 and is a separate mechanism with a separate fix in `add-gmg-newton-plasticity`.
- EC2 regression CI replay. Defence §8.4 #4, belongs to `add-gmg-perf-regression-ci`, naturally sequenced *after* this change so the wall-time numbers reflect the fixed preconditioner.
- Changing the transfer-operator algebra (full-weighting restriction, bilinear prolongation). These were measured correct on 41² SolCx; the bug is either a missing scaling factor or a boundary-row edge case, not a wrong operator.
- 3D extension.
- Any change to `lin_solver ∈ {0, 1, 2, -1}` paths.

## Decisions

### D1: Localisation protocol — isolate the offender by instrumented bisection, not guesswork

**Choice**: before writing any fix, do three scripted probes in order:

1. **Level-count bisection at 201²** — run the same SolViPerf at 201² with `gmg_levels ∈ {2, 3, 4, 5, 6}` and record per-restart residual at each. The step where ρ goes from ~0.1 to ~1.0 pinpoints the level at which the upleg breaks.
2. **Per-stage residual snapshots** — with `gmg_dump_vcycle = 1` on the smallest stalling configuration from step 1, dump the residual field *before and after each upleg stage* (post-smooth, prolongate, restrict on the way in, coarse solve, prolongate on the way out, post-smooth). A Python script computes `‖A·z − r‖ / ‖r‖` per stage. The stage whose ratio explodes is the offender.
3. **Constant-viscosity control** — replay step 2 on a constant-viscosity 201² SolCx fixture. If the amplification persists, the bug is operator-intrinsic (scaling / BC handling). If it disappears, the bug is viscosity-coarsening-dependent (harmonic averaging or coarse-operator consistency).

Only after these three probes commit to a fix. This is the same "diagnose before fixing" discipline that produced the Option-A stencil-bridge decision in the previous change — agent and human both have better judgement when they're deciding *where* the fix goes before writing code.

**Why**: the defence already pinned "upleg" but that's still a multi-hundred-line region of `MultigridStokes.c`. The three-probe protocol narrows it to a single function call, sometimes a single statement. It takes ~4 hours of wall-clock work and saves a week of wrong-direction fixes.

**Alternatives considered**
- *Read the code cold and guess*: rejected — the defence already spent 2h+ of agent time on the wrong-direction "fix" that way.
- *Write a full ω-sweep and level-count-sweep matrix first*: rejected — that's defence §8.4 #2's job, and confounds smoother and upleg variables.

### D2: The three most likely fix candidates (chosen after D1 pins the stage)

Given the signature (fine at one level, broken at many levels, grid-scale-dependent) the shortlist is:

**Candidate A — missing `dx²` or `1/dx²` factor in prolongation.** In staggered-grid multigrid the prolongation of a velocity correction from coarse grid to fine grid must be scaled by the ratio of cell volumes (or an equivalent geometric factor, depending on whether one defines the problem in `∇·σ = f` or `∫∇·σ = ∫f` form). If `MultigridStokes.c::ProlongateVelocity` or its continuity counterpart is missing such a factor, every prolongation under-corrects by `2²` or over-corrects by `2²` depending on the direction — exactly matching a 20–40× compounded amplification across 4–5 levels.

**Candidate B — post-smoother operating in the wrong residual frame.** After the prolongate, the fine-level state has been updated by the coarse correction. The post-smoother has to compute `r = b − A·x` on the *current* state, not on the pre-prolongate state. If the residual array isn't refreshed between prolongate and post-smooth, the smoother "corrects" nothing meaningful and amplifies numerical noise.

**Candidate C — boundary-row contributions dropped on the way up.** The fine-level boundary rows are identity-row (tag `0`, `11`, `30`, etc.) in the decoupled assembly. If the upleg transfer operator is discarding rather than passing-through the boundary-row residual, the fine level accumulates a persistent boundary-row imbalance that the smoother can't heal because it never visits boundary rows.

**Why three candidates and not one**: the signature fits all three. D1's probes narrow to one of the three in deterministic order; the fix is the same-shape ≤ 50 LOC edit in any case (a single scaling factor, a missing `ComputeResidual` call, or a missing boundary-row pass-through).

**Alternatives considered**
- *Full Galerkin RAP at the offending level*: heavyweight, defeats the D4 re-discretisation principle from the previous change. Keep as escape hatch for the interim-workaround branch.

### D3: Amplification-bound contract — ≤ 2× per level transfer, tested by regression fixture

**Choice**: the spec's new requirement says per-level upleg residual amplification SHALL stay ≤ 2× on a constant-viscosity 81² / 161² / 201² SolCx fixture. The regression test (`TESTS/DisabledBenchmarks/UplegAmplificationBound.cpp`) extracts this ratio from the `gmg_dump_vcycle = 1` HDF5 snapshots and asserts per-stage.

**Why 2×, not 1×**: perfect unit amplification across every stage is not achievable — prolongation inherently introduces some high-frequency error that the post-smoother then absorbs, so the residual transiently grows slightly between the two stages. 2× is a generous but finite bound: it allows physical intermediate growth, it catches the 20–40× bug immediately, and it gives margin for small operator tweaks to ship without regressing the test.

**Alternatives considered**
- *1.5×*: tight, risks flaking on anisotropy-enabled fixtures.
- *Geometric-mean bound across stages*: harder to debug when it fails; a per-stage bound localises the regression to one stage.

### D4: Convergence-predicate and log-message alignment

**Choice**: the predicate in `SolveStokesGMG` compares `‖r_k‖ / ‖r_0‖` against `gmg_fgmres_tol`. The log message on convergence reports `‖r_k‖` unadorned. These are different quantities — a relative-reduction tolerance of `1e-6` on an initial residual of 266 declares convergence at `‖r_k‖ ≈ 2.66e-4`, not at `2.66e-4 × 1e-6 = 2.66e-10`, so the log line `final_res = 3.04e+0 ... converged` with `tol = 1e-6` is *mathematically consistent* but *misleadingly worded*.

Align by changing the log message to report the same quantity the predicate tests: `final_res_rel = 1.14e-02, tol = 1e-6` (not converged) instead of `final_res = 3.04e+0 ... converged`. Unit test: a synthetic SolViPerf run with known `r_0` and `gmg_fgmres_tol`; verify the logged `final_res_rel` equals `‖r_k‖ / ‖r_0‖` to 1e-12 and that "converged" only appears when `final_res_rel ≤ tol`.

**Why align rather than reinterpret**: the predicate is doing what `gmg_fgmres_tol` was documented to do in the `add-gmg-stokes-solver` spec — relative reduction. The bug is purely in the log string. Changing the predicate to test absolute residual would be a behavioural change affecting every existing caller; changing the message is purely cosmetic.

**Alternatives considered**
- *Change the predicate to absolute*: breaks `gmg_fgmres_tol = 1e-6` expectations for any user already relying on the documented relative-reduction semantic.

### D5: `gmg_fgmres_max_restarts` parameter — default 20, valid range [1, 1000]

**Choice**: add to `params` struct, parse in `InputOutput.c` with `ReadInt2(fin, "gmg_fgmres_max_restarts", 20)`, validate ∈ [1, 1000] with `LOG_ERR + abort` on out-of-range. Thread through to `FGMRES_GMG` replacing the literal `20` at `MultigridStokes.c:2033`. Same lifecycle as the other `gmg_*` parameters (log at startup, default safe, optional override).

**Why expose this**: diagnostic work on the 201² stall benefits from runs with `max_restarts = 50` or `100` to distinguish "FGMRES needs more Krylov" from "V-cycle is broken". Today this requires a rebuild. Making it a parameter costs ~10 LOC and enables faster iteration.

**Alternatives considered**
- *Hide it behind a `-DDEV` build flag*: rejected — it's a legitimate user-facing tuning knob, same shape as `gmg_fgmres_restart`.

### D6: Interim-workaround trigger — 50-LOC hard threshold, with explicit user-visible note

**Choice**: if the localisation from D1 and the fix from D2 exceed 50 LOC of core-code change (the `add-gmg-stokes-solver` D8 boundary), the fallback is:

- Clamp `gmg_levels` to `min(auto_count, 3)` when any dimension `≥ 161` in `MultigridLevels.c::MultigridComputeDefaultLevelCount`.
- Emit a `LOG_INFO` line at startup naming the clamp: `"GMG: grid ≥ 161 with auto level-count, clamping levels to 3 pending add-gmg-upleg-fix follow-up"`.
- Add a STATUS.md entry to this change naming the clamp as known limitation, and add a defence-update note (separate docs change) calling it out in §8.
- Open a targeted follow-up change named `add-gmg-upleg-deep-fix` whose scope is exactly the root cause left unaddressed here.

The regression test still lands, but with `gmg_levels = 3` asserted against instead of the auto count.

**Why an explicit trigger**: the previous change's D8 was enforced case-by-case; that worked, but left ambiguity about "what does 'half a day' mean?". This change sets a numerical threshold so the fork decision is objective. If the fix is 30 LOC we ship the root-cause fix; if it's 300 LOC we ship the clamp and defer, both times without judgement calls mid-implementation.

**Alternatives considered**
- *Ship nothing until the root-cause fix lands*: rejected — leaves production users without a usable GMG path for an indefinite time.
- *Always ship the clamp and never fix the root cause*: rejected — makes the deep hierarchy worthless forever.

### D7: Test placement — both new fixtures under `TESTS/DisabledBenchmarks/`, neither `DISABLED_`

**Choice**: `UplegAmplificationBound.cpp` and `GmgConvergesAtDefaultLevels.cpp` both live under `TESTS/DisabledBenchmarks/` (which is the existing "GMG-specific regression fixtures" directory per the previous change), both registered as live tests (not `DISABLED_` prefixed). They're expected to fail before the fix and pass after, so they double as "fix verified" receipts and as future regression guards.

The 201² convergence fixture uses `SolViPerf` as its harness (same as the defence measurements) and asserts FGMRES converges in `≤ 60` iterations at `gmg_fgmres_tol = 1e-6` with default `gmg_levels`. Positive control is the archived `benchmarks/ec2/results/vcycle_probe_201_lvl3/` receipt.

**Why both live, not `DISABLED_`**: if either is wrong on the final commit, the fix didn't actually work. `DISABLED_` would hide the failure.

**Alternatives considered**
- *Land the tests first as `DISABLED_` and un-disable with the fix commit*: valid but splits the fix across two commits without gain.

### D8: Scope boundary with `add-gmg-smoother-tuning`

**Choice**: this change does *not* touch smoother damping (`ω`), smoother sweep counts (`nu_pre`, `nu_post`), or the Vanka block assembly logic. If the D1 probes reveal the post-smoother is doing something wrong specifically on the upleg (candidate B in D2), the fix is scoped to "ensure the post-smoother is called with the correct residual state", not "tune the post-smoother's parameters". Parameter-tuning belongs to defence §8.4 #2 and is deferred.

**Why the strict boundary**: confounds kill diagnostics. If this change fixes both upleg wiring and smoother damping, the follow-up `add-gmg-smoother-tuning` has no clean baseline to measure against.

## Risks / Trade-offs

- **D1 probes don't converge on a single offender** → Mitigated by bisection over `gmg_levels ∈ {2…6}` which deterministically localises the failure to a specific level transition; at worst we find that two of candidates A / B / C contribute simultaneously and the fix is a pair of ≤ 25-LOC edits still within the D6 budget.
- **Fix succeeds at 201² but breaks 41² / 81² (tiny-grid regression)** → Mitigated by running the existing dual-solver fixtures (`TopoRelaxGmgEquivalence`, `SolKzGmgEquivalence`, `AnisotropyShearBandGmgEquivalence`) from `add-gmg-stokes-defence` on every commit of this change; any regression surfaces within one CI run.
- **Fix succeeds on constant-viscosity but breaks 10³-contrast** → Mitigated by the D1 step 3 constant-viscosity-vs-variable-viscosity comparison: if the bug is present in both, fix is scale-invariant; if only in contrast-variable, fix is coarsening-specific and the test fixture covers both regimes.
- **`gmg_fgmres_max_restarts` default change breaks someone who was relying on max=20** → The current behaviour hits max=20 on the stalling problem. Since the old behaviour was observably broken, no one is relying on it. New default is still 20, preserving the shipped value.
- **Convergence-logging fix breaks downstream log-parsing** → Low likelihood (log lines are not parsed by any production tool; `benchmarks/ec2/analyse.py` reads the PERF_JSON sentinel, not convergence lines). Mitigated by keeping the log level and the message structure; only the reported quantity changes.
- **Interim workaround ships, root cause never gets fixed** → Mitigated by an explicit commitment in STATUS.md to open `add-gmg-upleg-deep-fix` with a scoped proposal within 30 days of this change landing.
- **New fixture `UplegAmplificationBound` is flaky across anisotropy or free-surface configurations** → Mitigated by making it constant-viscosity SolCx only; anisotropy and free-surface get their own dual-solver-equivalence fixtures from the previous change.

## Migration Plan

1. Land the `gmg_fgmres_max_restarts` parameter + convergence-logging alignment + `UplegAmplificationBound` regression fixture (as a currently-failing live test) — small, zero-risk prerequisite commit.
2. Run the D1 localisation protocol, write up findings under `openspec/changes/add-gmg-upleg-fix/STATUS.md` with the offender name + evidence.
3. Land the root-cause fix per D2 (one of candidates A / B / C). If ≤ 50 LOC per D6, proceed; if > 50 LOC, pivot to D6 workaround.
4. Verify: `UplegAmplificationBound` now green, `GmgConvergesAtDefaultLevels` now green, all existing dual-solver fixtures still green.
5. Re-run the 201² SolViPerf with `gmg_levels` default and confirm FGMRES converges at ≤ 60 iterations.
6. Update STATUS.md with final state, cross-reference the defence's §8.4 entry for this change.

**Rollback**: purely additive to the solver path; reverting is a clean `git revert`. The new parameter has a default = 20 matching the old hardcode, the logging fix is in the message string only, and the fix itself is ≤ 50 LOC within `MultigridStokes.c` or `MultigridLevels.c`. Existing `lin_solver ∈ {0, 1, 2, -1}` paths are not touched.

## Open Questions

- **Exact location in `MultigridStokes.c` of the offending scaling factor or missing residual refresh** — deliberately unresolved until D1 runs. Expected to close within the first session of implementation work.
- **Whether candidate C (boundary-row pass-through) also affects the downleg, not just the upleg** — if yes, the fix has a symmetric counterpart on the restrict side. The constant-viscosity probe from D1 step 3 will indicate.
- **Whether the 201² convergence fixture should require FGMRES iterations ≤ 30 (one restart, healthy) or ≤ 60 (two restarts, generous)** — pick after the fix lands based on the observed per-restart ρ. `≤ 60` is the proposal default; can tighten post-fix if numbers support.
- **Whether to also land a 801² convergence fixture** — out of scope for this change (cost is 50+ minutes wall-clock per run, does not fit CI). Defer to `add-gmg-perf-regression-ci`.
