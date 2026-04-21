## 1. Prerequisite scaffolding — parameter, logging fix, failing regression tests

Done first as one zero-risk commit: no numerical behaviour changes, only adds a new CLI knob, realigns a log-message string, and adds two regression tests that are *expected to fail* before the root-cause fix.

- [x] 1.1 Add `gmg_fgmres_max_restarts` field to `params` struct in `MDLIB/include/mdoodz.h` (int, default 20).
- [x] 1.2 Parse `gmg_fgmres_max_restarts` in `MDLIB/InputOutput.c` with `ReadInt2(fin, "gmg_fgmres_max_restarts", 20)`; validate `[1, 1000]` with `LOG_ERR + abort` on out-of-range.
- [x] 1.3 Log the active value at startup next to the other `gmg_*` parameters. (Implemented as a single `LOG_INFO("GMG-FGMRES config: …")` in `SolveStokesGMG` right after hierarchy allocation — echoes levels/nu_pre/nu_post/restart/max_restarts/tol/standalone/dump once per Stokes solve.)
- [x] 1.4 Thread `model.gmg_fgmres_max_restarts` through to `FGMRES_GMG` in `MDLIB/MultigridStokes.c` — replace the literal `20` at the `max_restarts = 20` hardcode (approx line 2033).
- [x] 1.5 In `SolveStokesGMG` / `FGMRES_GMG`, locate the convergence predicate. Confirm it tests `‖r_k‖ / ‖r_0‖ ≤ gmg_fgmres_tol` (relative reduction). (Old code tested `‖r_k‖ ≤ tol·‖b‖`; with `x_0 = 0` zeroed at entry these coincide. Replaced with explicit `r0_norm` capture on the first restart so the predicate now reads `‖r_k‖/‖r_0‖ ≤ tol` literally.)
- [x] 1.6 Rewrite the "converged" `LOG_INFO` message to report `final_res_rel = <‖r_k‖/‖r_0‖>, tol = <gmg_fgmres_tol>` instead of the raw `final_res = <‖r_k‖>`. Verify on a local SolViPerf run that `final_res_rel ≤ tol` holds whenever "converged" is printed. (Verified via `MultigridStokesTests --gtest_filter='*FGMRES*'`: reports `final_res_rel 4.319e-11 (tol 1.000e-10)`, converges at iter 9.)
- [x] 1.7 Rewrite the non-convergence `LOG_WARN` message similarly: report `final_res_rel` (predicate-comparable), not raw `‖r_k‖`. Ensure "converged" cannot appear in the non-convergence branch. (Non-converged message now starts `"GMG-FGMRES did not converge: …"` and echoes `final_res_rel=…, tol=…, ‖r_0‖=…, ‖r_k‖=…`.)
- [x] 1.8 Create `TESTS/DisabledBenchmarks/UplegAmplificationBound.cpp` — see §3.1. Register as a live test in `TESTS/CMakeLists.txt`. Expected to fail until the root-cause fix lands. (Landed with three TEST cases at 81²/161²/201², three new `.txt` fixtures `SolViUpleg{81,161,201}.txt` under `TESTS/SolViBenchmark/`. 81² smoke run: all upleg per-stage ratios pass (max 1.237×), whole-V-cycle ratio = **8.72× > 4× bound** — test fails as spec'd, confirming the pre-fix V-cycle amplification.)
- [x] 1.9 Create `TESTS/DisabledBenchmarks/GmgConvergesAtDefaultLevels.cpp` — see §3.2. Register as a live test. Expected to fail until the fix lands. (Landed with `SolViGmgDefaultLevels201.txt` fixture; test captures stdout via `freopen`, regex-greps for `"GMG-FGMRES converged: …"`, asserts iters ≤ 60 + no CHOLMOD fallback.)
- [ ] 1.10 Commit this bundle as one reviewable PR segment: new parameter + logging alignment + failing fixtures. Confirm the existing CHOLMOD / `lin_solver ∈ {0, 1, 2, -1}` paths are still green, and the previous change's dual-solver fixtures (`TopoRelaxGmgEquivalence`, `SolKzGmgEquivalence`, `AnisotropyShearBandGmgEquivalence`) still pass. (Partial: `MultigridStokesTests` 16/16 green post-change; the dual-solver equivalence binaries build but were not runtime-re-verified in this session — defer to §8.1 CI pass or to the commit-prep step.)

## 2. Convergence-predicate + log-message unit test

- [ ] 2.1 Write GTest `ConvergenceLoggingAlignment.LoggedQuantityMatchesPredicate` — drive a tiny `SolViPerf`-like problem with known `‖r_0‖` and capture the log stream; assert `final_res_rel` in the log equals `‖r_k‖ / ‖r_0‖` to 1e-12, and "converged" appears iff `final_res_rel ≤ gmg_fgmres_tol`.
- [ ] 2.2 Include a negative case: force non-convergence with `gmg_fgmres_tol = 1e-30` and `gmg_fgmres_max_restarts = 2`; assert log reports `final_res_rel > tol` and "converged" is absent.

## 3. Regression fixtures — content

### 3.1 `UplegAmplificationBound.cpp`

- [ ] 3.1.1 Drive a constant-viscosity 81² SolCx-analog fixture with `lin_solver = 3`, `gmg_dump_vcycle = 1`, `gmg_levels = 4`.
- [ ] 3.1.2 Parse the HDF5 snapshots from `${writer_subfolder}/vcycle_dump/` to reconstruct per-level pre-op and post-op residual fields.
- [ ] 3.1.3 For each upleg stage (post-smooth, prolongate, coarse-operator apply on the way up), compute `‖r_post‖₂ / ‖r_pre‖₂` and assert ≤ 2.0.
- [ ] 3.1.4 Repeat at 161² (`gmg_levels = 5`) and 201² (`gmg_levels = 5`) — three-resolution coverage per the spec.
- [ ] 3.1.5 Also check whole-V-cycle amplification `‖A·z − r‖ / ‖r‖ ≤ 4.0` per V-cycle on the same runs.

### 3.2 `GmgConvergesAtDefaultLevels.cpp`

- [ ] 3.2.1 Run `SolViPerf` at 201² with `lin_solver = 3`, `gmg_fgmres_tol = 1e-6`, `gmg_fgmres_restart = 30`, `gmg_fgmres_max_restarts = 20`, default auto-selected `gmg_levels`.
- [ ] 3.2.2 Assert FGMRES converges (predicate `‖r_k‖ / ‖r_0‖ ≤ 1e-6` reached) in ≤ 60 total inner iterations (at most 2 restarts).
- [ ] 3.2.3 Assert the dispatch layer did NOT fall back to CHOLMOD (no `LOG_WARN` non-convergence event on the run).

### 3.3 `GmgLevels3PositiveControl.cpp`

- [ ] 3.3.1 Run the same 201² `SolViPerf` but force `gmg_levels = 3` (the pre-fix positive control).
- [ ] 3.3.2 Assert convergence iteration count is ≤ the archived receipt at `benchmarks/ec2/results/vcycle_probe_201_lvl3/mdoodz.log` (within a 20% margin — the fix must not regress the shallow-hierarchy path).

## 4. D1 — localisation protocol (diagnose before fixing)

Runs on the MacBook M1 per STATUS.md §4 of the previous change; EC2 replay is out of scope here.

### 4.1 Probe A — level-count bisection at 201²

- [ ] 4.1.1 Run the same SolViPerf 201² fixture with `gmg_levels ∈ {2, 3, 4, 5, 6}`, record per-restart `‖r_k‖ / ‖r_0‖` from the log.
- [ ] 4.1.2 Record the step where per-V-cycle ρ transitions from ≈ 0.1 (healthy) to ≈ 1.0 (stalled). Document the offending level transition in `openspec/changes/add-gmg-upleg-fix/STATUS.md`.

### 4.2 Probe B — per-stage residual snapshots

- [ ] 4.2.1 Using the smallest stalling configuration from §4.1, run with `gmg_dump_vcycle = 1` and capture residual snapshots before and after each upleg stage.
- [ ] 4.2.2 Write `benchmarks/diagnostics/analyse_vcycle_upleg.py` (or extend the existing `make_gif.py`) that computes `‖r_post‖ / ‖r_pre‖` per stage per level from the HDF5 dumps.
- [ ] 4.2.3 Identify the single stage whose ratio explodes (> 2×). Record in STATUS.md: the level, the stage, the measured ratio, the snapshot filenames.

### 4.3 Probe C — constant-viscosity vs variable-viscosity discriminator

- [ ] 4.3.1 Replay §4.2 on a constant-viscosity 201² SolCx fixture. Record whether amplification persists.
- [ ] 4.3.2 If amplification persists → bug is operator-intrinsic (scaling or BC) → candidate A or candidate C in design D2.
- [ ] 4.3.3 If amplification disappears → bug is coarsening-dependent (harmonic averaging or coarse-operator consistency) → different investigation, possibly a different fix scope altogether. If this branch is hit, reassess the D6 workaround-vs-fix decision in STATUS.md before continuing.

### 4.4 Decision gate — pick the fix path

- [ ] 4.4.1 Based on §4.1–§4.3, pick one of design D2's candidates (A prolongation scaling, B post-smoother residual refresh, C boundary-row pass-through) — or flag a fourth unforeseen cause and write a short design amendment to `design.md` before coding.
- [ ] 4.4.2 Estimate LOC for the root-cause fix. Record in STATUS.md.
- [ ] 4.4.3 If estimated LOC ≤ 50 (D6 threshold): proceed to §5. Otherwise: pivot to §6 interim-workaround path and open the deferred-follow-up.

## 5. Root-cause fix (if ≤ 50 LOC per D6)

### 5.1 Implement the fix — one of candidates A / B / C

- [ ] 5.1.1 Apply the fix in the appropriate file (`MDLIB/MultigridStokes.c` or `MDLIB/MultigridLevels.c`) per §4.4.1.
- [ ] 5.1.2 Keep the edit surgical — only the identified statement/function. Document the one-line before/after in STATUS.md (quote the line + commit hash).
- [ ] 5.1.3 Rebuild; re-run the failing fixtures from §3.

### 5.2 Verify — the regression fixtures go green

- [ ] 5.2.1 `UplegAmplificationBound.cpp` — per-stage ≤ 2×, whole-V-cycle ≤ 4× at 81²/161²/201². Green.
- [ ] 5.2.2 `GmgConvergesAtDefaultLevels.cpp` — 201² converges in ≤ 60 inner iterations at default `gmg_levels`. Green.
- [ ] 5.2.3 `GmgLevels3PositiveControl.cpp` — shallow-hierarchy path is not regressed. Green.

### 5.3 Verify — no regressions on the previous change's fixtures

- [ ] 5.3.1 `TopoRelaxGmgEquivalence` — dVx, dVz, dP bounds unchanged. Green.
- [ ] 5.3.2 `SolKzGmgEquivalence` — 10⁷-contrast SolKz still matches CHOLMOD. Green.
- [ ] 5.3.3 `AnisotropyShearBandGmgEquivalence` — anisotropic fixture still passes its 5e-8 bound. Green.
- [ ] 5.3.4 All existing `lin_solver ∈ {0, 1, 2, -1}` tests — still green (sanity check, no core-path changes).

### 5.4 Symmetry check — does the fix affect the downleg?

- [ ] 5.4.1 If the D1 probes revealed the bug is symmetric (candidate C boundary-row pass-through, which lives on both restrict and prolongate), confirm the restrict-side counterpart is either already correct or needs the mirror fix.
- [ ] 5.4.2 Re-run §5.2 after any mirror fix.

## 6. Interim workaround (if > 50 LOC per D6)

Only executed if §4.4.3 triggers the D6 fallback.

- [ ] 6.1 In `MDLIB/MultigridLevels.c::MultigridComputeDefaultLevelCount`, clamp the auto-selected level count to `min(auto_count, 3)` when any dimension is ≥ 161.
- [ ] 6.2 Emit a `LOG_INFO` line at startup identifying the clamp and naming the deferred follow-up change: `"GMG: grid ≥ 161² with auto level-count, clamping levels to 3 — deep-hierarchy fix deferred to add-gmg-upleg-deep-fix"`.
- [ ] 6.3 Update `UplegAmplificationBound.cpp` (§3.1) to force `gmg_levels = 3` on the 161² / 201² variants; document this as a spec-local relaxation with a cross-reference to the follow-up change.
- [ ] 6.4 Still run `GmgConvergesAtDefaultLevels.cpp` — with the clamp, default levels = 3 on 201², so this fixture is effectively testing the positive-control path and should pass.
- [ ] 6.5 Write a proposal stub for `add-gmg-upleg-deep-fix` under `openspec/changes/add-gmg-upleg-deep-fix/proposal.md` (scaffolded but not yet registered) naming the offending stage + the fix approach that exceeded the budget here.
- [ ] 6.6 Add a STATUS.md entry naming the clamp as a known limitation with explicit scope + pointer to the follow-up change.

## 7. STATUS.md — write as results land

- [ ] 7.1 Create `openspec/changes/add-gmg-upleg-fix/STATUS.md` with the usual structure (executive summary, probes executed, offender identified, fix applied, receipts).
- [ ] 7.2 Entries to record: the D1 probe results (§4), the fix path chosen (§4.4), the measured per-stage amplification before and after (§5.2), the before/after FGMRES iteration count on the 201² convergence fixture (§3.2 assert), any regressions (should be none), the commit hash of the fix.
- [ ] 7.3 If §6 workaround path was taken: document the clamp, the measured before/after iteration counts at `gmg_levels = 3`, and the pointer to `add-gmg-upleg-deep-fix`.
- [ ] 7.4 Write a one-paragraph update for the defence document's §8.4 — does the upleg-fix follow-up close? Or does `add-gmg-upleg-deep-fix` now sit where this change was listed?

## 8. Final validation

- [ ] 8.1 Full CI suite green locally.
- [ ] 8.2 `SolViPerf` at 41² / 81² / 161² / 201² under `lin_solver = 3` at default `gmg_levels` — all converge, FGMRES iteration counts within expected bounds.
- [ ] 8.3 `openspec validate add-gmg-upleg-fix --strict` clean.
- [ ] 8.4 `gmg_fgmres_max_restarts` tunable at runtime — smoke test at `gmg_fgmres_max_restarts = 5` (forces early non-convergence → CHOLMOD fallback; LOG_WARN is emitted with `final_res_rel`, not raw residual).
- [ ] 8.5 Convergence-logging smoke test — `SolViPerf` at 41² with `gmg_fgmres_tol = 1e-6` — confirm `final_res_rel ≤ 1e-6` is printed verbatim next to "converged".
- [ ] 8.6 Update STATUS.md top-line from "implementation in progress" → "ready for archive" once §5 + §7 complete (or §6 + §7 if workaround path).
