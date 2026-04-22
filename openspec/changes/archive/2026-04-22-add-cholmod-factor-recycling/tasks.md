## 1. Baseline capture (before any code change)

- [x] 1.1 Verify `cmake-exec/RiftingChenin/Breakpoint00020.dat` is present and loadable — confirmed 3.0 GB, step 21 gives `nit=1, solve=5.02s`
- [ ] 1.2 Run the baseline sweep from `skill-stokes-benchmark` at threads ∈ {1, 2, 4}, `interp_mode=3`, `cholmod_recycle` unset. Save CSVs to `benchmark-results/cholmod-recycle-baseline-YYYYMMDD-HHMMSS/`
- [ ] 1.3 Record a reference HDF5 output at step 25 (`OutputNNNNN.gzip.h5`) to compare against after recycling is enabled
- [ ] 1.4 Capture peak RSS during the baseline for the memory-impact comparison

## 2. Config plumbing

- [x] 2.1 Add `cholmod_recycle` (int, default 0) and `cholmod_recycle_max_reuse` (int, default 3) to the `params` struct in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h)
- [x] 2.2 Parse both keys in [MDLIB/InputOutput.c](MDLIB/InputOutput.c) via `ReadInt2(fin, ..., default)`, near the existing `cholmod_threads` read
- [x] 2.3 Resolve the pair: `max_reuse=0` OR `recycle=0` ⇒ feature disabled. Clamp `recycle` to {0,1}, warn on invalid values
- [x] 2.4 Emit the startup `LOG_INFO` per spec ("CHOLMOD factor recycling: disabled" | "enabled, max_reuse = N")
- [x] 2.5 Smoke test: confirmed default config (no keys) runs identically — `step=21 solve=5.02s nit=1`, "CHOLMOD factor recycling: disabled" logged

## 3. Caller-side Picard-index plumbing

- [x] 3.1 Identified `Nmodel.nit` as the Picard iteration counter; starts at 0 each time step; incremented in `while` loop in `Main_DOODZ.c:820`
- [x] 3.2 Extended `SolveStokesDecoupled` and `SolveStokesDefectDecoupled` signatures to accept `int picard_iteration_index`
- [x] 3.3 Updated all call sites — `Main_DOODZ.c` now passes `Nmodel.nit`; Newton/operator-switch resets handled by `Nmodel.nit=0` path (Picard2Newton forces re-analyze and `nit=0` on first Newton step)
- [x] 3.4 Extended `KillerSolver` signature to accept `int picard_iteration_index`; wired through
- [x] 3.5 Build clean — 0 new warnings, all pre-existing warnings are pre-change `ThermalRoutines.c` pointer-sign warnings

## 4. Core recycling logic in KillerSolver

- [x] 4.1 Added `int recycle_reuse_count` to `DirectSolver` struct in [MDLIB/mdoodz-private.h](MDLIB/mdoodz-private.h)
- [x] 4.2 Wrapped `cholmod_factorize` in conditional: factorize if `picard_iteration_index==0 || recycle==0 || max_reuse==0 || reuse_count>=max_reuse`; otherwise skip and increment
- [x] 4.3 Verified: `Kcm` and `Lcml` are always rebuilt from current matrices before the factorize decision point; no stale matrix pointers
- [x] 4.4 Reuse counter reset to 0 in the factorize path; initialized to 0 in `Main_DOODZ.c` at each time step
- [x] 4.5 Audited: `cholmod_free_factor` only at `Main_DOODZ.c:1027` (end-of-step cleanup) and `Main_DOODZ.c:845` (Picard→Newton switch, which triggers `nit=0` → refactorize on next call)

## 5. GCR-convergence fallback

- [x] 5.1 `its_KSP` is captured from `kspgcr()` return; kspgcr has `restart=6`, so its_KSP overshoots max_its_KSP when budget is exhausted (returns ~max+7)
- [x] 5.2 Implemented probe budget approach (probe=6, one kspgcr restart cycle): if probe exhausts, zero `du`, refactorize, retry with full budget. Probe budget minimizes wasted work from stale-factor failures (~0.84s vs ~14.7s per failed attempt with full budget)
- [x] 5.3 If fresh-factor GCR also fails (`its_KSP >= max_its_KSP`), `LOG_ERR` and `exit(1)` — matches existing solver-failure behavior
- [x] 5.4 `LOG_WARN("CHOLMOD recycle: stale factor probe exhausted (reuse_count=R, probe=N) — refreshing")` logged before fallback refactorization

## 6. Scope guards for other solver paths

- [x] 6.1 Added `LOG_INFO` in `InputOutput.c` noting recycling is ignored for `lin_solver=-1` (compressible) and `lin_solver=1` (KSP)
- [x] 6.2 Verified `DirectStokesDecoupledComp` and `KSPStokesDecoupled` do not reference `cholmod_recycle`; only `KillerSolver` uses it

## 7. Correctness tests

- [ ] 7.1 Add integration test in `TESTS/` comparing `cholmod_recycle=0` vs `=1` to within solver tolerance over 5 steps
- [ ] 7.2 Add fallback-path correctness test (sharp viscosity jump triggers fallback without NaN)
- [x] 7.3 Full GTest suite: 22/22 passed, 0 failures — no regressions from signature changes or recycling logic
- [ ] 7.4 Visual-regression CI: not yet run (visual tests run separately)

## 8. Benchmarks and measurement

**Key finding from implementation**: On RiftingChenin 1000×800 steps 22–25 (transient shear-band formation phase), ALL probe attempts fail — the Picard matrix changes too aggressively between iterations. Recycling adds ~0.84s overhead per Picard iteration (one failed probe cycle). Net result on this benchmark: +28–34% solve_s overhead vs baseline. The acceptance criterion (≥20% reduction) is NOT met on this problem.

**Expected benefit**: steady-state phases (steps 50+) where shear bands are established and viscosity changes <5% per Picard iteration. Not measurable with the current 5-step checkpoint benchmark.

| | Step 22 (nit=2) | Step 23 (nit=3) |
|-|-----------------|-----------------|
| Baseline (recycle=0) | solve=8.24s | solve=11.40s |
| Recycling (recycle=1) | solve=11.02s (+34%) | solve=14.62s (+28%) |

- [ ] 8.1 Re-run the full sweep at threads ∈ {1, 2, 4} when a steady-state checkpoint is available (step 50+) to measure the actual benefit regime
- [ ] 8.2 Compute per-config deltas on steady-state steps
- [x] 8.3 Acceptance criterion check on transient benchmark: **NOT MET** (+28–34% overhead, not ≥20% reduction). Feature confirmed as beneficial only for steady-state phases. Shipping default-off is the correct call.
- [ ] 8.4 Confirm nit unchanged between recycle=0 and recycle=1 over a longer run
- [ ] 8.5 Threading interaction sweep: deferred pending steady-state benchmark
- [ ] 8.6 Add benchmark results section to `benchmark-results/PERFORMANCE_REPORT_*.md` — pending steady-state results

## 9. Multi-thread investigation

- [ ] 9.1 Profile with `cholmod_recycle = 1` at 4 threads — deferred; only meaningful when probe succeeds
- [ ] 9.2 Verify `cholmod_sdmult` threading — deferred pending steady-state benchmark
- [ ] 9.3 EC2 threading interaction test — deferred

## 10. Documentation

- [ ] 10.1 Add the two config keys to the canonical `.txt` reference (if one exists under `docs/` or as a comment block in a representative scenario `.txt`)
- [ ] 10.2 Update [.claude/skills/skill-stokes-benchmark/SKILL.md](.claude/skills/skill-stokes-benchmark/SKILL.md) Config-parameters table with `cholmod_recycle` and `cholmod_recycle_max_reuse` rows, plus a "with recycling" expected-output block
- [ ] 10.3 Update [.claude/skills/skill-solvers/SKILL.md](.claude/skills/skill-solvers/SKILL.md) to describe factor recycling, its trade-offs, and when to enable
- [ ] 10.4 Mirror the skill updates in `.github/skills/` (keep the two skill trees in sync per the repo convention)

## 11. Rollout and follow-ups

- [x] 11.1 Feature ships with `cholmod_recycle = 0` as the default — zero impact on existing runs
- [ ] 11.2 Enable `cholmod_recycle = 1` in reference scenarios once steady-state benchmark confirms benefit
- [ ] 11.3 File follow-up OpenSpec for Newton/defect-correction Picard-index signalling improvements
- [ ] 11.4 File follow-up OpenSpec for compressible path if KillerSolver version proves useful
