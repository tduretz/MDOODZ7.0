## Why

`lin_solver = 3` shipped in `add-gmg-stokes-solver` and was measured in `add-gmg-stokes-defence`. The measurement exposed that on production-sized grids (201²+) the GMG path is **60–120× slower than CHOLMOD** and FGMRES stalls at ~1% relative residual instead of reducing geometrically. The defence's V-cycle-dump analysis (`PERFORMANCE_REPORT.md` §7.5, defence §8.4 item 1) localised the stall precisely: each coarse-correction stage on the V-cycle **upleg** (restrict → deeper recursion → prolongate → post-smooth) amplifies the fine-level residual by 20–40× per level transfer; over 5 levels this compounds to `‖A·z − r‖ / ‖r‖ ≈ 20` per V-cycle, which saturates FGMRES's restart=30 Krylov basis. A positive control exists — `gmg_levels = 3` converges on 201² today — so the bug is specifically in deep-hierarchy upleg handling, not in the smoother, coarse solve, bridge, or assembly.

Until this is fixed, `lin_solver = 3` is not usable at production grids. Everything else in the GMG story (memory savings, correctness, test coverage) is already delivered and passing; the upleg fix is the single remaining blocker for a usable production GMG path.

## What Changes

- **NEW — localise the offending upleg operator**: extend the existing `gmg_dump_vcycle = 1` instrumentation to capture per-stage residual amplification (restrict, prolongate, post-smooth, coarse-operator apply) at every level; analyse the dumps to pin which operator is responsible for the 20–40× amplification factor.
- **FIX — bound the upleg amplification** once the offender is identified. Candidate fixes per defence §8.4 (picked after localisation):
  - missing / incorrect `dx²` (or `1/dx²`) scaling in the prolongation or coarse operator,
  - post-smoother doing the wrong thing with the prolongated correction,
  - boundary-row handling on the way back up discarding stencil contributions it should be adding.
  Expected fix size ≤ 50 LOC in `MDLIB/MultigridStokes.c` or `MDLIB/MultigridLevels.c` once located.
- **NEW — `gmg_fgmres_max_restarts` CLI parameter**: promote the `max_restarts = 20` hardcode at `MultigridStokes.c:2033` (approx) to a tunable parameter with default `20`, matching the pattern of the other `gmg_*` knobs from `add-gmg-stokes-solver`.
- **FIX — convergence-predicate logging bug**: the current `GMG-FGMRES converged: ..., final_res = 3.04e+0` log line fires at `final_res = 3.04` when `gmg_fgmres_tol = 1e-6`, because the predicate and the log message compare different quantities. Read through `SolveStokesGMG`'s convergence check; align predicate and log so the declared tolerance matches the one actually tested.
- **NEW — upleg regression test**: extend the existing V-cycle diagnostic machinery to assert that per-level upleg amplification stays below a contract bound (e.g. `≤ 2×` per level transfer on a constant-viscosity 81² / 161² / 201² SolCx fixture). Fails loud if a future refactor regresses the fix.
- **NEW — production-grid convergence fixture**: an integration test that drives `SolViPerf` at 201² with `lin_solver = 3` at default `gmg_levels` and asserts FGMRES convergence in `≤ 60` iterations (well above the healthy ~20 but well below the current stalling-and-maxing-out behaviour of 600+). Positive control: the same fixture at `gmg_levels = 3` is known to pass today.
- **INTERIM WORKAROUND (conditional)** — if the root-cause fix turns out to be larger than the 50-LOC budget per the `add-gmg-stokes-solver` D8 discipline: clamp `gmg_levels ≤ 3` for grids ≥ 161² in `MultigridLevels.c::MultigridComputeDefaultLevelCount`. This restores convergence at production grids at the cost of a less efficient hierarchy; the defence's §8 Limitations would carry the clamp as a known limitation with a pointer to a further follow-up change.

This change is **non-breaking**: `lin_solver ∈ {0, 1, 2, -1}` paths are untouched. The new parameter has a backward-compatible default. The logging fix changes a message string but not the numerical output. Production runs that previously set `lin_solver = 3` but hit the stall now converge; no production run that previously converged starts failing.

**Out of scope** (explicitly deferred to other follow-ups from defence §8.4):
- Fine-level smoother tuning (`add-gmg-smoother-tuning`, defence §8.4 #2) — should happen *after* this change lands so the measurements aren't polluted by the upleg stall.
- Drucker-Prager Newton-mode NaN (`add-gmg-newton-plasticity`, STATUS.md §1) — different mechanism, different fix.
- EC2 regression CI (`add-gmg-perf-regression-ci`, defence §8.4 #4) — replay the perf sweep on Ice-Lake silicon after this change.

## Capabilities

### New Capabilities

None. The fix refines existing behaviour; no new user-visible surface area beyond the single `gmg_fgmres_max_restarts` knob, which is a small extension of the existing `gmg_*` parameter family and belongs inside the existing capability.

### Modified Capabilities

- `gmg-stokes-solver`: gains three tightened requirements — (a) V-cycle upleg residual amplification SHALL stay bounded (contract to be defined precisely in the spec after localisation, nominal `≤ 2×` per level transfer on constant-viscosity fixtures), (b) the FGMRES convergence predicate and its log message SHALL report the same quantity relative to `gmg_fgmres_tol`, (c) the FGMRES max-restarts count SHALL be configurable via the new `gmg_fgmres_max_restarts` parameter with default `20`. No existing requirement is weakened; these either add new normative content or replace under-specified prior language.

## Impact

**New files**
- `TESTS/DisabledBenchmarks/UplegAmplificationBound.cpp` — the regression test asserting per-level upleg amplification stays below the contract bound.
- `TESTS/DisabledBenchmarks/GmgConvergesAtDefaultLevels.cpp` — 201² SolViPerf convergence fixture at default `gmg_levels`.
- Optional: `benchmarks/diagnostics/analyse_vcycle_upleg.py` if the analysis tooling grows beyond what the existing `make_gif.py` supports.

**Modified files**
- `MDLIB/MultigridStokes.c` — the actual upleg fix (location TBD by localisation), the `max_restarts` parameter plumbing, and the convergence-predicate / logging alignment.
- `MDLIB/MultigridLevels.c` — only if the interim workaround is needed (clamp `gmg_levels` at ≥ 161²).
- `MDLIB/include/mdoodz.h` — add `gmg_fgmres_max_restarts` field to `params`.
- `MDLIB/InputOutput.c` — parse + validate `gmg_fgmres_max_restarts`.
- `openspec/specs/gmg-stokes-solver/spec.md` — MODIFIED delta for the three tightened requirements above.

**APIs / user-facing**
- One new optional `.txt` parameter `gmg_fgmres_max_restarts`, default `20`. Backward-compatible.
- No change to HDF5 output schema, no change to existing parameter semantics.

**Dependencies**
- None added. Uses the existing `gmg_dump_vcycle` instrumentation and HDF5 infrastructure from `add-gmg-stokes-solver`.

**Risks**
- **Localisation takes more than one diagnostic pass**. Mitigated: we already have a positive control (`gmg_levels = 3` works) plus the preserved `perf_detail/` dumps from the defence's MacBook sweep, so the search space is narrow. If the first instrumentation pass doesn't pin the offender, a bisection across levels 3 → 4 → 5 at 201² is tractable.
- **Fix exceeds 50 LOC budget**. Mitigated by the interim-workaround clause: clamp `gmg_levels ≤ 3` for large grids, document in STATUS.md, continue with a targeted follow-up change. Ships a usable (if sub-optimal) production GMG path either way.
- **Fix regresses smaller-grid behaviour**. Mitigated: the `add-gmg-stokes-solver` dual-solver fixtures at 41², 81², 161² are re-run on every change per the previous change's D8 discipline; any regression surfaces immediately.
- **Convergence-logging fix changes observable log format in a way that breaks downstream tooling**. Low likelihood (the log lines are not parsed by anything production), mitigated by keeping the existing INFO-level message structure and only correcting the numerical quantity reported.
