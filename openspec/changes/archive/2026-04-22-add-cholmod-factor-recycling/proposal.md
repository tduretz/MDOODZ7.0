## Why

The Stokes solver in MDOODZ spends 60–70% of wall time inside `cholmod_factorize` on hard-solve steps (RiftingChenin 1000×800, 4 threads: `solve_s ≈ 12 s/iter`, `nit = 3–4`). Within a single time step, the Picard loop re-factorizes nearly identical matrices — same non-zero pattern, viscosity values change by a few percent per iteration. We currently throw away the previous iteration's factor and redo the full symbolic+numeric work every time. Reusing a stale factor as a preconditioner for iterative refinement should cut iterations 2+ of each step by roughly half, and the benchmark infrastructure (`Breakpoint00020.dat`, sweep script in `skill-stokes-benchmark`) is already in place to measure the gain.

## What Changes

- **New opt-in config flag** `cholmod_recycle` (int, default `0` = off). When `1`, KillerSolver skips numeric re-factorization on Picard iterations 2+ and solves the updated system by applying the stored factor as a preconditioner with iterative refinement (or a bounded FGMRES cycle). Symbolic analysis is already reused; this proposal extends reuse to the numeric factor.
- **Staleness policy** — the recycled factor is used for up to N iterations (configurable via `cholmod_recycle_max_reuse`, default 3) OR until iterative refinement fails to converge in a bounded budget, at which point we fall back to a fresh `cholmod_factorize` for the current iteration. This preserves correctness on steps where the Jacobian changes sharply.
- **Memory**: the KillerSolver factor (`cholmod_factor *L`) is kept alive across Picard iterations instead of being freed after each solve. On RiftingChenin 1000×800 this is ~1–2 GB of extra RSS, which the user has accepted as a worthwhile trade.
- **Threading hooks**: iterative refinement (SpMV + triangular solve + AXPY) is where the multi-thread gain lives — CHOLMOD's triangular solve is serial, but SpMV against the updated matrix and residual AXPYs can run under OpenMP. Design phase will evaluate whether refinement itself parallelizes or whether we should combine recycle + `cholmod_threads≠1` for the fallback factorization.
- **Benchmarks**: compare before/after `solve_s` on the 1000×800 RiftingChenin checkpoint at 1/2/4 threads with `cholmod_recycle = 0` and `= 1`. Acceptance: ≥20% reduction in mean `solve_s` over steps 22–25 at 4 threads, with identical `nit` and final residuals matching to round-off.
- **Backward compatibility**: default off. Every existing `.txt` file produces bit-identical behaviour. Documented in `skill-solvers` and `skill-stokes-benchmark` once the feature lands.

## Capabilities

### New Capabilities

- `cholmod-factor-recycling`: Opt-in reuse of the CHOLMOD numeric factor across Picard iterations within a single time step, with a staleness budget and correctness fallback. Covers the config surface (`cholmod_recycle`, `cholmod_recycle_max_reuse`), the factor lifecycle inside KillerSolver, and the iterative-refinement path used when the factor is reused.

### Modified Capabilities

None. The existing `cholmod-threading` spec governs `cholmod_threads` (BLAS-level threading inside factorization) and is orthogonal — recycling is a control-flow change around when we call `cholmod_factorize`, not how CHOLMOD uses threads internally.

## Impact

- **Code**
  - [MDLIB/Solvers.c](MDLIB/Solvers.c) around `KillerSolver` (~line 1054): factor lifecycle, iterative-refinement loop, fallback path
  - [MDLIB/InputOutput.c](MDLIB/InputOutput.c): `ReadInt2` entries for `cholmod_recycle` and `cholmod_recycle_max_reuse`
  - [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h): new fields on the `params` struct
  - [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c): ensure the recycled factor is freed at end-of-step and at shutdown
- **APIs / config** — `.txt` file surface only. No breaking changes; new keys are optional with sensible defaults.
- **Dependencies** — No new libraries. Uses existing CHOLMOD/SuiteSparse APIs already linked.
- **Memory** — +1–2 GB RSS on 1000×800 runs while the factor is cached. Acceptable per user direction (speed is priority, memory OK to grow).
- **Tests** — New integration test that enables `cholmod_recycle=1` on a small scenario (e.g. SolVi or RiftingChenin 150×100) and asserts velocity/pressure fields match the non-recycled baseline within solver tolerance across 5 time steps.
- **Skills** — `skill-stokes-benchmark` gets an example sweep that toggles the flag; `skill-solvers` documents the new knobs.
