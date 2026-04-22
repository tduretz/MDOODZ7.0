# GMG Stokes solver — retrospective

A standalone narrative of the GMG experiment in MDOODZ. It exists so a future reader (including future me) can reconstruct the story in one read without going spelunking through three archived OpenSpec changes, the archived performance report, and git logs.

**Scope of this write-up:** the arc from "let's try a multigrid alternative to CHOLMOD" to "we're removing it; CHOLMOD stays". Everything normative (requirements, specs, tasks, design decisions) lives in the archived changes; this file is the *retrospective*, not the record.

**Pointers:**

- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/` — build the thing
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/` — measure it
- `openspec/changes/archive/2026-04-18-add-gmg-upleg-fix/` — fix the big bug that surfaced
- `openspec/changes/archive/2026-04-21-remove-gmg-stokes-solver/` — remove it
- Measurement record (on disk where kept): `benchmarks/ec2/results/PERFORMANCE_REPORT.md` (gitignored)
- All implementation code remains in git history before the removal commit.

---

## 1. Motivation

CHOLMOD was the single Stokes backend in MDOODZ. Direct solvers are fast for 2D but:

- Factorisation cost is superlinear; memory is tight on big grids.
- Reusing factorisation across Picard/Newton iterations was not exploited.
- The code had no Krylov/iterative path to lean on when problems grew.

Geometric multigrid, as a preconditioner for FGMRES on the decoupled Stokes system, was the hypothesis: **near-linear scaling with problem size, lower peak memory, and a foundation for larger problems later**. The literature supports GMG for Stokes with Vanka smoothing on staggered grids, so the architectural bet was reasonable.

Plan: add a parallel `lin_solver = 3` path that uses a MAC-staggered GMG V-cycle as a right preconditioner for FGMRES, keeping CHOLMOD untouched as a fallback and baseline.

## 2. What was actually built

The implementation shipped in three archived changes and ran cleanly at the end:

1. **`add-gmg-stokes-solver`** (first experiment) — from zero to a working `lin_solver = 3` path.
   - `MDLIB/MultigridLevels.{c,h}` — hierarchy build, restrict/prolongate operators, active-mask propagation.
   - `MDLIB/MultigridStokes.{c,h}` — V-cycle driver, Vanka 5×5 block smoother, FGMRES outer loop, HDF5 `gmg_dump_vcycle` instrumentation.
   - `MDLIB/StokesAssemblyGMG.{c,h}` — a deliberately duplicated stencil-level matvec (D11 "Option A"), so the GMG fine-level operator could be bit-identical to MDOODZ's assembled CHOLMOD matrix without editing `StokesAssemblyDecoupled.c`.
   - Golden cross-check test: `BuildStokesOperatorDecoupled`'s assembled `A·x` compared to `ApplyStokesOperatorMDOODZ` to 1e-12 on random vectors, with a CI-only "drift canary" flag forcing a known perturbation to prove the test catches drift.
   - Tuning surface: `gmg_levels`, `gmg_nu_pre/post`, `gmg_fgmres_restart`, `gmg_fgmres_tol`, `gmg_standalone`, `gmg_dump_vcycle`.
   - Around ~4,900 LOC of solver code plus ~20 GMG-specific test fixtures.

2. **`add-gmg-stokes-defence`** — consolidated measurement + writeup.
   - `benchmarks/ec2/` harness (SolViPerf driver, yaml grid configs, CSV pipeline, provision/teardown for EC2, local-mode runner).
   - `benchmarks/vcycle_vis/` — V-cycle animation tooling, reading the HDF5 dumps.
   - `PERFORMANCE_REPORT.md` as the single source of numbers.

3. **`add-gmg-upleg-fix`** — **Fix F** (commit `ad39d7f`).
   - Split `use_mdoodz_matvec` into outer + Vanka flags; L0 Vanka now uses the textbook path. This unstalled the V-cycle upleg amplification.

At the moment of removal the path had: stencil-bridge matvec passing 1e-12 cross-check, Picard GMG ↔ CHOLMOD agreement to ~1e-10 on velocities / ~1e-7 on pressure, Newton GMG ↔ CHOLMOD agreement to ~1e-15 on velocities, end-to-end equivalence tests on SolVi / SolKz / TopoRelax / ShearBand / AnisotropicShearBand, and a V-cycle dumper producing per-stage HDF5 snapshots.

So **the code worked**. The question was always: does it actually go faster?

## 3. What happened when we measured

### 3.1 Pre-Fix F — the upleg amplification stall

On the first real perf sweep, GMG was **60–120× slower than CHOLMOD** on 201²+ SolVi. The V-cycle dumper pinned the cause precisely: each upleg stage (restrict → recurse → prolongate → post-smooth) amplified the fine-level residual by 20–40×, and across 5 levels this compounded to roughly 20× per V-cycle. FGMRES's restart=30 Krylov basis could not keep up. `gmg_levels = 3` worked fine as a positive control, so the bug was specific to deep-hierarchy upleg handling, not the smoother, coarse solve, bridge, or assembly.

### 3.2 Fix F

`add-gmg-upleg-fix` localised the offender and produced a one-flag split: the L0 (finest-level) Vanka pre-smoother was incorrectly going through the MDOODZ-bridge matvec path instead of the textbook multigrid path, so its block-matrix construction was inconsistent with the residual operator it was smoothing. Splitting `use_mdoodz_matvec` into an outer flag (used by the Krylov matvec) and a separate Vanka flag (which stays on the textbook path at L0) fixed the amplification without introducing anything exotic. ~50 LOC region.

Also folded in: `gmg_fgmres_max_restarts` CLI parameter (was a hardcode), alignment of the convergence predicate and its log message so they reported the same quantity, and an upleg amplification regression test asserting per-stage ratio ≤ 2.0 on the constant-viscosity 81² / 161² / 201² fixtures.

After Fix F, GMG **actually converged on production-sized grids** instead of stalling. That closed the infrastructure phase. Now the experiment was finally in a state where a fair apples-to-apples benchmark was possible.

### 3.3 Post-Fix F — the actual verdict

With Fix F in place, headline numbers on an Apple M1 MacBook Pro (4 P + 4 E cores, 16 GiB):

| Problem | Grid | Solver | Wall time | GMG/CHOL ratio |
|---|---|---|---:|---:|
| SolVi, κ=10³ | 81²–321² | CHOLMOD (1 t) | 0.3–4.4 s | — |
| SolVi, κ=10³ | 81²–321² | GMG (1 t) | 0.5–9.5 s | **1.8–2.5×** |
| SolVi, κ=10³ | 801² | CHOLMOD (4 t) | 17 s | — |
| SolVi, κ=10³ | 801² auto-level | GMG (4 t) | 497 s ⚠ fallback | 29× ⚠ stalled |
| SolVi, κ=10³ | 801², `gmg_levels = 4` | GMG (4 t) | 59 s | **3.4×** |
| Rifting, step-20 restart | 151×101 | CHOLMOD (4 t) | ≈ 1 s | — |
| Rifting, step-20 restart | 151×101 | GMG (4 t) | 67 s ⚠ fallback | ~67× ⚠ stalled |

Read row-by-row:

- **SolVi at 81²–321² (the "benchmark family" the solver was designed for):** GMG *converges*, but is **1.8–2.5× slower than CHOLMOD** across the range. Memory win is 7–18 %. 25–60× speedup vs. pre-Fix F, so Fix F worked, but not faster than CHOLMOD, ever.
- **SolVi at 801²:** the `MultigridComputeDefaultLevelCount` auto policy picks 5–6 levels, FGMRES stalls at ρ ≈ 0.96, falls back to CHOLMOD — 29× slower. Manually clamping `gmg_levels = 4` and tightening tol to `1e-7` got it down to 3.4× slower. Still never faster.
- **Rifting (lithospheric, power-law rheology, localised shear bands, viscosity contrasts 4+ orders of magnitude):** catastrophic. Restarting from post-warmup breakpoint, FGMRES drops the high-frequency modes once and then flat-lines at ρ ≈ 1.000 until the 600-iter budget is exhausted. Every time, every tested checkpoint (steps 10 / 15 / 20 / 25). Falls back to CHOLMOD every Picard iteration, wastes ~67 s per linear solve.

Root cause of the Rifting failure is textbook: rediscretised-coarse V-cycle GMG fails on strong, localised viscosity contrasts because the coarse operator can't represent the low-frequency error modes — harmonic viscosity restriction averages across η-jumps at crust/mantle/asthenosphere interfaces and shear bands, so the coarse grid has a *different problem* than the fine grid. This is exactly the regime Galerkin coarsening was invented for.

### 3.4 Auto-level policy as a footgun

Separate from the Rifting issue: the 801² SolVi stall showed the auto-level heuristic happily picks 6 levels on large grids, where the coarsest operator is still being chased by fine-scale modes. `gmg_levels = 4` worked. The heuristic needed a grid-aware ceiling — this would have been `~10 LOC` to add.

## 4. Attempts to make it better that we did *not* pursue

At the decision point we enumerated the known fixes for the Rifting failure mode and ruled each out:

| Candidate | Why not |
|---|---|
| **Galerkin coarsening** (`A_c = R·A_f·P`) | Widens the coarse-operator bandwidth, which the MDOODZ code path has been deliberately avoided. User ruled out up front. |
| **Algebraic multigrid (AMG)** | New external dependency, large integration surface, different license/maintenance profile. Out of scope. |
| **Rediscretised coarse with η-sensitive restriction** | Known in the literature to help marginally but not to rescue the Rifting-class regime; at best buys back a fraction of the stall, not a win. |
| **More smoother sweeps / tighter FGMRES** | Pre-Fix F analysis already showed the ceiling here: the coarse operator is wrong; smoothing harder at the fine level doesn't recover. |
| **Block-Schur / Uzawa preconditioner** | Different architecture, essentially a new change of comparable size to the original GMG work. Candidate for a later project, not a patch on this one. |

What we *did* try and land: Fix F itself (the big one), `gmg_fgmres_max_restarts` parameter, convergence-log alignment, the upleg amplification regression test, the 201² default-levels convergence fixture, the Picard/Newton equivalence tests, and the auto-level heuristic's interim clamp (inside `add-gmg-upleg-fix`).

## 5. Why we stopped

The user priority coming into the benchmarking phase was **speed, not memory**. CHOLMOD fits every current production grid in RAM on the target hardware, so the memory-saving argument for GMG evaporated. Without that, GMG had to beat CHOLMOD on wall time to justify its maintenance cost. It never did — at best 1.8× slower, at worst falling back to CHOLMOD and costing 67 s per solve on Rifting.

The maintenance tax was real: ~4,900 LOC of solver code, ~20 test fixtures, a separate stencil-bridge translation unit that had to stay in sync with `StokesAssemblyDecoupled.c` (enforced by the golden cross-check), an EC2/local benchmark harness, and V-cycle visualisation tooling. Every future change to the Stokes assembler would need to think about the bridge. Every solver parameter added to `params` was an extra field. Every .txt input with a `gmg_*` key was a potential point of user confusion.

The decision was to **close the experiment with a negative outcome, clean**. Keep the three archived OpenSpec changes as the permanent record (proposal + design + tasks + specs survive in the archive; the process-only files — status trackers, retrospectives, defence write-ups — were cleaned out). Remove the code in one surgical commit. Preserve the three non-GMG OMP-guard fixes that landed during the GMG era, because those had nothing to do with the experiment and we did not want to regress non-OMP builds on the way out.

## 6. What we learned

### What went well

- **Correctness-first infrastructure.** The stencil-bridge golden cross-check (1e-12 tolerance, drift canary, ran on every CI) meant we never had to wonder whether GMG was solving the *same problem* as CHOLMOD. When numbers disagreed, it was always at the linear-algebra stage, never at the assembly stage. That discipline is transferable and worth keeping in any future linear-algebra project.
- **Benchmark-driven closure.** `PERFORMANCE_REPORT.md` gave us unambiguous numbers to point at when deciding to stop. No "it's probably fast enough" discussions.
- **V-cycle dumper was the right diagnostic.** Without per-stage HDF5 snapshots of residual norms, Fix F would have taken weeks instead of days. For any future multigrid-adjacent work, build the dumper first.
- **Checkpoint-restart A/B benchmarking.** The Rifting comparison only became honest once we could restart both solvers from the same post-warmup breakpoint. Picking that strategy saved a lot of wasted running time.

### What went badly

- **We built the thing before benchmarking it.** The original plan assumed GMG would beat CHOLMOD on 2D because "that's what the literature says". It doesn't, for our matrix sizes and our viscosity regimes. A one-day spike on a toy problem before 4,900 LOC would have probably caught this.
- **Auto-level heuristic was too trusting.** `MultigridComputeDefaultLevelCount` was written to match the defence's small-grid fixtures and over-deepened at 801². The grid-aware ceiling ended up being a trivial fix (promoted to "follow-up #2" in the performance report before we closed the experiment entirely).
- **Pre-Fix F stall was a distraction, not a smoking gun.** The 60–120× pre-Fix F numbers made it feel like "just fix that bug and everything works". It actually took Fix F *plus* the level-clamp *plus* the tol tightening to even get on the same order of magnitude as CHOLMOD on SolVi 801². The Rifting failure was there the whole time; it just didn't surface clearly until Fix F unmasked it.
- **"Literature says it works" ≠ "works in MDOODZ".** Geomech rheologies with 4+ orders of magnitude of localised viscosity contrast are not the regime most GMG Stokes papers are benchmarking against. The textbook failure mode we hit (Rifting stall) is documented — we just didn't weight it heavily enough when we started.

### Intangibles

- Writing `StokesAssemblyGMG.c` forced a clean read-only restatement of the MDOODZ stencil, which is a useful pedagogical artifact even after the code leaves the tree. The golden-test pattern is also portable.
- The benchmark harness under `benchmarks/ec2/` (now removed) was well-designed in its own right — 5 repeats per config, median + IQR, throttle rejection, cost safeguard, provision/teardown. If we ever do a proper perf harness again, much of its structure is worth re-implementing.

## 7. Where to go next for Stokes-solver speed

Speed is still the bottleneck for Picard/Newton loops; GMG was just the wrong way to attack it. Ordered by likely payoff for the current 2D workloads:

1. **Reuse CHOLMOD factorisation across solves.** Symbolic analysis should *always* be reused; numeric factorisation can be reused conditionally when the matrix structure is unchanged and coefficients drift slowly. This is the recycled-CHOLMOD direction we discussed. Biggest likely win on existing workloads.
2. **Audit CHOLMOD/BLAS threading end-to-end.** Make sure `cholmod_threads`, BLAS backend threading, and `OMP_NUM_THREADS` aren't fighting each other. Oversubscription and wrong BLAS can look like "CHOLMOD doesn't scale" while actually being a config issue.
3. **Reduce number of Stokes solves per nonlinear step.** Better nonlinear tolerances / line search / Newton-switch heuristics can cut the count of direct solves, which usually dominates any per-solve micro-optimisation.
4. **Parallelise / optimise assembly.** Post-GMG we no longer have to think about the stencil-bridge, so assembly can be edited freely. If profiling shows assembly is a non-trivial share of the "Stokes phase", this is worth it.
5. **Block-Schur / Uzawa preconditioned Krylov, if and only if the above plateau.** A *new* project, sized similarly to the original GMG work, targeting the long-term thread-scalability question rather than a percentage improvement on current grids.

Explicitly out of scope forever for the current architecture (absent a reason to revisit):

- More GMG tuning. The ceiling is the coarse-operator representation; tuning does not lift it.
- Galerkin coarsening or AMG inside the current code path.

## 8. TL;DR

We hypothesised GMG-FGMRES would be faster than CHOLMOD for MDOODZ's 2D Stokes problem. We built it carefully (golden cross-check at 1e-12, full Picard/Newton agreement, end-to-end equivalence tests), benchmarked it rigorously, found a real upleg-amplification bug and fixed it, and the honest post-fix result is: **1.8–3.4× slower than CHOLMOD on SolVi, stalls catastrophically on Rifting**. The Rifting failure is a textbook rediscretised-coarse GMG failure mode that Galerkin coarsening would address, but Galerkin was ruled out for bandwidth reasons, and without that there is no path for current GMG to beat CHOLMOD on our workloads.

The implementation was removed cleanly in `remove-gmg-stokes-solver`. Three archived OpenSpec changes preserve the full experimental record. The next speed effort should be recycled CHOLMOD factorisation, not more multigrid.
