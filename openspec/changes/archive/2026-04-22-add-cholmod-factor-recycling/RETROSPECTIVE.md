# CHOLMOD factor recycling — retrospective

A standalone narrative of the recycled-CHOLMOD experiment. Exists so a future reader can understand why the approach was sound on paper but unworkable in practice, without re-running benchmarks or reading four archived spec files.

**Scope:** the arc from "reuse the factorisation we already computed" to "the matrix dimensions change every step; there is nothing to reuse". Normative material (requirements, design, tasks) lives in the archived change. This file is the *retrospective*, not the record.

**Pointers:**

- `openspec/changes/archive/YYYY-MM-DD-add-cholmod-factor-recycling/` — full spec-driven artifact set
- GMG retrospective (context for why this was next on the priority list): `openspec/changes/archive/2026-04-21-remove-gmg-stokes-solver/RETROSPECTIVE.md`
- All implementation code remains in git history before the revert commit on this branch.

---

## 1. Motivation

At the close of the GMG experiment, the single most promising remaining speed lever was identified as:

> **Reuse CHOLMOD factorisation across solves.** Symbolic analysis should always be reused; numeric factorisation can be reused conditionally when the matrix structure is unchanged and coefficients drift slowly.

The background: on a 1000×800 RiftingChenin step with `nit=2` Picard iterations, each iteration calls `cholmod_analyze` + `cholmod_factorize` + GCR solve — roughly 4 s per Picard step. The factorisation dominates. If the second Picard iteration could skip `cholmod_factorize` and use the previous factor as a GCR preconditioner, the theoretical saving was ~4 s per reuse (≈ 50 % on steps with two Picard iterations, more on steps with three).

The approach was to add a `cholmod_recycle` flag (default off) that controls whether `KillerSolver` skips the numeric factorisation on eligible iterations.

## 2. What was actually built

The implementation went through two distinct strategies before being retired.

### 2.1 Within-step Picard recycling (first strategy)

Recycle at `picard_iteration_index > 0`: the factor from Picard iteration 0 is kept alive and reused on iterations 1, 2, … up to `cholmod_recycle_max_reuse`.

The recycling decision in `KillerSolver`:

```c
int do_factorize = (picard_iteration_index == 0)
                || (recycle == 0)
                || (max_reuse == 0)
                || (reuse_count >= max_reuse);
```

A GCR **probe budget** (6 iterations, one restart cycle) was used to cheaply detect stale factors: if the probe exhausted, zero `du`, refactorize, and retry with the full budget. This minimised wasted work from stale-factor failures (~0.84 s per failed probe vs ~14.7 s per failed full-budget attempt).

Infrastructure added:

- `int cholmod_recycle` / `cholmod_recycle_max_reuse` in `params` struct and parsed in `InputOutput.c`
- `int recycle_reuse_count` in `DirectSolver` struct
- `int picard_iteration_index` parameter threaded through `SolveStokesDecoupled`, `SolveStokesDefectDecoupled`, `KillerSolver`
- `cholmod_start` / `cholmod_finish` moved to program scope (once before the step loop, once at teardown) — the standard CHOLMOD lifecycle pattern matching `ThermalSolver`
- Scope guards silencing the feature for `lin_solver=-1` (compressible) and `lin_solver=1` (KSP)
- Startup log line: `"CHOLMOD factor recycling: enabled, max_reuse = N"` / `"disabled"`

### 2.2 Cross-step recycling (pivot)

When within-step recycling showed ALL probe attempts failing (see §3.1), the feature was pivoted to recycle across time steps instead: keep the factor from step N's final solve and reuse it as the preconditioner for step N+1's first Picard iteration (`picard_iteration_index == 0`).

Rationale: consecutive steps are physically similar; the first Picard iteration of each step starts from the previous step's converged solution, so the factor should be a better preconditioner than within a single step's Picard sequence (which is precisely the sequence designed to move the solution).

A **dim-mismatch guard** was added before the recycle decision:

```c
int dim_changed = (pardi->Analyze == 0 && Lfact != NULL
                   && (size_t)Lcml->nrow != Lfact->n);
if (dim_changed) {
    LOG_INFO("CHOLMOD: matrix dims changed (%zu -> %zu); re-analyzing", ...);
    cholmod_free_factor(&Lfact, &c);
    pardi->Analyze = 1;
}
```

This prevents the `CHOLMOD error: A and L dimensions do not match` crash that would otherwise occur.

## 3. Why it failed

### 3.1 Within-step Picard recycling: the factor is wrong by design

On RiftingChenin steps 22–25 (transient shear-band formation), **every single probe attempt exhausted its 6-iteration budget**. The recycling fallback fired on 100 % of reuse attempts, adding ~0.84 s overhead per Picard iteration with zero benefit.

| | Step 22 (nit=2) | Step 23 (nit=3) |
|-|-----------------|-----------------|
| Baseline (recycle=0) | solve=8.24 s | solve=11.40 s |
| Within-step (recycle=1) | solve=11.02 s (+34 %) | solve=14.62 s (+28 %) |

Root cause: this is not a benchmark problem. Within a single Picard sequence, successive iterations are *designed* to change the matrix — yield activation, viscoplastic regularisation, power-law viscosity updates. A factor computed at iteration 0 is already stale by construction at iteration 1. The probe-budget fallback correctly detects this, but the detection itself costs ~0.84 s on top of the baseline factorisation. There is no configuration where within-step recycling helps on a problem with active plasticity.

### 3.2 Cross-step recycling: free-surface DOF instability

The pivot to cross-step recycling was structurally sound on periodic-BC or fixed-BC problems. On `free_surface=1` simulations, it fails at the symbolic-analysis level.

MDOODZ's free-surface advection (`free_surface=1`) moves the air/solid boundary each step via marker advection. Cells that were "solid" can become "air" and vice versa. This changes the number of *active* equations in the Stokes system: `neq_mom` decreases by a small but non-zero amount on almost every step (observed: 1519599 → 1519591 → 1519587 → 1519581 over steps 21–25 on the 1000×800 RiftingChenin grid).

The CHOLMOD factor `Lfact->n` encodes the symbolic structure for a specific matrix size. When `Lcml->nrow != Lfact->n`, the factor is invalid — CHOLMOD raises `A and L dimensions do not match` and the solve produces garbage or crashes. The dim-mismatch guard catches this and forces re-analysis, but that means:

- The guard fires on *every single step* (dims always change)
- Recycling never triggers (the only condition where recycling is attempted is `!dim_changed && picard_iteration_index == 0`)
- Net result: zero recycling events, pure overhead from the guard check

Benchmark confirmation (steps 21–25, recycle=1):

```
step 21 nit=1  dim_changed: 0 -> 1519599  [initial factorize]
step 22 nit=1  dim_changed: 1519599 -> 1519591  [guard fires, re-analyze]
step 23 nit=1  dim_changed: 1519591 -> 1519587  [guard fires, re-analyze]
step 24 nit=1  dim_changed: 1519587 -> 1519581  [guard fires, re-analyze]
step 25 nit=1  dim_changed: 1519581 -> 1519577  [guard fires, re-analyze]
```

Timing identical to baseline (within measurement noise). The feature does nothing.

### 3.3 The structural incompatibility

Free-surface simulations — the primary MDOODZ workload — are incompatible with CHOLMOD factor recycling because:

1. The **symbolic** factor encodes a specific sparsity pattern for a specific number of DOFs. Once the DOF count changes, the factor is invalid. You cannot recycle the symbolic analysis across steps with a moving free surface.
2. Even if you could fix the symbolic issue (e.g., by over-allocating and masking dead DOFs), the **numeric** factor would still be wrong because the preconditioner is built for a different active set.
3. The only remaining option would be to maintain a fixed-size super-matrix that includes "air" DOFs as identity rows, refactorize only the interior block, and use a permuted solve. This is a fundamentally different architecture from the current `KillerSolver` path and amounts to a partial rewrite of the Schur-complement solver.

Fixed-BC or no-free-surface simulations would work. But the benchmark problem (RiftingChenin) and all production geomechanics runs use `free_surface=1`.

## 4. What was not tried

| Candidate | Why not |
|---|---|
| **Super-matrix with masked air DOFs** | New architecture; amounts to rewriting `KillerSolver`. Complexity out of proportion to a "reuse a factor" optimisation. |
| **Recycle symbolic only, always refactorize numeric** | The numeric factorization is the expensive step (~4 s). Symbolic takes ~0.3 s. Saving 0.3 s at the cost of this work is not worth it. |
| **Detect DOF-stable windows** | Would work in principle (steps where no air/solid transitions occur), but the 1000×800 RiftingChenin grid had at least a few transitions every step in the measurement window (steps 21–25). Whether such windows exist at other stages would require a separate profiling effort. |
| **Factor recycling on no-free-surface runs** | Would work. But MDOODZ's target workload is always free-surface geology. Adding a feature that only helps the minority case — and adds code complexity and maintenance cost — was ruled out. |

## 5. What was correct

The core ideas were sound, and the infrastructure built during this experiment has value:

- **Program-scope `cholmod_start`/`cholmod_finish`** (matching `ThermalSolver`) is the right lifecycle. The per-step version that existed before was an unnecessary allocation/deallocation cycle. This change was reverted with the rest, but if any future work touches the CHOLMOD lifecycle, program-scope is the correct pattern.
- **GCR probe-budget pattern** (detect stale factors cheaply with 6-iteration trial) is a generally useful technique. The break-even analysis (0.84 s probe vs 4 s factorize) is worth keeping in mind.
- **Dim-mismatch guard** (`Lcml->nrow != Lfact->n` before recycling) is the correct safety check. Any future cross-step CHOLMOD reuse must include this guard.
- **22/22 GTests passed** throughout development. The signature changes to `KillerSolver`, `SolveStokesDecoupled`, and `SolveStokesDefectDecoupled` did not regress any tests.

## 6. Why we stopped

The experiment produced a working, correct, default-off implementation that provides zero benefit on the target workload and adds ~0.84 s overhead per iteration when the probe fires. There is no configuration tweak that addresses the structural incompatibility between free-surface DOF instability and CHOLMOD factor recycling.

The code was reverted cleanly. The openspec change was archived as the permanent record.

## 7. What we learned

### What went well

- **Probe-budget fallback** correctly contained the cost of bad recycling attempts. Without it, a naive `do_factorize = 0` path would have wasted the full GCR budget (100 iterations, ~14.7 s) on every failed reuse. The 6-iteration probe (0.84 s) was the right engineering choice.
- **Dim-mismatch guard** prevented silent wrong-answer behaviour. The guard turned a potential `CHOLMOD error: A and L dimensions do not match` crash into a controlled re-analyze path.
- **The benchmark revealed the root cause immediately.** Five steps of tracing showed `dim_changed` firing every step. This is exactly the kind of measurement that makes a clean kill possible — no ambiguity, no "maybe it works on longer runs".

### What went wrong

- **The free-surface DOF instability was not considered at design time.** The design document discussed within-step recycling and cross-step recycling but did not analyse whether the matrix dimensions were stable across steps. A one-step analysis of `neq_mom` across five consecutive steps would have caught the structural incompatibility before any code was written.
- **The within-step probe failure was treated as a signal to pivot.** On reflection, a Picard sequence with active plasticity always changes the matrix — within-step recycling was not salvageable by any probe strategy. The pivot to cross-step recycling was the right call, but the same analysis of DOF stability should have preceded it.

### Intangibles

- Moving `cholmod_start`/`cholmod_finish` to program scope (matching `ThermalSolver`) was a useful cleanup that surfaced the cross-step statefulness question explicitly. Reverted here, but the pattern is documented.
- The probe-budget technique and dim-mismatch guard are transferable to any future solver-reuse work.

## 8. Where to go next for Stokes-solver speed

The CHOLMOD factor recycling direction is closed for free-surface workloads. From the GMG retrospective's priority list:

1. **Audit CHOLMOD/BLAS threading end-to-end.** Ensure `cholmod_threads`, BLAS backend threading, and `OMP_NUM_THREADS` are not oversubscribed or fighting each other. Oversubscription on 4 P-cores + 4 E-cores looks like "CHOLMOD doesn't scale" while being a config issue.
2. **Reduce number of Stokes solves per nonlinear step.** Better nonlinear tolerances, line-search, or Newton-switch heuristics cut the *count* of direct solves, which dominates any per-solve micro-optimisation.
3. **Parallelise / optimise assembly.** Without the GMG stencil-bridge, `StokesAssemblyDecoupled.c` can be edited freely. If profiling shows assembly is a non-trivial fraction of the "Stokes phase", worth investigating.
4. **Block-Schur / Uzawa preconditioned Krylov, longer term.** A new project targeting thread-scalability on large grids.

Explicitly closed: factor recycling via CHOLMOD numeric reuse, for free-surface simulations. The dim-instability is fundamental to the moving free-surface algorithm, not a bug to be fixed.

## 9. TL;DR

We hypothesised that reusing the CHOLMOD numeric factor across Picard iterations (or time steps) would save the ~4 s factorisation cost per solve. We built a correct, default-off implementation with a GCR probe-budget fallback and a dim-mismatch guard. Measurement showed:

1. **Within-step**: 100 % probe failure rate — Picard matrices change too aggressively. Net: +28–34 % overhead, never a saving.
2. **Cross-step**: dim-mismatch guard fires every step — free-surface advection changes `neq_mom` by 4–20 DOFs per step, invalidating the symbolic factor. Net: zero recycling events, feature dead on arrival.

The structural incompatibility is fundamental: CHOLMOD factor recycling requires DOF-count stability across solves; free-surface MDOODZ simulations never have that. The code was reverted and the experiment was archived.
