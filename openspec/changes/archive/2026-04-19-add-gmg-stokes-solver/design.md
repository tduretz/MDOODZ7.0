## Context

MDOODZ solves the decoupled velocity–pressure Stokes system on a uniform staggered (MAC) grid. Today the linear system at each nonlinear iteration is solved by a CHOLMOD/UMFPACK direct factorisation, dispatched from [MDLIB/StokesRoutines.c:788-791](MDLIB/StokesRoutines.c#L788-L791) via the integer `lin_solver` parameter defined in [MDLIB/include/mdoodz.h:78](MDLIB/include/mdoodz.h#L78). Currently used values are `-1, 0, 1, 2` (augmented direct, direct, KSP, Powell–Hestenes augmented "killer solver"); value `3` is free for a new path. The assembly producing the matrices lives in [MDLIB/StokesAssemblyDecoupled.c](MDLIB/StokesAssemblyDecoupled.c), and the sparse CSC plumbing (fill-reducing ordering, factor/solve wrappers) is in [MDLIB/SparseTools.c](MDLIB/SparseTools.c).

The direct path costs roughly O(N^1.5) time and super-linear memory because of fill-in. Memory, not time, is the practical ceiling today — above ~1000×1000 the factor no longer fits on a single node. Because MDOODZ targets single-node + OpenMP and deliberately avoids MPI/PETSc, switching to geometric multigrid (GMG) is the highest-impact algorithmic lever available. The uniform staggered grid maps cleanly onto 2× coarsening in each direction, so GMG here is algorithmically lighter than for an FEM code with unstructured connectivity.

This design covers the GMG solver itself, its coupling to the existing Picard/Newton nonlinear loop in [MDLIB/StokesRoutines.c](MDLIB/StokesRoutines.c), and the test scaffolding that will protect future changes.

## Goals / Non-Goals

**Goals**
- Add a GMG V-cycle that operates on the existing staggered-grid Stokes stencils, with a Vanka block smoother and a CHOLMOD coarsest-level solve.
- Wrap GMG as a preconditioner for FGMRES so convergence remains robust at the large viscosity contrasts (crust/mantle/air, 1e6–1e10) typical of MDOODZ problems.
- Preserve the existing nonlinear loop (Picard + Newton + line search + safe_mode) unchanged — GMG must be a drop-in replacement at the linear-solve seam.
- Achieve grid-independent convergence (V-cycle reduction factor ρ < 0.15) at constant viscosity, and at least 3× wall-time speedup vs. CHOLMOD at 201×201.
- Make every GMG component independently unit-testable without running a simulation.
- Keep CHOLMOD the default; GMG is opt-in via `lin_solver = 3`.

**Non-Goals**
- 3D support. MDOODZ is currently 2D; a 3D extension is out of scope even though GMG generalises naturally.
- Algebraic multigrid (AMG). The uniform staggered grid makes GMG the natural fit and avoids the setup cost of AMG.
- Replacing CHOLMOD entirely. CHOLMOD is retained as the default solver and as the coarsest-level direct solve inside the GMG hierarchy.
- GPU or MPI parallelism. All GMG work stays single-node + OpenMP, matching the rest of MDOODZ.
- Rewriting `StokesAssemblyDecoupled.c`. The existing fine-level stencil assembly is reused; only the coefficients are exposed for Vanka block construction.
- Changing the Newton Jacobian formulation. GMG will be applied to the same Jacobian the outer Newton loop is already building.

## Decisions

### D1: Vanka block smoother as the relaxation scheme

**Choice**: Use a cell-based Vanka block smoother (4 velocities + 1 pressure per block, solved by direct 5×5 inversion per cell each sweep).

**Why**: Stokes is a saddle-point system, not SPD, so pointwise Gauss–Seidel diverges. Vanka is the standard smoother for staggered-grid GMG (StagYY, LaMEM's geometric path) and maps directly onto MDOODZ's stencils: the four velocities and central pressure around any cell are already what `StokesAssemblyDecoupled.c` couples together.

**Alternatives considered**
- *Uzawa / pressure-correction*: simpler but significantly slower under high-contrast viscosity and poor for the heavily variable-viscosity regime MDOODZ targets.
- *Braess–Sarazin*: a symmetric Uzawa variant. Better than plain Uzawa but still loses to Vanka when η contrasts span many decades.

Vanka is chosen for its robustness at the viscosity contrasts MDOODZ regularly handles; the cost of the 5×5 block solve is trivial relative to its convergence benefit.

### D2: GMG as an FGMRES preconditioner, not a standalone solver

**Choice**: Run the V-cycle inside an outer FGMRES iteration. `lin_solver = 3` invokes `FGMRES(Krylov) + GMG(precond)`, not a bare V-cycle.

**Why**: Constant-viscosity GMG converges in 5–10 V-cycles and doesn't need a Krylov outer loop. MDOODZ's real workloads have viscosity contrasts ≥ 1e6, where standard GMG degrades because the coarse grid cannot represent sharp viscosity jumps and restriction smears the viscosity field. Wrapping the imperfect V-cycle as a preconditioner for FGMRES lets the Krylov iteration "clean up" spectrally-poor modes and keeps convergence robust without requiring a heroic smoother or specialised transfer operators. This is the approach LaMEM uses, and is the most resilient option for a code that routinely sees 1e8+ contrasts.

**Alternatives considered**
- *Standalone V-cycle*: simplest, fastest per iteration, but fragile under large η contrasts.
- *Defect-correction iteration*: viable but typically slower than FGMRES+GMG and gives less observability (FGMRES residual history is a familiar diagnostic).

A `gmg_standalone` toggle will be supported for diagnostics only — useful for the `ConvergenceFactorGridIndependent` test which needs a clean V-cycle residual history.

### D3: Transfer operators — full-weighting restriction, bilinear prolongation, with staggering-aware variants

**Choice**:
- Velocity (Vx, Vz): full-weighting restriction (¼/8 + edges ⅛/4 + corners 1/16), bilinear prolongation.
- Pressure (cell-centred): full-weighting restriction, piecewise-constant (injection) prolongation.
- Viscosity field: harmonic-averaged restriction in regions of high contrast, arithmetic elsewhere.

**Why**: The staggered layout means Vx, Vz, and P live at different positions within a cell, so the stencils for restriction/prolongation differ per variable. Full-weighting is standard for GMG and satisfies the R = cP^T relationship required for symmetric V-cycles. Injection for pressure prolongation is standard for staggered-grid GMG since coarse pressure is natively defined at the coarse cell centre.

Harmonic-averaged viscosity restriction is critical: arithmetic averaging across a stiff/soft interface smears the contrast and wrecks coarse-grid operator quality. Harmonic averaging preserves the effective stiffness of the stiff material, which is the physically meaningful coarsening.

**Alternatives considered**
- *Injection for velocity restriction*: cheaper but loses smoothing information; rejected.
- *Operator-dependent transfer operators (Dendy-style)*: best for very high contrasts but far more complex to implement and debug; rejected for the first iteration — can be added later if needed.

### D4: Fine-level fidelity and coarse-level re-discretisation (two different rules)

The word "re-discretise" has to do different jobs at different levels. Splitting the rule:

**Fine level — identical to what CHOLMOD sees.** The fine-level GMG operator SHALL be bit-identical to the matrix MDOODZ currently assembles in [StokesAssemblyDecoupled.c](MDLIB/StokesAssemblyDecoupled.c) under the same `model` settings. That means: same Newton Jacobian coefficients (`D11`–`D34`) when Newton is active, same compressibility terms when `lin_solver = -1` semantics apply, same anisotropy-rotated coefficients, same BC-tag handling including tag `30` deactivation, same free-surface stencil adjustments, same penalty-augmentation under Powell–Hestenes. The way this rule is enforced concretely is covered by D11 below. A textbook Picard stencil *does not satisfy this rule* — it matches only on a narrow subset of configurations and drifts silently otherwise.

**Coarse levels — re-discretised Picard on restricted viscosity.** At every level below the fine one, the operator SHALL be a Picard Stokes stencil computed from the level-restricted viscosity field (harmonic where contrasts are high, arithmetic elsewhere, per D3). Coarse levels do not carry the Newton Jacobian, anisotropy rotation, compressibility, or full BC-tag logic. This is standard for staggered-grid GMG and is what StagYY and LaMEM's geometric path do.

**Why the asymmetry**: the fine-level operator defines the problem FGMRES is trying to solve, so any deviation there produces an answer to the wrong question. Coarse levels are preconditioners — their only job is to approximate the fine operator well enough that FGMRES converges quickly; small discrepancies there just cost extra FGMRES iterations rather than solving a different problem.

**Why not Galerkin for coarse levels**: `A_coarse = R A_fine P` is algorithmically cleaner and automatically consistent, but it widens the 5-point Stokes stencil to 13 points and breaks the Vanka 5×5 block structure. Memory cost grows 5–13×. Re-discretisation keeps the stencil banded and the Vanka block extraction simple.

**Alternatives considered**
- *Use re-discretised Picard at the fine level too* (the initial implementation): matches only on isotropic, Picard-mode, incompressible, closed-boundary problems. Fails on SolVi by ~15× L2 because it doesn't know about Newton, anisotropy, BC rewrites, or the penalty augmentation.
- *Galerkin for all levels*: see above — widened stencils, broken Vanka.
- *Galerkin for the top 1–2 levels only, re-discretisation below*: a reasonable hybrid if pure re-discretisation produces a poor coarse-grid preconditioner. Defer until we measure convergence degradation — the FGMRES outer iteration should absorb modest coarse-grid inconsistency.

### D5: Coarsest-level solve via the existing CHOLMOD path

**Choice**: Continue coarsening until the coarsest level has ≤ ~5000 DOFs (e.g. 63×42), then hand the level's `SparseMat` to the existing CHOLMOD direct solver.

**Why**: CHOLMOD is already in the build, already has the right input format (`SparseMat` CSC), and is extremely fast on small systems. No reason to reinvent this.

**Alternatives considered**
- *Recurse down to a single-cell system*: would require a specialised tiny-system solver and is pointless when CHOLMOD is already available.
- *Iterative solve at the coarsest level*: slower and defeats the O(N) property of the V-cycle.

### D6: Level count policy — geometric coarsening with a DOF floor

**Choice**: `n_levels = max(2, floor(log2(min(Nx,Nz) / 16)))`. Continue halving until either the minimum dimension falls below ~16 or the total DOF count drops below ~5000, whichever comes first. User can override via `gmg_levels`.

**Why**: ~16 cells per side is the smallest useful smoothing domain; below that the smoother degenerates into the coarse solve. The 5000-DOF floor matches the regime where CHOLMOD clearly beats multigrid and keeps the coarse solve under a few percent of total time.

**Alternatives considered**
- *Fixed `n_levels = 5`*: simple but wasteful at small grids and insufficient at huge grids.
- *Continue until the coarsest grid is 2×2*: smoothing at tiny grids is unnecessary overhead; better to bail out to CHOLMOD earlier.

### D7: Newton Jacobian — restrict the full Jacobian, smooth with the symmetric part

**Context — how Picard/Newton and GMG layer.** Picard and Newton are **not** alternative solvers to GMG; they are the outer nonlinear loop that produces the linear system Ax = b at every iteration. Picard produces a symmetric A (just current-iterate viscosity in the stencils). Newton produces a non-symmetric A that includes the full Jacobian coefficients D11–D34 assembled in [StokesAssemblyDecoupled.c](MDLIB/StokesAssemblyDecoupled.c). The linear solver (today CHOLMOD, after this change GMG-FGMRES when `lin_solver = 3`) is called from inside the outer loop at every Picard/Newton iteration and does not know which of the two produced A. The outer loop keeps running in either case.

**Target behaviour (eventual)**: When the outer loop is in Newton mode, the V-cycle SHALL use the full Jacobian when evaluating the residual that feeds FGMRES, and the symmetric (Picard) part of the operator when assembling each Vanka 5×5 block for the smoother. This is not a corner case — it is the always-on behaviour of GMG-FGMRES once Newton-mode GMG is fully implemented.

**Phased scope for this change**:

- **Phase 1 (Picard-mode GMG, first deliverable under this change)**: the full matvec accessor `ApplyStokesOperatorMDOODZ` (see D11) reproduces MDOODZ's Picard stencil exactly. When the outer loop is in Newton mode, the GMG dispatch SHALL detect this and return a "not applicable" status (RC `-3`) so that `StokesRoutines.c` falls back transparently to CHOLMOD. This interim behaviour is safe: no wrong answers, just no memory win under Newton.
- **Phase 2 (Newton-mode GMG, completion of this change)**: add the Jacobian variant of `ApplyStokesOperatorMDOODZ` that reproduces the full `D11`–`D34` stencil. The V-cycle then applies the symmetric-part-smoothing rule from the target behaviour paragraph above. At this point the Newton-mode guard is removed and GMG handles both Picard and Newton modes.

**Why phased**: Phase 1 is sufficient to deliver most of the memory-saving value (many MDOODZ runs spend most of their non-linear iterations in Picard before switching to Newton near convergence). Phase 2 is mechanically similar but requires auditing a larger set of Jacobian-assembly functions; splitting the work keeps each reviewable diff small.

**Why this layering in general**: Vanka block inversion is trivial for the symmetric Picard operator (5×5 SPD-ish block) but can be ill-conditioned for the non-symmetric Newton Jacobian, especially in strongly localising regions where the stress-rotation terms are large. Using the symmetric part inside the smoother keeps every Vanka sweep stable; using the full Jacobian for the outer FGMRES residual preserves Newton's quadratic convergence rate. This is standard in geodynamic GMG codes (LaMEM, StagYY).

**Alternatives considered**
- *Smooth with the full Jacobian*: risks smoother divergence at strong stress rotation.
- *Permanently fall back to Picard under Newton*: gives up Newton's convergence rate and the memory win on Newton iterations.
- *Apply D7 only when the Jacobian is "too asymmetric"*: no clean threshold, hard to tune, brittle. Always-on (in Phase 2) is simpler.
- *Skip Phase 1 and do everything at once*: possible but makes the reviewable diff much larger and couples correctness of Picard-mode GMG to correctness of Newton-mode Jacobian auditing, which are independent concerns.

### D8: Test taxonomy — components isolated from `RunMDOODZ`

**Choice**: Three test tiers:
1. **Component units** (no simulation, no HDF5): restriction full-weighting property, bilinear prolongation exactness on linear fields, PR=I identity on representable fields, Vanka single-block analytical solve, Vanka smoothing property (‖r₁‖ < ‖r₀‖), grid-hierarchy sizing for even/odd dimensions.
2. **Solver-level** (hand-built Ax=b, no simulation): constant-viscosity GMG-vs-CHOLMOD equivalence to 1e-10, grid-independent ρ < 0.15 across 21² → 161², and `ViscosityContrastScaling` documenting V-cycles needed at η contrasts 1, 1e2, 1e4, 1e6.
3. **Integration** (existing suite, parameterised over solver): run `BlankenBench.ConvectionDevelops`, `VelocityField.PureShear`, `TopoBench.Relaxation`, `NewtonIterationConvergence.*`, SolVi under both CHOLMOD and GMG and assert L2(Vx_gmg - Vx_chol)/L2(Vx_chol) < 1e-8 plus a similar bound for Vz and P.

**Why**: The component and solver tiers catch bugs during development without the cost (or opacity) of running a full simulation. The `ConvergenceFactorGridIndependent` test is the single most diagnostic check — a constant ρ across grid sizes is the mathematical signature of a correct implementation. The integration tier protects against subtle coupling bugs with Newton, free surface, and thermal feedback that the pure linear-algebra tests can't see.

### D9: Parameters exposed to `.txt` input

**Choice**: New parameters, all optional with defaults:
- `lin_solver = 3` → select the GMG-preconditioned FGMRES path.
- `gmg_levels` (default: auto from D6).
- `gmg_nu_pre`, `gmg_nu_post` (default: 2 pre, 2 post Vanka sweeps).
- `gmg_fgmres_restart` (default: 30), `gmg_fgmres_tol` (default: `nonlin_abs_mom` / 10).
- `gmg_standalone` (default: 0 = FGMRES-wrapped; 1 = bare V-cycle for diagnostics).

**Why**: Follows the existing pattern for `max_its_KSP`, `max_its_PH`, `diag_scaling`, `preconditioner` in [mdoodz.h:78](MDLIB/include/mdoodz.h#L78). All GMG knobs are optional with tuned defaults so users never need to touch them for ordinary runs.

### D10: OpenMP strategy — reuse the existing domain decomposition, red-black Vanka

**Context.** MDOODZ is single-node + OpenMP. The existing Stokes assembly uses a hand-rolled domain-decomposition pattern with `estart[ith]` / `eend[ith]` per-thread equation ranges, merged at the end via `MergeParallelMatrix`. See [StokesAssemblyDecoupled.c:1173, 1409, 1467, 1525](MDLIB/StokesAssemblyDecoupled.c#L1173). This pattern is already proven to scale and it avoids atomics by giving each thread a disjoint range of rows.

**Choice** — per GMG component:
- **Assembly at every level**: reuse the existing `AllocateDecompositionDecoupled` / `DomainDecompositionDecoupled` / `MergeParallelMatrix` pattern with per-level equation counts. No new parallelisation code.
- **Vanka smoother**: red–black cell colouring. Two passes per sweep; cells of the same colour share no velocity unknowns, so threads can update them in parallel without locks or atomics. Avoids ordering-dependence that plagues naive Gauss–Seidel-style parallelisation.
- **Matrix-vector product (for residual)**: `#pragma omp parallel for` over equation rows; no races because each row is written by exactly one thread.
- **Transfer operators (R, P, viscosity restriction)**: trivial parallel loops over coarse cells.
- **CHOLMOD at the coarsest level**: use CHOLMOD's own internal threading with whatever `OMP_NUM_THREADS` is set. At ≤ 5000 DOFs this is minor overall.
- **FGMRES**: Arnoldi orthogonalisation has a serial inner-product bottleneck. Parallelise vector updates (`x += α·p`) and matvecs; accept that orthogonalisation dominates at high thread counts.

**Expected scaling**: 16–32 threads effective on typical MDOODZ problem sizes, plateauing from there because FGMRES orthogonalisation is the bottleneck (Amdahl's law on a ~5% serial fraction gives a ~20× ceiling). This is similar to the current CHOLMOD path's scaling. No anti-scaling expected if red-black ordering is correctly applied.

**Alternatives considered**
- *Atomic updates over shared DOFs in Vanka*: simpler code, but modern x86 atomic contention scales poorly past 8 threads under heavy write traffic. Rejected.
- *Task-based parallelism (OpenMP tasks)*: good for unbalanced workloads, but GMG is highly regular; the overhead of task creation outweighs the benefit. Rejected.

**Non-goals**: OpenMP offload to GPU, explicit SIMD intrinsics in the smoother. Performance tuning of SIMD can come later; for now rely on the compiler's auto-vectorisation.

### D11: The stencil bridge — how the fine-level operator matches `StokesAssemblyDecoupled.c`

**Context.** D4 requires that the fine-level GMG operator be bit-identical to what MDOODZ's existing assembly produces. That rule needs a concrete mechanism. The question is how a matrix-vector product routine called from inside the V-cycle can share its coefficient formulas with the sparse-matrix assembly in [StokesAssemblyDecoupled.c](MDLIB/StokesAssemblyDecoupled.c) without either (a) rederiving the stencil from textbook and drifting from the MDOODZ version or (b) modifying hot-path assembly code that four other `lin_solver` paths rely on.

**Choice — Option A: additive duplication into a read-only matvec, with a golden cross-check test.**

A new translation unit `MDLIB/StokesAssemblyGMG.{c,h}` SHALL provide pure read-only accessors:

```c
void ApplyStokesOperatorMDOODZ(
    const grid *mesh, params model,
    const double *u, const double *v, const double *p,
    double *ru, double *rv, double *rp);

void StokesCellBlockMDOODZ(
    const grid *mesh, params model, int i, int j,
    double out_A[5][5], double out_b[5]);
```

These reproduce exactly the coefficient formulas in the existing `Xmomentum_InnerNodesDecoupled`, `Zmomentum_InnerNodesDecoupled`, `Continuity_InnerNodesDecoupled` (Picard, Phase 1) and `Xjacobian_InnerNodesDecoupled3`, `Zjacobian_InnerNodesDecoupled3` (Newton, Phase 2). `StokesAssemblyDecoupled.c` is **not modified**. `MultigridStokes.c::StokesApplyA` and `VankaBlockAssembleSolve` at the fine level SHALL call the new accessors. Coarse levels continue to re-discretise Picard Stokes on the restricted viscosity (D4).

Duplication is a real cost. It is paid for with an explicit correctness contract: a **golden cross-check test** that assembles the full sparse matrix `A` via the existing path, multiplies `A·x` for ≥ 10 random `x`, compares against `ApplyStokesOperatorMDOODZ(mesh, …, x)` to 1e-12, and fails loudly on any drift. This is the same pattern used in other scientific codebases (deal.II, PETSc) where hand-derived matvecs coexist with assembled-matrix paths. Any future change to a coefficient formula in `StokesAssemblyDecoupled.c` without the corresponding update to `StokesAssemblyGMG.c` is caught by this test before it can ship.

**Why Option A rather than B, C, or D (discussed fully in [STATUS.md §6](STATUS.md#6-possible-solutions))**:

- *Option B: refactor `StokesAssemblyDecoupled.c` to share coefficient code.* Single source of truth, no duplication. But: touches the hottest existing code path, diff is large and scattered across five `lin_solver` regression suites. Regression risk medium.
- *Option C: flag-driven matvec mode on existing builders.* Smaller diff than B, but every per-coefficient branch must be audited for the `matvec` vs `assemble` split. Partial duplication anyway.
- *Option D: freeze GMG as Picard-only, Newton stays on CHOLMOD forever.* Zero additional work but permanently weaker: loses Newton's quadratic convergence rate and the memory win on Newton iterations. Accuracy tests never enable.

Option A is picked because its regression risk to `lin_solver = 0, 1, 2, -1` is exactly zero (no existing function is touched), and the duplication cost is bounded by the golden cross-check test — which also serves as regression protection on the GMG side any time a coefficient formula evolves.

**Alternatives considered**
- See options B, C, D above. Also considered: *pre-assembling the full sparse matrix at every GMG level and doing a sparse matvec*. Rejected because it defeats the memory advantage of GMG — sparse storage is ~30× smaller than CHOLMOD fill-in but still dominates if kept per level.

**Mechanical invariant to preserve**: every coefficient formula in `StokesAssemblyDecoupled.c` that can be reached by any production `lin_solver` path must have a mirror in `StokesAssemblyGMG.c`, enforced by the golden test. The test SHALL exercise every BC tag that appears in production runs (0, −1, 2, 11, 13, 30, 31, −2, −12), every rheology branch (Picard, Newton), every anisotropy state (off, on), and a representative compressibility setting.

## Risks / Trade-offs

- **Convergence degradation at η contrasts > 1e6** → Mitigated by running GMG inside FGMRES (D2) and harmonic viscosity restriction (D3). The `ViscosityContrastScaling` test documents the envelope so regressions surface immediately.
- **Newton Jacobian asymmetry breaks the Vanka block solve** → Mitigated by smoothing with the symmetric (Picard) part while keeping the full Jacobian in the residual (D7). Validated by `NewtonIterationConvergence.*` integration tests running under both CHOLMOD and GMG.
- **Odd or prime `Nx`, `Nz` don't coarsen cleanly by factor 2** → Mitigated by per-level dimensions `ceil(N/2)` and handling of non-power-of-2 sizes in the grid hierarchy constructor, covered by `GridHierarchyOddSizes` test. Worst case the bottom level is slightly larger than the ideal target.
- **Irregular active domain under the free-surface marker chain** → MDOODZ does not use sticky air; it deactivates cells above the topographic marker chain by setting `BCu.type`, `BCv.type`, `BCg.type`, `BCp.type`, `BCt.type` to `30` in [FreeSurface.c:609-690](MDLIB/FreeSurface.c#L609-L690), and the assembly in [StokesAssemblyDecoupled.c:113-114, 610-611](MDLIB/StokesAssemblyDecoupled.c#L113-L114) zeroes stencil contributions from tag-`30` neighbours. Every GMG level must therefore carry a consistent active/inactive mask that agrees with the fine-level mask under coarsening. Mitigated by a dedicated `ActiveMaskRestrict` routine that marks a coarse cell inactive iff *all* of its four fine children are inactive, otherwise active; deactivated cells are skipped entirely by the Vanka sweep, the residual evaluation, and the transfer operators. Validated by `TopoBench.Relaxation` under both solvers. Risk is real but well-understood (embedded-boundary multigrid is a well-studied class).
- **Memory of the grid hierarchy is non-trivial at deep levels** → Each level roughly costs ¼ of the finer level, so the full hierarchy costs ~4/3 × fine-level arrays, far below the CHOLMOD factor's fill-in. Accepted.
- **Integration-test runtime doubles because each test runs under both solvers** → Mitigated by running the dual-solver check only on the representative subset in §D8, not the full 47-test suite. The rest continue to run under the default CHOLMOD path.
- **Vanka block inversion is `O(cells)` per sweep and could become a hotspot** → Mitigated by precomputing LU factors of the 5×5 blocks once per V-cycle (viscosity-dependent but not RHS-dependent). OpenMP-parallel over cells with red-black colouring to avoid race conditions.

## Migration Plan

1. Land the scaffolding under a disabled `lin_solver = 3` path first: new files compile, dispatcher recognises the value, but it currently returns a stub error. Existing solvers unchanged.
2. Land grid hierarchy + transfer operators + Vanka smoother with component tests. Still no integration with the nonlinear loop.
3. Land the V-cycle driver and the constant-viscosity equivalence test.
4. Land the FGMRES wrapper and the full `lin_solver = 3` dispatch.
5. Enable `lin_solver = 3` in the dual-solver integration tests; CI then protects both paths.
6. Document in `skill-solvers` via a separate docs change after merge.

**Rollback**: Purely additive — setting `lin_solver` back to `0, 1, 2` restores the previous behaviour bit-for-bit. No data migration, no schema changes, no dependency changes.

## Open Questions

- **Vanka block colouring vs. atomic updates under OpenMP**: red-black cell colouring is the planned default, but a comparison with atomic-updates-on-shared-dofs would tell us whether the extra colouring cost is worth it on modern CPUs. Defer to a benchmarking task during implementation.
- **Whether to ship `gmg_standalone` as a public parameter or as an internal testing hook only**. Leaning public for diagnostics but could hide it behind a `-DDEV` build flag.
- **Does the active-mask restriction lose accuracy when the marker chain is at a shallow angle?** A coarse cell straddling the topographic boundary may have 1, 2, or 3 of its four fine children active. The proposed rule "coarse active iff any fine child is active" is conservative (keeps the coarse grid larger than the fine domain), but it may produce coarse cells where the viscosity is a mix of true rock and "not really rock" cells whose values are stale. A more accurate alternative is to mark the coarse cell inactive if *fewer than half* of its children are active, but that shrinks the coarse domain and can lose low-frequency error near the surface. `TopoBench.Relaxation` under GMG will indicate which rule gives better convergence; we'll pick after measurement.
- **Coarse-level CHOLMOD factorisation caching across Newton iterations**: the coarsest viscosity field changes slowly; re-using a previous symbolic factorisation would save setup cost. Prototype during performance tuning.
