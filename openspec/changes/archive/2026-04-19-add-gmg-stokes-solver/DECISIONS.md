# `add-gmg-stokes-solver` — decision registry

Companion to `RETROSPECTIVE.md`. This file groups the decisions made
during the change by **level of scope** (strategic → architectural →
tactical → tooling), each with:

- **When** — which session / user message it was made in.
- **Alternatives** considered and why they lost.
- **Outcome** — measurable consequence (which test turned green, which
  bug surfaced, which follow-up got unblocked).

`design.md` is the canonical home for the architectural-level
decisions labelled **D1…D11**. This file cross-references those where
relevant and adds the strategic and tactical decisions that never
needed to land in the design doc.

---

## Strategic decisions (scope, phasing, review surface)

### S1. Land the GMG core behind `lin_solver = 3` with transparent CHOLMOD fall-back

- **When**: Session 1, derived from the attached `design.md` and
  confirmed in Session 2.
- **Alternative**: ship GMG as a replacement for the existing
  `lin_solver` paths.
- **Why chosen**: zero-regression guarantee for `lin_solver = 0, 1, 2,
  -1`; lets the GMG core be merged as its own reviewable diff before
  the stencil-bridge work. Users who don't opt in never see the new
  path.
- **Outcome**: every existing CI suite stayed green through every
  session; the GMG core and the stencil bridge were both reviewable as
  independent deltas.

### S2. Split Section 15 "Stencil bridge" into Phase 1 (Picard) and Phase 2 (Newton)

- **When**: user@391, when the user merged Section 15 into `tasks.md`.
- **Alternative**: land the full bit-identical bridge in one go.
- **Why chosen**: the Picard path is mathematically simpler and
  covers most MDOODZ workflows; it can be golden-tested independently
  of the Newton Jacobian. Phase 2 then reuses Phase 1's plumbing.
- **Outcome**: after Session 6b (Phase 1) the Picard SolVi
  dual-solver test hit `5.79e-11` — three decades under the spec bound
  — *before* any Newton code was written.

### S3. Defer performance / memory regression harness (§12) to a follow-up change

- **When**: Session 2, then re-confirmed in Session 7 (task 15.26).
- **Alternative**: land a 201×201 timing and memory test in this
  change.
- **Why chosen**: performance numbers are only meaningful after the
  full stencil bridge lands; otherwise the measurements would be
  CHOLMOD via the fall-back. Formal harness belongs to a dedicated
  "GMG performance and diagnostics" change with its own reviewable
  baseline.
- **Outcome**: `tasks.md` §12.1-§12.3 and §15.26 both carry an
  explicit deferral note with rationale, not an abandoned `[ ]`.

### S4. Defer integration-fixture wiring for `DISABLED_GmgStokesEquivalence.*` benchmarks

- **When**: Session 7 (task 15.25).
- **Alternative**: wire `PowerLawShearZoneL2Match`, `SolCxL2Match`,
  `FreeSurfaceTopography` into live integration tests before archiving.
- **Why chosen**: each needs a bespoke `.txt` twin + reference field;
  operator correctness for those benchmarks is *already* locked down
  at `1e-12` by the golden tests in `TESTS/StokesMatvecEquivalence.cpp`
  (constant-viscosity, D12/D21 coupling, fully anisotropic
  D13/D14/D31-D34, Newton Vx+Vz, free-surface stabilisation, D7
  cell-block invariant).
- **Outcome**: `tasks.md` §15.16 and §15.25 both cite the golden tests
  as sufficient for the correctness contract; fixture wiring tracked
  as a "GMG benchmark integration" follow-up.

### S5. Retrospective + decision registry as archive artefacts

- **When**: user@627 (this session).
- **Alternative**: archive with just `STATUS.md` describing the final
  state.
- **Why chosen**: user asked explicitly for a chronological
  retrospective capturing decisions, problems, and effort.
- **Outcome**: `RETROSPECTIVE.md` + this file (`DECISIONS.md`) both
  live under `openspec/changes/add-gmg-stokes-solver/` so they travel
  with the change to the archive.

---

## Architectural decisions (mirror of design.md — abbreviated)

These are the decisions formalised in `design.md` as **D1…D11**. Quick
reference only — see the design doc for the full rationale.

### D1 — UMFPACK (not CHOLMOD) at the coarse level

- **When**: Session 2, after CHOLMOD was tried first and rejected.
- **Why**: Stokes saddle-point is indefinite; CHOLMOD is SPD-only.
  UMFPACK is a general-sparse LU factorizer and lives in the same
  SuiteSparse bundle MDOODZ already links against.
- **Outcome**: `GMG_Coarse.UmfpackSolveMatchesExactSolution` at 1e-10.

### D2 — Vanka 5×5 block smoother (not symmetric Gauss-Seidel)

- **When**: Session 2, Section 4.
- **Why**: SGS fails on Stokes saddle-point because the momentum and
  continuity blocks don't commute; Vanka solves the local saddle-point
  exactly per cell.
- **Outcome**: `GMG_Vanka.ExactSolutionIsFixedPoint` and
  `SweepReducesResidual` green once ω=0.6 damping was added.

### D3 — FGMRES (not CG) as the outer driver

- **When**: Session 2, Section 8.
- **Why**: `A_Stokes` is indefinite; CG doesn't apply. FGMRES accepts
  a non-symmetric preconditioner — essential because our V-cycle is
  not symmetric.
- **Outcome**: `GMG_FGMRES.ConvergesOnConstantViscosity` at 1e-6 in 30
  iters on 33×33.

### D4 — Coarse-level *rediscretisation* on restricted viscosity (not Galerkin RAP)

- **When**: Session 5, formalised by user@386.
- **Alternative**: Galerkin RAP (coarse operator = R·A·P) for a
  theoretically cleaner coarse problem.
- **Why chosen**: RAP widens the 5-point Stokes stencil to 13-point
  and introduces cross-component coupling that breaks the 5×5 Vanka
  block structure, with 5×–13× memory cost. Rediscretisation is
  textbook for staggered-grid Stokes and used by LaMEM, StagYY,
  Aspect, and Underworld.
- **Outcome**: coarse levels share the Picard stencil code path with
  the fine level; Vanka block structure preserved everywhere.

### D5 — Pressure null-space projection (not a fixed-pressure gauge)

- **When**: Session 2 bug-hunt #5 in `RETROSPECTIVE.md`.
- **Why**: fine and coarse identity-row gauges live at different
  physical locations, breaking Galerkin consistency across levels.
  Null-space projection is gauge-invariant and compatible with
  multigrid.
- **Outcome**: pressure bound in dual-solver test is 1e-6 after mean
  subtraction (gauge residual); velocity bound is 1e-8 exactly.

### D6 — Self-consistent Picard matvec initially; full MDOODZ stencil via a later bridge

- **When**: Session 2, re-visited in user@381 / user@383.
- **Alternative**: implement the full stencil from day one.
- **Why chosen**: lets the GMG core ship and be reviewed independently
  of `StokesAssemblyDecoupled.c`.
- **Outcome**: GMG core merged and green (25 tests) before the bridge
  started.

### D7 — Symmetric-smoothing rule: Newton in matvec, *symmetric* Picard in the smoother

- **When**: user-supplied in the Phase 2 portion of `design.md` at
  user@391; formalised in task 15.21.
- **Why**: the full Newton Jacobian has off-diagonal cross terms
  (D13/D14 and D31/D32/D34) that destroy the SPD-like structure the
  Vanka block needs for stability. Splitting roles gives both
  quadratic nonlinear convergence (Newton in the matvec → FGMRES
  residual) *and* V-cycle stability (symmetric Picard in the block).
- **Outcome**: `CellBlockNewtonReducesToSymmetricPicardPart`
  at 1e-14 protects against future regressions; Newton SolVi
  dual-solver match at 1.89e-15 / 1.34e-15 / 7.52e-8.

### D8 — Newton-mode guard returning `-3` during Phase 1

- **When**: Session 4 (task 15.17).
- **Alternative**: let Newton silently solve with the Picard matvec.
- **Why chosen**: silently solving the wrong operator is user-hostile;
  a clean `-3` return lets the dispatch layer fall back to CHOLMOD
  transparently.
- **Outcome**: task 15.18 pinned the contract; **removed** in Phase 2
  task 15.23 once Newton-matvec landed.

### D9 — Deferred stencil-bundle accessor on `StokesAssemblyDecoupled.c`

- **When**: Session 3 (task 4.1), re-scoped in Session 5.
- **Alternative**: refactor `StokesAssemblyDecoupled.c` to expose a
  stencil accessor.
- **Why chosen**: Section 15's Option A (additive duplication) does
  *not* need this refactor. The two paths are no longer both needed.
- **Outcome**: §4.1, §4.5, §5.5 carry explicit deferral notes
  referencing design D11; no in-repo refactor of the hot assembly
  path.

### D10 — Active mask: only tag 30 ("air") is inactive

- **When**: Session 3 (task 5.1).
- **Alternative**: treat all non-interior tags as inactive.
- **Why chosen**: Dirichlet / Neumann / periodic / interior all need
  identity rows or boundary handling, not exclusion. Tag 30 is the
  only tag that corresponds to a physically-absent DOF.
- **Outcome**: `GMG_ActiveMask.FromBCTagsRespectsAirRegion` pins the
  behaviour; coarse masks propagate via "coarse active iff any fine
  child active".

### D11 — Stencil bridge via Option A: additive duplication, golden cross-check as contract

- **When**: user@383 (the proposal) → Session 5 (formalisation in
  `design.md`) → Sessions 6-7 (implementation).
- **Alternatives considered**: B (refactor existing assembly —
  touches hot path); C (flag-driven matvec mode inside the assembler
  — partial duplication + regression risk); D (freeze as Picard-only
  — abandons Newton/anisotropy users).
- **Why chosen**: perfect backward compatibility (no existing function
  modified); duplication risk mitigated by the golden cross-check at
  1e-12 and the `MDOODZ_GMG_DRIFT_CANARY` compile flag that proves
  the test actually catches drift.
- **Outcome**: 9 Picard + 3 Newton golden variants + 1 D7 cell-block
  invariant + 2 dual-solver integration tests (Picard and Newton
  SolVi), all green; 0 lines of `StokesAssemblyDecoupled.c` modified.

---

## Tactical decisions (implementation-level, in-session)

### T1. Under-relaxation ω = 0.6 for the Vanka sweep

- **When**: Session 2, bug-hunt #2.
- **Why**: classic Vanka-without-damping instability on saddle-point
  systems. ω = 0.6 is the value from Vanka 1986 and John & Tobiska
  2000.
- **Outcome**: 11×11 problem converges to 1e-7 in ≤ 200 sweeps; the
  V-cycle convergence factor dropped from unbounded to ≤ 0.6.

### T2. Zero boundary rows in the restricted residual

- **When**: Session 2, bug-hunt #3.
- **Why**: Vanka holds boundary DOFs as identity rows, so any non-zero
  coarse RHS there fossilises. Explicitly zero-out on restriction.
- **Outcome**: V-cycle residual went from exploding at 1e+136 to
  contracting at ≤ 0.6 per sweep.

### T3. Sign flip on momentum rows in the fine-level bridge

- **When**: Session 6b (task 15.8a in `tasks.md`).
- **Why**: MDOODZ assembler uses `-∇·σ + ∇p = f`; GMG textbook uses
  `+∇·σ - ∇p = f`. Mixing across fine/coarse makes the correction
  diverge. Negate momentum rows in
  `StokesApplyA_MDOODZ_bridge` / `BuildLevelRHSFromCompressed` /
  `VankaBlockAssembleSolve_MDOODZ_bridge`; continuity is
  sign-invariant so it stays untouched.
- **Outcome**: FGMRES converged in 35 iters to `1e-11` tolerance, no
  fall-back, Picard dual-solver 5.79e-11.

### T4. Compressed-equation-space plumbing for `SolveStokesGMG`

- **When**: Session 6b.
- **Why**: the caller (direct and defect-correction paths) passes
  `rhs_u` / `rhs_p` / `x` in compressed space. The original
  `SolveStokesGMG` was ignoring them and rebuilding from
  `mesh->roger_x`. Added `BuildLevelRHSFromCompressed` +
  `ExtractLevelToCompressed` so GMG is a true drop-in.
- **Outcome**: task 15.14's strict bound `L2(Vx) < 5e-2` flipped from
  red (0.94) to green (2.24e-2).

### T5. `comp = 1` override in Newton-mode `StokesCellBlockMDOODZ`

- **When**: Session 7 (task 15.21).
- **Why**: `ApplyStokesOperatorMDOODZ` forces `comp = 1` in its
  Newton branch (mirroring `BuildJacobianOperatorDecoupled`'s
  unconditional `comp = 1`). The cell block must use the same
  compressibility setting to stay spectrally compatible.
- **Outcome**: `CellBlockNewtonReducesToSymmetricPicardPart` passes at
  1e-14.

### T6. `lin_solver == 3` exemption in `InputOutput.c`'s Newton/anisotropy remap

- **When**: Session 7 (task 15.24).
- **Why**: pre-existing unconditional remap at `InputOutput.c:1277`
  was rewriting `lin_solver = 3 → 2` whenever `Newton == 1` or
  `anisotropy == 1`, silently bypassing GMG. Exempting `== 3` lets
  Phase 2 actually flow through `SolveStokesGMG`.
- **Outcome**: Newton dual-solver test went from *fake* 0.0 (both runs
  silently going through KillerSolver) to a *real* 1.89e-15 /
  1.34e-15 / 7.52e-8 comparison between GMG and KillerSolver.

### T7. Harmonic-mean viscosity restriction above contrast ratio 10

- **When**: Session 2 (task 3.3).
- **Why**: arithmetic mean of η with large contrast smooths out the
  physical interface; harmonic mean preserves it. Textbook threshold
  is ~10.
- **Outcome**: `GMG_ViscRestrict.HarmonicOnLargeContrast` green on a
  1:1e6 checkerboard.

### T8. Red-black Vanka sweep (not Jacobi, not in-order Gauss-Seidel)

- **When**: Session 2 (task 4.4).
- **Why**: in-order GS loses parallelism; Jacobi halves the
  effective convergence rate. Red-black is the standard compromise
  and is OpenMP-friendly for a future 14.* task.
- **Outcome**: `GMG_Vanka.ConvergesOnSmallProblem` at 1e-7 in ≤ 200
  sweeps; ordering is serial-correct and threadable.

---

## Tooling / workflow decisions

### W1. Stub `GmgStokesEquivalence.cpp` so CMake doesn't fail before the integration test is written

- **When**: Session 2.
- **Why**: the test file is registered in `TESTS/CMakeLists.txt` from
  the scaffolding phase; a stub keeps CI green until Section 11.1
  lands.

### W2. `MDOODZ_GMG_DRIFT_CANARY` compile flag as a self-test of the golden test

- **When**: Session 6a (task 15.13).
- **Why**: the whole point of the golden test is catching drift, so
  the test itself needs to be provably non-trivial. Wrapping one
  coefficient in a `1 + 1e-7` perturbation behind a compile flag proves
  the test turns red on a known-drift build.
- **Outcome**: default build green; `-DMDOODZ_GMG_DRIFT_CANARY` build
  red with a pointer to the divergent row. Documented inline.

### W3. Minimal `_OMP_` guards around pre-existing `omp_get_max_threads()` calls

- **When**: Session 1.
- **Why**: pre-existing bug (not introduced by this change) broke the
  `OMP=OFF` build chain the test rig uses. Minimal two-line
  `#ifdef _OMP_` guards at the call sites; *not* a refactor.
- **Outcome**: `OMP=OFF` build succeeds; no change to `OMP=ON`
  behaviour.

### W4. Add include guards to `mdoodz-private.h`

- **When**: Session 3.
- **Why**: missing include guards caused redefinition errors when
  `MultigridStokes.h` + a direct `#include "mdoodz-private.h"` both
  pulled it in. Agnostic, additive fix.

### W5. Mirror `SolViBenchmarkTests` fixture pattern exactly in `GmgStokesEquivalence.cpp`

- **When**: Session 4.
- **Why**: the first-attempt designated-initializer pattern for
  `SetParticles_ff` caused a segfault during material-phase read — the
  crash was *upstream* of the solver and entirely due to fixture
  construction order. Copying the working pattern unblocks in one
  edit.
- **Outcome**: `GmgSolViFixture.Res51GmgPipelineProducesBoundedSolution`
  runs end-to-end.

### W6. Twin `.txt` configs for dual-solver comparison (`SolViRes51_{gmg,chol,newton_gmg,newton_chol}.txt`)

- **When**: Session 6b (Picard) and Session 7 (Newton).
- **Why**: running the same benchmark through two different
  `lin_solver` values in one binary requires two configs and a
  side-by-side HDF5 diff (`computeL2Error` on the output fields). The
  apples-to-apples pairing is `lin_solver = 0` (not `-1` — the
  Powell-Hestenes penalty loop solves a *different* augmented
  system).
- **Outcome**: Picard fixture `SolViGmgMatchesCholmodWithin1e8`
  (5.79e-11 on Vx/Vz) and Newton fixture
  `NewtonSolViGmgMatchesCholmodWithin1e8` (1.89e-15 on Vx/Vz). Both
  now live in-repo as the canonical dual-solver contract tests.
