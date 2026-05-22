## Why

MDOODZ currently solves the saddle-point Stokes system exclusively with a CHOLMOD/UMFPACK direct factorisation. This scales as roughly O(N^1.5) in 2D and, more importantly, its memory footprint grows super-linearly because of fill-in. On the single-node + OpenMP platform MDOODZ targets, this is the dominant cost per nonlinear iteration and is the effective ceiling on resolution: simulations beyond ~1000×1000 hit a memory wall long before they hit a time-to-solution wall.

The MDOODZ discretisation — a uniform staggered (MAC) grid with constant dx/dz — is the ideal data structure for geometric multigrid (GMG). Coarsening by exactly 2× in each direction preserves the staggered variable layout, no unstructured mesh machinery is needed, and the momentum/continuity stencils already computed in `StokesAssemblyDecoupled.c` can be reused directly for a Vanka block smoother. Adding GMG therefore unlocks O(N) algorithmic scaling and drops memory consumption by a large factor, without introducing new external dependencies (PETSc, MPI, GPU).

## What Changes

- **NEW**: Geometric multigrid V-cycle for the decoupled Stokes system on the staggered grid (restriction, prolongation, Vanka block smoother, grid hierarchy).
- **NEW**: FGMRES Krylov outer iteration using the GMG V-cycle as a preconditioner, to keep convergence robust under large viscosity contrasts (>1e6) typical of crust/mantle/air problems.
- **NEW**: Coarsest-level direct solve via the existing CHOLMOD path, applied to a ~few-thousand-DOF problem at the bottom of the hierarchy.
- **NEW**: Unit tests for each GMG component (restriction full-weighting, bilinear prolongation, Vanka single-block solve, Vanka smoothing property, grid hierarchy sizing for even/odd `Nx`,`Nz`) that run without `RunMDOODZ` or HDF5 I/O.
- **NEW**: Solver-level tests: GMG-vs-CHOLMOD equivalence on a constant-viscosity manufactured system, grid-independent convergence factor (ρ < 0.15 across 21² → 161²), and documented convergence degradation under checkerboard viscosity contrasts up to 1e6.
- **NEW**: Performance regression tests: memory scaling exponent α < 1.2 and GMG ≥ 3× faster than CHOLMOD at 201×201.
- **MODIFIED**: `lin_solver` parameter gains a new value selecting the GMG-preconditioned FGMRES path; existing CHOLMOD/UMFPACK values are unchanged and remain the default.
- **MODIFIED**: Integration tests from the existing suite are parameterised to run under both CHOLMOD and GMG, asserting L2 equivalence of Vx/Vz/P within tolerance on representative benchmarks (SolVi, Blankenbach convection, pure shear, topographic relaxation, Newton-iteration convergence tests).

This change is **non-breaking**: GMG is opt-in via `lin_solver`, CHOLMOD remains the default, and the existing Stokes assembly code is reused rather than replaced.

## Capabilities

### New Capabilities

- `gmg-stokes-solver`: Geometric multigrid V-cycle with Vanka block smoother on the MDOODZ staggered grid, wrapped as an FGMRES preconditioner. Includes grid hierarchy construction with viscosity restriction, component operators (R, P, smoother, coarse solve), V-cycle driver, and integration with the existing `lin_solver` selection path.

### Modified Capabilities

None. The existing Stokes assembly, Newton–Picard iteration, penalty method, and CHOLMOD direct path retain their current requirements. GMG is an additive solver option selected at runtime, not a replacement. The `skill-solvers` documentation spec may be refreshed later via a separate docs change once the implementation lands.

## Impact

**New files**
- `MDLIB/MultigridStokes.c` / `MultigridStokes.h` — V-cycle driver, restriction, prolongation, Vanka smoother, FGMRES wrapper.
- `MDLIB/MultigridLevels.c` / `MultigridLevels.h` — grid hierarchy allocation, viscosity restriction, per-level stencil coefficients, memory management.
- `TESTS/` additions for GMG component, solver-level, and performance tests (GTest, no simulation driver required for most).

**Modified files**
- `MDLIB/Solvers.c` — dispatch new `lin_solver` value to the GMG-preconditioned FGMRES path.
- `MDLIB/StokesAssemblyDecoupled.c` — expose stencil coefficients (Vx/Vz momentum + continuity) per cell so they can be reused by the Vanka block assembly at every level without re-deriving them.
- `MDLIB/input.c` (or wherever `lin_solver` is parsed) — accept the new solver value, with per-level tuning parameters (`gmg_levels`, `gmg_nu_pre`, `gmg_nu_post`, FGMRES restart/tolerance).
- `TESTS/CMakeLists.txt` — register new test binaries; parameterise existing integration tests over solver choice.

**APIs / user-facing**
- One new `lin_solver` value and a small set of optional `gmg_*` parameters in the `.txt` input file.
- No change to HDF5 output schema, BC callbacks, or rheology interfaces.

**Dependencies**
- None added. The coarsest-level direct solve reuses the existing SuiteSparse / CHOLMOD build. No PETSc, no MPI, no GPU.

**Risks**
- Convergence degradation at extreme viscosity contrasts (>1e6). Mitigated by running GMG as an FGMRES preconditioner, not a standalone solver, and by viscosity-weighted restriction. The `ViscosityContrastScaling` test documents the degradation envelope so regressions are detectable.
- Interaction with the Newton Jacobian (`D11`–`D34` coefficients) at the fine level. Mitigated by restricting the full Jacobian stencil at each level and validating via the `NewtonIterationConvergence` integration tests running under both solvers.
