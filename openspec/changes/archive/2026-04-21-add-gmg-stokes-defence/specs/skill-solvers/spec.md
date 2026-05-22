## MODIFIED Requirements

### Requirement: Skill documents the sparse linear algebra backend
The skill SHALL explain that MDOODZ uses SuiteSparse (CHOLMOD + UMFPACK) for direct sparse matrix solves **and**, since `add-gmg-stokes-solver`, an additive geometric-multigrid-preconditioned FGMRES path under `lin_solver = 3`. The `DirectSolver` struct wraps CHOLMOD factorisation used by `lin_solver ∈ {0, 1, 2, -1}` and by the coarsest level of the GMG hierarchy (UMFPACK in practice — Stokes is saddle-point, SPD-only CHOLMOD is used only for the PH-augmented velocity block). The `lin_solver` parameter selects the solver variant. The skill SHALL mention that `diag_scaling` and `preconditioner` parameters can tune direct-solver performance, and that the `gmg_*` family (`gmg_levels`, `gmg_nu_pre`, `gmg_nu_post`, `gmg_fgmres_restart`, `gmg_fgmres_tol`, `gmg_standalone`, `gmg_dump_vcycle`) tunes the multigrid path. The skill SHALL note that `lin_solver = 3` falls back transparently to CHOLMOD when FGMRES fails to converge or when configurations outside its supported envelope (today: none — Phase 2 enabled Newton and anisotropy) are requested.

#### Scenario: User asks about solver performance
- **WHEN** the user asks about computational cost or solver tuning
- **THEN** the skill SHALL explain that the direct solver scales roughly as O(N^1.5) memory and wall-time for 2D problems, that `lin_solver = 3` scales as O(N) memory and O(N log N) wall-time with a grid-independent V-cycle convergence factor ρ ≤ 0.15 for well-posed problems, and that the regression bounds demonstrated on the `add-gmg-stokes-defence` pinned EC2 harness are α ≤ 1.2 for memory and ≥ 3× wall-time speedup vs CHOLMOD at 201×201. Suggest using OpenMP (`-DOMP=ON`) for parallel sparse factorisation and tuning `penalty` for condition number under direct paths; for GMG paths, defaults are tuned and rarely need to change.

#### Scenario: User asks which solver to pick for a new simulation
- **WHEN** a user asks which `lin_solver` value to use for a new problem
- **THEN** the skill SHALL explain the decision tree: start with `lin_solver = 2` (default, Killer / Powell–Hestenes, direct) for small-to-medium problems up to ~500×500 where memory is not the constraint; pick `lin_solver = 3` (GMG-preconditioned FGMRES) for larger grids or when the direct path runs out of memory during factorisation; and fall through to `lin_solver = 1` (Powell–Hestenes KSP) or `lin_solver = 0` (direct on the augmented velocity block) only for specific edge cases.

#### Scenario: User asks how to verify the GMG path is working
- **WHEN** a user running `lin_solver = 3` wants a sanity check that GMG is doing the work rather than falling back
- **THEN** the skill SHALL point the user at three diagnostics: (1) the per-run log lines reporting FGMRES convergence history and any fallback-to-CHOLMOD events, (2) the `gmg_dump_vcycle = 1` instrumentation that emits residual/error HDF5 snapshots at every V-cycle operator, (3) the dual-solver comparison pattern (run the same scenario under `lin_solver = 0` and `lin_solver = 3`, compare HDF5 outputs with a small Python script) which is the canonical GMG regression check in MDOODZ's CI suite.

