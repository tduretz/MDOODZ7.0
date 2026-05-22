## Purpose

Defines the behaviour contract for the `skill-solvers` agent skill: what the skill must teach a user about MDOODZ's numerical solvers (Stokes + thermal), its nonlinear iteration strategy, its sparse linear algebra backends (SuiteSparse direct and GMG-preconditioned FGMRES), convergence monitoring, and troubleshooting. The spec is consumed by `.github/skills/skill-solvers/SKILL.md`, which is the documentation users actually read — any drift between spec and skill must be reconciled before shipping.
## Requirements
### Requirement: Skill documents the main time-stepping loop
The skill SHALL explain the overall solver architecture in `Main_DOODZ.c` (`RunMDOODZ` function): read input → allocate grid and particles → initialise fields → enter time loop → (nonlinear Stokes solve → thermal solve → advect particles → update properties → write output) → repeat until Nt steps.

#### Scenario: User asks about the simulation flow
- **WHEN** the user asks what happens when a simulation runs
- **THEN** the skill SHALL describe the main loop sequence: (1) grid allocation and coordinate setup, (2) particle initialisation via callbacks, (3) time-stepping with mechanical solve, thermal solve, advection, reseeding, and output writing

### Requirement: Skill documents the staggered-grid finite difference discretisation
The skill SHALL explain that MDOODZ uses a staggered grid (Marker-And-Cell / MAC layout): pressure and normal stresses at cell centres, shear stress and strain rate at cell vertices, Vx on vertical cell edges, Vz on horizontal cell edges. Grid dimensions are defined by `Nx` and `Nz` in the `.txt` file.

#### Scenario: User asks about the numerical grid
- **WHEN** the user asks how the domain is discretised
- **THEN** the skill SHALL explain the staggered grid layout, the relationship between `Nx`/`Nz` and the number of cell centres (Nx-1, Nz-1) and vertices (Nx, Nz), and how `dx`/`dz` are computed from domain bounds

### Requirement: Skill documents the Stokes assembly
The skill SHALL explain the decoupled velocity–pressure Stokes formulation assembled in `StokesAssemblyDecoupled.c`. The momentum equations are discretised as a sparse matrix system, with the coefficient matrix incorporating variable viscosity (from rheology evaluation). The system is solved as a saddle-point problem with penalty method for pressure.

#### Scenario: User asks about the linear system
- **WHEN** the user asks how the Stokes equations are solved
- **THEN** the skill SHALL explain: momentum + continuity assembled into a sparse matrix (CSC format via `SparseMat` struct), solved with CHOLMOD/UMFPACK direct solver, with `penalty` parameter controlling pressure coupling strength

### Requirement: Skill documents Newton–Raphson nonlinear iteration
The skill SHALL explain the nonlinear iteration strategy:
- **Picard iteration**: initial iterations with simple substitution (viscosity from previous iteration)
- **Newton iteration**: switched on after `Picard2Newton_tol` residual threshold, using a Jacobian (`D11`–`D34` coefficients on the grid) for consistent linearisation
- Convergence criteria: `nonlin_abs_mom`, `nonlin_rel_mom`, `nonlin_abs_div`, `nonlin_rel_div` in the `Nparams` struct
- Line search for step-size control
- Maximum iterations `nit_max`, stagnation detection

#### Scenario: User encounters convergence issues
- **WHEN** the user reports that the nonlinear solver does not converge
- **THEN** the skill SHALL suggest checking: time step size (`dt`), penalty parameter, min/max viscosity bounds (`min_eta`, `max_eta`), mesh resolution, and whether the rheology is appropriate. Also suggest enabling `safe_mode` for automatic dt reduction on divergence.

### Requirement: Skill documents the thermal solver
The skill SHALL explain the thermal solver in `ThermalSolver.c` and `ThermalRoutines.c`: solves the heat equation with conduction (thermal conductivity `k`), shear heating (`shear_heating` switch), adiabatic heating (`adiab_heating`), radiogenic heat production (`Qr`), and subgrid diffusion for particle temperature (`subgrid_diffusion`). Temperature BCs are set via `SetBCT` callback.

#### Scenario: User asks about thermal evolution
- **WHEN** the user asks how temperature evolves
- **THEN** the skill SHALL explain: the energy equation is solved implicitly on the grid, then temperatures are interpolated back to particles with optional subgrid diffusion correction (modes 1 or 2)

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

### Requirement: Skill documents the GMG-preconditioned FGMRES path
The skill SHALL explain the GMG-preconditioned FGMRES Stokes solver introduced in `add-gmg-stokes-solver` and measured in `add-gmg-stokes-defence`: a geometric multigrid hierarchy of MAC-staggered grids with Vanka 5×5 block smoothing, coarse-level UMFPACK direct solve, FGMRES outer driver wrapping the V-cycle as a right preconditioner. The skill SHALL reference `MDLIB/MultigridLevels.{h,c}`, `MDLIB/MultigridStokes.{h,c}`, and `MDLIB/StokesAssemblyGMG.{h,c}` as the implementation files. The skill SHALL point at the reading library (`defence.pdf`, `stokes_internals.pdf`, `journey.pdf`) as the narrative sources.

#### Scenario: User asks how GMG works in MDOODZ
- **WHEN** a user asks to understand the GMG path architecture
- **THEN** the skill SHALL describe: V-cycle with restriction (full-weighting for velocity, injection for pressure) and prolongation (bilinear for velocity, piecewise-constant for pressure), harmonic-averaged viscosity restriction above contrast ratio 10, Vanka 5×5 block smoother with ω = 0.6 under-relaxation, UMFPACK direct solve at the coarsest level, FGMRES with restart 30 as outer driver; and SHALL direct the user to the primer PDF for the longer explanation.

#### Scenario: User asks why the fine-level operator is duplicated
- **WHEN** a user asks why `StokesAssemblyGMG.c` exists separately from `StokesAssemblyDecoupled.c`
- **THEN** the skill SHALL explain design decision D11 from `add-gmg-stokes-solver` (Option A — additive duplication protected by a 1e-12 golden cross-check test) and the rationale: backward compatibility for `lin_solver ∈ {0, 1, 2, -1}` is absolute; the duplication cost is bounded by the cross-check test which fails loud on any drift; the `MDOODZ_GMG_DRIFT_CANARY` compile flag proves the test catches drift.

### Requirement: Skill documents the performance characteristics of `lin_solver = 3`
The skill SHALL explain that `lin_solver = 3` has been performance-measured on the `add-gmg-stokes-defence` pinned EC2 harness (`c7i.8xlarge`, `eu-west-1`) and SHALL cite the measured memory-scaling exponent α and 201×201 wall-time ratio. If either measurement violates the regression bound the skill SHALL record the shortfall honestly rather than quote the bound as though it had been met.

#### Scenario: User asks about GMG's measured speedup
- **WHEN** the user asks how much faster or more memory-efficient GMG is in practice
- **THEN** the skill SHALL cite the measured α (e.g. `α = 1.12 ± 0.03`) and the 201² wall-time ratio (e.g. `3.4× faster than CHOLMOD`), both sourced from `defence.pdf` §6 Measurements, and SHALL link to the `benchmarks/ec2/results/summary.md` for per-grid-size detail.

### Requirement: Skill documents convergence monitoring
The skill SHALL explain the residual monitoring: momentum residuals (`resx`, `resz`), pressure/continuity residual (`resp`), relative and absolute tolerances, the `noisy` flag for verbose output, and the `gnuplot_log_res` option for live residual plotting.

#### Scenario: User wants to monitor solver convergence at runtime
- **WHEN** the user asks how to track convergence
- **THEN** the skill SHALL explain: set `noisy = 1` for iteration-by-iteration residual printout, and `gnuplot_log_res = 1` for gnuplot-based live residual plots. The `Nparams` struct stores per-iteration residual history in `rx_abs`, `rz_abs`, `rp_abs` arrays.

## MODIFIED Section: SolVi Convergence

### Requirement: Document L1 convergence findings
The SolVi convergence subsection of skill-solvers SHALL document that L1 norm at resolutions >80 cells recovers expected first-order pressure convergence (order ~1.0), while L2 norm at coarse resolutions (≤81) shows a misleadingly low ~0.75 due to the inclusion interface.

#### Scenario: Skills reflects L1 vs L2 insight
- **WHEN** a developer reads about SolVi convergence in skill-solvers
- **THEN** the section SHALL explain that L1 norm is recommended for evaluating P convergence and list measured orders at 101→201
