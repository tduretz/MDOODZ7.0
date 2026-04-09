## ADDED Requirements

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
The skill SHALL explain that MDOODZ uses SuiteSparse (CHOLMOD + UMFPACK) for direct sparse matrix solves. The `DirectSolver` struct wraps CHOLMOD factorisation. The `lin_solver` parameter selects the solver variant. The skill SHALL mention that `diag_scaling` and `preconditioner` parameters can tune solver performance.

#### Scenario: User asks about solver performance
- **WHEN** the user asks about computational cost or solver tuning
- **THEN** the skill SHALL explain that the direct solver is the main computational bottleneck, scaling roughly as O(N^1.5) for 2D problems. Suggest using OpenMP (`-DOMP=ON`) for parallel sparse factorisation and tuning `penalty` for condition number.

### Requirement: Skill documents convergence monitoring
The skill SHALL explain the residual monitoring: momentum residuals (`resx`, `resz`), pressure/continuity residual (`resp`), relative and absolute tolerances, the `noisy` flag for verbose output, and the `gnuplot_log_res` option for live residual plotting.

#### Scenario: User wants to monitor solver convergence at runtime
- **WHEN** the user asks how to track convergence
- **THEN** the skill SHALL explain: set `noisy = 1` for iteration-by-iteration residual printout, and `gnuplot_log_res = 1` for gnuplot-based live residual plots. The `Nparams` struct stores per-iteration residual history in `rx_abs`, `rz_abs`, `rp_abs` arrays.
