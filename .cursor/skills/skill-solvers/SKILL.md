---
name: skill-solvers
description: MDOODZ numerical solver architecture — Stokes solver, Newton-Raphson iteration, Picard iteration, convergence criteria, staggered grid, thermal solver, SuiteSparse/CHOLMOD direct solver, time stepping, penalty method, and troubleshooting convergence issues.
---

# Numerical Solvers

## Main Time Loop

The simulation is orchestrated by `RunMDOODZ` in `MDLIB/Main_DOODZ.c`:

```
Read input (.txt) → Allocate grid & particles → Initialise fields
→ TIME LOOP:
    1. Set boundary conditions (SetBCs)
    2. Interpolate particles → grid (P2Mastah)
    3. NONLINEAR ITERATION (Stokes solve):
       a. Evaluate rheology
       b. Assemble Stokes matrix
       c. Solve linear system (CHOLMOD/UMFPACK)
       d. Check convergence
       e. If not converged → update, repeat
    4. Solve thermal equation
    5. Advect particles (RK integration)
    6. Advect free surface
    7. Reseed particles if needed
    8. Interpolate grid → particles
    9. Write output (every writer_step)
→ REPEAT until Nt steps
```

## Staggered Grid (MAC Layout)

MDOODZ uses a Marker-And-Cell staggered grid:

```
    Vz[i,j+1]
       |
  ──Vx[i,j]──P[i,j]──Vx[i+1,j]──
       |      η_n      |
    Vz[i,j]  (centre)  
       |                |
  ──────────η_s─────────── (vertex)
```

| Quantity | Grid location | Array suffix |
|----------|---------------|-------------|
| Pressure, normal stress, T, ρ, η | Cell centres | `_n` |
| Shear stress, shear strain rate | Cell vertices | `_s` |
| Vx (horizontal velocity) | Vertical cell edges | `u` |
| Vz (vertical velocity) | Horizontal cell edges | `v` |

Grid dimensions: `Nx` × `Nz` nodes → (Nx-1) × (Nz-1) cell centres.

## Stokes Assembly

The momentum and continuity equations are assembled in `StokesAssemblyDecoupled.c`:

**Momentum** (x): ∂τ_xx/∂x + ∂τ_xz/∂z - ∂P/∂x + ρg_x = 0
**Momentum** (z): ∂τ_xz/∂x + ∂τ_zz/∂z - ∂P/∂z + ρg_z = 0
**Continuity**: ∂v_x/∂x + ∂v_z/∂z + P/(penalty·η) = 0

The system is stored as a sparse matrix in CSC format (`SparseMat` struct).

### Penalty Method

The pressure is coupled to the velocity via a penalty parameter:

| Parameter | `.txt` field | Typical value | Description |
|-----------|-------------|---------------|-------------|
| Penalty | `penalty` | 1e2 | Larger = stronger incompressibility enforcement |

## Nonlinear Iteration

### Picard-to-Newton Strategy

1. **Picard iterations** (initial): Simple substitution — viscosity from previous iteration, no Jacobian
2. **Newton iterations** (after tolerance reached): Full Jacobian linearisation using D-coefficients (`D11`–`D34` on grid)
3. Transition at: `Picard2Newton_tol` residual threshold

### Convergence Control (Nparams struct)

| Field | Description |
|-------|-------------|
| `nit_max` | Maximum nonlinear iterations per time step |
| `nonlin_abs_mom` | Absolute momentum residual tolerance |
| `nonlin_rel_mom` | Relative momentum residual tolerance |
| `nonlin_abs_div` | Absolute divergence residual tolerance |
| `nonlin_rel_div` | Relative divergence residual tolerance |
| `Picard2Newton` | Enable Picard→Newton switch (0/1) |
| `Picard2Newton_tol` | Residual threshold to switch |
| `let_res_grow` | Allow residual to grow temporarily (0/1) |
| `stagnated` | Stagnation detection flag |

### Residuals

| Residual | Code field | Physical meaning |
|----------|-----------|------------------|
| `resx` | Horizontal momentum | Force balance error (x) |
| `resz` | Vertical momentum | Force balance error (z) |
| `resp` | Continuity | Mass conservation error |

### Line Search

Step-size control for Newton updates. Parameters:
- `line_search`: Enable (0/1)
- `line_search_min`: Minimum step fraction

### Safe Mode

| Parameter | Description |
|-----------|-------------|
| `safe_mode` | Auto-reduce dt on divergence (0/1) |
| `safe_dt_div` | dt reduction factor |
| `max_num_stag` | Max stagnation count before dt reduction |

## Linear Solver

MDOODZ uses **SuiteSparse** direct solvers:

| Component | Purpose |
|-----------|---------|
| CHOLMOD | Sparse Cholesky factorisation |
| UMFPACK | Unsymmetric multifrontal LU |
| AMD/COLAMD | Matrix ordering for sparsity |

The `DirectSolver` struct wraps CHOLMOD:

| Field | Description |
|-------|-------------|
| `Lfact` | CHOLMOD factorisation object |
| `Analyze` | Whether symbolic analysis is done |
| `phase` | Solver phase (analyse/factorise/solve) |

### Tuning

| Parameter | `.txt` field | Description |
|-----------|-------------|-------------|
| Linear solver type | `lin_solver` | Solver variant selection |
| Diagonal scaling | `diag_scaling` | Pre-condition with diagonal scaling (0/1) |
| Preconditioner | `preconditioner` | Preconditioner type |
| Max KSP iterations | `max_its_KSP` | For iterative refinement |
| Relative KSP tolerance | `rel_tol_KSP` | Convergence for iterative refinement |

## Thermal Solver

Solves the energy equation (in `ThermalSolver.c` and `ThermalRoutines.c`):

ρCp · DT/Dt = ∇·(k∇T) + H_shear + H_adiab + Qr

| Source term | Switch | Description |
|-------------|--------|-------------|
| Shear heating | `shear_heating = 1` | H = τ_ij · ε̇_ij |
| Adiabatic heating | `adiab_heating = 1` | H = α·T·v_z·ρ·g |
| Radiogenic | Always (if Qr > 0) | H = Qr per phase |

**Subgrid diffusion** (`subgrid_diffusion`): Corrects particle temperatures for numerical diffusion. Mode 1 or 2.

## Convergence Monitoring

| Method | Parameter | Description |
|--------|-----------|-------------|
| Console output | `noisy = 1` | Print iteration-by-iteration residuals |
| Gnuplot plots | `gnuplot_log_res = 1` | Live residual evolution plots |
| Log arrays | `Nparams.rx_abs[]` etc. | Stored per-iteration residual history |

## Troubleshooting Convergence

| Symptom | Possible cause | Remedy |
|---------|---------------|--------|
| Residuals increase | Time step too large | Reduce `dt`, enable `safe_mode = 1` |
| Residuals stagnate | Picard not reaching Newton | Lower `Picard2Newton_tol` |
| NaN appears | Viscosity out of bounds | Check `min_eta`/`max_eta`, material params |
| Very slow convergence | Poor penalty value | Adjust `penalty` (try 1e1–1e4) |
| Oscillating residuals | Elastic step too large | Reduce `dt` relative to Maxwell time (η/G) |
| Newton diverges | Bad initial guess | Increase Picard iterations before Newton |

## Numerical Accuracy & Convergence Orders

Expected convergence orders for the staggered-grid marker-in-cell scheme (SolVi benchmark, circular inclusion η_c/η_m = 1000):

| Field | Observed order | Theoretical (smooth) | Limiting factor |
|-------|---------------|---------------------|------------------|
| Vx, Vz | ~1.0 | 2.0 | Viscosity discontinuity at inclusion boundary |
| P | ~0.75 | 1.0–2.0 | FD truncation error at material interface |

The error is dominated by the FD stencil at cells straddling the viscosity discontinuity, **not** by marker interpolation noise. Increasing markers per cell (from 4×4 to 16×16) gives negligible improvement (<12% for P, 0% for Vx) at 7× runtime cost.

### Viscosity averaging (`eta_average`)

The `eta_average` .txt parameter controls how marker viscosities are mixed to grid nodes within each cell (0 = arithmetic, 1 = harmonic, 2 = geometric). Arithmetic gives the best convergence orders; geometric gives the best absolute P error at fixed resolution. This is **not** cell-face averaging — the Stokes stencil in `StokesAssemblyDecoupled.c` uses pre-computed constitutive tensor components (`D11_n`, `D33_s`) without cell-face interpolation.

See `TESTS/SolViMarkerComparison.cpp` for full experimental data.

## Key Source Files

- Main loop: `MDLIB/Main_DOODZ.c` (`RunMDOODZ`)
- Stokes assembly: `MDLIB/StokesAssemblyDecoupled.c`
- Linear solvers: `MDLIB/Solvers.c`
- Sparse tools: `MDLIB/SparseTools.c`
- Thermal solver: `MDLIB/ThermalSolver.c`, `MDLIB/ThermalRoutines.c`
- Jacobian: `MDLIB/FD_Jacobian.c`
