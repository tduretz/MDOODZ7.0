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

MDOODZ exposes two linear-solver families for the Stokes system:

1. **SuiteSparse direct solvers** (`lin_solver ∈ {-1, 0, 1, 2}`): CHOLMOD + UMFPACK factorisation of the assembled decoupled Stokes matrix. The historical workhorse; O(N^1.5) memory and wall-time on 2D grids.
2. **GMG-preconditioned FGMRES** (`lin_solver = 3`): geometric-multigrid hierarchy of MAC-staggered grids + Vanka 5×5 block smoother + CHOLMOD (UMFPACK in practice) coarsest-level direct solve, wrapped in a flexible GMRES outer driver. O(N) memory and O(N log N) wall-time; added in `add-gmg-stokes-solver`, measured in `add-gmg-stokes-defence`.

Transparent fall-back: `lin_solver = 3` silently routes to CHOLMOD if FGMRES fails to converge or if a configuration outside its current envelope is requested. Default envelope (after Phase 2): Newton + anisotropy are both supported.

### Direct path — SuiteSparse

| Component | Purpose |
|-----------|---------|
| CHOLMOD | Sparse Cholesky factorisation (SPD blocks, e.g. PH-augmented velocity) |
| UMFPACK | Unsymmetric multifrontal LU (full saddle-point) |
| AMD/COLAMD | Matrix ordering for sparsity |

The `DirectSolver` struct wraps CHOLMOD:

| Field | Description |
|-------|-------------|
| `Lfact` | CHOLMOD factorisation object |
| `Analyze` | Whether symbolic analysis is done |
| `phase` | Solver phase (analyse/factorise/solve) |

### GMG path — `lin_solver = 3`

| File | Role |
|------|------|
| `MDLIB/MultigridLevels.{h,c}` | Hierarchy build, restrict/prolongate operators, active-mask propagation |
| `MDLIB/MultigridStokes.{h,c}` | V-cycle driver, Vanka smoother, FGMRES outer loop, one-shot `gmg_dump_vcycle` instrumentation |
| `MDLIB/StokesAssemblyGMG.{h,c}` | Fine-level matvec (stencil-bridge equivalent of `BuildStokesOperatorDecoupled`; guarded by a 1e-12 golden cross-check test — design D11) |

Operator choices:

| Operator | Velocity | Pressure |
|---------|----------|----------|
| Restriction | full-weighting | injection |
| Prolongation | bilinear | piecewise-constant |
| Viscosity restriction | harmonic-averaged above contrast ratio 10 | n/a |
| Smoother | Vanka 5×5 block, ω = 0.6 under-relaxation | coupled in Vanka block |
| Coarsest | UMFPACK direct solve on restricted Picard K | coupled |
| Outer driver | FGMRES, restart 30 | — |

Picard ↔ Newton phasing (design D7): the V-cycle's fine-level matvec tracks the current nonlinear mode, but the smoother and coarse-grid operator always use the Picard block — this keeps the preconditioner SPD and the coarse solve CHOLMOD-stable.

### Tuning

| Parameter | `.txt` field | Path | Description |
|-----------|--------------|------|-------------|
| Linear solver type | `lin_solver` | both | `{-1,0,1,2} = direct`, `3 = GMG-FGMRES` |
| Diagonal scaling | `diag_scaling` | direct | Pre-condition with diagonal scaling (0/1) |
| Preconditioner | `preconditioner` | direct | Preconditioner type |
| Max KSP iterations | `max_its_KSP` | direct | For iterative refinement |
| Relative KSP tolerance | `rel_tol_KSP` | direct | Convergence for iterative refinement |
| Multigrid levels | `gmg_levels` | GMG | 0 = auto (coarsen until ≤ 21²) |
| Pre-smooth sweeps | `gmg_nu_pre` | GMG | Typically 2 |
| Post-smooth sweeps | `gmg_nu_post` | GMG | Typically 2 |
| FGMRES restart | `gmg_fgmres_restart` | GMG | Typically 30 |
| FGMRES tolerance | `gmg_fgmres_tol` | GMG | Typically 1e-11 |
| Standalone (no fall-back) | `gmg_standalone` | GMG | 1 = hard-fail instead of routing to CHOLMOD |
| V-cycle snapshot dump | `gmg_dump_vcycle` | GMG | 1 = one-shot HDF5 dump per operator step |

### Which solver for what problem?

Quick decision tree:

* **Small to medium (≤ 500²), no memory pressure**: `lin_solver = 2` (Killer / Powell–Hestenes direct) — the historical default.
* **Larger grids or memory-bound**: `lin_solver = 3` (GMG-FGMRES) — O(N) memory, falls back to direct if it can't converge.
* **Edge cases** (e.g. specific compressibility, augmented-Lagrangian sweeps): `lin_solver ∈ {1, 0, -1}` as documented in `Solvers.c` and design D1 of `add-gmg-stokes-solver`.

### Verifying GMG is doing the work (not falling back)

Three diagnostics:

1. Per-run log lines report `FGMRES iter N | |r| = ...` and any `lin_solver = 3 → CHOLMOD fall-back` events; a solve that actually used GMG will have a block of FGMRES residual lines, a fallback will skip straight to the direct-solver header.
2. `gmg_dump_vcycle = 1` emits `${writer_subfolder}/vcycle_dump/level_{N}_{step}_{pre|post}_{seq}.h5` for the very first V-cycle; `benchmarks/vcycle_vis/make_gif.py` renders those into the animation + 6-frame stills used in the defence §10.
3. Dual-solver CI pattern: run the same `.txt` under `lin_solver = 0` (CHOLMOD twin) and `lin_solver = 3` (GMG twin), HDF5-compare the outputs. See `TESTS/GmgStokesEquivalence.cpp`, `TESTS/DisabledBenchmarks/TopoRelaxGmgEquivalence.cpp`, and `TESTS/DisabledBenchmarks/SolKzGmgEquivalence.cpp` for canonical templates (1e-8 relative L2 on velocities, 1e-6 on mean-subtracted pressure).

### Performance characteristics (measured)

`lin_solver = 3` has been measured on the `add-gmg-stokes-defence` pinned EC2 harness (`c7i.8xlarge`, `eu-west-1`). Headline regressions bounds:

* Memory-scaling exponent α ≤ 1.2 (theoretical bound for GMG: α = 1.0).
* 201×201 wall-time speedup ≥ 3× vs `lin_solver = 0` CHOLMOD baseline.
* V-cycle convergence factor ρ ≤ 0.15 on well-posed Stokes problems.

Point measurements and 95% confidence intervals live in `defence.pdf` §6 Measurements, with per-grid-size detail in `benchmarks/ec2/results/summary.md`. If a measurement falls short of a bound, the defence records it honestly in §8 Limitations rather than restating the bound as if met.

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
| Segfault in `cs_di_compress` or CHOLMOD | Memory corruption from buffer overflow | Build with ASAN (see skill-build-and-run) |

## Known Solver Bugs (Fixed)

### BCp.type guard in pressure-block setup

**Affected functions**: `DirectStokesDecoupledComp`, `KillerSolver`, `KSPStokesDecoupled` in `Solvers.c`.

The penalty/compressibility loop over pressure cells originally used `BCp.type[k] != 30 && != 31` as the guard condition. This is **incorrect** — `NumberStokes()` (in `StokesRoutines.c`) only assigns valid equation numbers (`eqn_p[k] >= 0`) for cells with `BCp.type[k] == -1`. All other types (0, 1, 30, 31) get `eqn_p[k] = -1`. The old guard let type=0 cells through, producing `i = -1 - matA->neq` (a wild pointer into `D1cm0->x` / `Dcm0->x`), causing heap-buffer-overflow → segfault.

**Fix**: Changed guard to `BCp.type[k] == -1`.

Additionally, the compressible branch `((double*)D1cm0->x)[k]` used cell index `k` instead of equation index `i`. Fixed to `D1cm0->x[i]`.

### Thermal loop bounds in EvaluateCourantCriterion

**Affected function**: `EvaluateCourantCriterion` in `AdvectionRoutines.c`.

The thermal CFL sub-loop used Vz-vertex-grid bounds `(k < Nx+1, l < Nz)` to iterate over centre-grid arrays `mesh->T` and `mesh->T0_n` of size `(Nx-1)*(Nz-1)`. This reads past the end of the array. **Fix**: bounds changed to `(k < Nx-1, l < Nz-1)`.

## Thermal Steady-State Timing

For convection benchmarks (e.g., Blankenbach), MDOODZ advects temperature via markers (not on the grid). This means **dt is CFL-limited by velocity**, not by the diffusion timescale. With adaptive dt (`constant_dt=0`), reaching thermal steady state (~2 diffusion times, td = H²/κ) can require 100k+ time steps. For a 41×41 grid this takes ~17 hours. Plan accordingly — either use long CI timeouts, coarser grids, or restructure such tests as smoke tests (verify physics start correctly, not final converged values).

## free_surface and Boundary Conditions

`free_surface=1` imposes **zero traction** (σ·n = 0) at the top boundary, not free-slip (v·n = 0). For closed-box benchmarks (Blankenbach, Rayleigh-Bénard) that require free-slip on all walls, use `free_surface=0`. Using `free_surface=1` for such benchmarks gives ~35% error in Vrms because material can flow through the top surface.

## Numerical Accuracy & Convergence Orders

Expected convergence orders for the staggered-grid marker-in-cell scheme (SolVi benchmark, circular inclusion η_c/η_m = 1000):

| Field | L2 order (41→81) | L1 order (101→201) | Theoretical (smooth) | Limiting factor |
|-------|------------------|---------------------|---------------------|------------------|
| Vx, Vz | ~1.0 | ~0.90 | 2.0 | Viscosity discontinuity at inclusion boundary |
| P | ~0.75 | ~0.85 | 1.0–2.0 | FD truncation error at material interface |

The L1 norm (`mean(|num - ana|)`) gives a cleaner convergence signal for pressure than relative L2, which is ill-conditioned near zero-crossings. At resolutions >80 cells the inclusion is well-resolved and L1 P order reaches ~0.85 — significantly better than L2 (~0.65). L2(P) order oscillates wildly between resolution pairs (including negative values) due to the relative norm amplifying errors where the analytical P is near zero.

The error is dominated by the FD stencil at cells straddling the viscosity discontinuity, **not** by marker interpolation noise. Increasing markers per cell (from 4×4 to 16×16) gives negligible improvement (<12% for P, 0% for Vx) at 7× runtime cost.

### Viscosity averaging (`eta_average`)

The `eta_average` .txt parameter controls how marker viscosities are mixed to grid nodes within each cell (0 = arithmetic, 1 = harmonic, 2 = geometric). Arithmetic gives the best convergence orders (Vx 1.02, P 0.75); harmonic gives the lowest absolute Vx error but worst convergence (0.45); geometric gives the best absolute P error but lower convergence. This is **not** cell-face averaging — the Stokes stencil in `StokesAssemblyDecoupled.c` uses pre-computed constitutive tensor components (`D11_n`, `D33_s`) without cell-face interpolation.

### Cell-face D-tensor averaging: tested and ruled out

Replacing D11E/D11W (etc.) in the stencil with their harmonic mean was attempted and **failed** — it produces a non-SPD stiffness matrix (CHOLMOD failure). In this staggered grid, D values at cell centres already represent the natural face values for the differencing scheme. Setting D11E = D11W = harmonic_mean destroys the viscosity contrast information and breaks the operator structure. This is fundamentally different from the "harmonic face averaging" in Deubelbeiss & Kaus (2008), which refers to marker-to-node interpolation (what `eta_average` already controls). Improving P convergence order beyond ~0.85 (L1) would require immersed-boundary or ghost-fluid methods.

See `TESTS/SolViMarkerComparison.cpp` for full experimental data.

## Key Source Files

- Main loop: `MDLIB/Main_DOODZ.c` (`RunMDOODZ`)
- Stokes assembly: `MDLIB/StokesAssemblyDecoupled.c`
- Linear solvers: `MDLIB/Solvers.c`
- Sparse tools: `MDLIB/SparseTools.c`
- Thermal solver: `MDLIB/ThermalSolver.c`, `MDLIB/ThermalRoutines.c`
- Jacobian: `MDLIB/FD_Jacobian.c`
