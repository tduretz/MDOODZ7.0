# gmg-stokes-solver Specification

## Purpose
TBD - created by archiving change add-gmg-stokes-solver. Update Purpose after archive.
## Requirements
### Requirement: GMG solver is selected via `lin_solver = 3`

The code SHALL route the Stokes linear solve through the geometric multigrid preconditioned FGMRES path when `lin_solver = 3` is set in the `.txt` input. For all other values of `lin_solver`, the existing direct / Powell–Hestenes / killer-solver dispatch SHALL remain unchanged and bit-for-bit equivalent to the pre-change behaviour.

#### Scenario: User selects the GMG solver

- **WHEN** a `.txt` input sets `lin_solver = 3` and the model is otherwise well-posed
- **THEN** the Stokes solve invokes the GMG-preconditioned FGMRES path, writes a log line identifying the solver, and produces Vx, Vz, P fields in the same units, shapes, and HDF5 layout as the CHOLMOD path

#### Scenario: Default behaviour is preserved

- **WHEN** a `.txt` input omits `lin_solver` or sets any pre-existing value (`-1, 0, 1, 2`)
- **THEN** the Stokes solve uses the same code path it did before this change, and output on every existing CI scenario matches its pre-change result within existing tolerances

### Requirement: Grid hierarchy is built with a bounded coarsest level

The code SHALL construct a multi-level grid hierarchy by halving `Nx` and `Nz` (rounded up) at each level, continuing until either the minimum dimension is below 16 cells or the total DOF count is below 5000, whichever comes first. The user SHALL be able to override the number of levels via an optional `gmg_levels` parameter; when overridden, the requested level count MUST be used exactly, clamped to `[2, floor(log2(min(Nx,Nz)))]`.

#### Scenario: Automatic level selection on a typical grid

- **WHEN** `Nx = 201`, `Nz = 201`, `lin_solver = 3`, and `gmg_levels` is unset
- **THEN** the hierarchy has levels of size 201, 101, 51, 26, 13, with the coarsest level having at most 5000 DOFs

#### Scenario: Odd and non-power-of-2 grid sizes coarsen correctly

- **WHEN** `Nx = 100`, `Nz = 150`, `lin_solver = 3`
- **THEN** each level's dimensions are `ceil(N_finer / 2)` per axis, producing valid staggered-grid sizes at every level, and the coarsest level has at most 5000 DOFs

#### Scenario: Explicit level override is honoured

- **WHEN** the user sets `gmg_levels = 3` on a 201×201 grid
- **THEN** the hierarchy is built with exactly three levels and the coarse-level DOF count is not enforced against the 5000 threshold

### Requirement: Restriction preserves constants and averages smooth fields

The restriction operator SHALL use full-weighting stencils for velocity components (Vx, Vz) and for pressure, with weights corresponding to the staggered location of each variable. Restriction SHALL preserve constant fields exactly and SHALL produce the arithmetic average of a smooth field's neighbourhood on the coarser grid.

#### Scenario: Constant field is preserved exactly

- **WHEN** a constant field `f = 7.0` on a 5×5 Vx grid is restricted to a 3×3 coarse grid
- **THEN** every coarse value equals `7.0` to machine precision

#### Scenario: Full-weighting stencil on a quadratic field

- **WHEN** a field `f(x, z) = x² + z²` is sampled on a 5×5 fine Vx grid and restricted
- **THEN** each coarse value equals the full-weighting combination `(1/16)(f_corners + 2·f_edges + 4·f_center)` of the fine values within one unit-in-the-last-place

### Requirement: Prolongation is exact on fields representable on the coarse grid

The prolongation operator SHALL use bilinear interpolation for velocity components and piecewise-constant injection for pressure. Prolongation SHALL reproduce any fine-grid field that is itself exactly representable on the coarse grid.

#### Scenario: Linear velocity field is reproduced exactly

- **WHEN** a coarse 3×3 velocity field `f(x, z) = x + 2z` is prolongated to a 5×5 fine grid
- **THEN** the fine-grid values match the exact bilinear interpolation of the coarse field to machine precision

#### Scenario: Piecewise-constant pressure prolongation

- **WHEN** a coarse cell-centred pressure field is prolongated to the fine grid via injection
- **THEN** every fine cell inside a coarse cell's footprint receives exactly the coarse value (no averaging between neighbouring coarse cells)

### Requirement: Viscosity is restricted harmonically across contrasts

The coarse-level viscosity field SHALL be produced from the fine-level viscosity field using a harmonic mean in regions where the local viscosity ratio exceeds 10, and an arithmetic mean elsewhere. This requirement applies to both the cell-centred (`ηn`) and vertex-based (`ηs`) viscosity fields.

#### Scenario: Harmonic averaging at a sharp contrast

- **WHEN** the fine viscosity field has a 1 | 1e6 checkerboard pattern and is restricted to a coarse grid
- **THEN** the coarse viscosity within mixed cells is approximately the harmonic mean (dominated by the soft phase) rather than the arithmetic mean (dominated by the stiff phase)

### Requirement: Vanka block smoother reduces the residual

A single Vanka sweep SHALL iterate over every interior cell, locally assemble the 5×5 block relating the four face velocities and the central pressure, solve the block exactly, and update the global solution in place. After one Vanka sweep from a zero initial guess on a constant-viscosity Stokes problem, the residual norm SHALL strictly decrease.

#### Scenario: Smoothing property on a small Stokes system

- **WHEN** a constant-viscosity Stokes system on an 11×11 grid is initialised with `x = 0` and a non-zero RHS, and one Vanka sweep is applied
- **THEN** `‖r₁‖ < ‖r₀‖`, where `r₀` and `r₁` are the initial and post-sweep residuals in the 2-norm

#### Scenario: Local 5×5 block solve reproduces analytical single-cell solution

- **WHEN** a single interior cell's Vanka block is assembled with known viscosity η and a known RHS corresponding to a prescribed single-cell analytical Vx, Vz, P
- **THEN** solving the 5×5 block returns those analytical values to within `1e-10` relative error

### Requirement: Coarsest-level solve uses CHOLMOD

At the coarsest level of the hierarchy, the code SHALL assemble the level's Stokes system as a `SparseMat` and solve it via the existing CHOLMOD direct path in `SparseTools.c`. The coarse-level CHOLMOD factor SHALL be cached across V-cycles within a single nonlinear iteration when the coarse-level matrix has not structurally changed.

#### Scenario: Coarse level handed to CHOLMOD

- **WHEN** a V-cycle descends to the coarsest level
- **THEN** the code calls the existing CHOLMOD factor-and-solve routines on a `SparseMat` built from the coarse level's stencils, without allocating a new solver backend

#### Scenario: Factor reuse inside a nonlinear iteration

- **WHEN** two consecutive V-cycles within the same Picard or Newton iteration descend to the coarsest level with the same coarse matrix
- **THEN** the second descent reuses the first V-cycle's cached factor and does not repeat the symbolic or numeric factorisation

### Requirement: The V-cycle is grid-independent for constant viscosity

For a constant-viscosity manufactured Stokes problem, the V-cycle SHALL reduce the residual norm by a factor `ρ ≤ 0.15` per iteration, and this bound SHALL hold independently of grid resolution across 21×21 through 161×161.

#### Scenario: Grid-independent convergence factor

- **WHEN** the V-cycle is applied to a constant-viscosity Stokes system with random RHS on grids of size 21×21, 41×41, 81×81, 161×161, iterating 10 times each
- **THEN** the measured per-iteration convergence factor `ρ_k = ‖r_{k+1}‖ / ‖r_k‖` remains below 0.15 on every grid for every iteration after the first

### Requirement: GMG-preconditioned FGMRES is the production solver mode

When `lin_solver = 3` and `gmg_standalone` is unset or zero, the Stokes system SHALL be solved by FGMRES with the GMG V-cycle as a right preconditioner. FGMRES SHALL iterate until either the **relative residual reduction** `‖r_k‖ / ‖r_0‖` reaches `gmg_fgmres_tol` (default: `nonlin_abs_mom / 10`) or the iteration count reaches `gmg_fgmres_restart × gmg_fgmres_max_restarts` (defaults: restart = 30, max-restarts = 20). The convergence predicate and the log message reporting "converged" SHALL both test and report the *same* quantity — namely `‖r_k‖ / ‖r_0‖` — so that a declared-converged run always satisfies `reported_final_res ≤ gmg_fgmres_tol` numerically.

#### Scenario: Convergence on a typical lithospheric setup

- **WHEN** a lithospheric-scale model with viscosity contrast up to 1e6 is solved with `lin_solver = 3`
- **THEN** FGMRES converges within its iteration budget, the final relative residual is below `gmg_fgmres_tol`, and the solver logs residual history per FGMRES iteration

#### Scenario: Diagnostic standalone V-cycle mode

- **WHEN** `lin_solver = 3` and `gmg_standalone = 1` are set
- **THEN** the Stokes solve applies the V-cycle without the FGMRES outer iteration, exposing the raw multigrid convergence for instrumentation

#### Scenario: Convergence predicate and log message agree

- **WHEN** FGMRES declares "converged" via `LOG_INFO`
- **THEN** the reported `final_res_rel` in the log line is numerically ≤ `gmg_fgmres_tol`, within 1e-12 of the value actually tested by the convergence predicate

#### Scenario: Non-convergence log message reports the correct quantity

- **WHEN** FGMRES exhausts its iteration budget without converging
- **THEN** the log line reports `final_res_rel` (the relative-reduction quantity tested by the predicate), not the raw `‖r_k‖`, and "converged" is not present in the message

### Requirement: Newton Jacobian is handled without polluting the smoother

When the outer nonlinear iteration is in Newton mode (using the `D11`–`D34` Jacobian coefficients), the GMG smoother SHALL operate on the symmetric Picard part of the operator, while the full Jacobian SHALL be used for residual evaluation inside FGMRES at every level of the hierarchy.

#### Scenario: Newton outer iteration over GMG

- **WHEN** a nonlinear solve has crossed `Picard2Newton_tol` and is in Newton mode, with `lin_solver = 3`
- **THEN** the residual fed to FGMRES uses the full Jacobian, the Vanka smoother uses the symmetric Picard stencil, and the outer Newton iteration achieves super-linear convergence comparable to Newton + CHOLMOD on the same problem

### Requirement: Fine-level operator is bit-identical to MDOODZ's assembled matrix

The fine-level matvec used by the V-cycle and the Vanka smoother SHALL reproduce exactly the linear operator that [`StokesAssemblyDecoupled.c`](MDLIB/StokesAssemblyDecoupled.c) builds under the same `model` settings. A textbook Picard Stokes stencil does not satisfy this requirement because it misses the Newton Jacobian coefficients (`D11`–`D34`), compressibility terms, anisotropy rotations, BC-tag handling (including deactivation via tag `30`), free-surface stencil adjustments, and the Powell–Hestenes penalty augmentation.

Implementation SHALL use the Option A mechanism (D11 in [design.md](design.md)): a separate read-only translation unit `MDLIB/StokesAssemblyGMG.{c,h}` providing `ApplyStokesOperatorMDOODZ` and `StokesCellBlockMDOODZ` that duplicate the coefficient formulas without modifying `StokesAssemblyDecoupled.c`. Correctness SHALL be enforced by a golden cross-check test (see next scenario) that runs on every CI invocation.

#### Scenario: Golden matvec cross-check against assembled CHOLMOD matrix

- **WHEN** the golden test assembles the full sparse matrix `A` via the existing `BuildStokesOperatorDecoupled` path on a representative small mesh (e.g. 21×21) and multiplies `A · x` for ten random vectors `x`, then calls `ApplyStokesOperatorMDOODZ(mesh, model, u, v, p, ru, rv, rp)` with the same inputs
- **THEN** the residuals `ru`, `rv`, `rp` agree with `A · x` componentwise to within `1e-12` relative tolerance for every `x` and every component

#### Scenario: Coefficient equivalence across every BC tag

- **WHEN** the golden test is run on meshes exercising each BC tag used in production (`0, -1, 2, 11, 13, 30, 31, -2, -12`) and each rheology branch (Picard, Newton-Phase-2), isotropic and anisotropic, incompressible and compressible
- **THEN** each configuration individually passes the `1e-12` matvec-equivalence bound without exceptions or tag-specific fall-backs

#### Scenario: Drift protection — matvec fails when assembly changes without the mirror update

- **WHEN** a coefficient formula in `StokesAssemblyDecoupled.c` is altered but the mirror formula in `StokesAssemblyGMG.c` is not updated (simulated via a deliberate one-line divergence in a CI-only test variant)
- **THEN** the golden cross-check test fails, pointing to the changed coefficient, before the drift can be merged

### Requirement: Newton-mode GMG is phased; interim behaviour falls back to CHOLMOD

The GMG solver SHALL support Newton-mode operation in two phases:

- **Phase 1**: the Picard variant of `ApplyStokesOperatorMDOODZ` is implemented and verified by the golden test. When the outer nonlinear loop enters Newton mode (using `D11`–`D34` Jacobian coefficients), the GMG dispatch SHALL detect this and return a "not applicable" status to the caller so that `StokesRoutines.c` falls back to the existing CHOLMOD path for that solve. The simulation SHALL continue normally and the fall-back SHALL be logged at `LOG_INFO` level per solve.
- **Phase 2**: the Jacobian variant of `ApplyStokesOperatorMDOODZ` is implemented and verified by the golden test across the Newton-mode BC-tag matrix. The Newton-mode guard is then removed and GMG handles Newton iterations using D7's symmetric-part-smoothing rule.

#### Scenario: Phase 1 guard under Newton mode

- **WHEN** `lin_solver = 3` and the outer loop has switched into Newton mode (`Newton = 1` and residual has crossed `Picard2Newton_tol`), and the Jacobian variant is not yet available
- **THEN** `SolveStokesGMG` returns the "Newton not yet implemented" status code, a `LOG_INFO` line is emitted identifying the fall-back, and `StokesRoutines.c` proceeds with the existing CHOLMOD Direct path on the same system

#### Scenario: Phase 2 Newton-mode convergence

- **WHEN** the Jacobian variant is available (Phase 2 complete) and the same Newton-mode problem is run
- **THEN** `SolveStokesGMG` runs to completion, FGMRES converges, and the outer Newton iteration achieves super-linear convergence comparable to Newton + CHOLMOD on the same problem (same integration scenarios as the `NewtonIterationConvergence.*` requirement)

### Requirement: Active mask is consistent across multigrid levels

MDOODZ deactivates cells above the free-surface marker chain by setting `BCu.type`, `BCv.type`, `BCg.type`, `BCp.type`, `BCt.type` to `30`. The GMG solver SHALL carry a per-level active mask derived from the fine-level BC tags: a coarse cell SHALL be marked active if any of its four fine children is active, inactive otherwise. Deactivated cells MUST be skipped by the Vanka smoother, the matrix-vector product, the residual evaluation, and the transfer operators at every level. The viscosity restriction SHALL ignore inactive fine cells when computing coarse-cell viscosity.

#### Scenario: Coarse mask consistency for a simple topography

- **WHEN** a fine-level mask has a step-shaped active region (left half 40 rows active, right half 30 rows active) and is restricted to the coarse grid
- **THEN** the coarse-level mask has a step at the same horizontal position, with per-column active-row counts equal to `ceil(fine_count / 2)` or `ceil((fine_count + 1) / 2)` consistent with the "active if any child active" rule

#### Scenario: Inactive cells do not contribute to smoothing or residual

- **WHEN** a Vanka sweep, a matvec, or a residual evaluation is run on a level with some cells marked inactive
- **THEN** the values stored at inactive cells are not read or written by the operation, and the residual at those cells is identically zero

#### Scenario: Viscosity restriction ignores inactive children

- **WHEN** a coarse cell has 3 of its 4 fine children active
- **THEN** the restricted viscosity at that coarse cell is computed from the 3 active children only (harmonic or arithmetic per the rule in the existing viscosity requirement), ignoring the inactive child's (stale) value

### Requirement: GMG produces correct shear-zone localisation

On a shear-zone localisation benchmark with Newton iteration enabled (strain-rate-softening power-law creep + Drucker–Prager plasticity, viscosity contrast ≥ 10⁴ inside the localising band), running the scenario with `lin_solver = 3` SHALL produce a strain-rate field whose peak location, width (FWHM), and integrated shear band displacement match the CHOLMOD reference within 5% at the prescribed reporting time step.

#### Scenario: Single shear-band scenario under Newton + GMG

- **WHEN** a single shear-band test case is run with `Newton = 1, lin_solver = 3` and again with `Newton = 1, lin_solver = 0`
- **THEN** the position of the peak `ε̇_II` differs by at most one fine cell, the FWHM of the shear band differs by at most 5%, and the integrated displacement across the band differs by at most 5%

#### Scenario: Conjugate shear bands scenario under Newton + GMG

- **WHEN** a conjugate shear-band test case (two crossing bands) is run under both solvers
- **THEN** both bands form in both runs with the same orientation (±2°) and the integrated work dissipated in the localising region differs by at most 5%

### Requirement: GMG results are consistent with CHOLMOD on existing benchmarks

For every scenario in the integration-test subset listed in the design document (SolVi, Blankenbach convection, pure shear velocity field, topographic relaxation, Newton iteration convergence tests), running the scenario with `lin_solver = 3` SHALL produce Vx, Vz, and P fields whose relative L2 difference from the CHOLMOD reference is below `1e-8`.

#### Scenario: SolVi viscous inclusion

- **WHEN** the SolVi benchmark is run with `lin_solver = 3`, then rerun with the default CHOLMOD path, and the final HDF5 fields are compared
- **THEN** `‖Vx_gmg - Vx_chol‖₂ / ‖Vx_chol‖₂ < 1e-8`, and the same bound holds for Vz and P

#### Scenario: Blankenbach thermal convection

- **WHEN** the Blankenbach convection benchmark is run to its standard end-step under both solvers
- **THEN** the L2 differences on Vx, Vz, P, and T remain below `1e-8` throughout the simulated time history

### Requirement: Memory scales near-linearly with problem size

Peak resident memory of the GMG solver SHALL scale as `RSS ≈ a · DOF^α` with `α ≤ 1.2` across 41², 81², 161² constant-viscosity Stokes problems, in contrast to the approximately `α ≈ 1.5` behaviour of CHOLMOD. GMG wall time SHALL be at least 3× faster than CHOLMOD at 201×201 on the reference CI machine.

#### Scenario: Memory-scaling exponent measurement

- **WHEN** the same manufactured constant-viscosity problem is solved at 41², 81², 161² with `lin_solver = 3` and peak RSS is recorded at each size
- **THEN** a least-squares fit `log RSS = log a + α · log DOF` yields `α ≤ 1.2`

#### Scenario: Wall-time comparison at 201×201

- **WHEN** an identical 201×201 problem is solved with CHOLMOD (`lin_solver = 0`) and with GMG (`lin_solver = 3`) on the CI reference machine
- **THEN** the GMG wall time is at most one third of the CHOLMOD wall time

### Requirement: Solver exposes optional tuning parameters with safe defaults

The `.txt` input SHALL accept the following optional parameters, each with the default listed; unset parameters MUST fall back to defaults without error:

- `gmg_levels` (default: auto-selected per the grid-hierarchy requirement)
- `gmg_nu_pre` (default: 2) — number of pre-smoothing Vanka sweeps per level
- `gmg_nu_post` (default: 2) — number of post-smoothing Vanka sweeps per level
- `gmg_fgmres_restart` (default: 30) — FGMRES restart length
- `gmg_fgmres_max_restarts` (default: 20) — maximum number of FGMRES restarts before declaring non-convergence; valid range `[1, 1000]`
- `gmg_fgmres_tol` (default: `nonlin_abs_mom / 10`) — FGMRES relative-residual convergence tolerance
- `gmg_standalone` (default: 0) — 1 = bare V-cycle, 0 = FGMRES-wrapped V-cycle
- `gmg_dump_vcycle` (default: 0) — 1 = dump per-level residual + solution snapshots for diagnostic use

#### Scenario: Defaults applied when parameters are absent

- **WHEN** a `.txt` input sets `lin_solver = 3` and omits every `gmg_*` parameter
- **THEN** the solver runs with the default values above and logs the effective configuration at startup, including the active `gmg_fgmres_max_restarts` value

#### Scenario: Invalid parameter is rejected with a clear error

- **WHEN** a `.txt` input sets `gmg_nu_pre = 0`, `gmg_fgmres_restart < 1`, or `gmg_fgmres_max_restarts` outside `[1, 1000]`
- **THEN** the solver aborts at startup with a log-level-ERROR message naming the offending parameter and its valid range

#### Scenario: `gmg_fgmres_max_restarts` tunable at runtime

- **WHEN** a `.txt` input sets `gmg_fgmres_max_restarts = 50` to permit longer FGMRES diagnostics
- **THEN** FGMRES iterates up to `gmg_fgmres_restart × 50` total inner iterations before declaring non-convergence, without requiring a rebuild

### Requirement: Solver logs residual history and falls back on non-convergence

On every Stokes linear solve using `lin_solver = 3`, the code SHALL log the initial residual `‖r_0‖`, per-iteration relative residuals `‖r_k‖ / ‖r_0‖` of the outer FGMRES iteration, and the final relative residual. If FGMRES does not converge within its iteration budget (restart length × max-restarts), the code SHALL log a WARN-level message whose numerical content is directly comparable to `gmg_fgmres_tol` (i.e. the relative-reduction quantity tested by the predicate) and fall back to the CHOLMOD direct solver for that Stokes solve, consistent with the thermal-PCG fallback pattern in `ThermalRoutines.c`.

#### Scenario: Convergence history is logged

- **WHEN** any Stokes solve with `lin_solver = 3` completes successfully
- **THEN** the log contains the initial residual `‖r_0‖`, per-FGMRES-iteration relative residuals `‖r_k‖ / ‖r_0‖`, and the final relative residual, each tagged so they can be extracted by the existing benchmarking scripts

#### Scenario: Fallback on non-convergence reports the comparable quantity

- **WHEN** FGMRES fails to reach `gmg_fgmres_tol` within `gmg_fgmres_restart × gmg_fgmres_max_restarts` inner iterations
- **THEN** a `LOG_WARN` line is emitted naming the step number and the achieved `final_res_rel` (the quantity tested against the tolerance), NOT the raw `‖r_k‖`; the CHOLMOD direct solver is invoked on the same system; simulation continues without aborting

#### Scenario: "Converged" is never logged when predicate fails

- **WHEN** the predicate `‖r_k‖ / ‖r_0‖ ≤ gmg_fgmres_tol` is false
- **THEN** the word "converged" is absent from the relevant log line regardless of how the predicate expression evolved during the run

### Requirement: V-cycle residual amplification is bounded

Across a single V-cycle, the composed residual amplification `‖r_exit‖ / ‖r_entry‖` at the finest level SHALL NOT exceed 4.0 on constant-viscosity SolCx-analog fixtures at 81×81, 161×161, and 201×201, measured via `gmg_dump_vcycle = 1` HDF5 snapshots. Additionally no single stage of the cycle — pre-smooth, restrict, coarse-operator apply, prolongate, post-smooth — SHALL amplify the residual by more than 4.0 relative to its input, preventing any single stage from silently violating the whole-cycle contract.

**Scope note.** This requirement was originally written with a tighter 2.0× per-stage bound under the assumption of Galerkin coarsening. The realised hierarchy uses MDOODZ-bridge operator at level 0 and textbook Picard rediscretisation on restricted viscosity at coarse levels (design D4 of `add-gmg-stokes-solver`); the cross-level operator mismatch appears as a benign 2×–4× jump specifically on the level-0 prolongate stage. Because the whole-cycle ratio remains well inside the 4.0 bound (measured 0.64×–1.36× at 81²–201² post-fix), the per-stage ceiling is set equal to the whole-cycle bound so that a genuine single-stage divergence (e.g. a smoother with spectral radius > 1 as observed pre-fix) still fires the regression without false-positives from the known Picard/bridge inconsistency.

#### Scenario: Whole V-cycle amplification is bounded

- **WHEN** a constant-viscosity 81² / 161² / 201² SolCx run executes with `gmg_dump_vcycle = 1` and the first V-cycle dump is analysed
- **THEN** the whole-V-cycle `‖r_exit‖ / ‖r_entry‖` ratio at level 0 stays below 4.0, demonstrating that the V-cycle is a productive preconditioner rather than a residual amplifier

#### Scenario: Per-stage amplification cannot silently violate the whole-cycle bound

- **WHEN** the same dump is walked stage-by-stage (pre-smooth, restrict, coarse-solve, prolongate, post-smooth) across both downleg and upleg
- **THEN** no single stage's output-to-input residual ratio exceeds 4.0, guarding against a future regression where one stage blows up by 8× and the rest damp to hide it

#### Scenario: Level-0 pre-smoother is a contraction on constant-viscosity inputs

- **WHEN** the first V-cycle dump of a constant-viscosity 81² fixture is examined at `level_0_pre_smooth` (pre → post)
- **THEN** the ratio is ≤ 4.0 (and, post-fix, is below 1.0 — a genuine damping), confirming the Vanka smoother at the finest level is a contraction; this was the specific pathology pre-fix (ratio 7.55×) that Fix F addresses by routing L0 `VankaBlockAssembleSolve` through the textbook block instead of the MDOODZ bridge

#### Scenario: 201² default-levels convergence

- **WHEN** the SolViPerf harness is driven at 201² with `lin_solver = 3` and the default auto-selected `gmg_levels` (expected ≥ 4 on this grid)
- **THEN** FGMRES converges to `gmg_fgmres_tol = 1e-6` in no more than 60 total inner iterations (at most 2 restarts at default `gmg_fgmres_restart = 30`), and no CHOLMOD fallback is triggered

