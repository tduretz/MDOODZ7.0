## Context

The thermal solver in MDOODZ uses a CHOLMOD Cholesky direct solver for the energy equation $\rho C_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + H$. The discretisation produces an SPD matrix with a 5-point stencil in CSR format. At 1001×801 resolution (~800K DOFs, ~4M non-zeros), the serial CHOLMOD factorisation costs 6.2–6.5s per step — 38% of wall time at 8 threads.

The matrix assembly already exists in `EnergyDirectSolve()` (ThermalRoutines.c) and produces CSR arrays `(A[], Ic[], J[])`. Currently these flow through a multi-stage conversion (CSR → triplet → CXSparse CSC → CHOLMOD sparse → CHOLMOD factor). The PCG solver will consume the CSR arrays directly, bypassing all format conversions.

## Goals / Non-Goals

**Goals:**

- Replace the serial CHOLMOD factorise+solve with an OpenMP-parallel PCG solver for the thermal system
- Achieve 10–20× speedup on the thermal portion at 8 threads (6.5s → 0.3–0.7s)
- Maintain bit-for-bit identical results when using the direct solver (backward compat)
- Zero new external dependencies — pure C + OpenMP
- Automatic fallback to CHOLMOD if PCG fails to converge

**Non-Goals:**

- Replacing the Stokes solver (different matrix, non-symmetric)
- Implementing multigrid or IC(0) preconditioner (future Phase 2/3)
- Supporting polar mode in the PCG path (rare; falls back to direct)
- Optimising the matrix assembly itself (already fast, ~1.7s)

## Decisions

### D1: Operate on CSR directly — skip all format conversions

**Decision**: The PCG solver reads the existing CSR arrays `(A[], Ic[], J[])` assembled by `EnergyDirectSolve`. No conversion to CXSparse or CHOLMOD format.

**Rationale**: CSR is the native SpMV format. The current conversion pipeline (CSR → triplet → CSC → CHOLMOD) exists only because CHOLMOD requires its own format. SpMV on CSR is cache-friendly (sequential row access) and trivially parallel. Eliminating the conversion also saves ~0.3s of allocation/copy overhead per step.

**Alternative considered**: Convert to CSC for column-oriented SpMV. Rejected — CSR row-wise SpMV is standard for CG and has better cache behaviour for the output vector.

### D2: Jacobi (diagonal) preconditioner

**Decision**: Use diagonal preconditioning $M^{-1} = \text{diag}(A)^{-1}$. Extract the diagonal during matrix assembly (it's the central stencil coefficient: $\rho C_p / \Delta t + \sum k_{neighbours} \cdot \text{FD\_coefs}$).

**Rationale**: For the thermal matrix with condition number $\kappa \approx 10^2$–$10^4$, Jacobi preconditioning reduces iterations from ~2000 (unpreconditioned) to ~100–300. Each iteration costs ~0.5ms at 800K DOFs on 8 threads, so 200 iterations ≈ 100ms total. The preconditioner is free to apply (one element-wise multiply) and free to construct (diagonal already computed during assembly).

**Alternative considered**: IC(0) incomplete Cholesky — stronger (50–100 iterations) but serial to apply and moderately complex to implement. Deferred to Phase 2 if Jacobi proves insufficient.

### D3: Dispatch inside EnergyDirectSolve, not in the caller

**Decision**: Add the solver selection branch inside `EnergyDirectSolve()` after matrix assembly. When `model.thermal_solver == 1`, call `SolveThermalPCG()` instead of `FactorEnergyCHOLMOD()` + `SolveEnergyCHOLMOD()`.

**Rationale**: The matrix assembly code is shared between both paths. The caller in `Main_DOODZ.c` does not need to know which solver is used — it already passes `&ThermalSolver` which carries the solver config. This keeps the API change minimal.

**Alternative considered**: Separate `EnergyIterativeSolve()` function. Rejected — would duplicate 300 lines of matrix assembly. The assembly is the expensive shared part; only the solve differs.

### D4: New function SolveThermalPCG in ThermalSolver.c

**Decision**: Add a single function:
```c
int SolveThermalPCG(double *A, int *Ic, int *J, double *b, double *x,
                    int neq, int max_its, double rel_tol);
```
Returns: iteration count (positive) on success, negative on failure (triggers fallback).

**Rationale**: Self-contained PCG implementation (~60 lines). Takes CSR matrix and RHS, writes solution to `x[]`. No external state, no CHOLMOD dependency. Easy to test in isolation.

### D5: Warm-start from previous timestep solution

**Decision**: Initialise $x_0$ from `mesh->T[]` (the temperature from the previous timestep) rather than $x_0 = 0$.

**Rationale**: Between timesteps, temperature changes by $O(\Delta t \cdot \dot{T})$ which is small. The initial residual $\|b - Ax_0\|$ starts near zero, reducing iteration count by ~30–50%. The current code already stores `mesh->T[]` from the previous step, so no extra storage needed.

### D6: Fallback to CHOLMOD on convergence failure

**Decision**: If PCG does not converge within `max_its_thermal` iterations, log a warning and call the existing CHOLMOD path (factorise + solve). The `DirectSolver` struct and CHOLMOD context are always initialised regardless of solver choice.

**Rationale**: Guarantees the simulation never crashes due to iterative solver failure. The CHOLMOD path is slow but unconditionally correct. This is critical for production runs where robustness is non-negotiable.

### D7: Skip PCG for polar mode

**Decision**: When `model.polar == 1`, always use the CHOLMOD path. The PCG path handles only the standard Cartesian case.

**Rationale**: Polar mode forms $A \cdot A^T$ (normal equations) which changes the matrix structure and is not SPD in the same sense. Supporting it in PCG would require forming the product matrix explicitly, negating the performance benefit. Polar mode is rarely used (only 2 of 79 scenarios).

### D8: Three new .txt parameters with safe defaults

| Parameter | Key | Type | Default | Meaning |
|-----------|-----|------|---------|---------|
| Solver mode | `thermal_solver` | int | `0` | 0 = direct (CHOLMOD), 1 = iterative (PCG) |
| Max iterations | `max_its_thermal` | int | `1000` | PCG iteration cap |
| Relative tolerance | `rel_tol_thermal` | double | `1e-8` | PCG stops when $\|r\| / \|b\| < \text{tol}$ |

Default `thermal_solver = 0` ensures all existing `.txt` files continue to use CHOLMOD without modification.

## Risks / Trade-offs

**[Risk] High conductivity contrast degrades convergence** → Mitigation: Jacobi preconditioning handles contrasts up to ~$10^3$. For extreme contrasts ($>10^4$, e.g. melt pockets), iteration count may exceed 500. The `max_its_thermal` cap + CHOLMOD fallback ensures the simulation continues. Log a warning when iteration count exceeds 50% of the cap so users know to consider direct mode.

**[Risk] PCG solution differs from CHOLMOD at machine precision** → Mitigation: This is expected and acceptable. PCG converges to a relative residual of $10^{-8}$ which gives ~8 digits of agreement. The thermal field has physical uncertainty much larger than this. For validation, run existing `ci-thermal-tests` with both modes and verify agreement to $10^{-6}$ relative.

**[Risk] Free-surface topology changes invalidate warm start** → Mitigation: When cells change from air (type 30) to active (or vice versa), `neq` changes and the dimension check in `EnergyDirectSolve` already resets the solver state. The PCG warm start uses `mesh->T[]` which is always the correct size for the current topology.

**[Risk] Memory for PCG work vectors** → Mitigation: PCG needs 5 vectors of size `neq` (r, z, p, Ap, diagonal). At 800K DOFs this is ~30 MB — negligible compared to the 12 GB RSS. All vectors are allocated once and reused across timesteps via the `DirectSolver` struct, or allocated/freed per solve call (simpler, still fast since it's 5 `malloc` calls).

**[Trade-off] Jacobi is not the optimal preconditioner** → Accepted. It's the simplest correct choice. If benchmarks show >300 iterations consistently, Phase 2 adds IC(0) as a drop-in replacement (same PCG outer loop, different preconditioner application).

## Resolved Questions

1. **PCG work vectors**: Persistent in `DirectSolver` — allocate once, reuse across timesteps. Avoids repeated malloc/free and keeps memory stable.
2. **Iteration count in perf.csv**: Yes — add a `pcg_its` column. Valuable for monitoring convergence stability across timesteps and detecting scenarios that approach the fallback threshold.
