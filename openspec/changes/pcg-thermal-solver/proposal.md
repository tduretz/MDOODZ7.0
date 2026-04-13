## Why

The CHOLMOD direct Cholesky solver for the thermal equation is the primary scaling bottleneck in MDOODZ. At 1001×801 resolution it consumes 6.2–6.5s per step (38% of wall time at 8 threads) and is entirely serial — adding threads cannot reduce it. This caps the overall parallel speedup at ~2.2× regardless of thread count. Replacing it with a Preconditioned Conjugate Gradient (PCG) solver converts the thermal solve from serial to parallel (SpMV is trivially OpenMP-parallel) and reduces the per-step cost by an estimated 10–20×.

## What Changes

- **New iterative thermal solver**: Implement PCG (Conjugate Gradient with Jacobi diagonal preconditioner) as an alternative to CHOLMOD for the energy equation.
- **Solver selection switch**: Add a `.txt` parameter (`thermal_solver`) to choose between `direct` (CHOLMOD, current default) and `iterative` (PCG). Default remains `direct` for backward compatibility.
- **Convergence monitoring**: PCG reports iteration count, residual norm, and elapsed time via the existing logging system. Falls back to direct CHOLMOD if PCG fails to converge within `max_its_thermal` iterations.
- **OpenMP parallelism**: The SpMV kernel and vector operations in PCG use `#pragma omp parallel for`, making the thermal solve scale with thread count for the first time.
- **No new dependencies**: PCG is implemented in pure C + OpenMP. No external libraries required.
- **CHOLMOD retained**: The direct solver remains available and is used as the fallback. The `DirectSolver` caching from the prior change continues to work when `thermal_solver = direct`.

## Capabilities

### New Capabilities

- `pcg-thermal-solve`: PCG iterative solver for the SPD thermal system — Jacobi preconditioner, OpenMP-parallel SpMV, convergence criteria, fallback to direct solver, `.txt` parameter integration.

### Modified Capabilities

_(none — the direct solver path is unchanged; the new solver is additive)_

## Impact

- **Code**: `MDLIB/ThermalSolver.c` (new PCG function), `MDLIB/ThermalRoutines.c` (solver dispatch), `MDLIB/mdoodz-private.h` (declarations), `MDLIB/InputOutput.c` (new `.txt` parameters).
- **Parameters**: New `.txt` keys: `thermal_solver` (0=direct, 1=iterative), `max_its_thermal` (default 1000), `rel_tol_thermal` (default 1e-8).
- **API**: `EnergyDirectSolve` signature gains a solver-type flag or the dispatch is internal. `FactorEnergyCHOLMOD` is unchanged.
- **Dependencies**: None added. SuiteSparse remains required (for Stokes + fallback).
- **Testing**: Existing `ci-thermal-tests` must pass with both solver modes. Add a PCG-specific convergence test.
- **Risk**: Low. The thermal matrix is SPD with typical condition number 10²–10⁴. PCG with Jacobi preconditioner converges reliably for standard geodynamic conductivity contrasts. Extreme k contrasts (>10⁴) may need more iterations but the fallback ensures correctness.
