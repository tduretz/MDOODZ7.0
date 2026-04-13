## 1. Parameter Reading

- [x] 1.1 Add `thermal_solver` (int, default 0), `max_its_thermal` (int, default 1000), `rel_tol_thermal` (double, default 1e-8) fields to the `params` struct in `mdoodz.h` (linear solver section)
- [x] 1.2 Read the three new parameters in `InputOutput.c` using `ReadInt2`/`ReadDou2` alongside existing solver parameters

## 2. PCG Solver Core

- [x] 2.1 Implement `SolveThermalPCG()` in `ThermalSolver.c` — PCG loop operating on CSR arrays `(A[], Ic[], J[])` with Jacobi diagonal preconditioner
- [x] 2.2 Implement CSR SpMV kernel (`csr_matvec`) with `#pragma omp parallel for` in `ThermalSolver.c`
- [x] 2.3 Implement OpenMP-parallel vector operations: dot product, axpy, element-wise multiply (preconditioner apply)
- [x] 2.4 Extract diagonal from CSR matrix for Jacobi preconditioner (invert to get `M_inv[]`)

## 3. Solver Dispatch

- [x] 3.1 Add solver dispatch branch in `EnergyDirectSolve()` (ThermalRoutines.c): PCG on `it==0` when `model.thermal_solver == 1 && model.polar == 0`; CHOLMOD for `it>=1` with cached factorization
- [x] 3.2 `model` struct is already in scope — params read directly from `model.thermal_solver`, `model.max_its_thermal`, `model.rel_tol_thermal`
- [x] 3.3 Implement warm-start: copy `mesh->T[]` (mapped through `eqn_t[]`) into `x[]` before calling PCG
- [x] 3.4 Implement CHOLMOD fallback: if `SolveThermalPCG()` returns negative, log `LOG_WARN` and call existing CHOLMOD factorise+solve path

## 4. Header Declarations

- [x] 4.1 Add `SolveThermalPCG()` declaration to `mdoodz-private.h`

## 5. Testing

- [x] 5.1 Build and run locally with `thermal_solver = 0` — verify identical results (backward compat)
- [x] 5.2 Build and run locally with `thermal_solver = 1` — PCG converges (310 its cold start, 13 its warm-start) on RiftingChenin 150×100
- [ ] 5.3 Run existing `ci-thermal-tests` with both solver modes — verify agreement to 1e-6 relative
- [x] 5.4 Test fallback: set `max_its_thermal = 1` and verify CHOLMOD fallback triggers with `LOG_WARN`

## 6. Benchmarking

- [ ] 6.1 Run RiftingComprehensive at 1001×801 with `thermal_solver = 1`, threads 1/2/4/6/8 — compare `thermal_s` against CHOLMOD baseline
- [ ] 6.2 Verify PCG iteration count is stable across timesteps (log or print per-step count)
