## Context

CHOLMOD is initialized in three places: Main_DOODZ.c (Stokes solver), ThermalRoutines.c (thermal solver), and ChemicalRoutines.c (chemical solver). All three call `cholmod_start(&c)` but none set `c.nthreads_max`. The `DirectSolver` struct has an unused `n_th` field. The `num_threads` variable is already captured from `omp_get_num_threads()` at startup (Main_DOODZ.c line 120) but never passed to CHOLMOD.

Benchmark data shows `thermal_s` constant at ~6.5s across 1-16 threads (1000×800 grid), which is 18% of wall time at optimal thread count (6t). The Stokes solve (`solve_s`) shows 0.0 in the current RiftingComprehensive scenario because there are no nonlinear iterations with solve calls, but other scenarios with actual solves would benefit.

## Goals / Non-Goals

**Goals:**
- Enable CHOLMOD multi-threaded sparse BLAS by setting `c.nthreads_max` after initialization
- Controllable via `.txt` parameter `cholmod_threads` (default 1 = single-threaded, backward compatible)
- Apply to all three CHOLMOD initialization points (Stokes, thermal, chemical)
- Benchmark on EC2 to measure thermal solver speedup
- Enable HDF5 writer in EC2 highres benchmarks to measure `output_s`

**Non-Goals:**
- Replacing CHOLMOD with MUMPS or another parallel solver (separate, larger change)
- Modifying CHOLMOD's ordering strategy (AMD is fine)
- Changing the `DirectSolver` struct's `n_th` field semantics (legacy field, leave as-is)

## Decisions

### D1: Thread count parameter design

**Decision:** New `.txt` parameter `cholmod_threads` with three modes:
- `1` (default): Single-threaded CHOLMOD — backward compatible, same as current behavior
- `-1`: Use all available OpenMP threads (`omp_get_num_threads()`)
- Any positive integer N: Use exactly N threads for CHOLMOD

**Rationale:** Default of 1 ensures no existing simulation changes behavior. Using `-1` as the "auto" value avoids collision with valid thread counts. Explicit count allows tuning (e.g., using 4 CHOLMOD threads while running 16 OpenMP threads, since sparse BLAS may not scale linearly).

**Alternatives considered:**
- *Boolean on/off:* Too coarse — doesn't allow tuning CHOLMOD threads independently of OpenMP threads.
- *Default to all threads:* Risky — could cause contention if CHOLMOD BLAS threads compete with OpenMP parallel-for threads in other subsystems within the same timestep.

### D2: Where to resolve thread count

**Decision:** Resolve the actual CHOLMOD thread count once at the parameter-read stage and store it in `model.cholmod_threads_resolved`. This avoids repeated `omp_get_num_threads()` calls at each solver invocation.

```c
// In InputOutput.c after ReadInt2:
if (model.cholmod_threads == -1) {
    model.cholmod_threads_resolved = num_threads; // from OMP
} else {
    model.cholmod_threads_resolved = model.cholmod_threads;
}
```

Then at each `cholmod_start()` call site: `c.nthreads_max = model.cholmod_threads_resolved;`

Actually, simpler: resolve at the `cholmod_start()` call sites since `num_threads` is available there. No need for a second struct field.

**Revised decision:** Keep it simple — one field `cholmod_threads` in the model struct. Resolve at each `cholmod_start()` site:
```c
cholmod_start(&c);
c.nthreads_max = (model.cholmod_threads == -1) ? num_threads : model.cholmod_threads;
```

### D3: Passing thread config to thermal/chemical solvers

**Decision:** The thermal solver (`EnergyDirectSolve` in ThermalRoutines.c) and chemical solver create their own local `cholmod_common c` inside the function. The `model` struct is already passed to these functions, so `model.cholmod_threads` is accessible. Set `c.nthreads_max` right after `cholmod_start(&c)`.

For the Stokes solver, `CholmodSolver.c` is a member of the `DirectSolver` struct initialized in Main_DOODZ.c where the model is available.

### D4: Benchmark configuration

**Decision:** Modify `benchmark-ec2-setup.sh` to set `writer = 1` and `writer_step = 1` in the highres `.txt` file. Run two sweeps:
1. `cholmod_threads = 1` (baseline, should match REPORT-20260413-113812.md)
2. `cholmod_threads = -1` (all threads)

Compare `thermal_s` and `solve_s` columns between sweeps.

## Risks / Trade-offs

**[BLAS contention with OpenMP]** → If CHOLMOD uses N BLAS threads while OpenMP parallel-for loops are also active, thread oversubscription occurs. → **Mitigation:** CHOLMOD runs in serial sections of the timestep (thermal solver is called outside the NL iteration's parallel-for). No contention with current code flow. The persistent-omp-thread-pool change (future) will need to account for this.

**[SuiteSparse version compatibility]** → `nthreads_max` field exists in all SuiteSparse versions (5.x through 7.x). → **Mitigation:** The field is part of `cholmod_common` in all versions. No version check needed.

**[Minimal gain if BLAS is already single-threaded]** → If the installed BLAS library (OpenBLAS, MKL, reference BLAS) doesn't support threading, setting `nthreads_max` has no effect. → **Mitigation:** EC2 Ubuntu has OpenBLAS which supports multi-threading via `OPENBLAS_NUM_THREADS`. CHOLMOD's `nthreads_max` controls its own threading, which leverages BLAS threading internally. Log the configured value at startup for diagnostics.

**[Potential numerical non-reproducibility]** → Multi-threaded BLAS may use different reduction ordering, producing slightly different results. → **Mitigation:** This is a known property of parallel BLAS. Differences are at machine epsilon level. The `.txt` parameter allows reverting to single-threaded for reproducibility.
