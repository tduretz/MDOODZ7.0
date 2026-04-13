## Why

`P2Mastah` (particle-to-grid interpolation) is the largest remaining bottleneck after the PCG thermal solver work. At 1000×800 on an 8-core machine it takes 3.5–7.1 s per step (20–50 % of wall time) and **anti-scales** beyond 4 threads (3.49 s at 4 threads → 7.06 s at 16 threads). Two root causes compound:

1. **Allocation churn**: Every call allocates ~51 MB of thread-local buffers (`Wm`, `BmWm`, `Wm_ph`) via `DoodzCalloc`, uses them in one of five separate fork-join regions, then immediately frees them — roughly 20 calls per timestep, totalling 1 GB+ of transient allocations and 100+ fork-join barriers per step.
2. **Thread-local buffer explosion**: The scatter loop maintains `nthreads` separate copies of each grid-sized array (`Wm[nthreads][Nx*Nz]`, `BmWm[nthreads][Nx*Nz]`, `Wm_ph[nthreads][Nb_phases][Nx*Nz]`), then runs a serial-style reduction loop to sum them. At 16 threads/8 phases this is 16 × (2 + 8) × 800 K doubles ≈ 100 MB — poor cache utilisation and O(nthreads × Nx × Nz) reduction work that grows linearly with thread count.

This change addresses both issues: persistent buffers eliminate the malloc/free overhead, and atomic scatter replaces the thread-local copies + reduction with direct accumulation into single global arrays, fixing the memory scaling and removing the reduction loop entirely.

## What Changes

- **Runtime switch** (`interp_mode` .txt parameter, int, default 0): Controls which interpolation strategy `P2Mastah` uses.
  - `0` — Legacy behaviour (thread-local buffers, per-call alloc/free). Fully backward compatible, bit-identical to current code.
  - `1` — Persistent buffers only. Same thread-local scatter + reduction logic, but buffers are allocated once and zeroed with `memset` before each reuse.
  - `2` — Persistent buffers + atomic scatter. Single global arrays with `#pragma omp atomic` accumulation; reduction loop removed.
- **Persistent buffers** (modes 1 and 2): Pre-allocate interpolation work arrays (`WM`, `BMWM`, and when `prop==1` the phase-percentage arrays) once at simulation startup. Replace per-call `DoodzCalloc`/`DoodzFree` with `memset`-zero + reuse.
- **Atomic scatter** (mode 2 only): Replace the `nthreads` thread-local `Wm[thread_num][]` / `BmWm[thread_num][]` / `Wm_ph[thread_num][][]` arrays with single global `WM[]` / `BMWM[]` arrays. Use `#pragma omp atomic` on the `+=` accumulations in the particle scatter loop. Remove the separate reduction loop that sums thread-local arrays into the global arrays.
- **Lifecycle API**: Add an `InterpBufPool` struct with init/free functions, called during simulation startup/cleanup. The pool holds pre-allocated buffers sized to `Nx*Nz`. Allocated when `interp_mode >= 1`; skipped entirely for mode 0.
- **Phase-percentage handling** (mode 2): For `prop==1` calls, `mesh->phase_perc_s/n` arrays are already global — apply atomics directly to them instead of accumulating into thread-local `Wm_ph` copies.
- **CI test coverage**: Add `interp_mode = 1` and `interp_mode = 2` CI test variants for all existing test suites except BlankenBench (too long-running). This covers: AnisotropyBenchmark, BoundaryCondition, Compressibility, ConvergenceRate, Density, FiniteStrain, FreeSurface, NeumannBC, ShearTemplate, Plasticity, RheologyCreep, ShearHeating, SolViBenchmark, SolViMarkerComparison, SolverMode, Thermal, TopoBench, VelocityField, ViscoElastic. Each variant duplicates the existing `.txt` parameter file with `interp_mode` set and asserts results match the mode-0 baseline within tolerance (bit-identical for mode 1, machine-epsilon for mode 2).

## Capabilities

### New Capabilities
- `interp-buffer-pool`: Persistent buffer pool and atomic-scatter interpolation for `P2Mastah` — allocation, zeroing, reuse, atomic accumulation, and cleanup API.

### Modified Capabilities
<!-- No existing spec-level behaviour changes. The interpolation results remain bit-identical;
     only the memory management strategy changes (persistent allocation + atomic scatter
     instead of thread-local copies + manual reduction). -->

## Impact

- **Code**: `MDLIB/ParticleRoutines.c` (`P2Mastah` — mode dispatch, scatter loop, reduction loop, cleanup section), new `InterpBufPool` struct + init/free functions, `InputOutput.c` (read `interp_mode`), `mdoodz.h` (new param field), callers in `Main_DOODZ.c` for pool lifecycle.
- **CI**: New GTest variants for `interp_mode = 1` and `interp_mode = 2` across 19 test suites (all except BlankenBench). Duplicate `.txt` param files with mode override. ~38 new test cases total.
- **Backward compatibility**: Default `interp_mode = 0` preserves the exact existing code path — no allocation changes, no atomics, bit-identical output. Existing `.txt` files require no modification.
- **Memory**: Mode 0 unchanged. Mode 1 keeps thread-local buffers alive across calls (same peak RSS, fewer allocator calls). Mode 2 **decreases** peak RSS substantially — from `nthreads × (2 + Nb_phases) × Nx*Nz` doubles (~100 MB at 16 threads) to 2 × `Nx*Nz` doubles persistent (~12.8 MB at 1000×800).
- **Performance**: Mode 1 — expected 20–30 % interp reduction from eliminating malloc/free churn. Mode 2 — additional improvement from removing the O(nthreads × Nx × Nz) reduction loop; should fix anti-scaling past 4 threads since memory footprint no longer grows with thread count.
- **Correctness**: Mode 0 and 1 are bit-identical to current code. Mode 2 matches to machine epsilon (atomic `+=` order may differ from thread-local reduction order).
- **Dependencies**: None. No new external libraries. Requires OpenMP 3.1+ (already required by the codebase).
