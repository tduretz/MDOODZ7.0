## Context

`P2Mastah` in `MDLIB/ParticleRoutines.c` (lines 1984–2335) is the sole particle-to-grid interpolation function. It is called 20–30+ times per timestep from `Main_DOODZ.c` and `CountPartCell`, operating on four different centroid types (`cent=1`: Nx-1 × Nz-1, `vert=0`: Nx × Nz, `vxnodes=-1`: Nx × Nz+1, `vznodes=-2`: Nx+1 × Nz). Each call currently:

1. Allocates `nthreads` copies of `Wm[Nx*Nz]` and `BmWm[Nx*Nz]` (+ `Wm_ph[Nb_phases][Nx*Nz]` when `prop==1`)
2. Runs a parallel scatter loop writing to `Wm[thread_num][kp]`, `BmWm[thread_num][kp]`
3. Runs a parallel reduction loop summing thread-local arrays into global `WM[]`, `BMWM[]`
4. Runs a parallel normalization loop computing final node values
5. Frees all buffers

At 1000×800 / 8 threads / 8 phases, each call allocates ~51 MB and the function anti-scales beyond 4 threads due to buffer explosion and reduction overhead.

The `interp_mode` parameter (default 0) selects between three strategies with increasing optimization and decreasing strict bit-reproducibility.

## Goals / Non-Goals

**Goals:**
- Eliminate per-call malloc/free overhead via persistent buffer pool (modes 1, 2)
- Eliminate thread-local buffer explosion and O(nthreads × N) reduction via atomic scatter (mode 2)
- Fix anti-scaling of `interp` time beyond 4 threads
- Maintain full backward compatibility with `interp_mode = 0` (default)
- CI coverage for modes 1 and 2 across all test suites except BlankenBench
- EC2 benchmark verification: mode 2 first (expect largest gains), then mode 1

**Non-Goals:**
- Cell-sorted particle reordering (future optimization, separate change)
- Merging multiple P2Mastah calls into batched multi-field interpolation
- Changing the interpolation stencil (1-cell vs 9-cell) or weighting scheme
- GPU offloading

## Decisions

### D1: Three-level `interp_mode` switch

| Mode | Strategy | Buffers | Scatter | Reduction | Bit-identical? |
|------|----------|---------|---------|-----------|----------------|
| 0 | Legacy | Per-call alloc/free | Thread-local | Explicit loop | Yes (unchanged path) |
| 1 | Persistent | Pool, memset-zero | Thread-local | Explicit loop | Yes |
| 2 | Persistent + atomics | Pool, memset-zero | `omp atomic` to global | Eliminated | Machine-epsilon |

**Rationale**: Mode 0 preserves the exact existing code path for users who need bit-reproducibility or want zero risk. Mode 1 isolates the allocation benefit for measurement. Mode 2 combines both optimizations. This mirrors the `thermal_solver` pattern already established in the codebase.

**Alternative considered**: Single switch (0/1 only) combining both strategies. Rejected because isolating the persistent-buffer benefit from the atomics benefit is important for benchmarking and diagnosing any issues. The persistent-buffer-only path (mode 1) is also bit-identical, making it a safer default upgrade path.

### D2: InterpBufPool struct with 4-centroid pre-allocation

The pool pre-allocates buffers for all four centroid grid sizes at init time:

```c
typedef struct {
    // Per-centroid-type persistent arrays (indexed 0..3 for cent, vert, vx, vz)
    double *WM[4];     // weight accumulator [Nx_c * Nz_c]
    double *BMWM[4];   // weighted-value accumulator [Nx_c * Nz_c]
    int     sizes[4];  // Nx_c * Nz_c for each centroid type
    // Mode 1 only: thread-local arrays
    double **Wm[4];    // [nthreads][Nx_c * Nz_c]  (NULL for mode 2)
    double **BmWm[4];  // [nthreads][Nx_c * Nz_c]  (NULL for mode 2)
    int     nthreads;
    int     interp_mode;
} InterpBufPool;
```

**Rationale**: P2Mastah calls use 4 distinct grid sizes depending on `centroid` (-2, -1, 0, 1). Pre-allocating all 4 avoids any per-call size lookup or conditional allocation. The sizes are known at grid init and never change during a simulation.

**Alternative considered**: Single max-size allocation reused for all centroid types. Rejected because the size differences are small (Nx±1, Nz±1) and using exact sizes avoids any out-of-bounds risk.

**Alternative considered**: Allocating `Wm_ph` thread-local arrays in the pool for mode 1 `prop==1` calls. Rejected — `prop==1` is only called from `CountPartCell` (once per step), making its buffer cost negligible. For mode 1, `Wm_ph` is allocated/freed per-call as today. For mode 2, `Wm_ph` is eliminated entirely (atomics write directly to `mesh->phase_perc_s/n`).

### D3: Centroid index mapping

Map `centroid` parameter to pool index at the top of `P2Mastah`:

```c
int pool_idx;
switch (centroid) {
    case  1: pool_idx = 0; break;  // cell centres
    case  0: pool_idx = 1; break;  // vertices
    case -1: pool_idx = 2; break;  // Vx nodes
    case -2: pool_idx = 3; break;  // Vz nodes
}
```

The pool is passed as `InterpBufPool *pool` — a new parameter added to the `P2Mastah` function signature. When `interp_mode == 0`, `pool` is `NULL` and the function uses the legacy code path unchanged.

### D4: Atomic scatter implementation (mode 2)

Replace thread-local accumulation:
```c
// Current (modes 0, 1):
Wm[thread_num][kp]   += weight;
BmWm[thread_num][kp] += mark_val * weight;

// Mode 2:
#pragma omp atomic
WM[kp] += weight;
#pragma omp atomic
BMWM[kp] += mark_val * weight;
```

For `prop==1` (phase percentages), mode 2 applies atomics directly to `mesh->phase_perc_s[p][kp]` or `mesh->phase_perc_n[p][kp]`, and the `Wm_ph` arrays are not allocated at all.

The reduction loop (lines 2258–2271) is skipped entirely in mode 2 since `WM`/`BMWM` already contain the final accumulated values.

**Rationale**: Particles from different threads rarely map to the same grid cell simultaneously (particles are spatially distributed), so atomic contention is low in practice. The cost of the atomics is far less than the cost of allocating and reducing `nthreads` full-grid copies.

**Risk**: High particle density in single cells (e.g., early simulation with clustered markers) could increase contention. The 9-cell stencil spreads writes across 9 nodes, further reducing contention.

### D5: Pool lifecycle in Main_DOODZ.c

```c
// After grid init, before first P2Mastah call:
InterpBufPool *pool = NULL;
if (input.model.interp_mode >= 1) {
    pool = InterpBufPoolInit(&mesh, input.model.interp_mode, omp_get_max_threads());
}

// ... all P2Mastah calls pass pool ...

// At simulation cleanup:
InterpBufPoolFree(pool);
```

The pool pointer is passed through to `P2Mastah` and `CountPartCell`. When `interp_mode == 0`, `pool` remains `NULL` and `P2Mastah` follows the legacy code path with no overhead.

### D6: P2Mastah signature change

```c
// Before:
void P2Mastah(params *model, markers particles, DoodzFP *mat_prop, grid *mesh,
              double *NodeField, char *NodeType, int flag, int avg,
              int prop, int centroid, int interp_stencil);

// After:
void P2Mastah(params *model, markers particles, DoodzFP *mat_prop, grid *mesh,
              double *NodeField, char *NodeType, int flag, int avg,
              int prop, int centroid, int interp_stencil,
              InterpBufPool *pool);
```

All ~30 call sites in `Main_DOODZ.c` and `CountPartCell` are updated to pass `pool`. The `pool` parameter is the last argument to minimize diff noise.

### D7: CI test strategy

For each of the 19 test suites (all except BlankenBench), add two new test cases:
- `<TestName>InterpMode1` — run with `interp_mode = 1`, assert results match mode-0 baseline **exactly** (bit-identical)
- `<TestName>InterpMode2` — run with `interp_mode = 1`, assert results match mode-0 baseline within relaxed tolerance (machine-epsilon level, e.g. `1e-10` relative)

Implementation: duplicate the existing `.txt` parameter file for each test with `interp_mode` appended (e.g., `GaussianDiffusionInterpMode1.txt`). The GTest code reuses the same setup/assertion pattern as the existing test, just pointing to the variant `.txt` file.

### D8: EC2 benchmark order

Run mode 2 first on EC2 (the full persistent+atomic optimization — expected to show the largest improvement and best thread scaling), then mode 1 (persistent buffers only). This follows the user's preference to see the "good news" first and is also practically useful: if mode 2 shows strong gains, mode 1 results provide the decomposition of how much comes from allocation elimination vs. atomics.

Benchmark configuration (same as Experiment 4 for comparability):
- Scenario: RiftingComprehensive, 1001×801, 10 steps
- Threads: 1, 2, 4, 6, 8, 12, 16
- `--sed 's/^interp_mode.*/interp_mode = 2/'` (then repeat with `= 1`)

## Risks / Trade-offs

**[Atomic contention under extreme particle clustering]** → In typical MDOODZ simulations, particles are well-distributed across cells. If a pathological setup concentrates many particles in a few cells, atomic contention in mode 2 could degrade performance. Mitigation: mode 1 (persistent buffers without atomics) remains available as a fallback; benchmarking on realistic scenarios will validate this.

**[Floating-point non-reproducibility in mode 2]** → Atomic `+=` accumulates in arrival order, which differs across runs and thread counts. Results match to machine epsilon but not bit-for-bit. Mitigation: mode 0 and 1 preserve exact reproducibility; mode 2 is opt-in.

**[Function signature change touches ~30 call sites]** → Adding the `pool` parameter requires updating every P2Mastah call. Mitigation: mechanical change (append `pool` to each call), low risk of logic errors, CI tests catch any missed sites at compile time.

**[Pool memory persists for simulation lifetime]** → Modes 1 and 2 keep buffers allocated throughout the simulation. For mode 1, this is `nthreads × 2 × max_grid_size × 8` bytes (~25 MB at 8 threads / 1000×800). For mode 2, it's `2 × 4 × max_grid_size × 8` bytes (~51 MB for all 4 centroid types). Mitigation: this memory was already being allocated and freed 20+ times per step, so peak RSS is actually unchanged (mode 1) or reduced (mode 2).

## Open Questions

- **Q1**: Should `interp_mode` be added to all SETS/*.txt files (like `thermal_solver = 0`)? Or rely on the `ReadInt2` default? For consistency with thermal_solver precedent, we should add it.
- **Q2**: Should CountPartCell also receive the pool? It calls P2Mastah with `prop==1`, which has different buffer needs (Wm_ph). For mode 2, atomics go directly to mesh arrays, so the pool is still useful. Decision: yes, pass pool to CountPartCell.
