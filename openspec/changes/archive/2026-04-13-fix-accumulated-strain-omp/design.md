## Context

`AccumulatedStrainII` in `RheologyParticles.c` interpolates strain-rate increments from cell centres to ~12M particles using bilinear weights. It runs once per timestep inside the `post_solve` phase. The function has two loops:

1. **Grid loop** (lines 213–230): Computes `strain_inc[c1]` for each cell centre. Already parallelized with `#pragma omp parallel for`.
2. **Particle loop** (lines 240–530): For each particle, reads 4 neighbouring `strain_inc` values, computes weighted sum, writes to `particles->strain*[k]`. Currently **single-threaded** — the `#pragma omp parallel for` is commented out.

The particle loop is embarrassingly parallel: each iteration reads shared grid arrays (read-only) and writes to `particles->strain*[k]` where `k` is the loop index (no overlap between iterations).

The original commented pragma is missing `dE_pl_vol` and `k` from its `private` clause, and the loop body contains `exit(0)` which is unsafe in a parallel region.

## Goals / Non-Goals

**Goals:**
- Re-enable OpenMP on the particle loop to parallelize `AccumulatedStrainII`
- Fix the `private` clause to include all local variables
- Replace the `exit(0)` with a thread-safe error handling approach
- Maintain bit-identical results (no floating-point reordering — each particle's computation is independent)

**Non-Goals:**
- Optimizing other functions in the `post_solve` timer (e.g., `MassSourceTerm`, `UpdateDensity`)
- Changing the interpolation algorithm or stencil
- Adding new parameters or configuration options

## Decisions

### 1. Uncomment and fix the existing pragma

**Choice**: Uncomment the existing `#pragma omp parallel for` and fix its variable lists rather than rewriting the loop.

**Rationale**: The original pragma was correct in structure — it just needs `dE_pl_vol`, `dE_el`, `k`, and `l` added to the `private` clause. The `shared`/`firstprivate` clauses are already correct. Minimal change = minimal risk.

**Fixed pragma**:
```c
#pragma omp parallel for shared(particles, strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pl_vol, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs, tag) \
    private(k, l, dE_tot, dE_el, dE_pl, dE_pl_vol, dE_pwl, dE_exp, dE_lin, dE_gbs, sumW, dst, dxm, dzm, iSW, iSE, iNW, iNE, i_part, j_part) \
    firstprivate(dx, dz, Nx, Nz)
```

### 2. Replace `exit(0)` with thread-safe error flag

**Choice**: Replace `exit(0)` on negative strain increment with `LOG_ERR` + continue (skip the particle). Do not abort.

**Rationale**: `exit()` in a parallel region causes undefined behaviour (other threads may be mid-write). The negative strain check is a diagnostic guard, not a correctness invariant — negative strain increments indicate a solver issue upstream, not a bug in this function. Logging the error and continuing is the safest approach. The first grid loop already validates `strain_inc[c1]` with the same check, so the particle loop check is redundant but harmless to keep as a warning.

**Alternative considered**: Use an atomic flag + `break` after the parallel region. Rejected — adds complexity for a condition that should never occur in practice.

### 3. Add `strain_inc_pl_vol` to shared clause

**Choice**: Add `strain_inc_pl_vol` to the `shared` clause — it was allocated before the loop but missing from the original pragma.

**Rationale**: The original commented pragma omitted `strain_inc_pl_vol` entirely (both from `shared` and from the interpolation in the original code, though it's used in the current code). This was likely the bug that caused the pragma to be commented out.

## Risks / Trade-offs

- **[Risk] Hidden reason for commenting out the pragma** → Mitigated by: (a) the loop is provably race-free (each particle writes to distinct `k` index), (b) all CI tests will validate correctness, (c) the pragma variable lists had genuine bugs (`dE_pl_vol` missing) which is the most likely reason it was disabled.
- **[Risk] `exit(0)` removal changes error behaviour** → Mitigated by: the condition should never trigger in a correct simulation. If it does, `LOG_ERR` will report it and the simulation will continue (possibly with wrong results), which is better than `exit()` in a parallel region causing memory corruption.
- **[Risk] No performance gain if post_solve is dominated by other functions** → Mitigated by: `AccumulatedStrainII` iterates over all ~12M particles with bilinear interpolation — it's the heaviest function in the post_solve phase by operation count.
