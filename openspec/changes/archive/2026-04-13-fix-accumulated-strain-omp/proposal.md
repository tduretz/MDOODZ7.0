## Why

`AccumulatedStrainII` in `RheologyParticles.c` has a commented-out `#pragma omp parallel for` on its main particle loop (line 240), making it single-threaded. This function runs inside the `post_solve` timer, which at 3.21s (35% of wall time at 16 threads on 1001×801) is now the dominant bottleneck after the PCG thermal solver and persistent interpolation buffer optimizations. The particle loop is embarrassingly parallel — each particle reads from shared grid arrays and writes only to its own `particles->strain[k]` — so re-enabling OpenMP should be straightforward and yield a significant `post_solve` reduction.

## What Changes

- Uncomment the `#pragma omp parallel for` on the particle loop in `AccumulatedStrainII` (line ~240 of `RheologyParticles.c`)
- Add `dE_pl_vol` to the `private` clause (it was missing from the original commented pragma)
- Verify thread-safety: all writes are to `particles->strain*[k]` (distinct per-particle, no races)
- Remove the `exit(0)` call inside the negative-strain-increment check (incompatible with parallel execution — replace with a thread-safe logging approach)

## Capabilities

### New Capabilities

_(none)_

### Modified Capabilities

_(none — this is a pure implementation fix with no spec-level behavior change; the function already exists and its contract is unchanged)_

## Impact

- **Code**: `MDLIB/RheologyParticles.c` — `AccumulatedStrainII` function only
- **Performance**: Expected ~50% reduction in `post_solve` time at ≥8 threads (from ~3.2s to ~1.6s), translating to ~15-17% overall wall time improvement
- **Correctness risk**: Low — the loop body only reads shared arrays and writes to per-particle arrays. The `exit(0)` on negative strain needs replacement since it's unsafe in a parallel region.
- **CI**: Existing tests cover this function implicitly through all scenarios that accumulate strain. No new test fixtures needed.
