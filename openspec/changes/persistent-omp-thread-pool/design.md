## Context

MDOODZ uses OpenMP across 76 pragma directives in 6 core files. The current pattern is per-function fork-join: each function (P2Mastah, Interp_Grid2P, NonNewtonianViscosityGrid, BuildStokesOperator, etc.) creates its own parallel region. Benchmarks on c5ad.4xlarge (1000×800) show scaling regression past 6-8 threads — wall time increases from 17.0s at 6t to 20.5s at 16t, with `interp_s` (P2Mastah) responsible for most of the degradation.

The 76 pragmas break down as:
- **49 simple `parallel for`** loops (RheologyParticles: 33, ThermalSolver: 6, ParticleRoutines: 10) — trivially convertible
- **11 `parallel` regions** that create new thread teams — 8 in StokesAssemblyDecoupled (manual work distribution via `omp_get_thread_num()`), 2 in ParticleRoutines (P2Mastah), 1 in Main_DOODZ.c
- **2 `barrier`** directives in StokesAssemblyDecoupled

No `omp_in_parallel()` usage exists in the current codebase. Thread count is set entirely via `OMP_NUM_THREADS` environment variable.

## Goals / Non-Goals

**Goals:**
- Eliminate fork-join overhead by using `omp_in_parallel()` checks — functions detect whether they are already inside an active parallel region and skip their own `#pragma omp parallel` if so
- New `.txt` parameter `omp_schedule` (default 0 = legacy, 1 = persistent region) controls behavior
- All functions remain callable in both persistent and legacy contexts
- Benchmark on EC2 (1000×800, 1-16 threads) to quantify improvement
- Enable HDF5 writer in benchmark to measure `output_s`
- No changes to simulation results — performance-only optimization

**Non-Goals:**
- Restructuring StokesAssemblyDecoupled's manual thread decomposition (too risky, separate change)
- Fusing P2Mastah's 3 internal parallel regions into 1 (optimization #3, separate change)
- Parallelizing the P2Mastah reduction loop (optimization #4, separate change)
- Particle sorting for cache locality (optimization #5, separate change)
- Replacing CHOLMOD with a parallel solver (optimization #2, separate change)
- Converting every function's internal parallel loops — only the top-level entry points need `omp_in_parallel()` guards

## Decisions

### D1: `omp_in_parallel()` guard approach vs. full refactor

**Decision:** Use `omp_in_parallel()` guards at function entry, not a full refactor of all 76 pragmas.

**Rationale:** A full refactor (converting all `#pragma omp parallel for` to `#pragma omp for`) would require changes to every function signature and would break the "callable standalone" contract. Instead, each function checks `omp_in_parallel()` at the top:
- If already in a parallel region → use `#pragma omp for` (worksharing within existing team)
- If not → use `#pragma omp parallel for` (create team as before)

This is safe, backward compatible, and incrementally testable. A function can be converted one at a time.

**Alternatives considered:**
- *Full refactor to worksharing only:* Would require all callers to provide a parallel context. Breaks modularity, not backward compatible.
- *Nested parallelism (`omp_set_nested(1)`):* Creates threads-of-threads, massive overhead, counterproductive.

### D2: Persistent region placement — around the timestep body in Main_DOODZ.c

**Decision:** Place a single `#pragma omp parallel` around the timestep body (from after `t_omp_step` to before perf.csv output), gated by `omp_schedule == 1`.

**Structure:**
```c
if (input.model.omp_schedule == 1) {
    #pragma omp parallel
    {
        // All timestep work happens here
        // Functions detect omp_in_parallel() and skip their own #pragma omp parallel
        #pragma omp single
        {
            // Serial orchestration: function calls, timers, logging
            // Each function internally uses #pragma omp for when it detects active region
        }
    }
} else {
    // Legacy code path — unchanged
}
```

Wait — this is wrong. `#pragma omp single` would serialize everything inside it. The correct pattern is:

**Revised structure:**
```c
if (input.model.omp_schedule == 1) {
    #pragma omp parallel
    {
        // Functions called from here will detect omp_in_parallel() == true
        // They will use #pragma omp for instead of #pragma omp parallel for
        // Serial sections (timers, logging, orchestration) use #pragma omp single/master
    }
} else {
    // Legacy: unchanged
}
```

**Problem:** The timestep body has complex serial logic between parallel loops (conditionals, timer bookkeeping, function dispatch). We cannot simply wrap everything in `#pragma omp parallel` because the serial orchestration needs `#pragma omp single` blocks, which adds complexity and risk.

**Revised Decision:** Instead of wrapping the entire timestep, make each function self-adaptive using `omp_in_parallel()`. The persistent region wraps only the heavy compute blocks where fork-join overhead matters most:

1. **Interpolation block** (P2Mastah calls) — wrap the full interpolation sequence
2. **Nonlinear iteration loop body** — wrap rheology + assembly per iteration
3. **Post-solve particle updates** — wrap the update sequence

This avoids restructuring the serial control flow while eliminating fork-join overhead in the 3 hottest blocks.

**Alternatives considered:**
- *Wrap entire timestep body:* Too much serial logic to convert to `#pragma omp single`. High risk of implicit barriers and deadlocks.
- *Single persistent region for the whole run:* The `for` timestep loop has complex adaptive logic (variable step count, I/O decisions) that cannot easily be parallelized.

### D3: Which functions get `omp_in_parallel()` guards

**Decision:** Add guards to functions in the three hot blocks, prioritized by benchmark impact:

**Priority 1 — interpolation (interp_s, 6.5s at 16t):**
- `P2Mastah` (ParticleRoutines.c) — 3 parallel regions → guarded
- `Interp_Grid2P` variants (ParticleRoutines.c) — simple `parallel for` → guarded

**Priority 2 — nonlinear iteration (rheology + assembly):**
- `NonNewtonianViscosityGrid` (RheologyParticles.c) — 33 `parallel for` → guarded
- `BuildStokesOperator` (StokesAssemblyDecoupled.c) — 8 manual parallel regions → guarded (complex; these use `omp_get_thread_num()` internally which still works inside a persistent region)

**Priority 3 — post-solve updates:**
- `UpdateNonLinearity`, `UpdateParticleEnergy`, `UpdateDensity`, etc. — simple `parallel for` → guarded

**Guard pattern:**
```c
void SomeFunction(...) {
    int in_parallel = omp_in_parallel();
    // ...
    if (in_parallel) {
        #pragma omp for schedule(static) private(k)
        for (k = 0; k < N; k++) { ... }
    } else {
        #pragma omp parallel for schedule(static) private(k)
        for (k = 0; k < N; k++) { ... }
    }
}
```

For P2Mastah specifically, the thread-local buffer allocation (`Wm[nthreads]`) is already done outside the parallel loops using `omp_get_num_threads()`. Inside a persistent region, `omp_get_num_threads()` returns the team size correctly, so no change needed to the allocation logic — only the pragma directives around the loops.

### D4: StokesAssemblyDecoupled — manual decomposition preserved

**Decision:** Keep StokesAssemblyDecoupled's manual `#pragma omp parallel` blocks with `omp_get_thread_num()` as-is. Add `omp_in_parallel()` guard: if already in a parallel region, skip the outer `#pragma omp parallel` and just execute the manual work distribution loop directly (thread IDs are still valid inside a persistent region).

```c
// Current:
#pragma omp parallel shared(estart, eend, ...) private(ith)
{
    ith = omp_get_thread_num();
    for (c = estart[ith]; c < eend[ith]+1; c++) { ... }
}

// With guard:
if (!omp_in_parallel()) {
    #pragma omp parallel shared(estart, eend, ...) private(ith)
    {
        ith = omp_get_thread_num();
        for (c = estart[ith]; c < eend[ith]+1; c++) { ... }
    }
} else {
    int ith = omp_get_thread_num();
    for (c = estart[ith]; c < eend[ith]+1; c++) { ... }
}
```

This works because inside the persistent region, each thread already has a valid `omp_get_thread_num()` and the `estart/eend` arrays are pre-computed for the correct number of threads.

**Risk:** The `estart/eend` decomposition must be computed from the same thread count as the persistent pool. This is safe because both use `omp_get_num_threads()` which is consistent within a team.

### D5: Parameter and control flow

**Decision:** Add `omp_schedule` as an `int` field in the `model` struct, read via `ReadInt2(fin, "omp_schedule", 0)` in InputOutput.c. Default = 0 (legacy). Value 1 enables persistent regions.

The persistent regions are placed in Main_DOODZ.c using a wrapper pattern:
```c
if (input.model.omp_schedule == 1) {
    #pragma omp parallel
    {
        // Interpolation block: multiple P2Mastah calls
        // Each P2Mastah detects omp_in_parallel() and uses #pragma omp for
    }
} else {
    // Same calls, but each function uses its own #pragma omp parallel for
}
```

Three such wrappers are added: around interpolation, around the NL iteration body, and around post-solve updates. The serial control flow (conditionals, timer reads, logging) stays outside the wrappers.

### D6: Benchmark configuration — enable HDF5 output

**Decision:** Modify `benchmark-ec2-setup.sh` to enable `writer = 1` and `writer_step = 1` in the highres `.txt` file before running. This populates the `output_s` column in perf.csv. The lowres/medres runs keep `writer = 0` for faster iteration.

## Risks / Trade-offs

**[Implicit barriers from `#pragma omp for`]** → Every `#pragma omp for` has an implicit barrier at the end. If two consecutive `for` loops operate on dependent data, this is correct behavior. If they don't, the barrier wastes time. → **Mitigation:** Use `#pragma omp for nowait` where data dependencies allow it. Start with barriers (safe default), profile later.

**[Thread count mismatch in P2Mastah buffers]** → P2Mastah allocates thread-local arrays based on `omp_get_num_threads()`. If the persistent region uses a different thread count than expected, buffer indexing would be wrong. → **Mitigation:** Inside a persistent region, `omp_get_num_threads()` returns the team size, which is the same value used for `omp_get_thread_num()`. No mismatch possible within a single team.

**[StokesAssemblyDecoupled `estart/eend` consistency]** → The element decomposition arrays are computed using a thread count query. If called from a persistent region, the thread count is the pool size, so decomposition is consistent. → **Mitigation:** Verify in testing that `estart/eend` allocation matches the persistent region's thread count.

**[Code duplication from if/else guard blocks]** → Every guarded function has an if/else with near-identical loop bodies (one `omp for`, one `omp parallel for`). → **Mitigation:** Use a macro `MDOODZ_OMP_FOR(clauses)` that expands to the correct pragma based on `omp_in_parallel()`. Or accept the duplication for clarity.

**[Correctness risk — silent numerical differences]** → If barrier placement differs between legacy and persistent mode, iteration order or accumulation order could change, producing bit-level numerical differences. → **Mitigation:** Run both modes on the same problem and compare all output fields to machine precision. No algorithmic changes means results should match exactly.

**[Serial sections inside persistent region]** → Timer reads (`omp_get_wtime()`), `LOG_TIME`, and control flow between parallel blocks must execute on a single thread. → **Mitigation:** Use `#pragma omp single` or `#pragma omp master` for serial sections within the persistent region. Keep persistent regions tightly scoped around compute blocks only.

## Open Questions

1. **Macro vs. explicit if/else for guard pattern?** A macro like `MDOODZ_PARALLEL_FOR(schedule, private_vars)` could reduce duplication, but macros hiding pragmas can be hard to debug. Start with explicit if/else, consider macro after validation.

2. **`nowait` opportunities?** Adjacent `#pragma omp for` loops on independent arrays could use `nowait` to avoid barrier overhead. Needs case-by-case analysis of data dependencies.

3. **Should the persistent region also wrap the thermal solver?** The thermal solver (6.5s, constant) is dominated by CHOLMOD (serial). Wrapping it in the persistent region would only help the 6 trivial `parallel for` loops in ThermalSolver.c (~0.01s total). Low priority.
