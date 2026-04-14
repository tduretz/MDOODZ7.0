## Context

The thermal solver calls `EnergyDirectSolve()` once per timestep. Inside that function, a local `cholmod_common c` is initialized via `cholmod_start()`, the thermal matrix is assembled, `cholmod_analyze()` computes the symbolic factorization (fill-reducing ordering + symbolic Cholesky tree), `cholmod_factorize()` performs numeric factorization, `cholmod_solve()` does back-substitution, and then everything is freed via `cholmod_free_factor()` + `cholmod_finish()`.

The sparsity pattern of the thermal matrix is determined by the 5-point stencil on cell centres and is **constant** across timesteps (grid topology never changes during a run). Only the numeric values change (due to updated conductivity, density×Cp, etc.). This means `cholmod_analyze()` — which computes AMD ordering, elimination tree, and column counts — produces identical output every timestep and is pure redundant work.

The Stokes solver already follows a persistent pattern: `DirectSolver` struct holds `cholmod_common c`, `cholmod_factor *Lfact`, and an `Analyze` flag. The thermal solver should adopt the same approach.

## Goals / Non-Goals

**Goals:**
- Eliminate redundant `cholmod_analyze()` calls by caching the symbolic factorization across timesteps
- Eliminate redundant `cholmod_start()`/`cholmod_finish()` calls by persisting the `cholmod_common` context
- Maintain identical numerical results (bit-for-bit)
- Preserve existing phase-timing instrumentation accuracy

**Non-Goals:**
- Parallelizing the thermal CHOLMOD factorization (prior experiment showed CHOLMOD threading has zero effect)
- Caching the numeric factorization (`cholmod_factorize` must run every timestep since coefficients change)
- Changing the thermal matrix assembly or stencil
- Supporting dynamic grid refinement / remeshing (not implemented in MDOODZ)

## Decisions

### Decision 1: Reuse the existing `DirectSolver` struct

**Choice:** Add a second `DirectSolver` instance (`ThermalSolver`) in `Main_DOODZ.c` alongside the existing `CholmodSolver`, rather than introducing a new struct type.

**Rationale:** The `DirectSolver` struct already contains exactly the fields needed: `cholmod_common c`, `cholmod_factor *Lfact`, and `int Analyze`. Using the same struct avoids code duplication and is consistent with established patterns. No new types in the header.

**Alternative considered:** A lightweight `ThermalCholmodCache` struct with only the three needed fields. Rejected because `DirectSolver` already exists and adding another struct would be gratuitous.

### Decision 2: Pass `DirectSolver*` into `EnergyDirectSolve`

**Choice:** Change the signature of `EnergyDirectSolve()` to accept a `DirectSolver*` parameter instead of creating a local `cholmod_common`.

**Rationale:** This is the minimal change that threads the persistent state through. The caller (`Main_DOODZ.c`) manages the lifecycle (init before timestep loop, finish after), matching the Stokes solver pattern exactly.

**Alternative considered:** Making the `DirectSolver` a global variable. Rejected — the codebase avoids globals for solver state; the Stokes solver passes it explicitly.

### Decision 3: Split `FactorEnergyCHOLMOD` into analyze + factorize paths

**Choice:** Add a flag check inside `FactorEnergyCHOLMOD`: when `Analyze == 1`, run `cholmod_analyze()` then `cholmod_factorize()`; when `Analyze == 0`, call only `cholmod_factorize()` on the existing cached factor.

**Rationale:** `cholmod_factorize()` accepts an existing `cholmod_factor*` and re-uses its symbolic structure. This is exactly how the SuiteSparse API is designed to be used — analyze once, factorize many times. The Stokes solver already uses this pattern.

**Alternative considered:** Calling `cholmod_analyze()` every time but caching a separate symbolic factor and copying it before numeric factorization. Rejected — unnecessary complexity; CHOLMOD's `cholmod_factorize()` already handles refactorization in-place.

### Decision 4: Lifecycle matches Stokes solver pattern

**Choice:** `cholmod_start()` before the timestep loop, `Analyze = 1` on first thermal solve, `Analyze = 0` thereafter, `cholmod_free_factor()` + `cholmod_finish()` after the loop.

**Rationale:** Proven pattern already in the codebase. The `Analyze` flag is set to 1 only once (before first thermal step). Since the grid never changes mid-run, no re-analysis is needed.

### Decision 5: No invalidation mechanism needed

**Choice:** Do not implement a cache invalidation / re-analyze trigger.

**Rationale:** MDOODZ does not support grid topology changes during a run. The grid is allocated once at startup. Adding invalidation logic would be dead code. If remeshing is ever added in the future, it can be added then (set `Analyze = 1`).

## Risks / Trade-offs

**[Risk] Memory footprint increase** → The cached `cholmod_factor` persists for the entire simulation instead of being freed each timestep. Mitigation: The factor for a thermal 5-point stencil is small relative to the Stokes system (which is already cached). At 1001×801 the thermal system is ~800k DOFs with ~4M non-zeros; the factor memory is a few hundred MB at most, already within the working set.

**[Risk] Stale factor if code is ever changed to modify grid** → Mitigation: Non-goal for now. A future change adding remeshing would need to set `ThermalSolver.Analyze = 1` after any grid change. This is documented as a known constraint.

**[Risk] `ThermalSteps()` also calls `EnergyDirectSolve`** → The function `ThermalSteps()` in ThermalRoutines.c calls `EnergyDirectSolve()` in a sub-loop for thermal-only stepping. It will also need to be updated to accept and pass through the `DirectSolver*`. Mitigation: It is called from `Main_DOODZ.c` in the same scope where `ThermalSolver` lives, so threading it through is straightforward.
