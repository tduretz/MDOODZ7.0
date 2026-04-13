## Why

The thermal solver (`EnergyDirectSolve`) recreates the entire CHOLMOD context (`cholmod_start`/`cholmod_finish`) and re-runs symbolic analysis (`cholmod_analyze`) every timestep, even though the sparsity pattern is constant (5-point stencil on a fixed grid). Benchmarks show `thermal_s` ≈ 7 s/step at 1001×801 (38% of wall time), entirely serial. Caching the symbolic factorization and CHOLMOD context across timesteps eliminates redundant work and may yield significant speedup.

## What Changes

- **Persist `cholmod_common` for the thermal solver** across timesteps instead of calling `cholmod_start`/`cholmod_finish` inside every `EnergyDirectSolve` call. Follow the existing Stokes solver pattern (`CholmodSolver` struct with persistent `.c` field).
- **Cache `cholmod_factor *Lfact` symbolic analysis** returned by `cholmod_analyze`. Reuse it across timesteps so only `cholmod_factorize` (numeric) and `cholmod_solve` (back-substitution) run each step.
- **Add an invalidation flag** so the cached factor is discarded and re-analyzed if the grid changes (e.g., remeshing — currently not used but future-proof).
- **Lifecycle management**: Initialize the thermal CHOLMOD context once before the timestep loop and finalize it after the loop exits, mirroring how the Stokes solver manages `CholmodSolver`.

## Capabilities

### New Capabilities
- `thermal-factorization-cache`: Caching of CHOLMOD symbolic analysis and persistent context for the thermal direct solver, with invalidation on grid topology change.

### Modified Capabilities
- `phase-timing`: The thermal solver timing instrumentation must continue to work correctly with the restructured CHOLMOD lifecycle (context no longer created/destroyed per timestep).

## Impact

- **MDLIB/ThermalRoutines.c**: `EnergyDirectSolve` refactored — remove per-call `cholmod_start`/`cholmod_finish`, accept persistent context + cached factor, skip `cholmod_analyze` when factor is valid.
- **MDLIB/ThermalSolver.c**: `FactorEnergyCHOLMOD` split into symbolic-only and numeric-only paths, or parameterized to skip analysis when a cached factor is provided.
- **MDLIB/Main_DOODZ.c**: Thermal CHOLMOD context and factor initialized before timestep loop, finalized after. Similar pattern to existing `CholmodSolver` struct.
- **MDLIB/include/mdoodz.h**: New struct or fields to hold persistent thermal CHOLMOD state.
