## ADDED Requirements

### Requirement: Persistent CHOLMOD context for thermal solver
The system SHALL initialize a single `cholmod_common` context for the thermal solver before the timestep loop and finalize it after the loop exits, rather than creating and destroying it on every call to `EnergyDirectSolve`.

#### Scenario: Thermal CHOLMOD context initialized once
- **WHEN** the simulation starts the timestep loop
- **THEN** `cholmod_start()` SHALL be called exactly once for the thermal solver context

#### Scenario: Thermal CHOLMOD context finalized once
- **WHEN** the simulation exits the timestep loop
- **THEN** `cholmod_finish()` SHALL be called exactly once for the thermal solver context

### Requirement: Cached symbolic factorization across timesteps
The system SHALL perform `cholmod_analyze()` only on the first timestep that invokes the thermal solver. On subsequent timesteps, the system SHALL reuse the cached `cholmod_factor` and call only `cholmod_factorize()` (numeric) and `cholmod_solve()` (back-substitution).

#### Scenario: First thermal solve performs full analysis
- **WHEN** `EnergyDirectSolve` is called for the first time in a simulation
- **THEN** `cholmod_analyze()` SHALL be called to compute the symbolic factorization
- **AND** the resulting `cholmod_factor` SHALL be stored for reuse

#### Scenario: Subsequent thermal solves skip analysis
- **WHEN** `EnergyDirectSolve` is called on the second or later timestep
- **THEN** `cholmod_analyze()` SHALL NOT be called
- **AND** `cholmod_factorize()` SHALL be called using the previously cached `cholmod_factor`

### Requirement: Bit-for-bit identical results
The system SHALL produce identical numerical output (temperature field, energy balance) whether or not the factorization cache is active, since `cholmod_analyze` with the same sparsity pattern produces the same ordering every time.

#### Scenario: Cached solve matches uncached solve
- **WHEN** a simulation runs with cached thermal factorization
- **THEN** the temperature field at each timestep SHALL be bit-for-bit identical to a run without caching (same ordering, same numeric factorization, same solve)

### Requirement: DirectSolver struct used for thermal state
The system SHALL use the existing `DirectSolver` struct (containing `cholmod_common c`, `cholmod_factor *Lfact`, `int Analyze`) to hold the persistent thermal solver state, matching the Stokes solver pattern.

#### Scenario: Thermal solver uses DirectSolver struct
- **WHEN** the thermal solver state is initialized
- **THEN** a `DirectSolver` instance SHALL be used with `Analyze` set to 1

#### Scenario: Analyze flag cleared after first solve
- **WHEN** the first thermal solve completes its `cholmod_analyze` call
- **THEN** the `Analyze` flag SHALL be set to 0

### Requirement: EnergyDirectSolve accepts external solver state
The function `EnergyDirectSolve` SHALL accept a `DirectSolver*` parameter for the CHOLMOD context and cached factor, instead of creating a local `cholmod_common`.

#### Scenario: EnergyDirectSolve signature includes DirectSolver
- **WHEN** `EnergyDirectSolve` is called
- **THEN** it SHALL receive a `DirectSolver*` parameter that provides the `cholmod_common` context and cached `cholmod_factor`

### Requirement: ThermalSteps passes through solver state
The function `ThermalSteps` SHALL accept and forward a `DirectSolver*` to its internal `EnergyDirectSolve` calls.

#### Scenario: ThermalSteps propagates DirectSolver
- **WHEN** `ThermalSteps` calls `EnergyDirectSolve`
- **THEN** it SHALL pass the same `DirectSolver*` received from its caller
