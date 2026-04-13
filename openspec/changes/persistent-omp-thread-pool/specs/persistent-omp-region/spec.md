## ADDED Requirements

### Requirement: Persistent OpenMP region parameter
The system SHALL support an integer parameter `omp_schedule` in the `.txt` configuration file, read via `ReadInt2(fin, "omp_schedule", 0)`. When set to `0` (default), the system SHALL use the legacy per-function fork-join parallelization. When set to `1`, the system SHALL use persistent `#pragma omp parallel` regions around compute-intensive blocks in the main timestep loop.

#### Scenario: Parameter absent from .txt file
- **WHEN** `omp_schedule` is not specified in the `.txt` configuration file
- **THEN** the system SHALL default to `0` (legacy fork-join behavior) and produce identical results to the current codebase

#### Scenario: Parameter set to 0
- **WHEN** `omp_schedule = 0` is specified in the `.txt` configuration file
- **THEN** the system SHALL use per-function `#pragma omp parallel for` regions (legacy behavior) with no change to runtime behavior or results

#### Scenario: Parameter set to 1
- **WHEN** `omp_schedule = 1` is specified in the `.txt` configuration file
- **THEN** the system SHALL wrap three compute blocks in persistent `#pragma omp parallel` regions: (1) the interpolation block (P2Mastah calls), (2) the nonlinear iteration body (rheology + assembly), and (3) the post-solve particle update sequence

### Requirement: Function-level parallel region detection
All MDLIB functions that contain `#pragma omp parallel for` directives and are called from within a persistent parallel region SHALL detect the active parallel context using `omp_in_parallel()` and use `#pragma omp for` (worksharing) instead of `#pragma omp parallel for` (fork-join) when already inside a parallel region.

#### Scenario: Function called from persistent region
- **WHEN** a guarded function is called while `omp_in_parallel()` returns true
- **THEN** the function SHALL use `#pragma omp for` worksharing directives instead of creating a new parallel region

#### Scenario: Function called standalone (legacy mode)
- **WHEN** a guarded function is called while `omp_in_parallel()` returns false
- **THEN** the function SHALL use `#pragma omp parallel for` to create its own parallel region, preserving legacy behavior

### Requirement: Manual thread decomposition compatibility
Functions that use manual work distribution via `omp_get_thread_num()` (e.g., StokesAssemblyDecoupled) SHALL detect an active parallel context using `omp_in_parallel()` and skip creating a new `#pragma omp parallel` region when already inside one, while preserving the manual `omp_get_thread_num()` based work distribution.

#### Scenario: Manual decomposition inside persistent region
- **WHEN** a function with manual `omp_get_thread_num()` work distribution is called while `omp_in_parallel()` returns true
- **THEN** the function SHALL skip the `#pragma omp parallel` directive and execute the manual work distribution loop directly using the thread ID from the enclosing persistent region

#### Scenario: Manual decomposition standalone
- **WHEN** a function with manual `omp_get_thread_num()` work distribution is called while `omp_in_parallel()` returns false
- **THEN** the function SHALL create its own `#pragma omp parallel` region as before

### Requirement: Numerical equivalence between modes
The system SHALL produce bit-identical simulation results regardless of whether `omp_schedule` is set to `0` or `1`. The persistent region optimization SHALL NOT alter the order of floating-point operations, reduction accumulation, or memory access patterns relative to legacy mode.

#### Scenario: Output comparison between modes
- **WHEN** the same simulation is run with `omp_schedule = 0` and `omp_schedule = 1` using the same thread count and input parameters
- **THEN** the HDF5 output fields SHALL be identical to machine precision (relative difference < 1e-10)

### Requirement: Persistent region scoping
Each persistent `#pragma omp parallel` region SHALL be tightly scoped around a compute-intensive block. Serial orchestration code (timer reads, LOG_TIME calls, conditionals, function dispatch) SHALL execute outside the persistent region or within `#pragma omp single` / `#pragma omp master` blocks.

#### Scenario: Timer accuracy preserved
- **WHEN** `omp_schedule = 1` is active
- **THEN** all `omp_get_wtime()` timer measurements SHALL be taken on a single thread (outside the persistent region or inside `#pragma omp master`) and the perf.csv values SHALL be consistent with legacy mode (within 10% for equivalent workloads)
