## ADDED Requirements

### Requirement: CHOLMOD thread count parameter
The system SHALL support an integer parameter `cholmod_threads` in the `.txt` configuration file, read via `ReadInt2(fin, "cholmod_threads", 1)`. The value `1` (default) SHALL configure CHOLMOD for single-threaded sparse BLAS. The value `-1` SHALL configure CHOLMOD to use all available OpenMP threads. Any positive integer N SHALL configure CHOLMOD to use exactly N threads.

#### Scenario: Parameter absent from .txt file
- **WHEN** `cholmod_threads` is not specified in the `.txt` configuration file
- **THEN** the system SHALL default to `1` and CHOLMOD SHALL operate in single-threaded mode, producing identical behavior to the current codebase

#### Scenario: Parameter set to -1
- **WHEN** `cholmod_threads = -1` is specified in the `.txt` configuration file
- **THEN** the system SHALL set `cholmod_common.nthreads_max` to the value of `omp_get_num_threads()` for all CHOLMOD instances (Stokes, thermal, chemical solvers)

#### Scenario: Parameter set to explicit count
- **WHEN** `cholmod_threads = N` (N > 0) is specified in the `.txt` configuration file
- **THEN** the system SHALL set `cholmod_common.nthreads_max = N` for all CHOLMOD instances

### Requirement: CHOLMOD threading applied to all solver instances
The system SHALL set `c.nthreads_max` after every `cholmod_start()` call in the codebase: the Stokes solver (Main_DOODZ.c), the thermal solver (ThermalRoutines.c), and the chemical solver (ChemicalRoutines.c).

#### Scenario: Stokes solver CHOLMOD threading
- **WHEN** the Stokes solver initializes CHOLMOD via `cholmod_start(&CholmodSolver.c)`
- **THEN** `CholmodSolver.c.nthreads_max` SHALL be set to the resolved thread count before any `cholmod_analyze` or `cholmod_factorize` calls

#### Scenario: Thermal solver CHOLMOD threading
- **WHEN** the thermal solver initializes CHOLMOD via `cholmod_start(&c)`
- **THEN** `c.nthreads_max` SHALL be set to the resolved thread count before any `cholmod_analyze` or `cholmod_factorize` calls

#### Scenario: Chemical solver CHOLMOD threading
- **WHEN** the chemical solver initializes CHOLMOD via `cholmod_start(&c)`
- **THEN** `c.nthreads_max` SHALL be set to the resolved thread count before any `cholmod_analyze` or `cholmod_factorize` calls

### Requirement: CHOLMOD thread count logging
The system SHALL log the configured CHOLMOD thread count at startup via `LOG_INFO` to aid diagnostics and benchmark analysis.

#### Scenario: Thread count logged
- **WHEN** CHOLMOD threading is configured
- **THEN** a `LOG_INFO` message SHALL be emitted with the format "CHOLMOD threads: <N>" showing the resolved thread count
