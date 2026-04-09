## ADDED Requirements

### Requirement: Skill documents CMake build configuration
The skill SHALL explain the CMake build system: minimum version 3.16, C11 standard, and the available CMake options:
- `-DJULIA=ON`: Build with Julia integration support
- `-DOMP=ON`: Enable OpenMP parallelisation
- `-DOPT=ON`: Enable optimised build (-O3)
- `-DTEST=ON`: Build test suite
- Default mode: debug (-g -Wall -O0 -fno-inline)

#### Scenario: User wants to build with optimisations
- **WHEN** the user asks how to build a fast version
- **THEN** the skill SHALL recommend `make build-dev OMP=ON OPT=ON SET=<Name>` or the equivalent CMake invocation with `-DOMP=ON -DOPT=ON`

### Requirement: Skill documents make targets
The skill SHALL list all make targets available in the top-level `makefile`:
- `make build SET=<Name>`: Build with default settings for a specific scenario
- `make build SET=<Name> TXT=<Custom.txt>`: Build with a custom parameter file
- `make build-dev OMP=ON OPT=ON SET=<Name>`: Development build with OpenMP and optimisations
- `make run SET=<Name>`: Execute a compiled simulation
- `make clean`: Remove build artifacts

#### Scenario: User wants to run a simulation end-to-end
- **WHEN** the user asks how to run a simulation
- **THEN** the skill SHALL provide the full workflow: `make build SET=RiftingChenin` → `make run SET=RiftingChenin`, explaining that the SET variable selects the `.c` file from `SETS/` and uses the matching `.txt` file by default

### Requirement: Skill documents dependency installation
The skill SHALL explain all required dependencies and how to install them:
- **SuiteSparse**: sparse linear algebra (CHOLMOD, UMFPACK, AMD, COLAMD, CCOLAMD, CAMD, SuiteSparse_config)
- **HDF5**: output file format
- **BLAS/LAPACK**: basic linear algebra
- **zlib**: compression for HDF5
- **pcg-c-basic**: random number generation (automatically fetched via CMake FetchContent)
- Optional: **OpenMP** runtime for parallelisation

#### Scenario: User encounters a missing dependency error
- **WHEN** the user gets a CMake or link error about SuiteSparse, HDF5, or BLAS
- **THEN** the skill SHALL provide platform-specific installation instructions (e.g., `brew install suite-sparse hdf5` on macOS, `apt install libsuitesparse-dev libhdf5-dev` on Ubuntu/Debian)

### Requirement: Skill documents env.cmake configuration
The skill SHALL explain the `env.cmake` file (created from `env.cmake.example`) for specifying custom paths to dependencies when they are not in standard system locations. The skill SHALL note that this file is gitignored and user-specific.

#### Scenario: User has dependencies in non-standard paths
- **WHEN** the user installs SuiteSparse or HDF5 in a custom location
- **THEN** the skill SHALL explain how to copy `env.cmake.example` to `env.cmake` and set `CMAKE_PREFIX_PATH` or individual `*_DIR` variables

### Requirement: Skill documents OpenMP parallelisation
The skill SHALL explain that MDOODZ supports shared-memory parallelisation via OpenMP: enabled at compile time with `-DOMP=ON`, controlled at runtime via `OMP_NUM_THREADS` environment variable. Key parallelised loops are in particle operations, rheology evaluation, and matrix assembly.

#### Scenario: User wants to use multiple cores
- **WHEN** the user asks about parallel execution
- **THEN** the skill SHALL explain: build with `OMP=ON`, set `export OMP_NUM_THREADS=<N>` before running, and note that the direct solver (CHOLMOD) also benefits from multiple threads

### Requirement: Skill documents common troubleshooting
The skill SHALL provide solutions for common build and runtime issues:
- SuiteSparse not found: check `FindSuiteSparse.cmake` paths
- HDF5 version mismatch: ensure consistent HDF5 version across library and headers
- Segfault on startup: check that `.txt` file exists and `Nb_phases` matches defined phases
- NaN in output: check time step size, viscosity bounds, and material parameters

#### Scenario: User gets a segfault when running
- **WHEN** the user reports a crash during simulation startup
- **THEN** the skill SHALL suggest: verify the `.txt` file is readable, check `Nb_phases` matches the number of `ID` blocks, ensure `Nx`/`Nz` are reasonable values, and try a debug build (`make build SET=<Name>`) for more informative error messages
