---
name: skill-build-and-run
description: Build, compile, and run MDOODZ simulations — CMake configuration, make targets (build, build-dev, run, clean), dependency installation (SuiteSparse, HDF5, BLAS/LAPACK), OpenMP parallelisation, env.cmake setup, and common troubleshooting.
---

# Building and Running MDOODZ7.0

## Build System Overview

MDOODZ uses **CMake** (minimum 3.16) with a convenience **makefile** wrapper. The C standard is C11.

## Make Targets

| Target | Command | Description |
|--------|---------|-------------|
| Build (optimised) | `make build SET=<Name>` | Builds with `-O3` and OpenMP enabled |
| Build (custom txt) | `make build SET=<Name> TXT=Custom.txt` | Uses a non-default parameter file |
| Build (dev) | `make build-dev OMP=ON OPT=ON SET=<Name>` | Explicit control over flags |
| Run | `make run SET=<Name>` | Execute compiled simulation |
| Run (copy txt) | `make run-copy SET=<Name>` | Copy `.txt` from SETS/ then run |
| Clean (one set) | `make clean SET=<Name>` | Remove build artifacts for one scenario |
| Clean (all) | `make clean-all` | Remove all build and exec directories |
| Install deps | `make deps` | Clone and build HDF5 + SuiteSparse from helper repo |
| Run tests | `make run-tests` | Run CTest suite |

### Typical Workflow

```bash
# 1. Build a scenario (optimised + OpenMP)
make build SET=RiftingChenin

# 2. Run it
make run SET=RiftingChenin

# Output files appear in cmake-exec/RiftingChenin/
```

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `-DOPT=ON` | OFF | Enable optimisation (`-O3 -ftree-vectorize -funroll-loops`) |
| `-DOMP=ON` | OFF | Enable OpenMP parallelisation |
| `-DJULIA=ON` | OFF | BinaryBuilder mode for Julia integration |
| `-DTEST=ON` | OFF | Build test suite (enables CTest) |
| `-DVIS=ON` | OFF | Build visual tests |

**Debug mode** (default when OPT=OFF): `-g -Wall -O0 -fno-inline -fno-omit-frame-pointer`

## Dependencies

| Library | Purpose | Required |
|---------|---------|----------|
| **SuiteSparse** | Sparse direct solver (CHOLMOD, UMFPACK, AMD, COLAMD, CCOLAMD, CAMD) | Yes |
| **HDF5** | Output file I/O | Yes |
| **BLAS/LAPACK** | Basic linear algebra (used by SuiteSparse) | Yes |
| **zlib** | Compression for HDF5 output | Yes |
| **pcg-c-basic** | Random number generation | Auto-fetched via CMake FetchContent |
| **OpenMP** | Shared-memory parallelisation | Optional (recommended) |

### Installing Dependencies

**Using the built-in helper:**
```bash
make deps  # Clones and builds HDF5 + SuiteSparse
```

**macOS (Homebrew):**
```bash
brew install suite-sparse hdf5 open-blas lapack
```

**Ubuntu/Debian:**
```bash
sudo apt install libsuitesparse-dev libhdf5-dev libblas-dev liblapack-dev zlib1g-dev
```

**Fedora/RHEL:**
```bash
sudo dnf install suitesparse-devel hdf5-devel blas-devel lapack-devel zlib-devel
```

## Environment Configuration (env.cmake)

If dependencies are in non-standard locations, create `env.cmake` from the template:

```bash
cp env.cmake.example env.cmake
```

Edit `env.cmake` to set custom paths:
```cmake
set(CMAKE_PREFIX_PATH "/path/to/suitesparse;/path/to/hdf5")
set(C_COMPILER gcc)  # or specify a full path
```

This file is gitignored and user-specific.

## OpenMP Parallelisation

Build with OpenMP:
```bash
make build-dev OMP=ON OPT=ON SET=RiftingChenin
```

Control thread count at runtime:
```bash
export OMP_NUM_THREADS=8
make run SET=RiftingChenin
```

Key parallelised operations: particle interpolation, rheology evaluation, matrix assembly, CHOLMOD factorisation.

## Common Troubleshooting

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| `SuiteSparse not found` | Missing or wrong path | Install SuiteSparse or set path in `env.cmake`. Check `FindSuiteSparse.cmake`. |
| `HDF5 not found` | Missing HDF5 | Install HDF5 dev package or run `make deps` |
| `HDF5 version mismatch` | Mixed HDF5 installations | Ensure a single consistent HDF5 version; set `HDF5_DIR` in `env.cmake` |
| Segfault on startup | Bad `.txt` file or missing file | Verify `.txt` exists, `Nb_phases` matches `ID` blocks, `Nx`/`Nz` > 0 |
| NaN in output | Unstable parameters | Check `dt`, `min_eta`/`max_eta`, material properties, try `safe_mode = 1` |
| Compilation warnings | Debug mode is verbose | Normal with `-Wall`; use `-DOPT=ON` for clean optimised build |
| Slow execution | No OpenMP or optimisation | Rebuild with `OMP=ON OPT=ON` |

## Debugging with AddressSanitizer (ASAN)

When a simulation segfaults or produces corrupted results, **always try ASAN first** before GDB. ASAN detects heap-buffer-overflows, use-after-free, and stack overflows that cause silent memory corruption in Release mode — often crashing many steps later in unrelated code (e.g., inside `cs_di_compress` or CHOLMOD).

### Building with ASAN

```bash
cmake -DCMAKE_C_FLAGS='-fsanitize=address -fno-omit-frame-pointer -g' \
      -DCMAKE_CXX_FLAGS='-fsanitize=address -fno-omit-frame-pointer -g' \
      -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . -j4
```

### Running with ASAN

```bash
ASAN_OPTIONS='detect_leaks=0' ./MyTest    # disable leak check for faster runs
```

### Known bugs found by ASAN

| Location | Bug | Fix |
|----------|-----|-----|
| `AdvectionRoutines.c` `EvaluateCourantCriterion` | Thermal loop used Vz-grid bounds `(Nx+1, Nz)` to iterate centre-grid arrays `mesh->T` of size `(Nx-1)*(Nz-1)` | Changed bounds to `(Nx-1, Nz-1)` |
| `Solvers.c` `DirectStokesDecoupledComp` | `BCp.type[k] != 30 && != 31` let type=0 cells through with `eqn_p=-1` → wild pointer `i = -1 - neq` | Changed guard to `BCp.type[k] == -1` |
| `Solvers.c` `KillerSolver`, `KSPStokesDecoupled` | Same guard bug as above; also `D1cm0->x[k]` used cell index `k` instead of equation index `i` | Same fix + `[k]` → `[i]` |

## Key Source Files

- Build configuration: `CMakeLists.txt`, `Deps.cmake`, `FindSuiteSparse.cmake`
- Make targets: `makefile`
- Environment template: `env.cmake.example`
- Library build: `MDLIB/CMakeLists.txt`
- Scenario build: `SETS/AddSet.cmake`
