# Design: blankenbach-steady-state

## Overview

Add a long-running (non-CI) Blankenbach Case 1a steady-state test to the existing BlankenBench suite. Two deliverables: a parameter file (`BlankenBenchSteady.txt`) and a new GTest (`NusseltAndVrms`) that computes Nusselt number and Vrms from the final HDF5 output and compares against published reference values.

The parameter file is tuned for faster convergence to steady state (higher Courant number → larger dt → fewer steps). The test must be run with 8 OpenMP threads.

## Architecture

### Reuse of existing fixture

The existing `BlankenBench` test fixture in `BlankenBenchTests.cpp` already defines all the physics callbacks (SetPhase, SetTemperature, SetDensity, SetBCVx, SetBCVz, SetBCT, SetBCPType). The new test is another `TEST_F(BlankenBench, ...)` in the same file — no new fixture or file needed.

### Single combined test

A single test `BlankenBench.NusseltAndVrms` both runs the simulation and checks results. This avoids running a multi-hour simulation twice. The test:

1. Calls `RunMDOODZ("BlankenBench/BlankenBenchSteady.txt", &setup)`
2. Opens the final HDF5 output (filename depends on Nt)
3. Computes Nu and Vrms
4. Asserts both against published values

### OpenMP parallelism

MDOODZ picks up thread count from `OMP_NUM_THREADS` at runtime. The test must be run with 8 threads:

```bash
BLANKENBACH_STEADY=1 OMP_NUM_THREADS=8 ./BlankenBenchTests --gtest_filter=BlankenBench.NusseltAndVrms
```

The test **verifies** the thread count before launching the simulation. When compiled with `-DOMP=ON`, the test calls `omp_get_max_threads()` inside an `#ifdef _OPENMP` guard and asserts it is at least 4. If the build does not have OpenMP, the test skips. The `BLANKENBACH_STEADY=1` environment variable ensures the test never runs accidentally during `ctest`.

```cpp
  // Guard: only run when explicitly requested
  const char *run_steady = std::getenv("BLANKENBACH_STEADY");
  if (!run_steady || std::string(run_steady) != "1") {
    GTEST_SKIP() << "Set BLANKENBACH_STEADY=1 to run this long test (~2-4 hours)";
  }

#ifdef _OPENMP
#include <omp.h>
  int nthreads = omp_get_max_threads();
  printf("BlankenBench: Running with %d OpenMP threads\n", nthreads);
  ASSERT_GE(nthreads, 4) << "Need at least 4 threads for reasonable runtime";
#else
  GTEST_SKIP() << "Build without OpenMP (-DOMP=ON). Skipping long-running test.";
#endif
```

This is documented in README.md.

### No CI registration

The test binary `BlankenBenchTests` is already registered with `add_test()` for the existing smoke tests. The new `NusseltAndVrms` test is only run manually via `--gtest_filter=BlankenBench.NusseltAndVrms`. No changes to CMakeLists.txt are required (the executable already builds with all TEST_F cases and `file(COPY ...)` already copies the whole `BlankenBench/` directory).

## Detailed Design

### 1. Parameter file — `BlankenBenchSteady.txt`

Copy of `BlankenBench.txt` with changes tuned for fastest convergence to steady state:

| Parameter | Smoke test | Steady state | Rationale |
|---|---|---|---|
| `Nt` | 500 | 100000 | Enough steps to reach ~3.6 diffusion times |
| `Courant` | 0.25 | 0.5 | Larger dt → fewer steps to cover same physical time |
| `writer_step` | 500 | 100000 | Write only the final output |
| `istep` | 001000 | 200000 | Prevent restart-file interference |

Raising Courant from 0.25 to 0.5 roughly doubles dt per step, halving the wall-clock time. For isoviscous Stokes (linear problem, Picard converges in few iterations) this is safe. The CFL for thermal diffusion at 41×41 is $\Delta t_{\text{diff}} = \Delta x^2 / (2\kappa) \approx 4 \times 10^{12}$ s, well above the advective dt.

All physics parameters remain identical to the smoke test. These are the **initial** settings — they may be refined after experimentation (see §7).

### 2. Nusselt number computation

**Definition**: $Nu = -\frac{H}{\Delta T} \langle \frac{\partial T}{\partial z} \rangle_{z=\text{top}}$

**Implementation**:

```
Read T field: readFieldAsArray(outFile, "Centers", "T")
  → ncx × ncz array, stored in Kelvin (SI), row-major, bottom-to-top

Top boundary heat flux via one-sided finite difference:
  For each column i in [0, ncx):
    T_centre = T[(ncz-1) * ncx + i]     // top cell-centre row
    dTdz_i   = (T_top_K - T_centre) / (dz / 2)

  Nu = -(H / DeltaT) * mean(dTdz)
```

Where:
- `T_top_K = 273.15` (0 °C in Kelvin)
- `DeltaT = 1000.0` K
- `H = 1.0e5` m
- `dz = H / ncz = 1e5 / 40 = 2500.0` m

This matches the existing informational Nu computation in `ConvectionDevelops`.

**Reference**: Nu = 4.884409 ± 5%

### 3. Vrms computation

**Definition**: $V_\text{rms} = \sqrt{\frac{1}{A} \int_\Omega (v_x^2 + v_z^2) \, dA}$

**Non-dimensionalisation**: HDF5 stores velocities in m/s. Blankenbach non-dimensional Vrms uses the thermal diffusivity scale:

$$\kappa = \frac{k}{\rho C_p} = \frac{3.3}{3300 \times 1250} = 8.0 \times 10^{-7} \; \text{m}^2/\text{s}$$

$$V_\text{rms}^{(\text{nd})} = V_\text{rms}^{(\text{SI})} \cdot \frac{H}{\kappa}$$

**Implementation**:

```
Read Vx: readFieldAsArray(outFile, "VxNodes", "Vx")  → Nx × (Nz+1) = 41 × 42
Read Vz: readFieldAsArray(outFile, "VzNodes", "Vz")  → (Nx+1) × Nz = 42 × 41

Interpolate to cell centres (ncx × ncz = 40 × 40):
  For each cell (i, j):
    vx_c = 0.5 * (Vx[j * Nx + i] + Vx[j * Nx + (i+1)])
    vz_c = 0.5 * (Vz[j * (Nx+1) + i] + Vz[(j+1) * (Nx+1) + i])

Uniform grid → equal cell areas:
  Vrms_SI = sqrt( sum(vx_c^2 + vz_c^2) / (ncx * ncz) )

Non-dimensionalise:
  Vrms_nd = Vrms_SI * H / kappa
```

Note on Vx grid layout: Vx is defined on a grid with Nx points in x and (Nz+1) points in z. The row index j runs over z-faces (Nz+1 = 42 rows), and column index i runs over cell centres in x (Nx = 41 columns). Averaging adjacent columns gives the cell-centre value.

Similarly, Vz has (Nx+1) columns and Nz rows.

**Reference**: Vrms = 42.864947 ± 5%

### 4. Tolerances

At 41×41 resolution, the Blankenbach reference values are obtained with much finer grids. A 5% relative tolerance is appropriate:

| Quantity | Published | Tolerance |
|---|---|---|
| Nu | 4.884409 | ±5% → [4.640, 5.129] |
| Vrms | 42.864947 | ±5% → [40.722, 45.008] |

Use `EXPECT_NEAR` with absolute delta = published × 0.05.

### 5. Test output

The test prints computed values even on pass, for manual inspection:

```
BlankenBench Steady State:
  Nu     = X.XXXX  (published: 4.884409)
  Vrms   = XX.XXXX (published: 42.864947)
```

### 6. README documentation

The Suite 18 section of README.md must clearly state:

1. **This test is not part of CI** — it takes hours, not seconds
2. **It is fully reproducible** — same code, same parameter file, deterministic result
3. **Exact run command** with `OMP_NUM_THREADS=8` and `--gtest_filter`
4. **Expected output** — what Nu/Vrms values to expect
5. **How to re-run** if parameters are tuned

Example section to add:

```markdown
### Steady-State Benchmark (manual, not in CI)

A long-running test that reaches thermal steady state and compares
Nu and Vrms against published Blankenbach et al. (1989) Case 1a values.

**Run with 8 threads:**

    cd build/TESTS
    BLANKENBACH_STEADY=1 OMP_NUM_THREADS=8 ./BlankenBenchTests --gtest_filter=BlankenBench.NusseltAndVrms

Expected runtime: ~2–4 hours at 41×41 with Courant=0.5.
Published references: Nu = 4.884409, Vrms = 42.864947.
```

### 7. Experimentation strategy

The parameter file is a starting point. After the first run, we can tune:

| Knob | Effect | Risk |
|---|---|---|
| **Courant** 0.5 → 0.75 | ~33% fewer steps | May cause advection overshooting |
| **Courant** 0.5 → 1.0 | ~50% fewer steps | Likely too aggressive for markers |
| **nit_max** 10 → 5 | Fewer Picard iterations per step | Isoviscous converges fast; likely safe |
| **Nx/Nz** 41 → 31 | ~40% fewer cells, faster per step | Coarser grid → larger Nu/Vrms error |
| **Nx/Nz** 41 → 61 | Better accuracy | More expensive per step and needs more steps |

The test tolerances (±5%) accommodate resolution-level discretisation error. After an initial successful run, we can tighten or loosen as needed.

### 8. File changes summary

| File | Action |
|---|---|
| `TESTS/BlankenBench/BlankenBenchSteady.txt` | **New** — parameter file |
| `TESTS/BlankenBenchTests.cpp` | **Edit** — add `TEST_F(BlankenBench, NusseltAndVrms)` |
| `TESTS/CMakeLists.txt` | **Edit** — link `OpenMP::OpenMP_CXX` to BlankenBenchTests when OMP=ON |
| `TESTS/README.md` | **Edit** — document steady-state benchmark in Suite 18 section |
| `MDLIB/InputOutput.c` | **Edit** — increase `Nb_part_max` multiplier from 4.1 to 10.0 |
