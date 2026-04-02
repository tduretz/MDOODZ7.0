---
description: "MDOODZ test suite — CI unit tests (GTest), visual regression tests (Eigen/gnuplot/HDF5), test authoring patterns, CMakeLists integration, build commands, CI workflow, and test coverage map."
---

# Skill: Testing Guide

## Overview

MDOODZ has two test tiers:

| Tier | Location | Framework | Build Flag | Run Command |
|------|----------|-----------|-----------|-------------|
| CI unit tests | `TESTS/` | GTest (FetchContent 1.12.1) | `-DTEST=ON` | `make run-tests` |
| Visual regression tests | `VISUAL_TESTS/` | C++17, Eigen3, HDF5pp, gnuplot | `-DVIS=ON` | `make run-vis` |

---

## CI Unit Tests (`TESTS/`)

### File Structure

Each test group has:
- **`TESTS/<Group>Tests.cpp`** — GTest fixture + test cases
- **`TESTS/<Group>/`** — directory of `.txt` parameter files (one per test case)
- **`TESTS/TestHelpers.h`** — shared HDF5 reading helpers

### Fixture Pattern

```cpp
#include "gtest/gtest.h"
extern "C" {
  #include "mdoodz.h"
}
#include "TestHelpers.h"

class MyTest : public ::testing::Test {
protected:
  // Setup callbacks as static members
  static int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    const double radius = input->model.user1 / input->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius)
      return 1;
    return 0;
  }
  static double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
  }

  void RunModel(const char *paramFile) {
    MdoodzSetup setup = {
      .SetParticles = &(SetParticles_ff){
        .SetPhase   = SetPhase,
        .SetDensity = SetDensity,
      },
      .SetBCs = &(SetBCs_ff){
        .SetBCVx = SetPureOrSimpleShearBCVx,
        .SetBCVz = SetPureOrSimpleShearBCVz,
      },
    };
    RunMDOODZ(paramFile, &setup);
  }
};
```

### Assertion Helpers (`TestHelpers.h`)

| Function | Returns | Purpose |
|----------|---------|---------|
| `getStepsCount(file)` | `int` | Number of completed time steps |
| `getFinalResidual(file, name)` | `double` | Named residual (e.g. `"rx_abs"`) |
| `getMaxFieldValue(file, group, dataset)` | `double` | Max value in HDF5 dataset |
| `getMinFieldValue(file, group, dataset)` | `double` | Min value in HDF5 dataset |
| `getMeanFieldValue(file, group, dataset)` | `double` | Mean value in HDF5 dataset |
| `getModelParam(file, index)` | `double` | Model parameter by index |
| `getFieldValueAt(file, group, dataset, index)` | `double` | Single value at 1D index |
| `getMaxAbsFieldValue(file, group, dataset)` | `double` | Max absolute value in dataset |

### Parameter File Conventions

- Use `user0`–`user4` for test-specific parameters (temperatures, radii, etc.)
- Grid: 21×21–31×31 for fast CI tests, 41×41 for free surface tests, 51×51 for thermal/solver tests
- Steps: 1–5 for basic convergence, up to 20 for evolution checks
- Scaling: use realistic values (`eta=1e22`, `L=1e4`) when using real flow laws (e.g. `linv=40`, `expv=40`)
- Unit scaling (`eta=1`, `L=1`) only for constant-viscosity or user-defined rheology (`pwlv=1`)

### Required Parameter File Settings for Tests

Every test `.txt` file **must** include:

| Parameter | Value | Why |
|-----------|-------|-----|
| `writer_subfolder` | `= <TestCaseName>` | Tests expect output at `<TestName>/Output*.gzip.h5` |
| `coh_soft` | `= 0` (unless testing softening) | Default is 1; if `C == Ce` MDOODZ calls `exit(122)` |

### Common Pitfalls

- **`Slim >= C`**: MDOODZ checks `Slim >= C` at startup. If violated, it prints "Upper stress limiter of phase X is lower than cohesion" and exits. Use `Slim = 1e91` when `C = 1e90` (disabled plasticity), or `Slim = 1e100` for absolute safety.
- **Softened end-value parameter names**: In `.txt` files, use `Ce` (not `C_end`), `phie` (not `phi_end`), `psie` (not `psi_end`). The code reads these via `ReadMatProps(fin, "Ce", ...)` (see `InputOutput.c:1318`).
- **`coh_soft` defaults to 1**: If not set in `.txt`, cohesion softening is enabled. When `C == Ce` (both defaulting to the same value), MDOODZ prints "Please set a difference in cohesion" and calls `exit(122)`. Always set `coh_soft = 0` explicitly unless the test needs softening.
- **GBS flow law (`gbsv=40`) + `elastic=0`**: When `elastic=0`, the code internally sets `G = 1e1` (line 468, `RheologyDensity.c`), which can produce NaN viscosities with GBS. Use `elastic = 1` for GBS tests.
- **Pure shear convergence**: With pure shear BCs and uniform/homogeneous viscosity, the solver may converge in **0 Newton iterations** because the initial guess is analytically exact. Use `ASSERT_GE(stepsCount, 0)` (not `ASSERT_GT`) and validate physics via field value checks instead.
- **`max_eta` caps viscosity contrast**: The non-dimensional `max_eta` cutoff (e.g. `1e6`) is applied after viscosity computation. If both phases' physical viscosities exceed `max_eta * eta_ref`, they are both clamped to the same value, erasing any grain-size or flow-law contrast. Do not assert `maxEta > minEta` when both phases may hit the cap.
- **Plasticity requires sufficient trial stress**: With `elastic=1`, the Maxwell visco-elastic trial stress is approximately `σ_trial ≈ 2·G·Δt·ε̇`. If `σ_trial < C·cos(φ)`, plasticity never activates and `eII_pl = 0`. Ensure the strain rate is high enough that the trial stress exceeds the yield stress — this depends on the shear modulus G and time step dt.
- **Strain softening requires yielding first**: Cohesion softening (`coh_soft=1`) only activates when plastic strain accumulates. If the trial stress is below the yield surface, no plastic strain → no softening at all. Same strain-rate requirement as DruckerPrager yield tests.
- **No-slip BCs + pure_shear_ALE = singular system**: No-slip enforces `V=0` on all boundaries, but `pure_shear_ALE=1` tries to apply kinematic BCs (`Vx=ε̇·x`). The conflicting conditions produce a singular stiffness matrix → segfault in Cholesky. For no-slip tests, set `pure_shear_ALE=0` and use gravity to drive flow.
- **Picard residual tolerance**: Picard (fixed-point) iteration converges linearly and may plateau at residuals slightly above the nonlinear tolerance. Use `EXPECT_LT(finalRes, 1e-8)` rather than `1e-9` for Picard tests to avoid flaky assertions.
- **Free surface tests need ≥ 41×41 grid**: With `free_surface=1`, coarser grids (21×21, 31×31) can cause CHOLMOD "not positive definite" failures due to the modified boundary stiffness. 41×41 is the minimum reliable resolution.
- **CRLF line endings (WSL/Windows)**: When `git config core.autocrlf=true` or when using `rsync` from Windows NTFS, C source files can get CRLF endings. The parameter file parser (`InputOutput.c`) fails to parse `\r\n` correctly, producing NULL pointers → SIGSEGV. Fix: `git checkout -- MDLIB/` and `sed -i 's/\r$//'` on test files.
- **Gravity + compressibility + dilatancy = CHOLMOD failure**: The combination of `gz < 0`, `compressible=1`, and `psi > 0` produces an ill-conditioned Stokes matrix. Avoid gravity in compressibility tests, or set `psi=0`.
- **Gravity + near-zero bkg_strain_rate = CHOLMOD failure**: With gravity enabled and very low strain rate (or `bkg_strain_rate=0`), the stiffness matrix can become ill-conditioned. Use dimensional scaling (`eta=1e22`, `L=1e4`) and standard `bkg_strain_rate=1e-14` for gravity tests.
- **Finite strain F stays at identity in ALE mode**: In `pure_shear_ALE=1`, the grid deforms with the material, so markers don't move relative to grid → F accumulates slowly. Don't expect analytical F values; check `F > 0` and `det(F) ≈ 1` instead.
- **type=13 is shear stress BC, not normal stress**: The `constant_shear_stress` BC (type=13) on E/W boundaries controls the tangential (shear) component, not the normal traction. You cannot apply normal stress BCs this way. Use `SetPureShearBCVx` which correctly uses type=13 on N/S (free-slip) and type=0 on E/W (velocity Dirichlet).
- **Thermal tests require SetBCT callback**: Tests with `thermal=1` and shear heating need explicit thermal boundary conditions via `SetBCT` in the `SetBCs_ff` struct. Without it, temperature may not evolve as expected.

### Valid Flow Law Indices

| Abbreviation | Index | Mineral / Law |
|-------------|-------|---------------|
| `pwlv=1` | User-defined (eta0, npwl, Qpwl) | Works with any scaling |
| `pwlv=40` | Olivine Dry Dislocation | Requires realistic scaling |
| `linv=40` | Olivine Dry Diffusion | Requires realistic scaling |
| `expv=40` | Olivine Peierls (Evans & Goetze) | Requires realistic scaling |
| `gbsv=40` | Olivine Low-T GBS | Requires realistic scaling |

> **Warning**: `linv=1`, `expv=1`, `gbsv=1` do NOT exist in `FlowLaws.c` and will call `exit(12)`.

### Adding a New CI Test

1. Create `TESTS/<Group>/<TestName>.txt` parameter file
2. Add test case in `TESTS/<Group>Tests.cpp`
3. In `TESTS/CMakeLists.txt`:
   ```cmake
   file(COPY ${PROJECT_SOURCE_DIR}/TESTS/<Group> DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
   add_executable(<Group>Tests <Group>Tests.cpp)
   target_link_libraries(<Group>Tests LINK_PUBLIC mdoodz GTest::gtest_main)
   add_test(NAME <Group>Tests COMMAND <Group>Tests)
   ```

### Existing CI Tests (15 suites, 39 test cases)

| File | Tests | What It Covers |
|------|-------|----------------|
| `NewtonIterationConvergence.cpp` | 8 tests (linear/nonlinear × pure/simple × iso/aniso) | Newton/Picard convergence |
| `RheologyCreepTests.cpp` | DiffusionCreepLinear, DiffusionCreepGrainSize, PeierlsCreep, GBSCreep, CompositePwlLin, **PowerLawCreep** | Flow law mechanisms |
| `PlasticityTests.cpp` | DruckerPragerYield, DruckerPragerNoYield, StrainSoftening, StressLimiter | Yield criteria |
| `ThermalTests.cpp` | GaussianDiffusion, SteadyStateGeotherm, RadiogenicHeat, DirichletBC | Energy equation |
| `BoundaryConditionTests.cpp` | PeriodicSimpleShear, FreeSlipPureShear, NoSlip | BC types |
| `SolverModeTests.cpp` | PicardConvergence, AdaptiveDt | Solver modes |
| `ViscoElasticTests.cpp` | StressAccumulation | Maxwell visco-elasticity |
| `DensityTests.cpp` | ThermalExpansion, HydrostaticPressure | Density EoS, hydrostatics |
| `ShearHeatingTests.cpp` | ViscousDissipation | Shear heating H=2ηε̇² |
| `FreeSurfaceTests.cpp` | SinkingBlock | Free surface, topography |
| `VelocityFieldTests.cpp` | PureShearVelocity | Velocity symmetry check |
| `CompressibilityTests.cpp` | DilatantFlow | Elastic compressibility |
| `FiniteStrainTests.cpp` | PureShearStrain | F tensor, det(F)≈1 |
| `NeumannBCTests.cpp` | StressBCPureShear | Stress BC type=13 |
| `ConvergenceRateTests.cpp` | NewtonPwl, PicardPwl, NewtonFasterThanPicard | Newton vs Picard comparison |

---

## Visual Regression Tests (`VISUAL_TESTS/`)

### File Structure

Each visual test has:
- **`.cpp`** — Plot functions reading HDF5, writing `.dat`, calling gnuplot
- **`.txt`** — MDOODZ parameter file
- **`.gnu`** — gnuplot script
- **`*Reference.h5`** or **`*_REF/`** — reference HDF5 output

### C++ Class Pattern (in `main.cpp`)

Scenario classes define static callbacks and a `run()` method:

```cpp
class MyScenario {
  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) { ... }
  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) { ... }

  public: void run() {
    MdoodzSetup setup = {
      .SetParticles = new SetParticles_ff{ .SetPhase = SetPhase, .SetDensity = SetDensity },
      .SetBCs = new SetBCs_ff{ .SetBCVx = SetPureOrSimpleShearBCVx, .SetBCVz = SetPureOrSimpleShearBCVz },
    };
    RunMDOODZ("MyScenario.txt", &setup);
  }
};
```

### Plot Function Pattern (in `<Name>.cpp`)

```cpp
#include "HDF5pp.h"
#include <Eigen/Core>
#include "visual-tests.h"

void PlotMyScenario() {
  H5p::File file = H5p::File("MyScenario.gzip.h5", "r");
  std::vector<double> params = file.read<std::vector<double>>("/Model/Params");
  const int nx = (int) params[3], nz = (int) params[4];
  // Read field, create Map<MatrixXf>, write .dat, call gnuplot
}
void PlotMyScenarioReference() { /* same, reading *Reference.h5 */ }
```

### L2 Norm Regression Checks

`visual-tests.h` provides:
- `relativeL2Error(computed, reference)` — returns `||A-B||/||B||`
- `VISUAL_ASSERT(condition, label)` — prints PASS/FAIL, sets `g_visualTestFailed`
- `VISUAL_ASSERT_L2(computed, reference, tol, label)` — convenience macro
- `CompareH5Field(compFile, refFile, fieldPath, tol, label)` — loads both files, extracts field, runs L2 check

### Adding a New Visual Test

1. **Create files**: `<Name>.cpp`, `<Name>.txt`, `<Name>.gnu`
2. **Add class** to `main.cpp` with callbacks and `run()` method
3. **Add to `RunTestCases()`**: call `run()`, then `rename("Output*.gzip.h5", "<Name>.gzip.h5")`
4. **Declare** `Plot<Name>()` and `Plot<Name>Reference()` in `visual-tests.h`
5. **Call** both plot functions from `main()`
6. **Add L2 check** using `CompareH5Field()` in `main()`
7. **Update CMakeLists.txt**: add `.cpp` to `add_executable`, `configure_file` for `.txt`, `.gnu`, `.h5`
8. **Generate reference**: `make generate-refs` (or `make build-dev VIS=ON && make run-vis`)

### Build & Run Commands

```bash
# CI tests
cmake -B ./cmake-build -DTEST=ON && cmake --build ./cmake-build
make run-tests

# Visual tests
cmake -B ./cmake-build -DVIS=ON && cmake --build ./cmake-build
make run-vis

# Generate reference HDF5 files
make generate-refs
```

---

## Analytical Benchmark Authoring

### L2 Error Framework

`TestHelpers.h` provides helpers for quantitative error measurement:

| Function | Purpose |
|----------|---------|
| `readFieldAsArray()` | Read HDF5 float field → `std::vector<double>` |
| `readCoordArray()` | Read 1D coordinate array from `Model/` group |
| `computeL2Error()` | Relative L2 norm with fallback to absolute when analytical ≈ 0 |

### Writing an L2 Benchmark Test

1. Run the simulation with `RunMDOODZ()`
2. Read numerical field from HDF5 with `readFieldAsArray()`
3. Construct analytical solution at grid coordinates
4. Call `computeL2Error(numerical, analytical)` and assert `EXPECT_LT(L2, threshold)`

### Grid-Convergence Order Testing

Run the same scenario at 3+ resolutions, compute L2 at each, then:
```cpp
double order = log(L2_coarse / L2_fine) / log(h_coarse / h_fine);
EXPECT_GE(order, expected_order);
```

### Threshold Calibration

1. Measure the actual L2 error by running the test
2. Set threshold to 2–5× above measured value
3. For convergence order: set threshold ~30–50% below measured order

### Common Pitfalls

- **Staggered grid sizes**: Vx is Nx×(Nz-1), Vz is (Nx-1)×Nz, P/T/σ are (Nx-1)×(Nz-1)
- **HDF5 scaling**: All fields are written in SI units (re-scaled from non-dimensional)
- **Marker-in-cell noise**: Reduces effective convergence order near material boundaries
- **ALE mode**: `pure_shear_ALE=1` deforms the mesh; avoid with large strain rates in tests
- **Heterogeneous domains**: Mean-field comparisons cancel out in domains with inclusions; use L2 or max

### Marker Density vs Accuracy (Experimental Finding)

Tested on SolVi benchmark (51×51 grid, inclusion η_c=1000, matrix η_m=1):

| Nx_part × Nz_part | Markers/cell | L2(Vx) | L2(P) | Runtime vs baseline |
|--------------------|-------------|--------|-------|-------------------|
| 4×4 (default) | 16 | 4.5e-2 | 7.5e-1 | 1.0× |
| 8×8 | 64 | 4.6e-2 | 6.6e-1 | 2.1× |
| 16×16 | 256 | 4.6e-2 | 6.6e-1 | 6.8× |

**Key finding: increasing markers per cell gives negligible accuracy improvement** (<12% for P, 0% for Vx) at massive runtime cost (up to 6.8×). The P convergence order actually worsens with more markers (0.75 → 0.59), meaning the small L2 reduction at fixed resolution does not carry over to finer grids.

The L2 error is dominated by the FD truncation error at the viscosity discontinuity, not by marker interpolation noise. More markers produce a better statistical average of the same inaccurate stencil at boundary cells.

**Recommendation**: Keep `Nx_part=4, Nz_part=4` for all tests. To improve accuracy, use cell-face harmonic viscosity averaging in the FD stencil (not yet implemented) rather than increasing marker count.

See `TESTS/SolViMarkerComparison.cpp` for the full experimental code and data.

### Viscosity Averaging — `eta_average` (Experimental Finding)

The `eta_average` parameter controls how marker viscosities are averaged to grid nodes **within each cell**: `0` = arithmetic (default), `1` = harmonic, `2` = geometric. This is phase-mixture averaging, **not** cell-face averaging in the FD stencil.

Tested on SolVi benchmark (4×4 markers/cell, 51×51 grid):

| `eta_average` | Method | L2(Vx) | L2(P) | Vx order (41→81) | P order (41→81) |
|--------------|--------|--------|-------|-------------------|------------------|
| 0 | Arithmetic | 4.5e-2 | 7.5e-1 | **0.97** | **0.75** |
| 1 | Harmonic | 3.9e-2 | 9.2e-1 | 0.93 | 0.44 |
| 2 | Geometric | 4.1e-2 | **6.1e-1** | 0.94 | 0.72 |

**Key findings**:
- **Geometric** gives the best absolute P error (18% lower than arithmetic) but slightly lower Vx convergence order
- **Harmonic** is worst for P (both absolute error and convergence order)
- **Arithmetic** gives the best convergence orders overall
- The `eta_average` parameter operates on marker-to-node interpolation within cells. True cell-face harmonic averaging in the stencil (which literature suggests adds ~0.3 convergence order for P) would require modifying `StokesAssemblyDecoupled.c`

**Recommendation**: Keep `eta_average=0` (arithmetic) for CI tests — it delivers the best convergence rates. Use `eta_average=2` (geometric) when absolute accuracy at a fixed resolution matters more than asymptotic convergence.
