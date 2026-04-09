## Context

MDOODZ 7.0 has 15 test suites (39 cases) that verify the solver runs correctly and produces qualitatively sensible output, but none measure quantitative error against known analytical solutions. Reviewers expect **L2 error norms** and **grid-convergence order** to demonstrate the code is mathematically verified, not just regression-tested.

The codebase already contains analytical solution code that has never been wired into the CI test suite:
- `eval_anal_Dani()` in `SETS/AnisotropyDabrowski.c` — full 2D Stokes analytical solution (Schmid & Podladchikov, 2003)
- Julia scripts in `misc/JuliaScripts/` — anisotropy director evolution and stress-vs-angle formulas

The test infrastructure uses GoogleTest 1.12.1, HDF5 1.10.7, and a shared `TestHelpers.h` header with helpers for reading fields from HDF5 output. All tests run inside CTest and the CI pipeline must stay under ~30 seconds total.

**Constraints:**
- C11 compilation (no C++17 features in MDLIB; tests are `.cpp` but link against C library)
- Staggered grid: Vx on `VxNodes` (Nx × (Nz-1)), Vz on `VzNodes` ((Nx-1) × Nz), P/T/σ on `Centers` ((Nx-1) × (Nz-1)), τ_xz on `Vertices` (Nx × Nz)
- HDF5 output is in SI units (already re-scaled from non-dimensional internal units)
- The SolVi analytical solution uses non-dimensional coordinates and a unit-viscosity matrix; the test must match the scaling used in the HDF5 output

## Goals / Non-Goals

**Goals:**
- Add a SolVi benchmark test that measures L2 error for Vx, Vz, and P against the analytical solution and verifies second-order spatial convergence
- Add an anisotropy benchmark that verifies director evolution and anisotropic stress computation against analytical formulas
- Upgrade 8 existing test suites with L2 error assertions against their respective analytical solutions
- Provide reusable L2 error computation helpers in `TestHelpers.h`
- Create a theory document (`TESTS/AnalyticalSolutions.md`) explaining each analytical solution
- Keep total CI runtime under 30 seconds

**Non-Goals:**
- Time-convergence (temporal order) testing — only spatial convergence is in scope
- Porting the Crameri topography benchmark — this requires free surface which adds complexity; defer to a future change
- Adaptive mesh refinement testing — MDOODZ uses a fixed grid
- Testing in 3D — MDOODZ is strictly 2D
- Performance benchmarking or profiling

## Decisions

### 1. L2 error computation approach

**Decision:** Implement `readFieldAsArray()` and `computeL2Error()` as static functions in `TestHelpers.h`, following the existing pattern of self-contained helpers.

`readFieldAsArray()` returns a `std::vector<double>` of the full field. `computeL2Error()` takes the numerical array, an analytical function pointer, coordinate arrays, and the field size, returning a normalised L2 norm:

$$L_2 = \sqrt{\frac{\sum_i (u_i^{num} - u_i^{ana})^2}{\sum_i (u_i^{ana})^2}}$$

If the analytical solution is identically zero (e.g. Vz in pure shear), use absolute L2 instead of relative.

**Alternatives considered:**
- *Separate L2 library (`.c`/`.h` pair)*: Would require CMake changes and a separate compilation unit. The existing helpers are all header-only statics; adding a library breaks the established pattern for little benefit.
- *Template-based generic L2*: Over-engineered for C-linked tests. A simple function pointer callback covers all cases.

### 2. SolVi benchmark architecture

**Decision:** Port `eval_anal_Dani()` directly into the SolVi test `.cpp` file as a static function. The test runs the `ViscousInclusion` scenario at a single resolution (e.g. 51×51), then reads HDF5 output parameters (Nx, Nz, coordinates) and computes L2 error for Vx, Vz, P, and σ'_xx.

For grid-convergence, use a second test case that runs 3 resolutions (21×21, 41×41, 81×81) with separate `.txt` parameter files. Compute L2 at each and verify convergence order ≥ 1.8 (allowing margin below the theoretical 2.0 for marker-in-cell noise).

The boundary conditions apply the analytical solution as Dirichlet on E/W faces and type=11 (analytical velocity) on N/S faces, matching the existing `AnisotropyDabrowski.c` pattern.

**Alternatives considered:**
- *Exposing `eval_anal_Dani` from MDLIB as a shared symbol*: Would require modifying the library's public API for test-only purposes. Duplicating the ~70-line function into the test file is simpler and keeps the library interface clean.
- *Running convergence at 4+ resolutions*: Diminishing returns — 3 points suffice to fit a slope and the 81×81 resolution already takes ~2 seconds. More resolutions would push CI time.

### 3. Coordinate reconstruction for L2 comparison

**Decision:** Read coordinate arrays from the HDF5 `Model/` group (`xc_coord`, `zc_coord` for centers; `xg_coord`, `zg_coord` for vertices). These are 1D arrays; reconstruct 2D (x, z) positions with `x = xcoord[i % Ncx]`, `z = zcoord[i / Ncx]` (row-major storage). This avoids recomputing grid geometry and ensures the test uses the exact same coordinates the solver used.

**Alternatives considered:**
- *Recomputing coordinates from Lx, Lz, Nx, Nz*: Fragile — rounding differences between test-reconstructed and solver-used coordinates would introduce spurious error.

### 4. Anisotropy benchmark design

**Decision:** Create two test cases within a single `AnisotropyBenchmarkTests.cpp` file:

1. **Director evolution under simple shear** — Run a homogeneous simple-shear setup with constant anisotropy factor δ for ~50 time steps. Read the director field (Nx, Nz components) from HDF5. Compute the analytical director at the final time using the Mühlhaus ODE integrated in the test, and compare L2 error. The analytical ODE: $\dot{\mathbf{n}} = \mathbf{W}\mathbf{n} - \mathbf{D}\mathbf{n}(\mathbf{n}^T\mathbf{n}) + (\mathbf{n}^T\mathbf{D}\mathbf{n})\mathbf{n}$.

2. **Stress invariant vs angle** — Run a single-step pure-shear setup with a prescribed director angle. Compare the computed stress invariant against the analytical formula: rotate strain rate to the principal anisotropy plane, apply $\mathbf{D} = 2\eta_N \text{diag}(1, 1, 1/\delta)$, rotate stress back, compute J2 invariant. Since this tests constitutive-level behaviour, a single coarse resolution (11×11) suffices.

**Alternatives considered:**
- *Testing anisotropy via SolVi with anisotropic inclusion*: The `AnisotropyDabrowski.c` scenario does this, but the analytical solution (`eval_anal_Dani`) is for isotropic inclusions — the anisotropy only affects the matrix. This would conflate two effects. Testing the anisotropy constitutive law in isolation is cleaner.

### 5. Upgrading existing tests with L2 assertions

**Decision:** Add L2 assertions as **additional** `EXPECT_*` statements in the existing test cases, not as separate test fixtures. This preserves backward compatibility: existing bounds checks remain, L2 checks are added alongside. The analytical formulas are simple closed-form expressions evaluated inline (no separate function needed for 1D cases like $T(z) = T_{top} + \Delta T \cdot z / H$).

For 1D analytical solutions applied to 2D fields (e.g. geotherm T(z), hydrostatic P(z)):
- Read the full 2D field and coordinate arrays from HDF5
- Evaluate the 1D formula at each grid point's z-coordinate
- Compute L2 over the full 2D grid (this works because the 1D solution is constant in x)

**Alternatives considered:**
- *Extracting a single column and comparing 1D*: Loses information — a bug affecting only certain x-positions would be missed. Full-grid L2 is more rigorous and only marginally more code.
- *Creating separate test fixtures for L2*: Would double the number of test files for no architectural benefit.

### 6. L2 threshold selection

**Decision:** Set thresholds empirically by running each test once and measuring the actual L2 error, then choosing a threshold 2–5× above the measured value. This gives headroom for platform variation while still catching genuine regressions. Document all threshold values and their measured baselines in `AnalyticalSolutions.md`.

Expected ranges:
- SolVi at 51×51: L2(Vx) ~ 1e-3 to 1e-2, L2(P) ~ 1e-2 to 1e-1 (pressure is noisier near the inclusion boundary)
- 1D analytical (geotherm, hydrostatic): L2 ~ 1e-6 to 1e-4 (these are trivial solutions the solver should nail)
- Anisotropy director: L2 ~ 1e-3 to 1e-2 (accumulated ODE integration error over many steps)

**Alternatives considered:**
- *Hard-coding theoretical convergence to set thresholds*: Theory gives the rate but not the constant. Empirical measurement is unavoidable.

### 7. Theory documentation structure

**Decision:** Create `TESTS/AnalyticalSolutions.md` with sections for each benchmark: problem statement, analytical formula, derivation sketch or reference, how it maps to the MDOODZ test, and the measured L2 baseline. Reference this from `TESTS/README.md` in a new "Analytical Benchmarks" section. Use LaTeX math notation (GitHub renders it).

## Risks / Trade-offs

**[SolVi L2 thresholds too tight or too loose]** → Measure actual error on the CI platform (Ubuntu, gcc) during implementation. Set thresholds as `measured * 3` initially. If CI flakes, widen by another factor. Document the measurement methodology so future maintainers can re-calibrate.

**[Marker-in-cell noise affects convergence order]** → The stochastic nature of particle positions means the measured convergence order may fluctuate between runs. Mitigate by accepting order ≥ 1.8 instead of exactly 2.0, and by using the same random seed (set via the `.txt` parameter file).

**[Anisotropy HDF5 output may not include director field]** → Verify that the director components (nx, nz) are written to HDF5 when anisotropy is enabled. If not, the director benchmark may need to check stress output instead. Fall back to stress-only comparison as a contingency.

**[CI runtime budget]** → SolVi at 3 resolutions (21+41+81) may take ~5 seconds. Combined with anisotropy benchmarks (~2 seconds) and existing 15 suites (~10 seconds), total should stay under 25 seconds. Monitor after implementation.

**[Coordinate system mismatch between analytical and numerical]** → The analytical solution uses mathematical convention (x right, z up) while HDF5 uses the same. Reading coordinates from HDF5 eliminates this risk entirely.

## Open Questions

1. **Does HDF5 output include the director field (nx, nz)?** If not, the anisotropy director-evolution benchmark must compare stress instead. This needs to be verified during implementation by checking `HDF5Output.c` for anisotropy-related writes.

2. **What is the exact staggered-grid layout for Vx and Vz nodes?** The L2 computation must use the correct coordinate arrays (`xvx_coord`/`zvx_coord` vs `xvz_coord`/`zvz_coord`). Verify array sizes match Nx×(Nz-1) for Vx and (Nx-1)×Nz for Vz during implementation.

3. **Should the convergence-order test be a separate CTest executable or a test case within SolViBenchmarkTests?** Keeping it in one executable simplifies CMake but means the multi-resolution run is all-or-nothing. A single executable with parameterised test cases (GTest `TEST_P` or multiple `TEST_F`) is the likely approach.
