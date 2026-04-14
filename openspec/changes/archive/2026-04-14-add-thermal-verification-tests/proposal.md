## Why

The existing CI thermal tests verify basic diffusion, geotherms, and radiogenic heating on simple setups. There are no tests for (1) oceanic cooling against the classical analytical erf solution, or (2) thermo-elastic coupling — the latter is critical for ongoing work by Lorenzo Candioti and Hannah. Adding these two verification tests closes gaps in thermal solver coverage. Crucially, each test must run with both `thermal_solver = 0` (CHOLMOD) and `thermal_solver = 1` (PCG) and cross-compare their results to verify that the PCG solver produces the same answer as the established CHOLMOD solver.

## What Changes

- Add a CI GTest for **oceanic cooling** (based on `SETS/ThermalDiffusion.c`/`.txt`) that runs the simulation with both `thermal_solver = 0` (CHOLMOD) and `thermal_solver = 1` (PCG), compares each vertical temperature profile against the analytical half-space cooling solution $T(z) = T_m \cdot \mathrm{erf}\!\bigl(\tfrac{-z}{\sqrt{4 t \kappa}}\bigr)$ with $T_m = 1400$ K, $\kappa = k / (\rho C_p)$, and cross-compares the two solver results to ensure they match within a tight tolerance.
- Add a CI GTest for **thermo-elastic coupling** (based on `SETS/ThermoElastic.c`/`.txt`) that runs the simulation with both solvers (`thermal_solver = 0` and `1`), verifies thermo-elastic stress/strain response, and cross-compares CHOLMOD vs PCG results.
- Add a C++ HDF5 data-extraction helper in the test harness to read 1D temperature profiles for comparison with analytical solutions.
- Add gnuplot scripts (not Julia) for visual verification of each test, producing PNGs that can be inspected during development without requiring a Julia environment.
- Add `.txt` parameter files under `TESTS/Thermal/` for both scenarios (CHOLMOD and PCG variants).

## Capabilities

### New Capabilities
- `ci-oceanic-cooling`: GTest verifying half-space oceanic cooling temperature profile against the analytical erf solution, running with `thermal_solver = 0` and `thermal_solver = 1`, cross-comparing CHOLMOD vs PCG results, with gnuplot visual output.
- `ci-thermoelastic-coupling`: GTest verifying thermo-elastic coupling with linearly increasing temperature, running with `thermal_solver = 0` and `thermal_solver = 1`, cross-comparing CHOLMOD vs PCG results, with gnuplot visual output.

### Modified Capabilities
- `ci-thermal-tests`: Add oceanic cooling and thermo-elastic test references to the existing thermal test spec.

## Impact

- **Code**: `TESTS/ThermalTests.cpp` gets new test cases; new `.txt` parameter files under `TESTS/Thermal/`; new gnuplot scripts under `TESTS/` or `VISUAL_TESTS/`.
- **Dependencies**: gnuplot (already available on CI and dev machines); HDF5 C++ API (already linked in test harness via `TestHelpers.h`).
- **CI**: Two new GTest cases added to `CTestTestfile.cmake`; test run time increases by ~10-20 seconds per scenario.
