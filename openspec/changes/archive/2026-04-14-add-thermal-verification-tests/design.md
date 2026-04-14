## Context

The thermal test suite (`TESTS/ThermalTests.cpp`) already covers Gaussian diffusion, steady-state geotherms, radiogenic heating, and Dirichlet BCs. Each test has CHOLMOD and PCG variants, with cross-comparison of L2 errors. The VISUAL_TESTS directory hosts C++ executables that extract HDF5 data and call gnuplot to produce PNGs.

Two new verification scenarios are needed: oceanic cooling (half-space erf analytical solution) and thermo-elastic coupling (linearly increasing temperature). Both must run with `thermal_solver = 0` (CHOLMOD) and `thermal_solver = 1` (PCG) and cross-compare.

## Goals / Non-Goals

**Goals:**
- Verify the thermal solver reproduces the classical oceanic cooling erf profile to within a quantified L2 tolerance
- Verify thermo-elastic coupling produces correct stress response to a linearly increasing temperature
- Cross-compare CHOLMOD vs PCG results for both scenarios to confirm the PCG solver is equivalent
- Produce gnuplot-based visual evidence (PNGs) for each test

**Non-Goals:**
- Full oceanic geodynamic scenario (no advection, no mechanical coupling for the cooling test)
- Testing Julia visualisation scripts — gnuplot is the verification tool
- Convergence rate studies with grid refinement
- Performance benchmarking of these tests

## Decisions

### 1. Add tests to existing `ThermalTests.cpp` rather than a new file

**Rationale:** All existing thermal verification tests live in a single `ThermalTests.cpp` with a shared `Thermal` test fixture. The fixture already provides `SetPhase`, `SetTemperature`, `SetDensity`, and `SetBCT` callbacks. Adding new `TEST_F` cases follows the established pattern and avoids CMake complexity.

**Alternative considered:** Separate `OceanicCoolingTests.cpp` and `ThermoElasticTests.cpp` — rejected because the fixture setup is identical and the tests are small.

### 2. Oceanic cooling: override callbacks in the test, don't link the SETS .c file

**Rationale:** The existing tests define their own inline callbacks (e.g. `SetTemperature`, `SetBCT`) rather than importing from `SETS/`. This keeps tests self-contained and avoids coupling to scenario files that may change. We'll replicate the essential physics from `SETS/ThermalDiffusion.c` in the test.

**Key physics to replicate:**
- Domain: 1D column (narrow in x), deep enough for the erf profile (e.g. 200 km)
- Initial T: mantle temperature `T_m = 1400 K` everywhere (or 1400 °C + 273.15 K per the Julia script)
- BCs: T = 0 °C at top (surface), T = T_m at bottom, Neumann (zero flux) on sides
- Run for `Nt` steps with a large enough `dt` to develop a visible cooling profile
- No mechanical coupling (`mechanical = 0`), advection off

**Analytical solution:** $T(z, t) = T_m \cdot \mathrm{erf}\!\left(\frac{-z}{\sqrt{4 t \kappa}}\right)$ where $\kappa = k / (\rho C_p)$

### 3. Thermo-elastic: use `FixTemperature` callback with linear ramp

**Rationale:** The test must exercise the `fix_temperature` code path (which overrides particle temperatures at each step). This is exactly what `SETS/ThermoElastic.c` does. We replicate the `FixTemperature` callback inline.

**Key physics:**
- Circular inclusion in a box, pure shear BCs, elastic + compressible
- `FixTemperature` returns `273.15/scaling.T + 0.1 * time` (linearly increasing)
- Verify: pressure change is consistent with $\Delta P \approx -K \alpha \Delta T$ (thermal pressure from thermal expansion)
- Cross-compare CHOLMOD vs PCG temperature and pressure fields

### 4. Parameter file convention: `<TestName>.txt` and `<TestName>PCG.txt`

**Rationale:** This is the existing pattern (see `GaussianDiffusion.txt` / `GaussianDiffusionPCG.txt`). The PCG variant differs only by adding `thermal_solver = 1` and changing `writer_subfolder`.

Files to create:
- `TESTS/Thermal/OceanicCooling.txt` (CHOLMOD)
- `TESTS/Thermal/OceanicCoolingPCG.txt` (PCG)
- `TESTS/Thermal/ThermoElastic.txt` (CHOLMOD)
- `TESTS/Thermal/ThermoElasticPCG.txt` (PCG)

### 5. Gnuplot scripts for visual output in `VISUAL_TESTS/`

**Rationale:** The user explicitly requested gnuplot over Julia. The `VISUAL_TESTS/` directory already contains `.gnu` scripts and a C++ `main.cpp` / `HDF5pp.h` infrastructure. We'll add:

- A C++ extraction program (or extend `main.cpp`) that reads the final HDF5 output, extracts a vertical T profile (column at `ix = ncx/2`), computes the analytical erf solution, and writes both to a CSV.
- A gnuplot script `OceanicCooling.gnu` that plots numerical vs analytical T(z) for both CHOLMOD and PCG.
- A C++ extraction for thermo-elastic: reads final T and P fields, computes expected ΔP from ΔT, writes to CSV.
- A gnuplot script `ThermoElastic.gnu` that plots T evolution and ΔP comparison.

### 6. Profile extraction helper: add `readColumnProfile` to `TestHelpers.h`

**Rationale:** The oceanic cooling test needs a 1D vertical column of T values at a fixed x-index. The existing helpers (`readFieldAsArray`, `readCoordArray`) read flat arrays. We need a helper that extracts `T[ix_mid, :]` from the flat `ncx × ncz` array. This is a small utility that benefits other future tests.

```cpp
static std::vector<double> readColumnProfile(
    const char *hdf5FileName, const char *groupName,
    const char *datasetName, int ix, int ncx, int ncz);
```

### 7. L2 tolerance values

- **Oceanic cooling vs analytical:** `L2 < 5e-2` (FD discretization on a coarse grid; the existing Gaussian L2 test uses `2e-2` on 51×51)
- **CHOLMOD vs PCG cross-comparison:** `|L2_cholmod - L2_pcg| < 1e-6` (same as existing tests)
- **Thermo-elastic:** compare max pressure with analytical $\Delta P = -K \alpha \Delta T$ within 5% relative error

## Risks / Trade-offs

- **Oceanic cooling domain size vs runtime:** A tall, narrow domain (e.g. 200 km deep) with enough resolution to resolve the erf boundary layer needs at least 51 cells vertically. Using a 11×51 grid keeps runtime under 5 seconds. → Risk: low resolution in x may trigger edge effects. Mitigation: Neumann BCs on sides and a domain wider than the boundary layer.

- **ThermoElastic analytical validation is approximate:** The thermo-elastic pressure response $\Delta P = -K \alpha \Delta T$ is exact only for a uniform, unconstrained body. With a circular inclusion and BCs, stress concentrations make the comparison approximate. → Mitigation: check mean pressure change, not pointwise, and use a generous tolerance (5%).

- **gnuplot availability on CI:** The CI runner may not have gnuplot installed. → Mitigation: gnuplot scripts are for local visual inspection only; the GTest assertions don't depend on gnuplot. The C++ extractor writes CSVs that gnuplot reads, but the CI test only checks L2 norms.
