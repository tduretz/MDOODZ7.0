## Why

The test suite has zero thermo-mechanical coupling tests. All 47 existing tests exercise either the thermal solver or the Stokes solver in isolation. Blankenbach et al. (1989) Case 1a is the community-standard benchmark for isoviscous Rayleigh-Bénard convection (Ra = 10⁴) with published Nusselt number and RMS velocity values. Adding it would validate the Boussinesq buoyancy coupling (ρ = ρ₀(1 − αΔT)), temperature advection by the velocity field, and long-integration stability over hundreds of time steps — all currently untested.

## What Changes

- Add a new scenario pair `SETS/BlankenBench.c` + `SETS/BlankenBench.txt` implementing the Blankenbach Case 1a setup: unit square, constant viscosity, T=1 bottom / T=0 top, free-slip walls, Ra = 10⁴
- Add a GTest file `TESTS/BlankenBenchTests.cpp` with two test cases:
  - **NusseltAndVrms**: run to steady state, compute Nusselt number (top heat flux) and RMS velocity, assert within ~1–2% of published values (Nu = 4.884, Vrms = 42.865)
  - **TemperatureProfile**: validate mid-depth maximum temperature against published value (T_max = 0.4266)
- Add a gnuplot script `TESTS/BlankenBench/plot_convection.gp` that produces a presentable temperature-field image with a heat colormap and velocity-vector overlay, readable from HDF5/ASCII output
- Add parameter file `TESTS/BlankenBench/BlankenBench.txt` for the CI test (moderate resolution for ~30s CI runtime)
- Register the new test executable in `TESTS/CMakeLists.txt`

## Capabilities

### New Capabilities
- `ci-blankenbach-convection`: GTest CI benchmark validating Nusselt number, RMS velocity, and mid-depth temperature for Blankenbach Case 1a (Ra = 10⁴) against published values
- `scenario-blankenbach`: New SETS/ scenario pair for isoviscous thermal convection in a square box (reusable beyond CI)
- `visual-blankenbach`: Gnuplot visualization script producing a temperature field plot with velocity vectors and a heat-palette colormap

### Modified Capabilities
_(none — no existing spec requirements change)_

## Impact

- **Code:** New `SETS/BlankenBench.c`, `TESTS/BlankenBenchTests.cpp`, `TESTS/BlankenBench/*.txt`, `TESTS/BlankenBench/plot_convection.gp`
- **Build:** `TESTS/CMakeLists.txt` — register new test executable
- **Docs:** `TESTS/AnalyticalSolutions.md` — new §5 documenting the Blankenbach analytical/published values
- **Skill:** `skill-testing-guide` and `skill-scenario-gallery` — add entries
- **CI time:** ~30–60s additional (the coupled simulation needs ~200+ steps to reach steady state at moderate resolution)
- **Dependencies:** None new — uses existing HDF5, SuiteSparse, GTest infrastructure
