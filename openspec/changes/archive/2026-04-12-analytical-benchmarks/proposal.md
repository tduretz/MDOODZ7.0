## Why

The existing CI test suite (15 suites, 39 cases) catches regressions by checking qualitative bounds ("viscosity > 0", "temperature increases"), but never compares results against known analytical solutions. For a submission, reviewers expect **quantitative verification** — measured L2 error norms and grid-convergence order — proving the solver is mathematically correct, not just crash-free. Several existing tests already solve problems with exact analytical solutions but don't use them. The codebase already contains analytical solution code (`eval_anal_Dani` in `SETS/AnisotropyDabrowski.c`, Julia scripts in `misc/JuliaScripts/`) that has never been wired into the automated test suite.

## What Changes

### A. New 2D Stokes benchmark — SolVi (Schmid & Podladchikov, 2003)
- Port the existing `eval_anal_Dani()` function from `SETS/AnisotropyDabrowski.c` into a new CI test. This gives the exact 2D analytical velocity ($V_x, V_z$), pressure ($P$), and deviatoric stress ($\sigma'_{xx}, \sigma'_{zz}$) around a circular viscous inclusion under pure shear or simple shear.
- Add a **grid-convergence test** that runs SolVi at multiple resolutions (e.g. 21→41→81) and verifies $O(h^2)$ error reduction, proving second-order spatial accuracy.

### B. New anisotropy benchmark — Schmalholz simple shear director evolution
- Port the analytical solution from `misc/JuliaScripts/SimpleShear_Anisotropy_fixed.jl` (comparison against Stefan Schmalholz's reference data). Under simple shear with constant anisotropy factor $\delta$, the director $(N_x, N_y)$ evolves analytically via the Mühlhaus director equation, and the anisotropic shear stress $\tau_{xy}$ follows the constitutive tensor $\mathbf{S} = \mathbf{C}_{ISO} + 2(\mu_N - \mu_S)\mathbf{C}_{ANI}(\theta)$.
- This tests the full anisotropy module: director rotation, anisotropic viscosity tensor assembly, and stress computation.

### C. New anisotropy benchmark — Viscous anisotropy stress vs orientation angle
- Port the analytical per-angle stress formula from `misc/JuliaScripts/Julia_anisotropy/ViscousAnisotropySimple.jl`. For each director angle $\theta$, the stress is computed by rotating strain rate to the principal anisotropy plane and applying $\mathbf{D} = 2\eta_N \text{diag}(1, 1, 1/\delta)$.
- This verifies stress invariant orientation dependence for the anisotropic constitutive law.

### D. Upgrade existing tests with L2 error assertions
- **SteadyStateGeotherm**: $T(z) = T_{top} + (T_{bot} - T_{top})(z_{top} - z)/H$
- **HydrostaticPressure**: $P(z) = \rho g |z|$
- **PureShearVelocity**: $V_x = \dot\varepsilon x$, $V_z = -\dot\varepsilon z$
- **StressBCPureShear**: $\sigma'_{xx} = 2\eta\dot\varepsilon$, $\sigma_{xz} = 0$
- **StressAccumulation**: Maxwell $\sigma(t) = 2\eta\dot\varepsilon(1 - e^{-Gt/\eta})$
- **ViscousDissipation**: $\bar{T}(t) = T_0 + 2\eta\dot\varepsilon^2 t / (\rho C_p)$
- **RadiogenicHeat**: $\bar{T}(t) = T_0 + Q_r t / (\rho C_p)$
- **ThermalExpansion**: $\rho(T) = \rho_0(1 - \alpha \Delta T)$

### E. Infrastructure & documentation
- Add `computeL2Error` and `computeFieldL2` helpers to `TestHelpers.h`.
- Create a **theory document** (`TESTS/AnalyticalSolutions.md`) referenced from `TESTS/README.md` that explains: the analytical solutions, how 1D analytical models are applied to 2D numerical grids, the SolVi derivation, the anisotropy director evolution, and full references.

## Existing Codebase Assets

These files already contain the analytical solutions we need — they only need to be ported into GTest:

| Asset | Location | What it provides |
|---|---|---|
| `eval_anal_Dani()` | `SETS/AnisotropyDabrowski.c` | Full 2D Vx, Vz, P, σxx, σzz for circular inclusion (Schmid 2002) |
| Anisotropy simple shear | `misc/JuliaScripts/SimpleShear_Anisotropy_fixed.jl` | Director evolution + τxy under simple shear vs Schmalholz reference |
| Anisotropy stress vs angle | `misc/JuliaScripts/Julia_anisotropy/ViscousAnisotropySimple.jl` | Stress invariant as function of director angle |
| Crameri topography | `misc/Python_visualisation/Read_Doodz_TopoBenchCase1.py` | Exponential topography decay: $h(t) = h_{max} e^{-\dot{h} t}$ |
| `ViscousInclusion.c` | `SETS/ViscousInclusion.c` | Same `eval_anal_Dani` for isotropic inclusion |
| `ViscousInclusionTestSuite.c` | `SETS/ViscousInclusionTestSuite.c` | Multi-resolution inclusion setup |

## Capabilities

### New Capabilities
- `benchmark-solvi`: SolVi 2D Stokes benchmark — port `eval_anal_Dani()` from SETS/, analytical velocity/pressure for circular inclusion, L2 error measurement, grid-convergence order verification
- `benchmark-anisotropy`: Anisotropy benchmarks — director evolution under simple shear (Schmalholz reference), stress invariant vs orientation angle, anisotropic constitutive law verification
- `l2-error-framework`: Shared L2 error computation helpers in TestHelpers.h, methodology for comparing 1D analytical solutions against 2D numerical fields (central column, full-grid), theory documentation

### Modified Capabilities
- `ci-thermal-tests`: Add L2 assertions to SteadyStateGeotherm and RadiogenicHeat against analytical solutions
- `ci-density-tests`: Add L2 assertion to HydrostaticPressure ($P = \rho g|z|$) and ThermalExpansion ($\rho = \rho_0(1-\alpha\Delta T)$)
- `ci-velocity-field-tests`: Add L2 assertion to PureShearVelocity ($V_x = \dot\varepsilon x$, $V_z = -\dot\varepsilon z$)
- `ci-neumann-bc-tests`: Add L2 assertion to StressBCPureShear ($\sigma'_{xx} = 2\eta\dot\varepsilon$)
- `ci-viscoelastic-tests`: Add L2 assertion to StressAccumulation (Maxwell analytical curve)
- `ci-shear-heating-tests`: Add L2 assertion to ViscousDissipation ($\Delta T = 2\eta\dot\varepsilon^2 t / \rho C_p$)
- `skill-testing-guide`: Document analytical benchmark patterns, L2 error thresholds, convergence-order testing

## Impact

- **`TESTS/TestHelpers.h`**: New functions `computeL2Error`, `computeFieldL2`, `readFieldAsVector`, and ported `eval_anal_Dani()`
- **`TESTS/SolViBenchmark/`**: New test suite directory with parameter files at multiple resolutions
- **`TESTS/SolViBenchmarkTests.cpp`**: New GTest file — SolVi L2 error on Vx/Vz/P, grid-convergence order check
- **`TESTS/AnisotropyBenchmark/`**: New test suite directory for anisotropy benchmarks
- **`TESTS/AnisotropyBenchmarkTests.cpp`**: New GTest file — director evolution, stress vs angle
- **6–8 existing `.cpp` test files**: Modified with additional `EXPECT_LT(L2, threshold)` assertions
- **`TESTS/CMakeLists.txt`**: New SolViBenchmarkTests and AnisotropyBenchmarkTests executables
- **`TESTS/AnalyticalSolutions.md`**: New theory document with derivations and references
- **`TESTS/README.md`**: Updated to reference the theory document and document new benchmark suites
- **CI runtime**: SolVi at 81×81 + anisotropy adds ~5–10 seconds; total suite stays under 30 seconds
