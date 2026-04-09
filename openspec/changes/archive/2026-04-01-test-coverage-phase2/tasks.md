## 1. TestHelpers Extensions

- [ ] 1.1 Add `getFieldValueAt(filename, group, dataset, index)` to `TESTS/TestHelpers.h` for reading single values at specific grid indices
- [ ] 1.2 Add `getTopoFieldMax(filename)` and `getTopoFieldMin(filename)` helpers for reading `Topo/z_grid` (1D array, different from 2D Centers fields)

## 2. ViscoElastic Tests

- [ ] 2.1 Create `TESTS/ViscoElastic/StressAccumulation.txt` — pure shear, elastic=1, G=1e10, cstv=1, eta0=1e22, 5 steps, 31×31
- [ ] 2.2 Create `TESTS/ViscoElasticTests.cpp` — fixture with ShearTemplate-like callbacks, tests: StressAccumulation (sxxd step5 > sxxd step1), ElasticStrainRate (eII_el > 0)
- [ ] 2.3 Add ViscoElasticTests to `TESTS/CMakeLists.txt`

## 3. Density Tests

- [ ] 3.1 Create `TESTS/Density/ThermalExpansion.txt` — two phases with different T, α=3e-5, no deformation (bkg_strain_rate=0), gz=-10, 1 step
- [ ] 3.2 Create `TESTS/Density/HydrostaticPressure.txt` — uniform density, no-slip BCs, pure_shear_ALE=0, bkg_strain_rate=0, gz=-10, 1 step
- [ ] 3.3 Create `TESTS/DensityTests.cpp` — ThermalExpansion (min_rho < max_rho), HydrostaticPressure (maxP > 0 and within 50% of ρgH)
- [ ] 3.4 Add DensityTests to `TESTS/CMakeLists.txt`

## 4. ShearHeating Tests

- [ ] 4.1 Create `TESTS/ShearHeating/ViscousDissipation.txt` — cstv=1, eta0=1e22, thermal=1, shear_heating=1, bkg_strain_rate=1e-14, 3 steps, 31×31
- [ ] 4.2 Create `TESTS/ShearHeatingTests.cpp` — ViscousDissipation (maxT_final > maxT_init, ΔT within factor-of-2 of 2ηε̇²Δt/ρCp)
- [ ] 4.3 Add ShearHeatingTests to `TESTS/CMakeLists.txt`

## 5. FreeSurface Tests

- [ ] 5.1 Create `TESTS/FreeSurface/SinkingBlock.txt` — free_surface=1, free_surface_stab=1, dense inclusion (ρ=3300 vs ρ=3000), gz=-10, 41×41, 3 steps
- [ ] 5.2 Create `TESTS/FreeSurfaceTests.cpp` — SinkingBlock (Topo/z_grid changes from flat, values are finite)
- [ ] 5.3 Add FreeSurfaceTests to `TESTS/CMakeLists.txt`

## 6. VelocityField Tests

- [ ] 6.1 Create `TESTS/VelocityField/PureShearVelocity.txt` — same as FreeSlipPureShear but re-used for Vx/Vz validation
- [ ] 6.2 Create `TESTS/VelocityFieldTests.cpp` — PureShearVelocity (max Vx ≈ ε̇·0.5 within 20%), NoSlipVelocity (read from existing NoSlip output, verify bounded V)
- [ ] 6.3 Add VelocityFieldTests to `TESTS/CMakeLists.txt`

## 7. Compressibility Tests

- [ ] 7.1 Create `TESTS/Compressibility/DilatantFlow.txt` — compressible=1, plast=1, psi=10°, high strain rate, 1 step
- [ ] 7.2 Create `TESTS/Compressibility/IncompressibleRef.txt` — compressible=0, same otherwise
- [ ] 7.3 Create `TESTS/CompressibilityTests.cpp` — DilatantFlow (max|divu| > 0), IncompressibleRef (max|divu| < 1e-6·ε̇)
- [ ] 7.4 Add CompressibilityTests to `TESTS/CMakeLists.txt`

## 8. FiniteStrain Tests

- [ ] 8.1 Create `TESTS/FiniteStrain/PureShearStrain.txt` — finite_strain=1, pure shear, 5 steps, 31×31
- [ ] 8.2 Create `TESTS/FiniteStrainTests.cpp` — PureShearStrain (Fxx > 1.0, Fzz < 1.0, all F finite)
- [ ] 8.3 Add FiniteStrainTests to `TESTS/CMakeLists.txt`

## 9. NeumannBC Tests

- [ ] 9.1 Create `TESTS/NeumannBC/StressBCPureShear.txt` — constant viscosity, stress BCs (type 13) on east/west, 31×31, 1 step
- [ ] 9.2 Create `TESTS/NeumannBCTests.cpp` with SetBCVx/SetBCVz callbacks setting type=13 with known traction — StressBCPureShear (sxxd non-zero, solver converges)
- [ ] 9.3 Add NeumannBCTests to `TESTS/CMakeLists.txt`

## 10. ConvergenceRate Tests

- [ ] 10.1 Create `TESTS/ConvergenceRate/NewtonPwl.txt` — Newton=1, line_search=1, pwlv=40, nit_max=50
- [ ] 10.2 Create `TESTS/ConvergenceRate/PicardPwl.txt` — Newton=0, same problem otherwise
- [ ] 10.3 Create `TESTS/ConvergenceRateTests.cpp` — NewtonVsPicard (Newton iterations ≤ Picard iterations), NewtonFlag (LogIsNewtonStep has non-zero entries)
- [ ] 10.4 Add ConvergenceRateTests to `TESTS/CMakeLists.txt`

## 11. Strengthen Existing Tests — RheologyCreep

- [ ] 11.1 Add PowerLawCreep test case to `TESTS/RheologyCreepTests.cpp` — reads eII_pwl, asserts > 0
- [ ] 11.2 Create `TESTS/RheologyCreep/PowerLawCreep.txt` — pwlv=40, linv=0, cstv=0
- [ ] 11.3 Add eII_lin > 0 assertion to existing DiffusionCreepLinear test
- [ ] 11.4 Add StrainRatePartitioning test case to `TESTS/RheologyCreepTests.cpp` — reads eII_pwl and eII_lin from CompositePwlLin, asserts both > 0

## 12. Strengthen Existing Tests — BoundaryCondition

- [ ] 12.1 Add Vx field check to FreeSlipPureShear test in `TESTS/BoundaryConditionTests.cpp` — maxVx ≈ ε̇·0.5
- [ ] 12.2 Add Vx field check to NoSlip test — verify velocity is bounded

## 13. Strengthen Existing Tests — Thermal

- [ ] 13.1 Add meanT ≈ (T_top + T_bot)/2 check to SteadyStateGeotherm test in `TESTS/ThermalTests.cpp`
- [ ] 13.2 Add quantitative ΔT check to RadiogenicHeat test — within factor of 5 of Qr·Nt·Δt/(ρCp)

## 14. Build & Validate

- [ ] 14.1 Copy updated repo to WSL native filesystem, build with `cmake -DTEST=ON`, verify compilation succeeds
- [ ] 14.2 Run `ctest --output-on-failure` in WSL — all tests must pass (0 failures)
- [ ] 14.3 Fix any failing tests (parameter adjustments, assertion relaxation) and re-run until 100% pass

## 15. Documentation & Skills

- [ ] 15.1 Add all new test suites to `TESTS/README.md` with governing equations, cause→effect, parameter tables, and code assertions
- [ ] 15.2 Update `skill-testing-guide` with new test patterns and pitfalls
- [ ] 15.3 Update `skill-physics-theory` with shear heating, hydrostatic pressure equations
- [ ] 15.4 Update `skill-rheology` with power-law test pattern, strain partitioning
- [ ] 15.5 Update `skill-model-setup` with free surface requirements, Neumann BC setup
- [ ] 15.6 Update `skill-parameter-validation` with new switch compatibility rows
