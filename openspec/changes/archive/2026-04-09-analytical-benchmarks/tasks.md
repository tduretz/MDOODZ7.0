## 1. L2 Error Framework (TestHelpers.h)

- [x] 1.1 Add `readFieldAsArray()` to `TESTS/TestHelpers.h` — reads a float field from HDF5 group/dataset into `std::vector<double>`
- [x] 1.2 Add `readCoordArray()` to `TESTS/TestHelpers.h` — reads a 1D coordinate array from HDF5 `Model/` group into `std::vector<double>`
- [x] 1.3 Add `computeL2Error()` to `TESTS/TestHelpers.h` — relative L2 norm with fallback to absolute L2 when analytical solution is zero
- [x] 1.4 Verify HDF5 field sizes match staggered grid expectations (Vx: Nx×(Nz-1), Vz: (Nx-1)×Nz, Centers: (Nx-1)×(Nz-1)) by adding a smoke test or assertion in one existing test

## 2. SolVi Benchmark — New Test Suite

- [x] 2.1 Create `TESTS/SolViBenchmark/` directory with `.txt` parameter files for resolutions 21, 41, 51, and 81 (single time step, no thermal, mc=1e3, inclusion radius via `user1`)
- [x] 2.2 Create `TESTS/SolViBenchmarkTests.cpp` with ported `eval_anal_Dani()` function and GTest fixture class
- [x] 2.3 Implement SetBCVx/SetBCVz callbacks applying analytical BCs (Dirichlet type=0 on E/W, type=11 on N/S)
- [x] 2.4 Implement L2 error test case at 51×51: read Vx, Vz, P, sxxd from HDF5, compute L2 against analytical, assert thresholds
- [x] 2.5 Implement grid-convergence test case: run at 21, 41, 81; compute L2(Vx) and L2(P) at each; compute convergence order; assert order ≥ 1.8 for Vx and ≥ 1.5 for P
- [x] 2.6 Add SolViBenchmarkTests to `TESTS/CMakeLists.txt` (file copy + executable + link + add_test)
- [x] 2.7 Build and run SolVi tests in WSL; record measured L2 values; calibrate thresholds to measured × 3

## 3. Anisotropy Benchmark — New Test Suite

- [x] 3.1 Check `MDLIB/HDF5Output.c` for director field output (nx, nz datasets) — determine if stress-only fallback is needed
- [x] 3.2 Create `TESTS/AnisotropyBenchmark/` directory with `.txt` parameter files for simple shear (director evolution) and pure shear (stress-vs-angle)
- [x] 3.3 Create `TESTS/AnisotropyBenchmarkTests.cpp` with GTest fixture, simple-shear and pure-shear callbacks
- [x] 3.4 Implement director evolution test: integrate Mühlhaus ODE analytically in the test, compare against numerical director field L2
- [x] 3.5 Implement stress-vs-angle test: compute analytical stress invariant at prescribed θ, compare against numerical τ_II
- [x] 3.6 Add AnisotropyBenchmarkTests to `TESTS/CMakeLists.txt`
- [x] 3.7 Build and run anisotropy tests in WSL; calibrate thresholds

## 4. Upgrade Existing Tests with L2 Assertions

- [x] 4.1 `ThermalTests.cpp` — SteadyStateGeotherm: add L2 of T field vs T(z) = T_top + (T_bot − T_top)(z_top − z)/H
- [x] 4.2 `ThermalTests.cpp` — RadiogenicHeat: add L2 of mean T vs T̄(t) = T₀ + Qr·t/(ρ·Cp)
- [x] 4.3 `DensityTests.cpp` — HydrostaticPressure: add L2 of P field vs P(z) = ρ·|g|·|z|
- [x] 4.4 `DensityTests.cpp` — ThermalExpansion: add L2 of ρ field vs ρ(T) = ρ₀(1 − α(T − T_ref))
- [x] 4.5 `VelocityFieldTests.cpp` — PureShearVelocity: add L2 of Vx vs ε̇·x and Vz vs −ε̇·z
- [x] 4.6 `NeumannBCTests.cpp` — StressBCPureShear: add L2 of sxxd vs 2η·ε̇ and absolute L2 of sxz ≈ 0
- [x] 4.7 `ViscoElasticTests.cpp` — StressAccumulation: add L2 of sxxd vs Maxwell σ(t) = 2η·ε̇·(1 − exp(−G·t/η))
- [x] 4.8 `ShearHeatingTests.cpp` — ViscousDissipation: add L2 of mean T vs T̄(t) = T₀ + 2η·ε̇2·t/(ρ·Cp)
- [x] 4.9 Build and run all 15+ suites in WSL; verify existing tests still pass; calibrate new L2 thresholds

## 5. Documentation

- [x] 5.1 Create `TESTS/AnalyticalSolutions.md` — SolVi section: problem statement, complex-variable formulas, Schmid & Podladchikov reference, measured L2 at each resolution
- [x] 5.2 `TESTS/AnalyticalSolutions.md` — 1D analytical solutions section: geotherm, hydrostatic P, pure shear V, Maxwell stress, viscous dissipation, radiogenic heat, thermal expansion
- [x] 5.3 `TESTS/AnalyticalSolutions.md` — Anisotropy section: Mühlhaus director ODE, anisotropic constitutive tensor, stress-vs-angle formula
- [x] 5.4 Update `TESTS/README.md` — add "Analytical Benchmarks" section linking to `AnalyticalSolutions.md` and summarising the new benchmark suites

## 6. Skill Update & Final Validation

- [x] 6.1 Update `skill-testing-guide` — add sections for analytical benchmark authoring, convergence-order testing pattern, L2 threshold calibration, and analytical benchmark pitfalls
- [x] 6.2 Full CI validation: build with `cmake -DTEST=ON`, run `ctest`, verify all suites pass and total runtime < 30 seconds
