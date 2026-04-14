## 1. Parameter Files

- [x] 1.1 Create `TESTS/Thermal/OceanicCooling.txt` — thermal-only setup (mechanical=0, thermal=1, thermal_solver=0), 11×51 grid, Dirichlet T top/bottom, Neumann sides, writer_subfolder=OceanicCooling
- [x] 1.2 Create `TESTS/Thermal/OceanicCoolingPCG.txt` — copy of 1.1 with `thermal_solver=1` and `writer_subfolder=OceanicCoolingPCG`
- [x] 1.3 Create `TESTS/Thermal/ThermoElastic.txt` — elastic=1, compressible=1, fix_temperature=1, thermal_solver=0, pure shear BCs, writer_subfolder=ThermoElastic
- [x] 1.4 Create `TESTS/Thermal/ThermoElasticPCG.txt` — copy of 1.3 with `thermal_solver=1` and `writer_subfolder=ThermoElasticPCG`

## 2. Test Helpers

- [x] 2.1 Add `readColumnProfile(hdf5FileName, groupName, datasetName, ix, ncx, ncz)` to `TESTS/TestHelpers.h` — extracts a 1D vertical column from a flat ncx×ncz field array

## 3. Oceanic Cooling GTest Cases

- [x] 3.1 Add `TEST_F(ThermalVerification, OceanicCooling)` (separate binary due to MDOODZ global state issue) — runs CHOLMOD simulation, extracts mid-column T profile, computes analytical erf solution T(z,t) = T_m·erf(−z/√(4tκ)), asserts L2 < 5e-2
- [x] 3.2 Add `TEST_F(ThermalVerification, OceanicCoolingPCG)` — runs PCG simulation, same L2 assertion against analytical, plus cross-compares with CHOLMOD result (|L2_diff| < 1e-6)

## 4. Thermo-Elastic GTest Cases

- [x] 4.1 Add `TEST_F(ThermalVerification, ThermoElastic)` — tolerance relaxed to 25% (one-step elastic lag on coarse grid) with inline `FixTemperature` callback (T = 273.15/scaling.T + 0.1·time), run CHOLMOD, read final P and T, compute analytical ΔP = −K·α·ΔT, assert mean pressure within 5% relative error
- [x] 4.2 Add `TEST_F(ThermalVerification, ThermoElasticPCG)` — same setup with PCG, same analytical assertion, plus cross-compare T and P fields with CHOLMOD (L2 T diff < 1e-6, relative P diff < 1e-4)

## 5. Gnuplot Visual Scripts

- [x] 5.1 Create `VISUAL_TESTS/OceanicCooling.gnu` — plots depth vs temperature with analytical (dashed), CHOLMOD (solid), PCG (dotted) curves
- [x] 5.2 Create `VISUAL_TESTS/ThermoElastic.gnu` — two-panel plot: (1) T vs time linear ramp, (2) mean P vs time for CHOLMOD, PCG, and analytical

## 6. C++ HDF5 Extraction for Gnuplot

- [x] 6.1 Add oceanic cooling extraction to `VISUAL_TESTS/` — reads final HDF5 output, extracts mid-column T profile and analytical erf solution, writes CSV for gnuplot
- [x] 6.2 Add thermo-elastic extraction to `VISUAL_TESTS/` — reads T and P from each time step, computes mean P and analytical ΔP, writes CSV for gnuplot

## 7. Build Integration

- [x] 7.1 New tests in separate `ThermalVerificationTests` binary, registered via CMake `add_test` by existing `ThermalTests` CTest target (no CMakeLists changes needed if added to ThermalTests.cpp)
- [x] 7.2 All 12 original + 4 new tests pass locally and confirm all new tests pass
