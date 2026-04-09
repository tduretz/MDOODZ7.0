## 1. CI Test Infrastructure

- [x] 1.1 Create `TESTS/RheologyCreepTests.cpp` with GTest fixture and setup callbacks (ShearTemplate-like: circular inclusion, SetPhase, SetDensity, SetBCVx/Vz)
- [x] 1.2 Create parameter files `TESTS/RheologyCreep/DiffusionCreepLinear.txt` (linv=40, 31×31, 1 step) and `TESTS/RheologyCreep/DiffusionCreepGrainSize.txt` (linv=40, gs_ref=1e-3)
- [x] 1.3 Create parameter files `TESTS/RheologyCreep/PeierlsCreep.txt` (expv=40, 31×31, 1 step) and `TESTS/RheologyCreep/GBSCreep.txt` (gbsv=40)
- [x] 1.4 Create parameter files `TESTS/RheologyCreep/CompositePwlLin.txt` (pwlv=40, linv=40)
- [x] 1.5 Implement test cases: `DiffusionCreepLinear` (assert 1 iter), `DiffusionCreepGrainSize` (assert eta bounds), `PeierlsCreep` (assert 1–20 iters), `GBSCreep` (assert <20 iters), `CompositePwlLin` (assert 1–15 iters + viscosity bound check)
- [x] 1.6 Add `RheologyCreepTests` executable to `TESTS/CMakeLists.txt` and register with CTest

## 2. CI Plasticity Tests

- [x] 2.1 Create `TESTS/PlasticityTests.cpp` with GTest fixture and setup callbacks
- [x] 2.2 Create parameter files `TESTS/Plasticity/DruckerPragerYield.txt` (plast=1, C=20e6, phi=30°, high strain rate) and `TESTS/Plasticity/DruckerPragerNoYield.txt` (same but strain rate 1e-18)
- [x] 2.3 Create parameter file `TESTS/Plasticity/StrainSoftening.txt` (coh_soft=1, C_end=5e6, 3 steps)
- [x] 2.4 Create parameter file `TESTS/Plasticity/StressLimiter.txt` (Slim=100e6)
- [x] 2.5 Implement test cases: `DruckerPragerYield` (assert eII_pl > 0), `DruckerPragerNoYield` (assert eII_pl ≈ 0), `StrainSoftening` (assert cohesion < C_init), `StressLimiter` (assert stress bounded)
- [x] 2.6 Add HDF5 field-reading helpers: `getMaxFieldValue(filename, datasetPath, ncx, ncz)` and `getMinFieldValue()` in `TestHelpers.h`
- [x] 2.7 Add `PlasticityTests` to `TESTS/CMakeLists.txt`

## 3. CI Thermal Tests

- [x] 3.1 Create `TESTS/ThermalTests.cpp` with GTest fixture and thermal-only setup (mechanical=1, thermal=1)
- [x] 3.2 Create parameter files `TESTS/Thermal/GaussianDiffusion.txt` (thermal perturbation, 5 steps) and `TESTS/Thermal/SteadyStateGeotherm.txt` (fixed T BCs, 20 steps)
- [x] 3.3 Create parameter file `TESTS/Thermal/RadiogenicHeat.txt` (Qr=1e-6, 5 steps)
- [x] 3.4 Create parameter file `TESTS/Thermal/DirichletBC.txt` (T_top=273.15, T_bottom=1600, 1 step)
- [x] 3.5 Implement test cases: `GaussianDiffusion` (assert peak T decreases), `SteadyStateGeotherm` (assert max T > min T), `RadiogenicHeat` (assert mean T increases), `DirichletBC` (assert T bounded)
- [x] 3.6 Add `ThermalTests` to `TESTS/CMakeLists.txt`

## 4. CI Boundary Condition Tests

- [x] 4.1 Create `TESTS/BoundaryConditionTests.cpp` with GTest fixture
- [x] 4.2 Create parameter files `TESTS/BoundaryCondition/PeriodicSimpleShear.txt` (periodic_x=1, shear_style=1) and `TESTS/BoundaryCondition/FreeSlipPureShear.txt`
- [x] 4.3 Create parameter file `TESTS/BoundaryCondition/NoSlip.txt`
- [x] 4.4 Implement test cases: `PeriodicSimpleShear` (assert P finite), `FreeSlipPureShear` (assert sxz bounded), `NoSlip` (assert convergence)
- [x] 4.5 Add `BoundaryConditionTests` to `TESTS/CMakeLists.txt`

## 5. CI Solver Mode Tests

- [x] 5.1 Create `TESTS/SolverModeTests.cpp` with GTest fixture
- [x] 5.2 Create parameter files `TESTS/SolverMode/PicardConvergence.txt` (Newton=0, nit=50, pwlv=1) and `TESTS/SolverMode/AdaptiveDt.txt` (constant_dt=0, Courant=0.5, 3 steps)
- [x] 5.3 Implement test cases: `PicardConvergence` (assert final residual + iteration count), `AdaptiveDt` (assert both steps complete)
- [x] 5.4 Add `SolverModeTests` to `TESTS/CMakeLists.txt`

## 6. Visual Test — Advanced Rheology

- [x] 6.1 Create `VISUAL_TESTS/PeirlsCreep.cpp` plot functions (with expv=40, ~5 steps, 51×51)
- [x] 6.2 Create `VISUAL_TESTS/PeirlsCreep.txt` parameter file and `PeirlsCreep.gnu` gnuplot script
- [x] 6.3 Generate reference HDF5 `PeirlsCreepReference.h5` from a known-good run
- [x] 6.4 Create `VISUAL_TESTS/AnisoVEP.cpp` plot functions (elastic=1, anisotropy=1, plasticity, ~10 steps)
- [x] 6.5 Create `VISUAL_TESTS/AnisoVEP.txt` parameter file and `AnisoVEP.gnu` gnuplot script
- [x] 6.6 Generate reference HDF5 `AnisoVEPReference.h5`
- [x] 6.7 Add both scenarios to `RunTestCases()` and plotting functions in `main.cpp`
- [x] 6.8 Add new files to `VISUAL_TESTS/CMakeLists.txt` (source files + copy txt/h5 to build dir)

## 7. Visual Test — Thermal and Phase Transition

- [x] 7.1 Create `VISUAL_TESTS/ThermalDiffusion.cpp` plot functions (mechanical=1, thermal=1, Gaussian perturbation, ~10 steps)
- [x] 7.2 Create `VISUAL_TESTS/ThermalDiffusion.txt` and `ThermalDiffusion.gnu`
- [x] 7.3 Generate reference HDF5 `ThermalDiffusionReference.h5`
- [x] 7.4 Create `VISUAL_TESTS/QuartzCoesite.cpp` plot functions (based on SETS/QuartzCoesite.c setup)
- [x] 7.5 Create `VISUAL_TESTS/QuartzCoesite.txt` and `QuartzCoesite.gnu`
- [x] 7.6 Generate reference HDF5 `QuartzCoesiteReference.h5`
- [x] 7.7 Add both scenarios to `main.cpp` and `VISUAL_TESTS/CMakeLists.txt`

## 8. Visual Test Numeric Assertions

- [x] 8.1 Add `relativeL2Error()` helper function to `visual-tests.h` using Eigen
- [x] 8.2 Add tolerance-based assertion macro/function (`VISUAL_ASSERT`, `VISUAL_ASSERT_L2`) and `g_visualTestFailed` flag
- [x] 8.3 Retrofit existing visual tests (ShearTemplate, ShearTemplate1, ShearTemplateAniso, ShearHeatingDuretz14, RiftingChenin) with `CompareH5Field()` L2 checks
- [x] 8.4 New visual tests use `CompareH5Field()` once reference HDF5 files are generated
- [x] 8.5 Added `CompareH5Field()` for field-level scalar/L2 assertions
- [x] 8.6 `main()` returns nonzero exit code if `g_visualTestFailed` is true

## 9. CI Integration

- [x] 9.1 Add `visual-tests` job to `.github/workflows/ci.yml`: install Eigen3, build with `VIS=ON`, run `./visualtests`, check exit code
- [x] 9.2 Gate visual test job with `if: github.event_name == 'push' || contains(github.event.pull_request.labels.*.name, 'run-visual-tests')`
- [x] 9.3 Add `generate-refs` make target to `makefile` that runs visual tests and copies output HDF5 to reference locations

## 10. Testing Skill

- [x] 10.1 Create `.github/skills/skill-testing-guide/SKILL.md` with YAML frontmatter
- [x] 10.2 Write CI test authoring section (fixture pattern, .txt conventions, assertion patterns, CMakeLists integration)
- [x] 10.3 Write visual test authoring section (C++ class pattern, Eigen extraction, L2 comparison, gnuplot, reference generation)
- [x] 10.4 Write test execution section (build commands, run commands, CI behaviour)
- [x] 10.5 Write test coverage map table (feature → test file → CI/visual)
- [x] 10.6 Add skill-testing-guide to `.github/copilot-instructions.md` routing table

## 11. Verification

- [x] 11.1 Build with `TEST=ON` and run `ctest` — all new CI tests pass
- [ ] 11.2 Build with `VIS=ON` and run `./visualtests` — exit code is 0
- [ ] 11.3 Verify CI workflow runs successfully on a test push
