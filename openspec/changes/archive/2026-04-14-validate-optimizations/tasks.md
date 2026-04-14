## 1. PCG Residual Export (MDLIB)

- [x] 1.1 Add `export_pcg_residuals` parameter read in `InputOutput.c` (default 0), store in `model`
- [x] 1.2 Pass `model.export_pcg_residuals` and step number to `SolveThermalPCG` (add to call in `ThermalRoutines.c`)
- [x] 1.3 Implement CSV writing in `SolveThermalPCG`: when enabled, write `pcg_residuals_step<NNNNN>.csv` with header + per-iteration rows (iteration, abs_residual, rel_residual)
- [x] 1.4 Verify residual export has no effect on solver results (run existing `GaussianDiffusionPCG` test with and without flag)

## 2. Thermal Accuracy Comparison Tests

- [x] 2.1 Create `TESTS/Thermal/GaussianDiffusionL2.txt` (51×51, thermal_solver=0, 5 steps, writer_subfolder=GaussianDiffusionL2)
- [x] 2.2 Create `TESTS/Thermal/GaussianDiffusionL2PCG.txt` (same but thermal_solver=1, writer_subfolder=GaussianDiffusionL2PCG)
- [x] 2.3 Add `GaussianDiffusionL2` test case in `ThermalTests.cpp`: compute analytic T(r,t) at final step, compare via `computeL2Error`, assert L2 < 2e-2
- [x] 2.4 Add `GaussianDiffusionL2PCG` test case: same L2 check, plus assert difference between CHOLMOD and PCG L2 errors < 1e-6
- [x] 2.5 Enhance existing `SteadyStateGeothermPCG` test to print L2 alongside CHOLMOD L2 for side-by-side comparison
- [x] 2.6 No CMakeLists change needed (auto-registered via existing ThermalTests target); all 12 tests pass

## 3. Interp Equivalence Test

- [x] 3.1 Create `TESTS/InterpEquiv/InterpEquiv_Mode0.txt` (51×51, thermal=1, interp_mode=0, 5 steps, Gaussian diffusion setup)
- [x] 3.2 Create `TESTS/InterpEquiv/InterpEquiv_Mode3.txt` (identical except interp_mode=3 and writer_subfolder)
- [x] 3.3 Create `TESTS/InterpEquivalenceTests.cpp`: run both configs, read T/P/eta_n from final HDF5, assert max absolute diff < 1e-6
- [x] 3.4 Register `InterpEquivalenceTests` in `TESTS/CMakeLists.txt` and verify it passes (all diffs = 0)

## 4. Visual Tests (C++ + gnuplot)

- [x] 4.1 Create `VISUAL_TESTS/Validation.cpp` with `PlotThermalAccuracy()`, `PlotPCGConvergence()`, `PlotInterpEquivalence()`
- [x] 4.2 Create `VISUAL_TESTS/ThermalAccuracy.gnu`: grouped bar chart of L2 errors
- [x] 4.3 Create `VISUAL_TESTS/PCGConvergence.gnu`: log-scale convergence plot
- [x] 4.4 Create `VISUAL_TESTS/InterpEquivalence.gnu`: 2D heatmap with colour bar
- [x] 4.5 Create `VISUAL_TESTS/PCGConvergence.txt` with `export_pcg_residuals=1`
- [x] 4.6 Add Validation class to `main.cpp` with simulation runs and plot calls
- [x] 4.7 Add `Validation.cpp` to `VISUAL_TESTS/CMakeLists.txt` with configure_file entries
- [x] 4.8 Fix outdated SetBCT/MutateInput signatures in existing visual tests; verify `make build` with `-DVIS=ON`

## 5. Documentation

- [x] 5.1 Create `OPTIMIZATIONS.md` at repo root: overview of PCG thermal solver and fused P2Mastah interpolation
- [x] 5.2 Add performance summary table (Experiments 1–8: wall-time reduction, per-phase breakdown, thread scaling)
- [x] 5.3 Add validation section with embedded figure references and interpretation
- [x] 5.4 Add usage instructions: `thermal_solver = 1`, `max_its_thermal`, `rel_tol_thermal`, `interp_mode = 3`
- [x] 5.5 Add "when NOT to use" section: small grids < 100×100, single-thread runs, already-fast CHOLMOD with factorisation cache

## 6. Integration and Verification

- [x] 6.1 Build with `-DTEST=ON -DVIS=ON` and run full test suite: all 12 thermal + 1 interp equiv pass
- [x] 6.2 Visual tests compile with VIS=ON (Eigen3 + HDF5pp); PNG generation deferred to `make run-vis`
- [x] 6.3 Review `OPTIMIZATIONS.md` for completeness and accuracy
