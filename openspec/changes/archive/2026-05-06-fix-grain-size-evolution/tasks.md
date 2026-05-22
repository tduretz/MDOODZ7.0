## 1. Cleanup: remove dead `gsel` switch

- [x] 1.1 Delete `int gsel` declaration at [MDLIB/InputOutput.c:1077](MDLIB/InputOutput.c#L1077)
- [x] 1.2 Delete `gsel = ReadMatProps(...)` line at [MDLIB/InputOutput.c:1419](MDLIB/InputOutput.c#L1419)
- [x] 1.3 Delete `if (abs(gsel) > 0) LOG_INFO("Grain size evolution activated")` line at [MDLIB/InputOutput.c:1487](MDLIB/InputOutput.c#L1487)
- [x] 1.4 Build and run all GTest binaries to confirm no test references `gsel` directly; existing `.txt` files setting `gsel = 0` are silently ignored with no functional change (8/8 RheologyCreep tests pass)
- [x] 1.5 Manually run `cmake-exec/PinchSwellGSE/PinchSwellGSE` and confirm "Grain size evolution activated" is no longer logged (still crashes at init pre-D1, as expected)

## 2. Fix: NaN guard becomes fatal with diagnostic

- [x] 2.1 In `LocalIterationViscoElasticGrainSize` at [MDLIB/RheologyDensity.c:412-414](MDLIB/RheologyDensity.c#L412-L414), replace the `LOG_INFO("Local iterations from grain size and viscoisity exploded ...")` with `LOG_ERR("GSE local iter produced non-finite eta_ve: phase=%d T=%2.2e Eii=%2.2e last eta_ve=%2.2e", params.phase, params.T, params.Eii, eta_ve)` followed by `exit(346)` (also covers `isinf`)
- [x] 2.2 Add `T` (temperature) to the `LocalIterationParams` struct in [MDLIB/RheologyDensity.h](MDLIB/RheologyDensity.h) and populate `.T = T` in the params initialiser at [MDLIB/RheologyDensity.c:686](MDLIB/RheologyDensity.c#L686)
- [x] 2.3 Verified: PinchSwellGSE now emits `ERROR | GSE local iter produced non-finite eta_ve: phase=1 T=1.56e+00 Eii=5.00e-01 last eta_ve=nan` and exits with status 90 (= 346 mod 256, matching `exit(346)`). Natural divergence triggers the path, no test injection needed.
- [x] 2.4 N/A — no temporary injection was needed (2.3 verified via natural divergence)

## 3. Fix: bisection-then-Newton in `LocalIterationViscoElasticGrainSize`

- [x] 3.1 Replaced the dead bisection comment block + active Newton-only loop with a unified loop in `LocalIterationViscoElasticGrainSize`. Bisection runs for first `bisection_steps = 20` iterations, Newton for remaining `nitmax = 50 - 20 = 30`.
- [x] 3.2 Newton step gated on `it >= bisection_steps` (sequential, not interleaved)
- [x] 3.3 Replaced the original sign-comparison logic (which used a stale `r_eta_lo` from before the loop) with monotonicity-based update: r decreases with eta_ve, so r > 0 → narrow from below (`eta_up = eta_ve`), r < 0 → narrow from above (`eta_lo = eta_ve`). Documented inline.
- [x] 3.4 Added two early-exits: `res < tol/100` breaks the loop entirely (back-compat with prior behavior); `it == bisection_steps - 1 && res < tol` breaks at the bisection→Newton handoff.
- [x] 3.5 Used **log-space (geometric-mean) bisection** instead of linear: `eta_ve = sqrt(eta_up * eta_lo)`. Required because the bracket spans 20+ orders of magnitude — linear bisection would need 70+ steps; log-space halves the dex per step and converges in ~6.
- [x] 3.6 Replaced the silent `LOG_INFO("Visco-Elastic iterations failed!"); exit(0)` failure mode with a `LOG_ERR(...)` that includes phase, T, Eii, residual, and exits 346.
- [x] 3.7 Bumped `nitmax` from 20 to 50 (20 bisection + 30 Newton) — bisection convergence is linear-in-bracket-log so needs more iterations than Newton.
- [x] 3.8 `RheologyCreepTests` (8 tests) all pass — no regression.
- [x] 3.9 PinchSwellGSE runs to completion: exit=0, 11 output files at `Nt=100`, `d_n` at step 100 has shape (10000,), min=1.00e-05, max=2.78e-04, mean=3.69e-05, all finite, all positive.

## 4. New CI test: `RheologyCreep.GrainSizeSteadyState`

- [x] 4.1 Chose `pwlv = 15` (Calcite Renner et al. 2002) — pairs cleanly with `gs = 10` (Calcite paleowattmeter Austin & Evans 2002), no need for `cstv` fallback. PinchSwellGSE already uses this combo.
- [x] 4.2 Created [TESTS/RheologyCreep/GrainSizeSteadyState.txt](TESTS/RheologyCreep/GrainSizeSteadyState.txt): 11×11 grid, `Nb_phases = 1`, `user1 = 0` (no inclusion), `elastic = 0`, `thermal = 0`, `Nt = 1`, `pwlv = 15`, `gs = 10`, `gs_ref = 1e-3`, `bkg_strain_rate = 1e-14`, `user0 = 350` (°C). Yields `d_ss ≈ 124 µm`.
- [x] 4.3 Inlined constants in the test ($K_g = 2.5\times 10^{-9}$, $Q_g = 175$ kJ/mol, $\gamma = 1$, $\lambda = 0.1$, $c_g = \pi$, $p = 3$). For Renner power-law: $A = 1.5849\times 10^{-25}$, $n = 4.7$, $Q = 297$ kJ/mol, $F = (1/6) \cdot 2^{1/n} \cdot 3^{(n-1)/(2n)}$ (axial-compression correction `tpwl=1`).
- [x] 4.4 Added `TEST_F(RheologyCreep, GrainSizeSteadyState)` in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp).
- [x] 4.5 Used existing helpers only: `readFieldAsArray`, `computeL2Error`. Added `<algorithm>` and `<cmath>` to includes; no new test-helper utilities.
- [x] 4.6 Calibrated threshold: observed L2 = 7.05e-4 (analytical 1.2435e-4 m vs measured 1.2444e-4 m, 0.07% systematic offset due to MDOODZ prefactor handling — Newton iteration itself converges to ~1e-13 internally). Set L2 threshold to 5e-3 (5× measured per AnalyticalSolutions.md §"Threshold Calibration Methodology"). All 9 RheologyCreep tests now pass.
- [x] 4.7 CI verification — user confirmed "CI works"

## 5. Visualisation: gnuplot script + PNG

- [x] 5.1 Test writes `grain_size_benchmark.dat` (single row with `Eii T tau_II d_ss_ana d_mean_mdoodz`) as a side-effect.
- [x] 5.2 Created [TESTS/RheologyCreep/plot_grain_size.gp](TESTS/RheologyCreep/plot_grain_size.gp) with inline `d_ss(tau)` function and the `.dat` overlay; log-log axes spanning $10^4$–$10^9$ Pa (x) and $10^{-6}$–$10^{-2}$ m (y).
- [x] 5.3 Legend "Analytical d_{ss}(τ_II)..." for the line, "MDOODZ measurement" for the point; LaTeX-style axis labels.
- [x] 5.4 Verified PNG renders correctly (~37 KB, MDOODZ point lands on the analytical curve at ~23 MPa, ~125 µm).
- [x] 5.5 Copied `grain_size_benchmark.png` to [TESTS/](TESTS/) for commit alongside other benchmark PNGs.

## 6. Documentation

- [x] 6.1 Added "## 8. Grain Size Evolution: Steady-State Paleowattmeter" section at the end of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) (after §7 SolCx). Structure mirrors §3.4: Source, Parameter file, Test, Reference, Problem, Analytical, Implementation, L2 metric, Measured accuracy, Code assertions, Visualisation, PNG embed.
- [x] 6.2 Added two rows to the summary table — Grain size L2 (5e-3) and Grain size mean (1%).
- [x] 6.3 PNG embedded via `![Grain size benchmark...](grain_size_benchmark.png)` in §8.

## 7. Verification: PinchSwellGSE no longer crashes

- [x] 7.1 PinchSwellGSE runs to `Nt=100`, exit=0, 11 output files (verified during Group 3, task 3.9).
- [x] 7.2 Spot-checked: `Output00100.gzip.h5` `/Centers/d` shape (10000,), min=1.00e-05, max=2.78e-04, all finite, all positive.
- [x] 7.3 VISUAL_TESTS GSE driver — equivalent verification done via post-fix `gse.gif` regeneration. Compared against `gseref.gif`: 0.08-0.22% pixel difference, max channel diff ≤27/255 (gnuplot dithering noise). The previously committed `gse.gif` differed by up to 10% pixels with channel diff 255 — was generated when GSE was broken. Post-fix output reproduces the working-state reference.

## 8. CI integration

- [x] 8.1 Confirmed: `add_test(RheologyCreepTests RheologyCreepTests)` in [TESTS/CMakeLists.txt](TESTS/CMakeLists.txt) registers ctest discovery for the whole binary. Adding `TEST_F` to the existing `.cpp` requires no CMakeLists change. The `file(COPY ${PROJECT_SOURCE_DIR}/TESTS/RheologyCreep DESTINATION ...)` line auto-picks up `GrainSizeSteadyState.txt`.
- [x] 8.2 N/A — no new file created (test added to existing `RheologyCreepTests.cpp`).
- [x] 8.3 Push branch + observe green CI — user confirmed "CI works".
