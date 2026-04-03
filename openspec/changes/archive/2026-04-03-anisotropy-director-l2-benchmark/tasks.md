## 1. Parameter Files

- [x] 1.1 Create `TESTS/AnisotropyBenchmark/DirectorEvolution_dt01.txt` — copy from `DirectorEvolution.txt`, set dt=0.1, Nt=5
- [x] 1.2 Create `TESTS/AnisotropyBenchmark/DirectorEvolution_dt025.txt` — copy from `DirectorEvolution.txt`, set dt=0.025, Nt=20
- [x] 1.3 Update `TESTS/CMakeLists.txt` to copy the new parameter files to the build directory

## 2. Analytical Solution Helpers

- [x] 2.1 Add `analyticalDirectorAngle(theta0, gamma_dot, t)` helper in `AnisotropyBenchmarkTests.cpp` — returns arctan(tan(θ₀) − γ̇·t)
- [x] 2.2 Add `analyticalStressTII(theta, delta, eta, exx)` helper — implements the rotation+scaling formula from ViscosityConciseAniso and returns τ_II

## 3. Upgrade DirectorEvolution Test

- [x] 3.1 Replace qualitative `EXPECT_LT(meanAngle, 44.0)` with L2 error computation: read nx/nz from HDF5, compute θ = atan2(nz, nx) per cell, compute L2 vs analytical θ_final
- [x] 3.2 Assert L2(θ) < 0.01 radians
- [x] 3.3 If 11×11 grid is too noisy, increase to 21×21 or 31×31 and update parameter file

## 4. Add Director Dt-Convergence Test

- [x] 4.1 Create `TEST_F(AnisotropyBenchmark, DirectorDtConvergence)` that runs at dt=0.1, 0.05, 0.025 and computes L2(θ) at each
- [x] 4.2 Compute convergence order from finest pair (dt=0.05 → dt=0.025), assert order ≥ 0.8
- [x] 4.3 Assert monotonic error decrease: L2(dt=0.025) < L2(dt=0.05) < L2(dt=0.1)

## 5. Upgrade StressAnisotropy Test

- [x] 5.1 Replace qualitative `EXPECT_GT(meanTauII, 0.1)` with analytical τ_II comparison using the rotation formula
- [x] 5.2 Assert relative error < 2% for τ_II
- [x] 5.3 Also check sxxd relative error < 2%

## 6. Documentation

- [x] 6.1 Add anisotropy section to `TESTS/AnalyticalSolutions.md` — director ODE, closed-form, stress formula, measured L2 errors, convergence orders
- [x] 6.2 Update `.github/skills/skill-testing-guide/SKILL.md` — add anisotropy benchmarks to test coverage map with L2 thresholds
- [x] 6.3 Update `.github/skills/skill-anisotropy/SKILL.md` — add reference to quantitative benchmarks and measured accuracy

## 7. Verify

- [x] 7.1 Build and run all anisotropy tests — confirm all pass
- [x] 7.2 Record measured L2 errors and convergence orders in AnalyticalSolutions.md
