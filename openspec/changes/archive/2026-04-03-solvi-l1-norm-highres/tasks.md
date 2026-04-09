## 1. L1 Error Norm Helper

- [x] 1.1 Add `computeL1Error()` static inline function in `TESTS/TestHelpers.h` right after `computeL2Error()`. Signature: `static inline double computeL1Error(const double* num, const double* ana, int N)`. Computes `(1.0/N) * sum(fabs(num[i] - ana[i]))`.

## 2. High-Resolution Parameter Files

- [x] 2.1 Create `TESTS/SolViBenchmark/SolViRes101.txt` — copy of `SolViRes41.txt` with `Nx = 101`, `Nz = 101`, `writer_subfolder = SolViRes101`
- [x] 2.2 Create `TESTS/SolViBenchmark/SolViRes151.txt` — copy of `SolViRes41.txt` with `Nx = 151`, `Nz = 151`, `writer_subfolder = SolViRes151`
- [x] 2.3 Create `TESTS/SolViBenchmark/SolViRes201.txt` — copy of `SolViRes41.txt` with `Nx = 201`, `Nz = 201`, `writer_subfolder = SolViRes201`

## 3. HighResL1Convergence Test

- [x] 3.1 Add `TEST(SolViBenchmarkTests, HighResL1Convergence)` in `TESTS/SolViBenchmarkTests.cpp` that runs SolVi at resolutions 41, 81, 101, 151, 201
- [x] 3.2 Compute both `computeL2Error()` and `computeL1Error()` for Vx and P at each resolution
- [x] 3.3 Print a formatted table: resolution, h, L2(Vx), L1(Vx), L2(P), L1(P) for each resolution
- [x] 3.4 Compute and print convergence orders for all consecutive pairs (41→81, 81→101, 101→151, 151→201) and the 101→201 pair
- [x] 3.5 Add assertions: `EXPECT_GE(order_P_L1_101_201, 0.9)` and `EXPECT_GE(order_Vx_L1_101_201, 0.8)`

## 4. Existing GridConvergence L1 Reporting

- [x] 4.1 In the existing `GridConvergence` test in `SolViBenchmarkTests.cpp`, compute and print L1 norms alongside L2 for each resolution (41, 81). No new assertions.

## 5. Build and Run

- [x] 5.1 Sync files to WSL, build, verify all existing tests still pass (no regression)
- [x] 5.2 Run `HighResL1Convergence` test, record actual L1 and L2 values and orders
- [x] 5.3 If P L1 order at 101→201 is < 0.9, adjust assertion threshold based on measured value and document finding

## 6. Documentation

- [x] 6.1 Update `TESTS/AnalyticalSolutions.md` with high-resolution L1 convergence table and findings
- [x] 6.2 Update `skill-testing-guide/SKILL.md` with L1 norm findings and HighResL1Convergence test listing
- [x] 6.3 Update `skill-solvers/SKILL.md` convergence section with L1 vs L2 insight and measured orders
