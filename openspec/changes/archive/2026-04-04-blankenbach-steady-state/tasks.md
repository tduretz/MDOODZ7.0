# Tasks: blankenbach-steady-state

## 1. Parameter file

- [x] Create `TESTS/BlankenBench/BlankenBenchSteady.txt` — copy `BlankenBench.txt`, change `Nt=100000`, `Courant=0.5`, `writer_step=100000`, `istep=200000`

## 2. Test implementation

- [x] Add `TEST_F(BlankenBench, NusseltAndVrms)` to `TESTS/BlankenBenchTests.cpp`:
  - [x] OpenMP guard: `#ifdef _OMP_` → assert `omp_get_max_threads() == 8`, else `GTEST_SKIP()`
  - [x] Run simulation: `RunMDOODZ("BlankenBench/BlankenBenchSteady.txt", &setup)`
  - [x] Compute Nusselt number from top-row `dT/dz` finite difference
  - [x] Compute Vrms: interpolate Vx/Vz to cell centres, RMS, non-dimensionalise with `H/κ`
  - [x] `EXPECT_NEAR(Nu, 4.884409, 0.05 * 4.884409)`
  - [x] `EXPECT_NEAR(Vrms, 42.864947, 0.05 * 42.864947)`
  - [x] Print computed values for manual inspection

## 3. Documentation

- [x] Update `TESTS/README.md` Suite 18 section: add "Steady-State Benchmark" subsection with run command (`OMP_NUM_THREADS=8 ./BlankenBenchTests --gtest_filter=BlankenBench.NusseltAndVrms`), expected runtime, published reference values, and note that it is not part of CI

## 4. Build and verify

- [x] Build with `make -j$(nproc)` — confirm `BlankenBenchTests` compiles with the new test
- [x] Run smoke tests (`ctest`) — confirm existing tests still pass
- [ ] Run `OMP_NUM_THREADS=8 ./BlankenBenchTests --gtest_filter=BlankenBench.NusseltAndVrms` — confirm thread check passes and simulation starts
