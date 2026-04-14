## 1. Setup

- [x] 1.1 Copy `TESTS/BlankenBench/BlankenBench.txt` to `SETS/BlankenBench.txt`
- [x] 1.2 Build BlankenBench scenario: `make build SET=BlankenBench`
- [x] 1.3 Verify BlankenBench runs locally: quick test with resolution mode (grid mode fails — BlankenBench physics tuned for 41×41)

## 2. Baseline Benchmark Sweep

- [x] 2.1 Run baseline sweep: `./misc/benchmark.sh --scenario BlankenBench --resolutions "default" --steps 1000 --writer --threads "1 2 4 6 8" --skip-build` → `benchmark-results/20260414-134249/`
- [x] 2.2 Verify results directory created with summary.csv and per-run perf.csv files
- [x] 2.3 Spot-check: confirm `avg_thermal_s` shows CHOLMOD times (should be >0) and `avg_output_s` shows HDF5 write times

## 3. Optimized Benchmark Sweep

- [x] 3.1 Run optimized sweep: `./misc/benchmark.sh --scenario BlankenBench --resolutions "default" --steps 1000 --writer --threads "1 2 4 6 8" --skip-build --sed 's/^thermal_solver.*/thermal_solver = 1/' --sed 's/^interp_mode.*/interp_mode = 3/'` → `benchmark-results/20260414-140451/`
- [x] 3.2 Verify results directory created with summary.csv
- [x] 3.3 Spot-check: PCG times identical to CHOLMOD at 41×41 — matrix too small for PCG to help

## 4. Analysis and Documentation

- [x] 4.1 Compare baseline vs optimized summary.csv: all differences within ±2% (noise)
- [x] 4.2 M1 observations: no scaling past 4t (E-cores add overhead), interp=78-81% of wall, HDF5=4-5%
- [x] 4.3 Documented as Experiment 8 in `.github/skills/skill-benchmarking/SKILL.md`

## 5. Cleanup

- [x] 5.1 Remove `SETS/BlankenBench.txt` (temporary copy)
- [x] 5.2 HDF5 outputs left in benchmark-results/ directories (self-contained)
