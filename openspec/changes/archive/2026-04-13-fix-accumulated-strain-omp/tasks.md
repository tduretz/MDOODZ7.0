## 1. Fix the OpenMP pragma

- [x] 1.1 Uncomment the `#pragma omp parallel for` on the particle loop in `AccumulatedStrainII` (line ~240 of `MDLIB/RheologyParticles.c`)
- [x] 1.2 Fix the `private` clause: add `dE_pl_vol`, `k`, `l` and `dE_el` (missing from the original commented pragma)
- [x] 1.3 Add `strain_inc_pl_vol` to the `shared` clause

## 2. Replace exit(0) with thread-safe error handling

- [x] 2.1 Replace `exit(0)` in the grid loop's negative-strain check (line ~222) with `LOG_ERR` + `exit(1)` (this is outside the particle parallel region, but `exit(0)` is wrong for an error)
- [x] 2.2 Verify no `exit()` or `abort()` calls exist inside the particle loop body

## 3. Validation

- [x] 3.1 Build locally and run CI test suite — all 19 tests pass
- [x] 3.2 Run a quick local benchmark (RiftingComprehensive default, 1 and 4 threads) to verify `post_solve` improvement

## 4. EC2 Benchmarking

- [x] 4.1 Run EC2 benchmark: RiftingComprehensive 1001×801, 10 steps, threads 1/2/4/6/8/12/16 with `thermal_solver = 1`, `interp_mode = 2` (same config as Experiment 5 mode 2)
- [x] 4.2 Download results, compare `post_solve_s` and `wall_s` against Experiment 5 mode 2 baseline
- [x] 4.3 Document as Experiment 6 in skill-benchmarking SKILL.md
