## 1. Parameter infrastructure

- [ ] 1.1 Add `int omp_schedule` field to the `model` struct in `MDLIB/include/mdoodz.h`
- [ ] 1.2 Add `ReadInt2(fin, "omp_schedule", 0)` in `MDLIB/InputOutput.c` alongside other model parameters
- [ ] 1.3 Add `LOG_INFO("OpenMP schedule mode: %d", model.omp_schedule)` at startup to confirm parameter is read
- [ ] 1.4 Build and verify parameter reads correctly (default 0, set to 1)

## 2. Guard ParticleRoutines.c (interpolation block — Priority 1)

- [ ] 2.1 Add `omp_in_parallel()` guards to all `#pragma omp parallel for` loops in `Interp_Grid2P` and its variants (~4 loops) — use `#pragma omp for` when already in parallel
- [ ] 2.2 Add `omp_in_parallel()` guard to P2Mastah's thread-count query region (line ~2046) — skip `#pragma omp parallel` when already in parallel, use `omp_get_num_threads()` directly
- [ ] 2.3 Add `omp_in_parallel()` guards to P2Mastah's main scatter loop (line ~2072) and reduction loops (lines ~2258, 2281) — convert `parallel for` to `for` when in parallel
- [ ] 2.4 Add `omp_in_parallel()` guards to remaining `parallel for` loops in ParticleRoutines.c (particle init, reseeding, etc.)

## 3. Guard RheologyParticles.c (NL iteration — Priority 2)

- [ ] 3.1 Add `omp_in_parallel()` guards to all 33 `#pragma omp parallel for` loops in `NonNewtonianViscosityGrid` and related functions — convert to `#pragma omp for` when in parallel

## 4. Guard StokesAssemblyDecoupled.c (NL iteration — Priority 2)

- [ ] 4.1 Add `omp_in_parallel()` guards to all 8 `#pragma omp parallel` regions that use manual `omp_get_thread_num()` decomposition — skip `#pragma omp parallel` when already in parallel, execute work distribution directly
- [ ] 4.2 Verify `estart/eend` arrays are computed using the same thread count as the persistent region (both use `omp_get_num_threads()`)

## 5. Guard ThermalSolver.c

- [ ] 5.1 Add `omp_in_parallel()` guards to all 6 `#pragma omp parallel for` loops in ThermalSolver.c

## 6. Guard post-solve functions

- [ ] 6.1 Identify all functions called in the post-solve block of Main_DOODZ.c that contain `#pragma omp parallel for` directives
- [ ] 6.2 Add `omp_in_parallel()` guards to those functions (UpdateNonLinearity, UpdateParticleEnergy, UpdateDensity, Advection helpers, etc.)

## 7. Persistent regions in Main_DOODZ.c

- [ ] 7.1 Add persistent `#pragma omp parallel` region around the interpolation block (P2Mastah calls), gated by `input.model.omp_schedule == 1`
- [ ] 7.2 Add persistent `#pragma omp parallel` region around the nonlinear iteration body (rheology + assembly calls), gated by `input.model.omp_schedule == 1`
- [ ] 7.3 Add persistent `#pragma omp parallel` region around the post-solve particle update sequence, gated by `input.model.omp_schedule == 1`
- [ ] 7.4 Ensure timer reads (`omp_get_wtime()`) and `LOG_TIME` calls are outside persistent regions or in `#pragma omp master` blocks
- [ ] 7.5 Verify the legacy code path (`omp_schedule == 0`) is completely unchanged

## 8. Local validation

- [ ] 8.1 Build with `make build`
- [ ] 8.2 Run RiftingComprehensive lowres (100×80, 10 steps) with `omp_schedule = 0` — verify identical to baseline
- [ ] 8.3 Run RiftingComprehensive lowres (100×80, 10 steps) with `omp_schedule = 1` — verify it runs without crashes
- [ ] 8.4 Compare HDF5 output between `omp_schedule = 0` and `omp_schedule = 1` runs — verify numerical equivalence (relative diff < 1e-10)
- [ ] 8.5 Verify perf.csv has 25 columns and coverage > 95% in both modes

## 9. Benchmark setup

- [ ] 9.1 Modify `misc/benchmark-ec2-setup.sh` to set `writer = 1` and `writer_step = 1` in the highres `.txt` configuration
- [ ] 9.2 Add a benchmark mode or parameter to run with `omp_schedule = 1` alongside the default `omp_schedule = 0` baseline

## 10. Commit and push

- [ ] 10.1 Commit all changes with message "feat: persistent OMP thread pool with omp_schedule parameter"
- [ ] 10.2 Push to branch

## 11. EC2 benchmark

- [ ] 11.1 Start c5ad.4xlarge instance
- [ ] 11.2 Sync code and rebuild on EC2
- [ ] 11.3 Run highres (1000×800, 10 steps) sweep with `omp_schedule = 0`: threads 1, 2, 4, 6, 8, 12, 16
- [ ] 11.4 Run highres (1000×800, 10 steps) sweep with `omp_schedule = 1`: threads 1, 2, 4, 6, 8, 12, 16
- [ ] 11.5 Download results to `benchmark-results/<timestamp>/`
- [ ] 11.6 Compare wall times between modes — quantify improvement at 8, 12, 16 threads
- [ ] 11.7 Verify `output_s` column is non-zero in writer-enabled runs
- [ ] 11.8 Stop EC2 instance

## 12. Report and skills

- [ ] 12.1 Generate benchmark report comparing `omp_schedule = 0` vs `omp_schedule = 1`
- [ ] 12.2 Update `skill-benchmarking` SKILL.md with `omp_schedule` parameter documentation
- [ ] 12.3 Update `skill-code-glossary` SKILL.md with `omp_schedule` field and `omp_in_parallel()` pattern
