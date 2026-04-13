## 1. Add new timer variables to Main_DOODZ.c

- [x] 1.1 Declare `dt_interp`, `dt_stokes_setup`, `dt_nl_overhead`, `dt_post_solve` timer variables (initialised to 0.0) alongside existing timers
- [x] 1.2 Wrap the particle-to-grid interpolation block (~lines 519–660) with `omp_get_wtime()` start/stop, accumulate into `dt_interp`
- [x] 1.3 Wrap the Stokes pre-loop setup block (~lines 663–725) with `omp_get_wtime()` start/stop, accumulate into `dt_stokes_setup`
- [x] 1.4 Add `t_iter_start = omp_get_wtime()` at the top of each nonlinear iteration and compute `dt_nl_overhead += (iter_total - dt_rheology - dt_assembly - dt_solve_iter)` at the end
- [x] 1.5 Wrap the post-solve block (~lines 920–1290, excluding already-timed melting/GSE) with `omp_get_wtime()` start/stop, accumulate into `dt_post_solve`

## 2. Update perf.csv output format

- [x] 2.1 Update the perf.csv header fprintf to 25 columns: insert `interp_s,stokes_setup_s,nl_overhead_s,post_solve_s` after `output_s` and before `nit`
- [x] 2.2 Update the per-step fprintf data row to include the 4 new timer values in the correct column positions

## 3. Update LOG_TIME messages

- [x] 3.1 Add `LOG_TIME("Stokes setup: %.3f sec", dt_stokes_setup)` after the setup block
- [x] 3.2 Update the per-iteration LOG_TIME to include overhead: `"Iteration %02d: rheology=%.3f assembly=%.3f solve=%.3f overhead=%.3f sec"`
- [x] 3.3 Add `LOG_TIME("Post-solve updates: %.3f sec", dt_post_solve)` after the post-solve block
- [x] 3.4 Update the per-timestep breakdown LOG_TIME to include all new timers and coverage percentage

## 4. Local validation

- [x] 4.1 Build locally with `make build`
- [x] 4.2 Run RiftingComprehensive lowres for 3 steps, verify perf.csv has 25 columns
- [x] 4.3 Verify coverage: `sum(tracked) / wall_s >= 0.95` for each step
- [x] 4.4 Run RiftingComprehensive default for 3 steps, verify coverage >= 95%

## 5. Update benchmark scripts

- [x] 5.1 Update `misc/benchmark.sh` summary.csv header to include `avg_interp_s,avg_stokes_setup_s,avg_nl_overhead_s,avg_post_solve_s`
- [x] 5.2 Update `misc/benchmark.sh` awk extraction to compute averages for columns 15–18 from the 25-column perf.csv
- [x] 5.3 Update `misc/benchmark-report.sh` Results Summary table to include the 4 new columns
- [x] 5.4 Update `misc/benchmark-report.sh` to show a Coverage column instead of Other/overhead

## 6. Commit and push

- [x] 6.1 Commit all changes with message "feat: close perf.csv instrumentation gap — 25-column format with 4 new timers"
- [x] 6.2 Push to `add-performance-metrics` branch

## 7. Re-run AWS benchmark

- [x] 7.1 Start the c5ad.4xlarge instance (i-0df5754e400cb3ec1)
- [x] 7.2 Sync code to EC2 and rebuild
- [x] 7.3 Run the full 28-run sweep: `--resolutions "lowres default medres highres" --threads "1 2 4 6 8 12 16" --steps 10`
- [x] 7.4 Download results to local `benchmark-results/<timestamp>/`
- [x] 7.5 Verify coverage >= 95% across all 28 runs
- [x] 7.6 Stop the EC2 instance

## 8. Updated benchmark report

- [x] 8.1 Generate the updated benchmark report with `benchmark-report.sh`
- [x] 8.2 Review the report: confirm "Other/overhead" is replaced by explicit coverage percentage
- [x] 8.3 Save the report to `benchmark-results/20260413-105848/REPORT-20260413-113812.md`

## 9. Update skills

- [x] 9.1 Update `skill-logging` SKILL.md: perf.csv column table expanded from 21 to 25
- [x] 9.2 Update `skill-code-glossary` SKILL.md: add `dt_interp`, `dt_stokes_setup`, `dt_nl_overhead`, `dt_post_solve` to timing variables table
