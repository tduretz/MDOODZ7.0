## Why

MDOODZ currently uses dozens of independent `#pragma omp parallel for` regions per timestep — each one forks and joins a thread team. Benchmark data (c5ad.4xlarge, 1000×800, 16 threads) shows wall time *increasing* from 17.0s at 6 threads to 20.5s at 16 threads, with `interp_s` (P2Mastah) jumping from 3.1s to 6.5s. The fork-join overhead from ~45+ parallel region entries per step is the primary cause of this scaling regression. A persistent thread pool eliminates this overhead and is the prerequisite for all subsequent parallelization improvements.

Additionally, the current benchmark suite runs with HDF5 output disabled (`writer = 0`), so the `output_s` column in perf.csv is always 0.0. We need at least one benchmark configuration with output enabled to quantify I/O cost.

## What Changes

- Add a new `.txt` parameter `omp_schedule` (default `0` = current behavior, `1` = persistent thread pool) that controls whether the main timestep loop uses a single persistent `#pragma omp parallel` region with `#pragma omp for` / `#pragma omp single` inside, or the legacy per-function fork-join pattern
- Refactor `Main_DOODZ.c` main timestep body to support a persistent parallel region when `omp_schedule = 1`: wrap the timestep body in `#pragma omp parallel` and convert inner parallel-for loops to worksharing `#pragma omp for`
- Ensure all functions called within the persistent region are safe for execution inside an already-active parallel context (no nested `#pragma omp parallel` when already in a team)
- Update `misc/benchmark-ec2-setup.sh` to enable HDF5 output (`writer = 1`, `writer_step = 1`) for the highres benchmark configuration
- Run highres-only benchmark sweep on EC2 (1000×800, threads 1/2/4/6/8/12/16) and compare with baseline report REPORT-20260413-113812.md

## Capabilities

### New Capabilities
- `persistent-omp-region`: Specification for the persistent OpenMP parallel region mode — parameter definition, activation logic, thread safety requirements for called functions, and fallback to legacy behavior when disabled

### Modified Capabilities
- `phase-timing`: The perf.csv `output_s` column must be populated when `writer = 1`; benchmark runs should include at least one output-enabled configuration to validate this column

## Impact

- **MDLIB/Main_DOODZ.c**: Primary change — persistent parallel region wrapping the timestep body
- **MDLIB/ParticleRoutines.c**: P2Mastah and Interp_Grid2P must detect whether they are already inside a parallel region (via `omp_in_parallel()`) and skip their own `#pragma omp parallel` if so
- **MDLIB/HeatSolve.c, StokesAssemblyDecoupled.c, Advection.c, RheologyParticles.c**: Same pattern — functions must be callable from both persistent and legacy contexts
- **SETS/*.txt**: New `omp_schedule` parameter (backward compatible — defaults to 0 / legacy behavior)
- **misc/benchmark-ec2-setup.sh**: Add output-enabled highres config
- **Backward compatibility**: Fully backward compatible. Default `omp_schedule = 0` preserves current behavior exactly. No changes to simulation results — this is a performance-only optimization.
