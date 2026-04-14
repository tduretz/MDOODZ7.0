## Why

CHOLMOD (the sparse Cholesky solver from SuiteSparse) is used for both the Stokes and thermal solvers, but it is never configured for multi-threaded BLAS. Benchmark data (c5ad.4xlarge, 1000×800) shows `thermal_s` flat at ~6.5s regardless of thread count (1-16), confirming CHOLMOD runs single-threaded. Setting `c.nthreads_max = num_threads` after `cholmod_start()` enables CHOLMOD's internal multi-threaded sparse BLAS operations with zero algorithmic changes.

Additionally, the benchmark suite runs with `writer = 0`, so `output_s` is always 0.0. Enabling HDF5 output during EC2 benchmarks provides real I/O timing data.

## What Changes

- Set `CholmodSolver.c.nthreads_max` to the OpenMP thread count after `cholmod_start()` in Main_DOODZ.c (Stokes solver)
- Set `c.nthreads_max` after `cholmod_start()` in ThermalRoutines.c (thermal solver) and ChemicalRoutines.c (chemical solver)
- Add a `.txt` parameter `cholmod_threads` (default `1` = current single-threaded behavior, `-1` = use all available OpenMP threads, or explicit count) to control CHOLMOD threading — backward compatible
- Update `misc/benchmark-ec2-setup.sh` to enable `writer = 1` and `writer_step = 1` for the highres configuration
- Run highres-only benchmark sweep on EC2 and compare with baseline REPORT-20260413-113812.md

## Capabilities

### New Capabilities
- `cholmod-threading`: Specification for CHOLMOD multi-threaded BLAS configuration — parameter definition, initialization points, thread count propagation

### Modified Capabilities
- `phase-timing`: Benchmark runs should include HDF5 output enabled configuration to validate the `output_s` column

## Impact

- **MDLIB/Main_DOODZ.c**: Set `CholmodSolver.c.nthreads_max` after `cholmod_start()` (line ~715)
- **MDLIB/ThermalRoutines.c**: Set `c.nthreads_max` after `cholmod_start()` (line ~110)
- **MDLIB/ChemicalRoutines.c**: Set `c.nthreads_max` after `cholmod_start()` (line ~110)
- **MDLIB/InputOutput.c**: Add `ReadInt2(fin, "cholmod_threads", 1)` parameter
- **MDLIB/include/mdoodz.h**: Add `int cholmod_threads` field to model struct
- **SETS/*.txt**: New parameter (backward compatible, default 1 preserves current behavior)
- **misc/benchmark-ec2-setup.sh**: Enable writer for highres benchmark
- **Backward compatibility**: Fully backward compatible. Default `cholmod_threads = 1` preserves single-threaded CHOLMOD exactly.
