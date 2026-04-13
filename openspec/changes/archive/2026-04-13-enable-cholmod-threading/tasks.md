## Tasks

### Group 1: Add parameter and struct field

- [x] T01: Add `int cholmod_threads` field to the `model` struct in `mdoodz.h`
- [x] T02: Add `ReadInt2(fin, "cholmod_threads", 1)` call in `InputOutput.c` to read the parameter from `.txt` files
- [x] T03: Add `LOG_INFO("CHOLMOD threads: %d", resolved)` after resolving the thread count

### Group 2: Configure CHOLMOD threading at all init sites

- [x] T04: After `cholmod_start(&CholmodSolver.c)` in `Main_DOODZ.c`, set `CholmodSolver.c.nthreads_max = cholmod_nthreads`
- [x] T05: After `cholmod_start(&c)` in `ThermalRoutines.c`, set `c.nthreads_max` using `omp_get_max_threads()`
- [x] T06: After `cholmod_start(&c)` in `ChemicalRoutines.c`, set `c.nthreads_max` using `omp_get_max_threads()`

### Group 3: Update scenario .txt files

- [x] T07: Add `cholmod_threads = 1` to `SETS/RiftingChenin.txt` (benchmark scenario) as documentation — default is already 1 so existing .txt files need no change for backward compatibility

### Group 4: Local validation

- [x] T08: Build with `make build` and verify compilation succeeds
- [x] T09: Run lowres RiftingChenin with default `cholmod_threads` (absent from .txt) — verify identical output to current baseline
- [x] T10: Run lowres RiftingChenin with `cholmod_threads = -1` — verify LOG_INFO shows all threads and simulation completes

### Group 5: Enable HDF5 writer in benchmark config

- [x] T11: Update `benchmark.sh` to support `--writer` flag and `--sed` for extra .txt patching

### Group 6: EC2 benchmark

- [x] T12: Start EC2 instance, pull changes, build
- [x] T13: Run highres benchmark sweep with `cholmod_threads = 1` (threads: 1,2,4,6,8,12,16) — baseline (writer disabled due to disk)
- [x] T14: Run highres benchmark sweep with `cholmod_threads = -1` (threads: 1,2,4,6,8,12,16) — multi-threaded CHOLMOD
- [x] T15: Download both perf.csv sets and generate comparison report

### Group 7: Analysis and documentation

- [x] T16: Compare `thermal_s` and `solve_s` columns between sweeps; quantify speedup
- [x] T17: Update skill-benchmarking if new parameter or workflow changes warrant it
- [x] T18: Commit and push all changes

## Lessons Learned

### CHOLMOD `nthreads_max` does NOT parallelize factorization
Setting `cholmod_common.nthreads_max` only controls threading for CHOLMOD's internal sparse BLAS calls (e.g., during supernodal factorization). For the thermal solver's sparsity pattern at 1001×801 (~800k cells, 5-point stencil → ~4M nonzeros), the BLAS calls are negligible. The dominant cost is the **symbolic analysis and numeric factorization**, which are serial in CHOLMOD.

**Evidence**: `thermal_s` = ~6.97s ± 0.05s across both sweeps (cholmod_threads=1 vs -1), all 7 thread counts. Zero measurable speedup.

### Parallelizing the thermal solver requires a different approach
Options for future work:
1. **Reuse symbolic factorization** — `cholmod_analyze` is called every timestep but the sparsity pattern is constant. Caching `cholmod_factor *Lfact` across timesteps could save ~30-50% of thermal time.
2. **Parallel direct solver** — MUMPS or PaStiX support distributed-memory parallelism for the factorization itself.
3. **Iterative solver** — CG/PCG with ICC preconditioner would parallelize naturally via SpMV.

### EC2 disk space is tight for writer-enabled benchmarks
With `writer=1, writer_step=1` at 1001×801, each step writes ~300MB of HDF5. A 10-step benchmark run produces ~3GB, exceeding the 3.1GB free on the 6.8GB root volume. Either resize the EBS volume or use a separate data volume for benchmark output.

### The `--writer` and `--sed` flags added to benchmark.sh are useful infrastructure
Even though CHOLMOD threading didn't help, the `--sed` flag enables arbitrary `.txt` parameter sweeps from the command line, which will be valuable for future A/B benchmarks (e.g., persistent OMP thread pool).
