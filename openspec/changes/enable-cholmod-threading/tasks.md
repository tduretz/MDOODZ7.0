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

- [ ] T12: Start EC2 instance, pull changes, build
- [ ] T13: Run highres benchmark sweep with `cholmod_threads = 1` (threads: 1,2,4,6,8,12,16) — baseline with writer enabled
- [ ] T14: Run highres benchmark sweep with `cholmod_threads = -1` (threads: 1,2,4,6,8,12,16) — multi-threaded CHOLMOD
- [ ] T15: Download both perf.csv sets and generate comparison report

### Group 7: Analysis and documentation

- [ ] T16: Compare `thermal_s` and `solve_s` columns between sweeps; quantify speedup
- [ ] T17: Update skill-benchmarking if new parameter or workflow changes warrant it
- [ ] T18: Commit and push all changes
