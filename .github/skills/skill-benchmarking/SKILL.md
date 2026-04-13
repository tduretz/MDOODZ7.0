---
name: skill-benchmarking
description: MDOODZ performance benchmarking — running benchmarks, interpreting perf.csv, grid scaling studies, thread scaling, comparing machines, and the benchmark.sh script.
---

# MDOODZ Performance Benchmarking

## Quick Start

```bash
# Build and run benchmarks (from repo root)
chmod +x misc/benchmark.sh
misc/benchmark.sh

# Quick mode (small grids, few steps)
misc/benchmark.sh --quick

# Custom grid-based configuration
misc/benchmark.sh --grids "301x201 501x501" --threads "1 4 8" --steps 5

# Resolution-based with a specific scenario
misc/benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default medres" --threads "1 4 8"

# Validate first (3 steps per resolution), then benchmark
misc/benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default" --validate
```

## Benchmark Script (`misc/benchmark.sh`)

The script automates grid-scaling and thread-scaling benchmarks:

1. Builds the scenario with `-O3` and OpenMP
2. Optionally validates each resolution (3 steps, 1 thread, bail on failure)
3. For each grid/resolution × thread count combination:
   - Patches or copies the `.txt` file (Nx, Nz, Nt, disables HDF5 output)
   - Runs with `OMP_NUM_THREADS` set
   - Collects `perf.csv` per run
4. Produces `summary.csv` with per-run averages

Two modes:
- **Grid-based** (default): Uses `--grids` to patch Nx/Nz into a single .txt file
- **Resolution-based**: Uses `--resolutions` to select pre-defined `<Scenario>_<res>.txt` files from `SETS/`

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--scenario NAME` | `RiftingBasic` | Scenario name (builds with `-DSET=NAME`) |
| `--grids "NxAxNzA ..."` | `150x100 301x201 501x501` | Grid sizes (grid mode) |
| `--resolutions "res1 res2"` | — | Resolution names (resolution mode). Use `default` for base .txt |
| `--threads "1 2 4 8"` | `1 2 4 8 12 16` | Thread counts |
| `--steps N` | `10` | Timesteps per run |
| `--validate` | — | Run 3-step validation before full benchmark |
| `--quick` | — | Shorthand: small grids, 3 steps, 1 and 4 threads |
| `--skip-build` | — | Skip the build step |

### Output Structure

```
benchmark-results/
  20260410-143000/
    system.txt              # hostname, OS, arch, RAM, CPU count
    summary.csv             # one row per (grid × threads) combination
    150x100_1t/
      bench.txt             # patched parameter file used
      perf.csv              # per-timestep metrics
      stdout.log            # full simulation output
    150x100_4t/
      ...
    301x201_1t/
      ...
```

### summary.csv Columns

```
host,os,arch,grid,nx,nz,threads,steps,avg_wall_s,total_wall_s,
avg_rheology_s,avg_assembly_s,avg_solve_s,avg_thermal_s,avg_advection_s,
avg_melting_s,avg_anisotropy_s,avg_gse_s,avg_output_s,
avg_interp_s,avg_stokes_setup_s,avg_nl_overhead_s,avg_post_solve_s,
avg_nit,peak_rss_mb
```

## perf.csv (Per-Timestep Metrics)

Every MDOODZ run writes `perf.csv` automatically — see the `skill-logging` skill for column definitions.

## What to Look For

### Grid Scaling (fix threads, vary grid)

Compare `avg_wall_s` across grid sizes at the same thread count:
- Linear scaling (2× grid → 2× time) = memory-bandwidth bound
- Quadratic scaling (2× grid → 4× time) = solver-dominated (expected for direct solve)

### Thread Scaling (fix grid, vary threads)

Compare `avg_wall_s` across thread counts at the same grid:
- Perfect scaling: 2× threads → 0.5× time
- If speedup plateaus early → memory bandwidth saturated
- If speedup degrades → OpenMP overhead or false sharing

### Memory (peak_rss_mb)

- Should grow roughly proportional to grid size × particles per cell
- Unexpected growth between timesteps → particle reseeding leak
- If approaching system RAM → expect swap slowdown or OOM

### Cross-Machine Comparison

Run on multiple machines, then compare `summary.csv` files:
```bash
# Merge results from two machines
head -1 machine_a/summary.csv > combined.csv
tail -n+2 machine_a/summary.csv >> combined.csv
tail -n+2 machine_b/summary.csv >> combined.csv
```

Key comparisons:
- ARM (M1) vs x86 (Zen 2) at same grid/threads
- Memory bandwidth effect: compare `avg_solve_s` (CHOLMOD is bandwidth-bound)
- Per-thread efficiency: `avg_wall_s × threads` (lower = better parallel efficiency)

## Platform Notes

### macOS (M1/Apple Silicon)

```bash
# Install dependencies
brew install suite-sparse hdf5 libomp

# May need to set OpenMP paths
export LDFLAGS="-L$(brew --prefix libomp)/lib"
export CPPFLAGS="-I$(brew --prefix libomp)/include"
```

- `ru_maxrss` reports bytes (auto-handled in code via `__APPLE__` ifdef)
- M1 has 4 performance + 4 efficiency cores — benchmark with `--threads "1 2 4 8"` to see the perf→efficiency core transition
- 16 GB RAM limits max grid to roughly 500×500 with particles

### Linux (Hetzner/AWS)

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt install build-essential cmake libsuitesparse-dev libhdf5-dev libblas-dev liblapack-dev
```

- Use `numactl --interleave=all` for NUMA machines (32+ cores)
- Close other processes for clean measurements

## Interpreting Results

A healthy benchmark should show:
- `peak_rss_mb` stable after first 2-3 timesteps (no memory leaks)
- `avg_wall_s` decreasing roughly linearly as threads increase (up to memory bandwidth limit)
- `avg_solve_s` dominating `avg_wall_s` for large grids (CHOLMOD is the bottleneck)
- `avg_rheology_s` small relative to solve (particle interpolation is cheap)
- Subsystem columns (`thermal_s`, `advection_s`, `melting_s`, `anisotropy_s`, `gse_s`) reveal where time is spent beyond the Stokes solver
- Aggregate timers (`interp_s`, `stokes_setup_s`, `nl_overhead_s`, `post_solve_s`) close the instrumentation gap — coverage typically exceeds 99%

## AWS Benchmarking

Run full benchmark sweeps on a dedicated EC2 instance for reproducible results.

### Prerequisites

1. **AWS CLI** installed and configured (`aws configure`)
2. **EC2 instance** — `c5ad.4xlarge` (16 vCPUs, 32 GB RAM) or similar compute-optimised instance
3. **SSH key** — `.pem` file for the instance (e.g. `~/.ssh/key-mdoodz.pem`, permissions `600`)
4. **S3 bucket** — for result uploads

### Required IAM Permissions

```json
{
  "Effect": "Allow",
  "Action": [
    "ec2:StartInstances",
    "ec2:StopInstances",
    "ec2:DescribeInstances",
    "s3:PutObject",
    "s3:GetObject",
    "s3:ListBucket",
    "s3:DeleteObject",
    "sts:GetCallerIdentity"
  ],
  "Resource": "*"
}
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BENCH_EC2_INSTANCE` | `i-0df5754e400cb3ec1` | EC2 instance ID |
| `BENCH_SSH_KEY` | `~/.ssh/key-mdoodz.pem` | Path to SSH private key |
| `BENCH_S3_BUCKET` | `mdoodz-bench-s3-bucket` | S3 bucket name |
| `BENCH_SSH_USER` | `ubuntu` | SSH username |
| `BENCH_REGION` | `eu-central-1` | AWS region |

### Quick Start (AWS)

```bash
# 1. Run preflight to verify connectivity
./misc/benchmark-preflight.sh

# 2. Launch full benchmark (starts instance, syncs, builds, benchmarks, stops)
./misc/benchmark-aws.sh --scenario RiftingComprehensive --resolutions "lowres default medres highres"

# 3. Results are downloaded locally and uploaded to S3
ls benchmark-results/
```

### Scripts

| Script | Runs where | Purpose |
|--------|-----------|--------|
| `misc/benchmark-preflight.sh` | Local | Check AWS CLI, EC2, SSH, S3 |
| `misc/benchmark-aws.sh` | Local | Orchestrate: start → rsync → benchmark → download → stop |
| `misc/benchmark-ec2-setup.sh` | EC2 | Install deps, build, run benchmark, upload to S3 |
| `misc/benchmark.sh` | Either | Core benchmark runner |
| `misc/benchmark-report.sh` | Either | Generate Markdown report from results |

### Safety

- `benchmark-aws.sh` has a **trap handler** that stops the instance on exit/interrupt/error
- Remote benchmark runs with a **6-hour timeout** (configurable via `--timeout`)
- Partial results are uploaded to S3 even on timeout
- Instance is **stopped** (not terminated) — preserves disk for reuse
- `*.pem` files are in `.gitignore`

### Cost Guidance

- `c5ad.4xlarge`: ~$0.69/hr (eu-central-1, on-demand)
- Full sweep (4 resolutions × 7 threads × 10 steps): typically 30–60 minutes ≈ $0.35–0.70
- Instance is stopped automatically — no ongoing charges after benchmark completes
- The attached nvme storage on c5ad instances is ephemeral — results are uploaded to S3

## Troubleshooting and Lessons Learned

### env.cmake breaks EC2 builds

The macOS `env.cmake` (pointing at `/opt/homebrew`) gets synced to EC2 and causes build failures (`No rule to make target libomp.dylib`). **Fix**: `benchmark-aws.sh` excludes `env.cmake` from rsync. The `benchmark-ec2-setup.sh` script creates its own minimal `env.cmake` on EC2. If building manually on EC2, remove or rename the macOS `env.cmake` first.

### SSH drops kill the benchmark

Long-running SSH sessions to EC2 can drop (exit 255), which triggers the EXIT trap and stops the instance mid-benchmark. **Fix**: `benchmark-aws.sh` uses `nohup` to launch the remote benchmark as a detached process, then polls via short SSH connections every 60 seconds. The EXIT trap only stops the instance when `STOP_ON_EXIT=1` (set explicitly after successful download, not automatically on disconnect).

### AWS CLI not installed on EC2

Fresh Ubuntu instances don't have `aws` CLI. The `benchmark-ec2-setup.sh` S3 upload will fail silently if `awscli` is not installed. **Fix**: Either install `awscli` in the EC2 setup, or download results locally and upload from the local machine (which has AWS CLI).

### Scenario .c files hardcode the .txt filename

`RiftingComprehensive.c` calls `RunMDOODZ("RiftingComprehensive.txt", ...)` — it ignores any CLI arguments. **Fix**: `benchmark.sh` copies the bench `.txt` file to `${EXEC_DIR}/${SCENARIO}.txt` so the exe finds it. The `.txt` is cleaned up after validation to prevent stale parameters from leaking into benchmark runs.

### Empirical Performance Characteristics (c5ad.4xlarge, RiftingComprehensive)

Based on the 28-run sweep (4 resolutions × 7 thread counts × 10 steps):

- **Optimal thread count: 6–8** on 16-vCPU machines. Beyond 8 threads, OpenMP overhead and serial bottlenecks degrade performance.
- **Thermal solver is the primary bottleneck** at medium-to-high resolution. CHOLMOD direct solve takes ~6.4s/step at 1000×800, regardless of thread count (~38% of wall time at 8 threads). This is inherently serial.
- **Rheology and advection scale well** with threads (near-linear to 4–6 threads).
- **Assembly does not scale** well — stays near-constant or increases with threads.
- **Memory scales linearly** at ~15.5 MB per 1000 grid cells, independent of thread count.
- **Grid scaling is linear** $O(n)$: doubling resolution in each dimension (4× cells) costs ~4× wall time.
- **10 steps is sufficient** for benchmarking — per-step cost is stable across steps.

## Optimization Experiments Log

### Experiment 1: CHOLMOD Multi-Threading (2026-04-13) — NO EFFECT

**Hypothesis**: Setting `cholmod_common.nthreads_max` to match OMP thread count would parallelize the thermal solver's CHOLMOD factorization, reducing `thermal_s`.

**Implementation**: Added `cholmod_threads` parameter to `.txt` config (default 1). Values: `1` = single-thread, `-1` = all OMP threads, `N` = explicit count. Applied to all 3 `cholmod_start` sites (Main_DOODZ.c, ChemicalRoutines.c, ThermalRoutines.c).

**Results** (c5ad.4xlarge, 1001×801, 10 steps, writer=0):

| Threads | thermal_s (threads=1) | thermal_s (threads=auto) | Delta |
|---------|----------------------|--------------------------|-------|
| 1 | 6.97s | 6.97s | 0.0% |
| 4 | 6.96s | 6.96s | 0.0% |
| 8 | 6.95s | 6.95s | 0.0% |
| 16 | 6.95s | 6.96s | 0.0% |

**Conclusion**: Zero measurable speedup. `nthreads_max` only controls CHOLMOD's internal sparse BLAS calls, which are negligible for the 5-point stencil sparsity pattern. The dominant cost is serial symbolic analysis + numeric factorization. Parallelizing thermal solve requires either (a) caching `cholmod_analyze` across timesteps (sparsity pattern is constant), (b) a parallel direct solver (MUMPS/PaStiX), or (c) an iterative solver (CG+ICC).

**Data**: `benchmark-results/cholmod-sweep1-threads1/`, `benchmark-results/cholmod-sweep2-threads-auto/`

### Experiment 2: Internal Persistent OMP Regions (2026-04-13) — NO EFFECT / SLIGHT REGRESSION

**Hypothesis**: Consolidating multiple `#pragma omp parallel for` (fork-join) loops into single `#pragma omp parallel` regions within P2Mastah and NonNewtonianViscosityGrid would reduce thread pool fork-join overhead, improving `interp_s` and `rheology_s` at high thread counts.

**Implementation**:
- P2Mastah: Merged 5 separate fork-join loops (init, scatter, reduction, node values) into 1 persistent `#pragma omp parallel` region using `#pragma omp for` worksharing. Changed `omp_get_num_threads()` to `omp_get_max_threads()` for pre-allocation.
- NonNewtonianViscosityGrid: Merged 2 fork-join loops (centroid + vertex viscosity) into 1 persistent region.
- Added `omp_schedule` parameter (placeholder, default 0).

**Results** (c5ad.4xlarge, 1001×801, 10 steps, writer=0):

| Threads | Baseline wall | New wall | Delta | Baseline interp | New interp | Delta |
|---------|--------------|----------|-------|-----------------|------------|-------|
| 1 | 36.44s | 36.75s | +0.8% | 7.45s | 7.49s | +0.5% |
| 2 | 25.56s | 25.74s | +0.7% | 4.63s | 4.37s | -5.6% |
| 4 | 18.38s | 19.68s | +7.1% | 3.16s | 3.17s | +0.3% |
| 6 | 16.98s | 18.24s | +7.4% | 3.08s | 3.24s | +5.2% |
| 8 | 17.24s | 18.41s | +6.8% | 3.62s | 4.09s | +13.0% |
| 12 | 20.68s | 20.53s | -0.7% | 6.25s | 5.75s | -8.0% |
| 16 | 20.48s | 20.88s | +2.0% | 6.50s | 6.47s | -0.5% |

**Conclusion**: No improvement; 4–8 threads show a consistent ~7% regression in wall time. The fork-join overhead at 1001×801 is negligible compared to actual computation — the OS thread pool reuses threads efficiently. The regression at 4–8 threads may be due to `omp_get_max_threads()` over-allocating thread-local buffers (16 instead of actual 4–8), causing cache pressure. Internal persistent regions are not worthwhile for this codebase at production grid sizes.

**Key lesson**: OpenMP fork-join overhead is only significant for very small loops or very high thread counts. At 1001×801 with ~12M particles, each loop body does enough work that fork-join cost is in the noise. Future optimization should focus on algorithmic improvements (thermal solver caching, assembly parallelism) rather than reducing fork-join overhead.

**Data**: Baseline `benchmark-results/20260413-105848/`, Persistent OMP `benchmark-results/20260413-151213/`
**Code**: Reverted — these changes were not kept.

### Experiment 3: Cache Thermal Factorization (2026-04-13) — MODEST GAIN

**Hypothesis**: Caching the `cholmod_factor` across timesteps (skip `cholmod_analyze` on steps 2+) would reduce `thermal_s` by eliminating repeated symbolic analysis. The sparsity pattern is constant when the free-surface topology doesn't change.

**Implementation**: Added persistent `DirectSolver ThermalSolver` to `Main_DOODZ.c`. Split `FactorEnergyCHOLMOD` into analyze+factorize paths: first call does `cholmod_analyze` + `cholmod_factorize`, subsequent calls reuse cached `cholmod_factor` with only `cholmod_factorize`. Added dimension invalidation check — when `neq` changes (free surface), cached factor is freed and re-analyzed.

**Results** (c5ad.4xlarge, 1001×801, 10 steps):

| Threads | Baseline thermal_s | Cached thermal_s | Delta | Wall baseline | Wall cached |
|---------|-------------------|-----------------|-------|--------------|-------------|
| 1 | 6.557s | 6.192s | −5.6% | 36.67s | 36.27s |
| 2 | 6.562s | 6.178s | −5.9% | 25.64s | 25.04s |
| 4 | 6.387s | 6.184s | −3.2% | 18.49s | 19.35s |
| 6 | 6.396s | 6.174s | −3.5% | 17.00s | 19.21s |
| 8 | 6.439s | 6.166s | −4.2% | 16.97s | 18.93s |
| 12 | 6.277s | 6.168s | −1.7% | 17.88s | 20.56s |
| 16 | 6.499s | 6.154s | −5.3% | 20.33s | 20.82s |

**Conclusion**: Caching saves ~0.2–0.4s/step (3–6%) on thermal_s by eliminating `cholmod_analyze`. The numeric factorization (`cholmod_factorize`) at ~6.2s/step remains the dominant serial cost. Overall wall time shows run-to-run variance in `interp_s` that masks the thermal gain at higher thread counts. The change is correct and low-risk but insufficient alone — an iterative solver (PCG) is needed for order-of-magnitude thermal speedup.

**Data**: Baseline `benchmark-results/20260413-105848/`, Cached `benchmark-results/20260413-170125/`
**Code**: Committed as `6285d3d` on `add-performance-metrics` branch. Kept.

### Experiment 4: PCG Thermal Solver (2026-04-13) — MAJOR GAIN

**Hypothesis**: Replacing serial CHOLMOD direct solve with a Preconditioned Conjugate Gradient (PCG) iterative solver using Jacobi preconditioner and warm-start from the previous timestep's temperature would dramatically reduce `thermal_s`. The thermal matrix is SPD with a 5-point stencil — ideal for CG. Warm-start should converge in very few iterations since temperature changes slowly between steps.

**Implementation**: Added `SolveThermalPCG()` in `ThermalSolver.c` — full CG with Jacobi (diagonal) preconditioner, OpenMP-parallel SpMV/dot/axpy. Dispatch in `EnergyDirectSolve()`: PCG on `it==0` (first non-linear iteration), CHOLMOD fallback on PCG failure or `it>=1`. Warm-start x from `mesh->T` via equation numbering. New params: `thermal_solver` (0=CHOLMOD, 1=PCG), `max_its_thermal` (default 1000), `rel_tol_thermal` (default 1e-8).

**Convergence pattern**: Step 1 (cold start, no warm-start) → PCG fails at 1000 iterations, CHOLMOD fallback. Steps 2+: PCG converges in **6–7 iterations** with warm-start.

**Results** (c5ad.4xlarge, 1000×800, 10 steps, comparison vs cached-factorization CHOLMOD baseline):

| Threads | CHOLMOD wall | PCG wall | Wall Δ | CHOLMOD thermal | PCG thermal | Thermal Δ |
|---------|-------------|---------|--------|----------------|------------|-----------|
| 1 | 36.27s | 30.24s | −16.6% | 6.192s | 0.159s | −97.4% |
| 2 | 25.04s | 18.98s | −24.2% | 6.178s | 0.115s | −98.1% |
| 4 | 19.35s | 13.31s | −31.2% | 6.184s | 0.101s | −98.4% |
| 6 | 19.21s | 13.24s | −31.1% | 6.174s | 0.098s | −98.4% |
| 8 | 18.93s | 12.83s | −32.2% | 6.166s | 0.096s | −98.4% |
| 12 | 20.56s | 14.36s | −30.2% | 6.168s | 0.099s | −98.4% |
| 16 | 20.82s | 14.67s | −29.5% | 6.154s | 0.099s | −98.4% |

**Conclusion**: PCG reduces thermal solver time by **~98%** (from ~6.2s to ~0.1s per step). Overall wall time drops 17–32% depending on thread count, with optimal improvement at 8 threads. Thermal went from **33% of wall time to under 1%**. The remaining bottleneck is now `interp_s` (P2Mastah particle interpolation), which anti-scales above 4 threads: 3.5s at 4 threads → 7.1s at 16 threads. The optimal thread count remains 6–8 due to `interp` anti-scaling.

**Key observation**: Wall time no longer scales beyond 8 threads — the `interp` anti-scaling dominates. Next optimization target should be P2Mastah (persistent buffers to eliminate malloc/free churn, and/or atomic scatter to eliminate thread-local buffer memory explosion).

**Data**: CHOLMOD baseline `benchmark-results/20260413-170125/`, PCG `benchmark-results/20260413-182325/`
**Code**: Committed on `add-performance-metrics` branch. Kept.
