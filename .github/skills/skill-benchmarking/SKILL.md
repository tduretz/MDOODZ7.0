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
