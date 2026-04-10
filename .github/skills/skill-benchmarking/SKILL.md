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

# Custom configuration
misc/benchmark.sh --grids "301x201 501x501" --threads "1 4 8" --steps 5
```

## Benchmark Script (`misc/benchmark.sh`)

The script automates grid-scaling and thread-scaling benchmarks:

1. Builds `RiftingBasic` with `-O3` and OpenMP
2. For each grid size × thread count combination:
   - Patches the `.txt` file (Nx, Nz, Nt, disables HDF5 output)
   - Runs with `OMP_NUM_THREADS` set
   - Collects `perf.csv` per run
3. Produces `summary.csv` with per-run averages

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--grids "NxA×NzA NxB×NzB"` | `150x100 301x201 501x501` | Grid sizes to test |
| `--threads "1 2 4 8"` | `1 2 4 8` | Thread counts |
| `--steps N` | `10` | Timesteps per run |
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

`host,os,arch,grid,nx,nz,threads,steps,avg_wall_s,total_wall_s,avg_rheology_s,avg_assembly_s,avg_solve_s,avg_nit,peak_rss_mb`

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
