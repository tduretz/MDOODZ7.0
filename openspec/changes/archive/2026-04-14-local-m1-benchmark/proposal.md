## Why

All Experiments 1–7 benchmarked on EC2 c5ad.4xlarge (x86_64, 16 vCPUs) using RiftingComprehensive at 1000×800 with 10 steps — a large-grid, few-step, rheology-heavy workload. We have no data on how the optimizations perform on Apple Silicon (M1, 8 cores, unified memory) with a fundamentally different workload: BlankenBench (41×41 grid, 1000 steps, thermal-convection dominated, single phase, isoviscous). This also lets us measure HDF5 write overhead (writer=1), which was disabled in all EC2 runs.

## What Changes

- Run BlankenBench locally on MacBook M1 (8 cores, 16 GB) with full 1000-step CI config
- Two configurations compared:
  - **Baseline**: default settings (thermal_solver=0 CHOLMOD, interp_mode=0 legacy)
  - **Optimized**: thermal_solver=1 (PCG), interp_mode=3 (fused P2Mastah)
- Thread scaling: 1, 2, 4, 6, 8 threads
- Writer ON (HDF5 output enabled, writer_step=500 → writes at steps 500 and 1000)
- Collect per-step perf.csv, analyze phase breakdown, identify M1-specific bottlenecks
- Document results as Experiment 8 in skill-benchmarking SKILL.md

## Capabilities

### New Capabilities

- `local-m1-benchmark`: Specification for running local M1 benchmarks — scenarios, thread counts, writer settings, comparison methodology between optimized and baseline configurations on Apple Silicon

### Modified Capabilities

_(none — no requirement-level changes to existing specs)_

## Impact

- **No code changes** — this is a benchmark-only change using existing interp_mode and thermal_solver features
- **skill-benchmarking/SKILL.md** will be updated with Experiment 8 results (M1 BlankenBench)
- **benchmark.sh** may need minor adjustments if BlankenBench scenario isn't directly supported in grid mode (it uses TESTS/ .txt files, not SETS/ .txt files)
- **Disk**: ~1000 steps × 2 writes = small HDF5 output, not a storage concern
