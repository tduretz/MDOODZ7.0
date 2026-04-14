## Context

Experiments 1–7 benchmarked on EC2 c5ad.4xlarge (x86_64, 16 vCPUs, 32 GB) using RiftingComprehensive at 1000×800 with 10 steps. That workload is large-grid, few-step, rheology-heavy with output disabled. We now want to profile on Apple Silicon (M1 Pro, 8 cores — 4P+4E, 16 GB unified memory) with a fundamentally different workload: BlankenBench thermal convection (41×41, 1000 steps, single phase, isoviscous, HDF5 output enabled).

Key differences from EC2 benchmarks:
- **Small grid, many steps**: 41×41 = 1,681 cells (vs 800,000 on EC2). Each step is cheap, so overhead (I/O, logging, loop bookkeeping) becomes visible.
- **Thermal-dominated**: Blankenbach is all about thermal convection — thermal solver runs every step. PCG benefit should be dramatic.
- **Writer ON**: HDF5 writes at steps 500 and 1000 (writer_step=500). First time measuring I/O overhead.
- **M1 unified memory**: No NUMA, shared L2, different atomics performance than x86.

## Goals / Non-Goals

**Goals:**
- Measure full 1000-step BlankenBench wall time and phase breakdown on M1 at 1, 2, 4, 6, 8 threads
- Compare baseline (thermal_solver=0, interp_mode=0) vs optimized (thermal_solver=1, interp_mode=3)
- Quantify HDF5 write overhead (visible in `output_s` column)
- Identify M1-specific scaling characteristics (P-core → E-core transition at 4→6 threads)
- Document as Experiment 8 in skill-benchmarking SKILL.md

**Non-Goals:**
- No code changes — pure benchmark study
- Not comparing M1 vs EC2 directly (different grid sizes, step counts)
- Not optimizing the BlankenBench scenario itself
- Not testing higher grid resolutions (41×41 is the CI config)

## Decisions

### D1: Use benchmark.sh with BlankenBench .txt symlinked into SETS/

**Problem**: `benchmark.sh` reads `.txt` from `SETS/${SCENARIO}.txt`, but BlankenBench's `.txt` lives in `TESTS/BlankenBench/BlankenBench.txt`, not `SETS/`.

**Decision**: Copy `TESTS/BlankenBench/BlankenBench.txt` to `SETS/BlankenBench.txt` before benchmarking, then remove it after. This is simpler than modifying benchmark.sh.

**Alternative considered**: Add a `--txt-dir` flag to benchmark.sh. Rejected — over-engineering for a one-off benchmark.

### D2: Run full 1000 steps, not reduced 10-step runs

**Problem**: EC2 benchmarks used 10 steps to keep runtime short. BlankenBench at 41×41 is fast enough to run all 1000 steps (~5 min for CI).

**Decision**: Use `--steps 1000` to match the CI config exactly. This gives us:
- Statistically meaningful per-step averages (1000 samples vs 10)
- Realistic HDF5 write pattern (writes at 500 and 1000)
- PCG warm-start settling behavior visible across many steps
- Total runtime ~5 min × 5 threads × 2 configs = ~50 min

### D3: Two benchmark sweeps — baseline then optimized

**Sweep 1 (baseline)**: `--steps 1000 --writer --threads "1 2 4 6 8"` with no `--sed` overrides (uses defaults: thermal_solver=0, interp_mode=0).

**Sweep 2 (optimized)**: Same but with `--sed 's/^thermal_solver.*/thermal_solver = 1/' --sed 's/^interp_mode.*/interp_mode = 3/'`.

Both with `--writer` to enable HDF5 output.

### D4: Thread counts 1, 2, 4, 6, 8

M1 has 4 performance + 4 efficiency cores. The 4→6 and 6→8 transitions cross from P-cores to E-cores, which is where we expect scaling anomalies. No point testing 12/16 (only 8 hardware threads).

### D5: Writer step = 1 (every step)

The benchmark.sh `--writer` flag sets `writer = 1` and `writer_step = 1`, meaning HDF5 output every step. We keep this default — writing every step makes `output_s` clearly visible in the phase breakdown and gives an accurate measurement of I/O overhead. At 41×41 each HDF5 file is small (~few KB), so 1000 files won't be a disk concern.

## Risks / Trade-offs

- **M1 P-core/E-core asymmetry**: macOS scheduler may place threads on E-cores unpredictably at 5–8 threads, causing run-to-run variance. Mitigation: accept this as a data point about M1 scheduling, not a code bug.
- **Small grid amplifies overhead**: At 41×41 (1,681 cells, ~27K particles), per-step work is tiny. Timer resolution, function call overhead, and memset costs become proportionally larger. Ratios won't match EC2 large-grid results.
- **Fused P2Mastah less impactful at small grids**: With only 1,681 grid nodes, the stencil weight computation is very cheap. The fusion benefit (avoiding redundant weight computation) might be noise-level. This is expected and informative.
- **Disk I/O**: HDF5 writes at 41×41 are small (~few KB per file). If `output_s` is noise, that confirms I/O isn't a bottleneck for small grids. The real question is whether it matters at large grids — but we can't test that locally (no large-grid .txt for BlankenBench).

## Open Questions

_(none)_
