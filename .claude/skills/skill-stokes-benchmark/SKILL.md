---
name: skill-stokes-benchmark
description: Benchmark the MDOODZ Stokes solver on a hard-solve state — build a shear-band checkpoint with a warmup run, restart from the .dat file, and sweep thread counts / interp_mode to measure solve_s under realistic high-viscosity-contrast conditions. Use when the user wants to measure solver cost (not wall time dominated by advection/interp) or compare optimizations (recycled factorization, CHOLMOD threading, Krylov, etc).
---

# Stokes Solver Benchmarking

## Why this skill exists

`skill-benchmarking` covers the general `benchmark.sh` workflow (cold-start, multi-scenario, grid sweeps). This skill is narrower: it focuses on **measuring solver cost in isolation** by running from a checkpoint where shear bands have already formed.

The problem with cold-start benchmarking: on a fresh model, most steps have `nit=0` (convergence at iteration 0, no solve). You measure interpolation and advection, not the solver. To benchmark the solver you need a model state where:

1. Shear bands exist (steep viscosity contrasts)
2. Each step does 1–4 Picard iterations
3. Each iteration's factorization is large enough to actually time (≥1s)

Pattern: run a warmup, save `Breakpoint00020.dat` (or similar), restart from it, sweep configurations.

## When to use this skill

- "Which solver optimization gives the biggest gain?"
- "How does CHOLMOD/KillerSolver scale with threads?"
- "Is recycled factorization worth implementing?"
- "Baseline the solver before/after a change"
- Any request where advection and interpolation should NOT dominate the measurement

**Do not use** for:
- Cold-start runtime comparison — use `skill-benchmarking`
- Grid-scaling studies — use `benchmark.sh --grids`
- Warm-restart production workloads where `nit=0` everywhere (solver never runs)

## The core workflow

### Phase 1 — Warmup (once per scenario+grid)

Run from `irestart=0` long enough for shear bands to form. Save periodic checkpoints with `writer=1`.

For **RiftingChenin 1000×800** on an M1 laptop this takes ~5 minutes and produces ~3 GB checkpoints at step 10 and 20. Step 19+ has `nit≥1`.

### Phase 2 — Benchmark (repeat per configuration)

Set `irestart=1`, `istep=<checkpoint step>`, run 5 steps, record `perf.csv`.

### Phase 3 — Sweep

Loop over `OMP_NUM_THREADS` and `interp_mode`. Average steps 2–5 of each run (step 1 after restart has cold-cache overhead).

## Executable ↔ config wiring

**Each scenario's executable has the config filename hardcoded** in its `main()`:

```c
// e.g. SETS/RiftingChenin.c
int main() { RunMDOODZ("RiftingChenin.txt", &setup); }
```

Passing a different filename as argument is **silently ignored**. To run a different config, **edit the `.txt` file in place** (or `cp` a variant on top of it). The `.txt` must live next to the executable:

```
cmake-exec/RiftingChenin/RiftingChenin        # executable
cmake-exec/RiftingChenin/RiftingChenin.txt    # config it reads
cmake-exec/RiftingChenin/perf.csv             # output (written after each step)
cmake-exec/RiftingChenin/Breakpoint00020.dat  # checkpoint
```

Always keep `<scenario>_orig.txt` as a backup before editing.

## Config parameters that matter

### Restart and checkpoints

| Param | Warmup value | Benchmark value | Notes |
|-------|--------------|-----------------|-------|
| `irestart` | 0 | 1 | 0 = cold start, 1 = load `.dat` |
| `istep` | 0 | `<checkpoint step>` | Which `BreakpointNNNNN.dat` to load |
| `writer` | 1 | 0 | 1 = save `.dat` checkpoints |
| `writer_step` | 10 | — | Period between checkpoint writes |
| `Nt` | warmup depth | `istep + 5` | Loop runs while `step <= Nt` |

**Gotcha on `Nt`**: the time loop is `for (; step <= Nt; step++)`. Restarting at step 20 with `Nt=5` does nothing (`21 > 5`). For 5 benchmark steps from step 20 you need `Nt=25`.

### Solver selection (after Apr 2026 rebuild)

`InputOutput.c:1327` force-redirects `lin_solver=0` → 2 (KillerSolver) whenever `lin_solver==0`, `Newton==1`, or `anisotropy==1`. Valid values are `-1, 0, 1, 2`. The three solver paths that actually run are:

| `lin_solver` | Code path (`StokesRoutines.c:788`) | Notes |
|--------------|------------------------------------|-------|
| `-1` | `DirectStokesDecoupledComp` | CHOLMOD direct, compressible (BlankenBench uses this) |
| `1` | `KSPStokesDecoupled` | Krylov |
| `2` | `KillerSolver` (`Solvers.c:1054`) | Block-preconditioned Schur-complement; the current default |
| `3` | — | Removed (was GMG-FGMRES) |

If you want CHOLMOD directly, set `lin_solver = -1` (only works for compressible setups).

### Threading

| Param | Purpose |
|-------|---------|
| `OMP_NUM_THREADS` (env) | OMP thread count. On M1: 4 = 4 P-cores (optimal), 8 mixes in E-cores and often regresses |
| `cholmod_threads = 1` | CHOLMOD internal threading; keep at 1 on M1 (measured optimal) |

### Physics knobs that change solve cost

| Param | Effect on solver load |
|-------|----------------------|
| `nit_max` | Cap on Picard iterations per step. Set ≥ 5 for benchmarking (don't artificially clip) |
| `Newton` | 0 = Picard, 1 = Newton. Newton forces `lin_solver=2` |
| `anisotropy` | 1 forces `lin_solver=2` |
| `interp_mode` | 0 = classic, 3 = fused P2Mastah (faster at 1t, does NOT scale with threads — see below) |
| `thermal_solver` | 1 = PCG (fast), 0 = direct. PCG saves ~1–2s/step at 1000×800 |

## Recommended benchmark setup

Target grid/scenario by machine:

| Machine | Scenario | Grid | Expected solve/iter |
|---------|----------|------|--------------------|
| M1 laptop (16 GB) | RiftingChenin | 1000×800 | ~11–13 s (at nit=3) |
| M1 laptop | RiftingChenin | 501×401 | ~1.7 s |
| M1 laptop | RiftingChenin | 150×100 | ~0.08 s (too fast to benchmark meaningfully) |
| EC2 / workstation | RiftingChenin | 1000×800 | ~3–5 s |
| High-mem / HPC | RiftingChenin | 2000×1600 | ~100 s+ |

**Rule of thumb**: `solve_s` scales as `O(N^1.5)` where N = DOF count ≈ 1.86·Nx·Nz. To double solve time, multiply each dimension by 1.26.

For "solver dominates wall time" on M1, **RiftingChenin at 1000×800 is the sweet spot** — step 25 gives `nit=4, solve=21s, wall=30s` (70% solver), memory fits in 16 GB (~3 GB peak RSS).

## Running a warmup — full recipe (RiftingChenin 1000×800)

```bash
cd /Users/romankulakov/MDOODZ7.0/cmake-exec/RiftingChenin

# 1. Backup original config
cp RiftingChenin.txt RiftingChenin_orig.txt

# 2. Patch for warmup
sed -i '' \
  -e 's/^istep = .*/istep = 0/' \
  -e 's/^irestart = .*/irestart = 0/' \
  -e 's/^writer = .*/writer = 1/' \
  -e 's/^writer_step = .*/writer_step = 10/' \
  -e 's/^Nx = .*/Nx = 1000/' \
  -e 's/^Nz = .*/Nz = 800/' \
  -e 's/^Nt = .*/Nt = 25/' \
  -e 's/^nit_max.*=.*/nit_max = 10/' \
  -e 's/^interp_mode = .*/interp_mode = 3/' \
  RiftingChenin.txt

# Add lin_solver and thermal_solver if missing
grep -q '^lin_solver' RiftingChenin.txt || \
  sed -i '' '/^cholmod_threads/a\'$'\n''lin_solver = 0\nthermal_solver = 1' RiftingChenin.txt

# 3. Run (4 threads, ~5 min; saves Breakpoint00010.dat, Breakpoint00020.dat)
rm -f perf.csv Breakpoint*.dat
OMP_NUM_THREADS=4 ./RiftingChenin > warmup.log 2>&1

# 4. Verify the shear band transition
awk -F',' 'NR>1 {printf "step=%s solve=%.2f nit=%s\n", $1, $6, $19}' perf.csv
# Look for: step 17+ with nit>=1, solve_s >= 4s
```

The expected output shows `nit=0, solve=0` for steps 1–16, then `nit=1, solve≈5s` at step 17, climbing to `nit=4, solve≈21s` by step 25.

## Running a benchmark from the checkpoint

```bash
cd /Users/romankulakov/MDOODZ7.0/cmake-exec/RiftingChenin

# 1. Patch for restart
sed -i '' \
  -e 's/^istep = .*/istep = 20/' \
  -e 's/^irestart = .*/irestart = 1/' \
  -e 's/^writer = .*/writer = 0/' \
  -e 's/^Nt = .*/Nt = 25/' \
  RiftingChenin.txt

# 2. Single-config smoke test
rm -f perf.csv
OMP_NUM_THREADS=4 ./RiftingChenin > /dev/null 2>&1
awk -F',' 'NR>1 {printf "step=%s wall=%.1f solve=%.2f(%d%%) nit=%s\n", $1,$2,$6, int($6/$2*100), $19}' perf.csv
# Expected: steps 21-25, wall 10-30s, solve 4-21s, nit 1-4
```

## Full sweep script (copy-paste)

```bash
#!/usr/bin/env bash
# Run from cmake-exec/RiftingChenin/ after warmup
set -e
OUTDIR="../../benchmark-results/stokes-bench-$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUTDIR"
cp RiftingChenin.txt "${OUTDIR}/config.txt"

echo "=== Stokes solver sweep ==="
for IMODE in 3 0; do
  sed -i '' "s/^interp_mode = .*/interp_mode = $IMODE/" RiftingChenin.txt
  for THREADS in 1 2 4 8; do
    rm -f perf.csv
    OMP_NUM_THREADS=$THREADS ./RiftingChenin > /dev/null 2>&1
    cp perf.csv "${OUTDIR}/im${IMODE}_t${THREADS}.csv"
    # Skip step 21 (cold cache after restart), average 22-25
    awk -F',' -v im=$IMODE -v t=$THREADS 'NR>=3 && NR<=6 {
      w+=$2; s+=$6; a+=$5; i+=$15; nit+=$19; n++
    } END {
      printf "  im=%s t=%s: wall=%.1fs solve=%.2fs(%d%%) assem=%.3fs interp=%.2fs nit=%.1f\n",
             im, t, w/n, s/n, int(s/w*100*n/n), a/n, i/n, nit/n
    }' perf.csv
  done
done

# Restore a sensible default
sed -i '' "s/^interp_mode = .*/interp_mode = 3/" RiftingChenin.txt
echo "Done → $OUTDIR"
```

Expected output on M1 with RiftingChenin 1000×800:

```
im=3 t=1: wall=27.2s solve=11.65s(42%) assem=0.852s interp=2.62s nit=3.0
im=3 t=2: wall=20.6s solve=11.15s(54%) assem=0.607s interp=1.87s nit=3.0
im=3 t=4: wall=18.9s solve=11.95s(63%) assem=0.567s interp=1.45s nit=3.0
im=3 t=8: wall=18.9s solve=12.02s(63%) assem=0.784s interp=1.52s nit=3.0
im=0 t=4: wall=19.0s solve=11.41s(59%) assem=0.510s interp=2.03s nit=3.0
```

## perf.csv column reference

Columns (1-indexed):

```
1  step          2  wall_s         3  time_ma        4  rheology_s
5  assembly_s    6  solve_s        7  thermal_s      8  advection_s
9  free_surface_s  10  reseeding_s  11  melting_s     12  anisotropy_s
13  gse_s         14  output_s       15  interp_s      16  stokes_setup_s
17  nl_overhead_s  18  post_solve_s  19  nit           20  n_particles
21  neq_mom       22  neq_cont      23  peak_rss_mb   24  user_cpu_s
25  sys_cpu_s
```

Known quirks:
- `post_solve_s` can be negative on hard steps — the phase accounting double-counts in places (sum of phases > wall_s). Interpret as a residual, not a real timing.
- Step 1 after a restart has cold-cache overhead; always drop it from averages.
- `n_particles` is the total (≈ 16 × Nx × Nz for 4×4 particles/cell).
- `nit=0` means the initial residual already satisfies the convergence check — no solve was called, `solve_s=0`.

## Where checkpoints live

`Breakpoint<NNNNN>.dat` is written to the **current working directory** (the directory you run the executable from — normally `cmake-exec/<Scenario>/`). Size is roughly proportional to `Nx * Nz * n_particles_per_cell`:

| Grid | Checkpoint size |
|------|-----------------|
| 150×100 | ~57 MB |
| 501×401 | ~920 MB–1 GB |
| 1000×800 | ~3 GB |

Restart requires `Breakpoint<istep>.dat` present next to the executable. If you move the executable or run from a different directory, symlink or copy the `.dat`.

**Clean up after benchmarking** — these files are large and don't compress well.

## Reading the solver cost

In the benchmark sweep above, `solve_s ≈ 11.65s` at 1 thread drops to ≈ 11.95s at 4 threads — roughly flat. That's expected: both CHOLMOD (via KillerSolver) and `cholmod_threads=1` are single-threaded on our builds. Everything else (interp, advection, reseeding, post-solve) scales with threads; the solver does not. As thread count rises, the solver fraction of wall time climbs (42% → 63% at 1t → 4t), which is the measurement you want.

**To baseline a solver change**: run the sweep above before and after. Focus on `solve_s` at 4 threads. Anything that reduces `solve_s` by >10% on this workload is a meaningful gain.

## Common mistakes

1. **Passing config filename as argument** — it's ignored. Always edit the `.txt` in place.
2. **Setting `Nt` as "number of steps to run"** — it's actually "max step number". Use `Nt = istep + N_steps`.
3. **Forgetting `writer=0` for benchmark runs** — writing a 3 GB `.dat` takes seconds and contaminates `wall_s`.
4. **Averaging step 1 of the restart** — cold cache inflates timings. Skip it.
5. **Comparing `solve_s` across different `nit` counts** — divide by `nit` to get per-iteration cost.
6. **Using `interp_mode=3` at 8 threads** — it scales poorly; at 4+ threads `interp_mode=0` is comparable or better.
7. **Using `lin_solver=0` and expecting CHOLMOD** — it's silently redirected to `KillerSolver` (type 2). Use `-1` for true CHOLMOD (compressible only).
8. **Running the warmup at the wrong grid** — checkpoints are grid-specific. A `.dat` from 1000×800 cannot restart a 501×401 run.

## Related skills

- `skill-benchmarking` — general-purpose `benchmark.sh` multi-scenario grid sweeps
- `skill-solvers` — what CHOLMOD, KillerSolver, Picard, Newton do internally
- `skill-build-and-run` — how the executables are built and where they live
- `skill-logging` — interpreting the per-step log output
