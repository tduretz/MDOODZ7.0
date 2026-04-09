# Spec: steady-state-parameter-file

## Overview

A separate MDOODZ parameter file `BlankenBench/BlankenBenchSteady.txt` that runs the Blankenbach Case 1a convection problem long enough to reach thermal steady state.

## Requirements

### R1: File location and naming

- File: `TESTS/BlankenBench/BlankenBenchSteady.txt`
- Must coexist with the existing `BlankenBench.txt` (500-step smoke test) without conflict

### R2: Step count sufficient for steady state

- `Nt = 100000` (100k steps)
- With `Courant = 0.5` (doubled from smoke test), dt is roughly doubled → ~2–4 hours wall clock with 8 threads
- The diffusion time is $t_d = H^2/\kappa \approx 5 \times 10^{15}$ s; 100k steps covers ~3.6 diffusion times at the original Courant, or ~7 diffusion times with Courant=0.5 — more than sufficient

### R3: Output only the final step

- `writer_step = 100000` — writes a single HDF5 file (`Output100000.gzip.h5`)
- `writer_subfolder = BlankenBench` — same subfolder as the smoke test
- No intermediate outputs to avoid filling disk (~40 MB per output × 100k steps would be catastrophic)

### R4: Identical physics to the smoke test

All physical parameters must be identical to `BlankenBench.txt`:

| Parameter | Value |
|---|---|
| Nx, Nz | 41, 41 |
| Domain | [-5e4, 5e4] × [-5e4, 5e4] (H = 1e5 m) |
| η | 1e20 Pa·s |
| ρ | 3300 kg/m³ |
| α | 2.42e-5 K⁻¹ (via density model 3) |
| k | 3.3 W/m/K |
| Cp | 1250 J/kg/K |
| g | -10 m/s² |
| T_top | 273.15 K, T_bot | 1273.15 K |
| BCs | Free-slip all walls, insulating sidewalls |
| constant_dt | 0 (adaptive) |
| Solver | Picard (Newton=0), nit_max=10 |

### R5: Differences from smoke test

Only these parameters change:

| Parameter | Smoke test | Steady state |
|---|---|---|
| `Nt` | 500 | 100000 |
| `Courant` | 0.25 | 0.5 |
| `writer_step` | 500 | 100000 |
| `istep` | 001000 | 200000 |

`Courant = 0.5` doubles dt per step, halving the number of steps needed.

### R6: No CI registration

- The parameter file is **not** referenced in any `add_test()` call in CMakeLists.txt
- It is only used when the user manually runs `./BlankenBenchTests --gtest_filter=BlankenBench.NusseltNumber` (or similar)
