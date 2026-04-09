# Benchmark Experiments Log

Systematic parameter experiments for the 5 non-SolVi analytical benchmarks.
Each section records raw numbers from every experiment run.

---

## 1. SteadyStateGeotherm

**Analytical**: T(z) = T_top + (T_bot - T_top)·(z_top - z)/H
**File**: ThermalTests.cpp, Thermal/SteadyStateGeotherm.txt

### Baseline

| Nx | Nt | dt | penalty | rel_tol_KSP | L2(T) | L1(T) | Runtime |
|----|----|----|---------|-------------|-------|-------|---------|
| 31 | 20 | 1e12 | 1e5 | 1e-4 | 6.82e-1 | 6.03e+2 | 0.78s |

### Time Stepping Experiments

| Nx | Nt | dt | L2(T) | L1(T) | Notes |
|----|----|----|-------|-------|-------|
| 31 | 20 | 1e13 | 5.22e-1 | 4.68e+2 | total=2e14, 1.7% of τ |
| 31 | 20 | 1e14 | 1.14e-1 | 1.04e+2 | total=2e15, 17% of τ |
| 31 | 20 | 1e15 | 2.42e-6 | 2.21e-3 | total=2e16, 1.7τ — converged! |
| 31 | 20 | 5e15 | 3.96e-8 | 3.13e-5 | total=1e17, 8.7τ — saturated |
| 31 | 20 | 1e16 | 3.96e-8 | 3.13e-5 | total=2e17, 17τ — same as 5e15 |

### Resolution Experiments

| Nx | Nt | dt | L2(T) | L1(T) | Notes |
|----|----|----|-------|-------|-------|
| 31 | 20 | 1e15 | 2.42e-6 | 2.21e-3 | baseline at best dt |
| 61 | 20 | 1e15 | 2.37e-6 | 2.16e-3 | negligible gain, 1889ms (3×) |

### Penalty Experiments

| Nx | Nt | penalty | L2(T) | L1(T) | Notes |
|----|----|----|-------|-------|-------|
| 31 | 20 | 1e5 | 2.42e-6 | 2.21e-3 | at dt=1e15 |
| 31 | 20 | 1e6 | 2.42e-6 | 2.21e-3 | identical |
| 31 | 20 | 1e7 | 2.42e-6 | 2.21e-3 | identical |

### Summary & Recommendation

**Root cause**: Baseline dt=1e12 gives total time = 2e13 s, only 0.17% of thermal diffusion
timescale τ = L²ρCp/k ≈ 1.15e16 s. The simulation never reaches steady state.

**Fix**: Change dt from 1e12 to 1e15. Total time becomes 1.7τ — fully converged.
L2 drops from 0.682 to 2.4e-6.  No resolution or penalty changes needed.
Runtime unchanged (~660ms).

**Recommended CI parameters**: dt=1e15, Nx=31, penalty=1e5 (unchanged), threshold L2 < 1e-4.

---

## 2. RadiogenicHeat

**Analytical**: T_mean(t) = T0 + Qr·t/(ρ·Cp)
**File**: ThermalTests.cpp, Thermal/RadiogenicHeat.txt

### Baseline

| Nx | Nt | dt | meanT_final | T_ana | Abs error | Rel error | Runtime |
|----|----|----|-------------|-------|-----------|-----------|---------|
| 31 | 5 | 1e12 | 774.587 | 774.593 | 0.006 K | 0.0008% | 0.22s |

### Time Stepping Experiments

| Nx | Nt | dt | meanT_final | T_ana | Abs error | Rel error | Notes |
|----|----|----|-------------|-------|-----------|-----------|-------|
| 31 | 5 | 1e13 | 786.944 | 787.580 | 0.636 K | 4.4% | error scales linearly with total time |

### Resolution & BC Experiments

_Skipped_ — error is already 0.008% at baseline. BC losses scale with time, not resolution.

### Summary & Recommendation

**Root cause**: ±50% tolerance is wildly loose. Actual error at baseline is 0.42% of ΔT.
Error grows linearly with total simulated time (boundary diffusion losses).
Current dt=1e12, Nt=5 is the sweet spot.

**Fix**: Tighten tolerance from ±50% to ±2% (5× margin over measured 0.42%).
No parameter file changes needed.

**Recommended CI parameters**: unchanged (Nx=31, Nt=5, dt=1e12).

---

## 3. HydrostaticPressure

**Analytical**: P(z) = ρ·g·(z_top - z)
**File**: DensityTests.cpp, Density/HydrostaticPressure.txt

### Baseline

| Nx | penalty | rel_tol_KSP | L2(P) | L1(P) | Runtime |
|----|---------|-------------|-------|-------|---------|
| 21 | 1e5 | 1e-4 | 8.60e-2 | 1.34e+7 | 0.05s |

### Resolution Experiments

| Nx | penalty | L2(P) | L1(P) | Notes |
|----|---------|-------|-------|-------|
| 21 | 1e5 | 8.60e-2 | 1.34e+7 | baseline |
| 41 | 1e5 | 4.26e-2 | 6.64e+6 | 101ms, first-order convergence |
| 61 | 1e5 | 2.81e-2 | 4.38e+6 | 216ms |

Convergence order 21→41: log(0.086/0.043)/log(41/21) = 1.05 (first-order)

### Penalty Experiments

| Nx | penalty | L2(P) | L1(P) | Notes |
|----|---------|-------|-------|-------|
| 21 | 1e5 | 8.60e-2 | 1.34e+7 | baseline |
| 21 | 1e6 | 8.60e-2 | 1.34e+7 | identical |
| 21 | 1e7 | 8.60e-2 | 1.34e+7 | identical |

### Linear Solver Tolerance Experiments

| Nx | rel_tol_KSP | L2(P) | L1(P) | Notes |
|----|-------------|-------|-------|-------|
| 21 | 1e-4 | 8.60e-2 | 1.34e+7 | baseline |
| 21 | 1e-6 | 8.60e-2 | 1.34e+7 | identical |

### Summary & Recommendation

**Root cause**: Error is purely spatial discretization — first-order convergence.
Penalty and solver tolerance have zero effect. Negative pressure at top is a
staggered-grid boundary artifact (not a bug).

**Fix**: Increase Nx from 21 to 41. L2 drops from 0.086 to 0.043.
Runtime increases from 51ms to 101ms (still fast).
Tighten threshold from 0.2 to 0.06.

**Recommended CI parameters**: Nx=41, penalty=1e5 (unchanged), threshold L2 < 0.06.

---

## 4. PureShearVelocity

**Analytical**: Vx = ε̇·x, Vz = -ε̇·z (perturbed by inclusion η_c/η_m=1000)
**File**: VelocityFieldTests.cpp, VelocityField/PureShearVelocity.txt

### Baseline (with inclusion)

| Nx | penalty | L2(Vx) | L1(Vx) | L2(Vz) | L1(Vz) | Runtime |
|----|---------|--------|--------|--------|--------|---------|
| 21 | 1e5 | 2.03 | 5.00e-1 | 1.98 | 5.00e-1 | 0.06s |

**BUG FOUND**: Analytical formula had wrong sign! MDOODZ pure_shear_ALE with positive
bkg_strain_rate = shortening in x. Correct: Vx = -ε̇·x, Vz = +ε̇·z (test had both flipped).
Also: Vx array has Nx×(Nz+1) entries, Vz has (Nx+1)×Nz — test used wrong dimensions.

### After fix (corrected analytical formula)

| Nx | inclusion | L2(Vx) | L1(Vx) | L2(Vz) | L1(Vz) | Runtime |
|----|-----------|--------|--------|--------|--------|---------|
| 21 | no | 2.04e-8 | 4.47e-9 | 2.04e-8 | 4.47e-9 | 43ms |
| 21 | yes (η=10:1) | 2.43e-2 | 4.49e-3 | 2.43e-2 | 4.49e-3 | 77ms |
| 41 | yes | 2.55e-2 | 4.87e-3 | 2.55e-2 | 4.87e-3 | 150ms |
| 61 | yes | 2.39e-2 | 4.52e-3 | 2.39e-2 | 4.52e-3 | 286ms |

### Resolution Experiments

See table above — Nx=21/41/61 with inclusion.
Error is flat at ~0.024 (dominated by inclusion perturbation, not spatial discretization).

### Without Inclusion (user1=0)

See table above — L2=2.04e-8 (machine precision).

### Penalty Experiments

_Skipped_ — error is dominated by inclusion perturbation, not solver parameters.

### Summary & Recommendation

**Root cause**: Analytical formula had both signs wrong AND Vx/Vz array dimensions
were incorrect (read 420 of 462 entries, out-of-bounds access). L2=2.0 was
100% due to sign flip, not actual numerical error.

**Fix**: Corrected signs (Vx=-ε̇·x, Vz=+ε̇·z) and array dimensions.
L2 drops from 2.0 to 0.024 (with inclusion) or 2e-8 (without).
Threshold tightened from 3.0 to 0.05.
No parameter file changes needed.

**Recommended CI parameters**: unchanged (Nx=21, with inclusion). Threshold L2 < 0.05.

---

## 5. ViscousDissipation

**Analytical**: ΔT = 4η·ε̇_II²·t/(ρ·Cp)  *(corrected: MDOODZ Wdiss = τ_II²/η = 4ηε̇²)*
**File**: ShearHeatingTests.cpp, ShearHeating/ViscousDissipation.txt

### Baseline

| Nx | Nt | dt | dT_num | dT_ana | Abs error | Rel error | Runtime |
|----|----|----|--------|--------|-----------|-----------|---------|
| 31 | 3 | 1e11 | 0.321 K | 0.173 K | 0.148 K | 85% | 0.14s |

**TWO BUGS FOUND**:
1. Analytical formula used 2ηε̇² but MDOODZ computes Wdiss = τ_II²/η = 4ηε̇_II² (factor 2 error)
2. Dirichlet T BCs on top/bottom caused ~7% heat leakage

### After fix (4ηε̇² formula + Neumann BCs)

| Nx | Nt | dt | BCs | dT_num | dT_ana | Rel error | Runtime |
|----|----|----|-----|--------|--------|-----------|---------|
| 31 | 3 | 1e11 | Neumann | 0.349 K | 0.346 K | 0.8% | 0.13s |
| 31 | 3 | 1e12 | Neumann | 3.490 K | 3.463 K | 0.8% | 0.13s |

### Time Stepping Experiments

See table above — error is stable at 0.8% across different total times (dt=1e11 and 1e12).

### Resolution Experiments

_Skipped_ — error is ~0.8% of analytical with Neumann BCs, already excellent.

### BC Sensitivity Experiments

| Nx | Nt | T BCs | dT_num | dT_ana | Rel error | Notes |
|----|----|----|--------|--------|-----------|---------|
| 31 | 3 | Dirichlet (N/S) | 0.321 | 0.346 | -7.2% | heat leaks through fixed-T boundaries |
| 31 | 3 | Neumann (all) | 0.349 | 0.346 | +0.8% | no heat loss, excellent match |

### Summary & Recommendation

**Root causes** (two bugs):
1. Analytical formula: MDOODZ Wdiss = τ_II²/η = (2ηε̇_II)²/η = 4ηε̇_II², not 2ηε̇_II².
2. Dirichlet T BCs at top/bottom drained ~7% of the shear heating signal.

**Fix**: Corrected factor (2→4 in dT_ana), switched T BCs to all-Neumann.
Error drops from 85% to 0.8%. Tolerance tightened from ±200% to ±5%.

**Recommended CI parameters**: unchanged (Nx=31, Nt=3, dt=1e11), Neumann T BCs.
