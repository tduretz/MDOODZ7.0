## Context

Five non-SolVi analytical benchmarks exist in the CI test suite, each comparing numerical output against a known formula. Current thresholds were set conservatively when the tests were first written — some allow errors of 200–300%. The tests run at low resolution (21–31 cells) with default parameters that may not be optimal for accuracy.

The SolVi L1 investigation showed that norm choice, resolution, and parameter tuning can dramatically change measured convergence quality. This design applies the same experimental approach to the remaining benchmarks.

Current test inventory with analytical comparisons:

| Test | File | Current threshold | Resolution | Physics |
|------|------|------------------|------------|---------|
| SteadyStateGeotherm | ThermalTests.cpp | L2(T) < 1.0 | 31×31, Nt=20 | Linear geotherm T(z) |
| RadiogenicHeat | ThermalTests.cpp | ±50% of ΔT | 31×31, Nt=5 | T_mean = T0 + Qr·t/(ρ·Cp) |
| HydrostaticPressure | DensityTests.cpp | L2(P) < 0.2 | 21×21 | P(z) = ρ·g·z |
| PureShearVelocity | VelocityFieldTests.cpp | L2(Vx,Vz) < 3.0 | 21×21 | Vx = ε̇·x (with inclusion) |
| ViscousDissipation | ShearHeatingTests.cpp | ±200% of dT_ana | 31×31, Nt=3 | ΔT = 2η·ε̇²·t/(ρ·Cp) |

## Goals / Non-Goals

**Goals:**
- Measure baseline L2 and L1 errors for all 5 tests with current parameters
- Systematically explore which parameters improve accuracy (resolution, Nt/dt, penalty, rel_tol_KSP, physics parameters)
- Document every experiment result in `TESTS/BenchmarkExperiments.md`
- Tighten thresholds based on measured data (target: threshold ≈ 2× measured error)
- Update parameter files where improvements are clear and CI runtime stays reasonable
- Compact findings into skill-testing-guide

**Non-Goals:**
- Modifying the solver or stencil code (C source files in MDLIB/)
- Adding new analytical benchmarks
- Achieving "perfect" errors — the goal is understanding what's achievable and setting honest thresholds
- Convergence-order analysis for these tests (most are single-resolution; SolVi already covers this)

## Decisions

### 1. Experiment-first workflow: measure before changing

**Decision**: Run every test at current parameters first (baseline), then vary one parameter at a time, recording results in `TESTS/BenchmarkExperiments.md`. Only update thresholds and .txt files after patterns are clear across multiple experiments.

**Rationale**: Previous SolVi work showed assumptions can be wrong (cell-face averaging seemed promising, then failed). Measuring first avoids wasted effort. The experiment log provides a permanent record of what was tried.

### 2. One experiment log file, structured by test

**Decision**: All experiment results go into `TESTS/BenchmarkExperiments.md`, organized by test name with tables for each parameter sweep. After all experiments, compacted findings go into `skill-testing-guide/SKILL.md`.

**Rationale**: One file is easier to maintain than five. The test-by-test structure makes it searchable. Skills get only the distilled patterns, not raw data.

### 3. Parameter exploration priority order

**Decision**: For each test, explore parameters in this order:
1. **Baseline** (current params, add L1 reporting)
2. **Resolution** (double grid, e.g. 21→41 or 31→61)
3. **Time integration** (increase Nt or dt for transient tests)
4. **Penalty** (1e5 → 1e6 → 1e7 for Stokes-coupled tests)
5. **Linear solver tolerance** (rel_tol_KSP: 1e-4 → 1e-6)
6. **Physics parameters** (k, Cp, scaling, viscosity contrast)

**Rationale**: Resolution and time stepping are the most likely to help and are cheapest to test. Physics parameters are explored last because they change the problem definition, not just the numerical accuracy.

### 4. Test-specific experiment plans

**SteadyStateGeotherm**: The analytical solution is the steady-state linear geotherm. The test runs 20 time steps with dt=1e12. If it hasn't reached steady state, the error is time-integration inadequacy, not spatial discretisation. Experiment: increase Nt to 50, 100. Also try larger dt. Also try turning off mechanical solver (if possible) to isolate thermal solve.

**RadiogenicHeat**: Analytical T_mean = T0 + Qr·t/(ρ·Cp). The 50% tolerance may be because diffusion spreads heat to boundaries where it's lost. Experiment: check if insulating (zero-flux) BCs on all sides gives better agreement. Also increase Nt.

**HydrostaticPressure**: P(z) = ρ·g·z is a static problem — should be exact up to discretisation error. L2 < 0.2 at 21×21 suggests ~15–20% error, which is high for a linear profile. Experiment: raise resolution to 41, 61. Try higher penalty (pressure accuracy is directly penalty-dependent). Try rel_tol_KSP=1e-6.

**PureShearVelocity**: Has a viscosity inclusion (η_c=1000× in phase 1) that perturbs the pure shear field. The L2 < 3.0 threshold is measuring the inclusion perturbation, not the background field accuracy. Experiment: (a) raise resolution to 41, 61; (b) try without inclusion (user1=0) to measure pure shear accuracy; (c) compare L1 vs L2 — L1 may be less dominated by the near-inclusion error.

**ViscousDissipation**: ΔT = 2η·ε̇²·t/(ρ·Cp) with shear_heating=1. The ±200% threshold suggests the signal is tiny or thermal BCs are leaking heat. Experiment: increase Nt from 3 to 10, 20. Check if Dirichlet T BCs on N/S are draining the heat signal. Try zero-flux BCs. Try higher resolution.

### 5. CI runtime budget

**Decision**: Updated parameter files should keep each individual test under ~30 seconds and total runtime for all 5 tests under ~2 minutes. There is no hard CI time limit, but tests should not take several minutes each.

**Rationale**: These are CI tests that run on every commit. The experiment log captures results at any resolution, but the committed .txt files should be reasonably fast. Higher resolution or more time steps are acceptable if they meaningfully improve accuracy.

**Alternative considered**: Moving slow configurations to an experimental test file (like SolViMarkerComparison.cpp). Rejected for now — keep it simple unless needed.

### 6. L1 reporting: add to all tests with pointwise comparisons

**Decision**: Add `computeL1Error()` calls alongside existing L2 in every test that builds an analytical array (SteadyStateGeotherm, HydrostaticPressure, PureShearVelocity). For tests using scalar comparisons (RadiogenicHeat, ViscousDissipation, StressAccumulation), L1 vs L2 doesn't apply — these compare single values.

**Rationale**: `computeL1Error()` is already in TestHelpers.h. Adding it to field-comparison tests costs one extra line and provides a direct comparison with L2 in the experiment log.

## Risks / Trade-offs

**[Risk] Some tests may already be at their accuracy limit with current methods**
→ Mitigation: Document this as a finding rather than forcing an improvement. An honest "this test achieves L2=0.5 at 31×31 and that's the best we can do" is a valid outcome.

**[Risk] Parameter changes may cause existing test assertions to fail**
→ Mitigation: Experiments run separately (manual, not CI). Only update .txt files and thresholds together in a coordinated commit.

**[Risk] Experiment log may grow large**
→ Mitigation: Structure with clear headers and tables. Each section is self-contained. Compacted findings go to skills; the log is reference material.

**[Risk] Some "improvements" may be parameter-specific and not generalisable**
→ Mitigation: Record the full parameter set for each experiment. If a finding only works with specific values, note that explicitly.
