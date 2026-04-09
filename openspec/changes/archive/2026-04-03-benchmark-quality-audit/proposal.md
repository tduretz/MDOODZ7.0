## Why

Several non-SolVi analytical benchmark tests have sub-optimal error thresholds and accuracy. The three loosest are: PureShearVelocity (L2 < 3.0), ViscousDissipation (±200% tolerance), and SteadyStateGeotherm (L2 < 1.0). These tests would pass even with grossly wrong results, undermining their value as regression guards.

The same three questions proven useful in the SolVi L1 work apply here: (1) does L1 norm give a more representative error than L2? (2) does increasing resolution improve errors? (3) can parameter tuning (more time steps, higher penalty, tighter tolerances) improve accuracy?

This is an **experimental change** — we run many parameter sweeps, document every result, iterate freely, and only commit to threshold/parameter changes after patterns emerge. Results are documented per-test in a dedicated experiment log (`TESTS/BenchmarkExperiments.md`) and compacted findings go into skills.

## What Changes

For each of the 5 tests with analytical comparisons (SteadyStateGeotherm, RadiogenicHeat, HydrostaticPressure, PureShearVelocity, ViscousDissipation), systematically:

- **Add L1 error reporting** alongside existing L2 where applicable
- **Run parameter experiments** to find if simple changes improve accuracy:
  - Resolution: try 41×41 and 61×61 alongside current 21–31
  - SteadyStateGeotherm: increase `Nt` (currently 20) or `dt` to reach true steady state
  - ViscousDissipation: increase `Nt` (currently 3), check thermal BC sensitivity
  - PureShearVelocity: increase resolution (currently 21×21), note this test has an inclusion that perturbs the field
  - HydrostaticPressure: increase resolution (currently 21×21)
  - Any other .txt parameter that might matter (penalty, solver tolerances, marker density)
- **Explore physics parameter sensitivity**:
  - **Penalty parameter** (`penalty = 1e5` everywhere): higher values enforce incompressibility more tightly, directly affects pressure accuracy. Try 1e6, 1e7
  - **Viscosity contrast** (PureShearVelocity has η_c/η_m = 1000 inclusion): does reducing contrast improve velocity L2?
  - **Thermal conductivity** `k` and **heat capacity** `Cp`: affect how fast SteadyStateGeotherm reaches equilibrium and ViscousDissipation temperature signal
  - **Time step `dt`**: too large may cause under-resolved thermal evolution; too small wastes steps. Both SteadyStateGeotherm (dt=1e12) and ViscousDissipation (dt=1e11) may benefit from tuning
  - **Gravity `gz`**: HydrostaticPressure uses gz=-10; verify the analytical formula matches
  - **Scaling parameters** (`eta`, `L`, `V`, `T`): non-dimensionalisation affects condition number and solver accuracy
  - **`rel_tol_KSP = 1e-4`**: loose linear solver tolerance may limit accuracy, especially for pressure. Try 1e-6
- **Document every experiment** in `TESTS/BenchmarkExperiments.md` — raw numbers, parameter combinations, what worked and what didn't
- **Tighten thresholds** only after patterns are clear — target ~2× measured error as the assertion bound
- **Compact findings into skills** — what parameters matter, what doesn't, L1 vs L2 recommendations

### Experimental Methodology

Each test goes through this cycle (potentially multiple iterations):

1. **Baseline**: Run with current parameters, record L2 and L1 errors
2. **Resolution sweep**: Re-run at 2–3 resolutions, record errors, compute convergence orders if applicable
3. **Parameter sweep**: Vary one parameter at a time (Nt, dt, penalty, etc.), record errors
4. **Document**: Write all raw numbers to `TESTS/BenchmarkExperiments.md`
5. **Decide**: Based on measured data, choose parameter updates and threshold tightening
6. **Implement**: Update .txt files and test assertions
7. **Verify**: Re-run to confirm

There is no fixed number of iterations — we experiment until each test has optimal (or at least well-understood) accuracy.

## Capabilities

### New Capabilities
- `benchmark-quality-audit`: The experimental investigation — parameter sweeps, L1 vs L2 comparison, threshold tightening for all 5 analytical tests. Includes the experiment log `TESTS/BenchmarkExperiments.md`.

### Modified Capabilities
- `ci-thermal-tests`: Threshold updates for SteadyStateGeotherm and RadiogenicHeat if parameters change
- `ci-density-tests`: Threshold updates for HydrostaticPressure if parameters change
- `ci-velocity-field-tests`: Threshold updates for PureShearVelocity if parameters change
- `ci-shear-heating-tests`: Threshold updates for ViscousDissipation if parameters change
- `skill-testing-guide`: Document compacted findings — which parameters matter, what doesn't, L1 vs L2 recommendations, per-test accuracy expectations

## Impact

- **Test files**: `ThermalTests.cpp`, `DensityTests.cpp`, `VelocityFieldTests.cpp`, `ShearHeatingTests.cpp` — add L1 reporting, tighten thresholds
- **Parameter files**: Potentially update `.txt` files under `TESTS/Thermal/`, `TESTS/Density/`, `TESTS/VelocityField/`, `TESTS/ShearHeating/` (resolution, Nt, dt)
- **TestHelpers.h**: `computeL1Error()` already available from the SolVi L1 change
- **New file**: `TESTS/BenchmarkExperiments.md` — detailed experiment log with raw numbers from every run
- **Documentation**: `TESTS/AnalyticalSolutions.md` (update measured errors), `skill-testing-guide/SKILL.md` (compacted patterns)
- **Runtime**: Higher resolutions or more time steps will increase CI time; budget ~30s total for all 5 tests
- **No breaking changes**: only tightening existing thresholds and adding informational output
