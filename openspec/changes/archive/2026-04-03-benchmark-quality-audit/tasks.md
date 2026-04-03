## 0. Setup

- [x] 0.1 Add L1 error reporting to SteadyStateGeotherm test (print L1 alongside L2)
- [x] 0.2 Add L1 error reporting to HydrostaticPressure test (print L1 alongside L2)
- [x] 0.3 Add L1 error reporting to PureShearVelocity test (print L1 alongside L2 for Vx and Vz)
- [x] 0.4 Create `TESTS/BenchmarkExperiments.md` with header and structure for all 5 tests
- [x] 0.5 Sync to WSL, build, run all 5 tests at current parameters, record baseline results in BenchmarkExperiments.md

## 1. SteadyStateGeotherm Experiments

- [x] 1.1 Baseline: record L2(T), L1(T) at Nx=31, Nt=20, dt=1e12
- [x] 1.2 Time stepping: try Nt=50, Nt=100 — record L2(T), L1(T) for each
- [x] 1.3 Time stepping: try larger dt (1e13, 1e14, 1e15, 5e15, 1e16) — record results
- [x] 1.4 Resolution: try Nx=61 at best dt — record results
- [x] 1.5 Penalty: try penalty=1e6, 1e7 — record results
- [x] 1.6 Summarise findings and choose optimal CI parameters
- [x] 1.7 Update SteadyStateGeotherm.txt and tighten L2 threshold in ThermalTests.cpp

## 2. RadiogenicHeat Experiments

- [x] 2.1 Baseline: record meanT_final, T_ana, absolute error, relative error at Nt=5
- [x] 2.2 Time stepping: try dt=1e13 — error scales linearly with time (boundary losses)
- [x] 2.3 Resolution: skipped — error already 0.008%, boundary-driven not resolution-driven
- [x] 2.4 BC sensitivity: error scaling confirms boundary heat leakage is the limiting factor
- [x] 2.5 Summarise findings and choose optimal CI parameters
- [x] 2.6 Update tolerance in ThermalTests.cpp from ±50% to ±2%

## 3. HydrostaticPressure Experiments

- [x] 3.1 Baseline: record L2(P), L1(P) at Nx=21
- [x] 3.2 Resolution: Nx=21,41,61 — first-order convergence (L2∝1/Nx)
- [x] 3.3 Penalty: 1e5,1e6,1e7 — zero effect
- [x] 3.4 Linear solver tolerance: 1e-4,1e-6 — zero effect
- [x] 3.5 Summarise findings: purely spatial discretization error, first-order
- [x] 3.6 Update HydrostaticPressure.txt (Nx=41) and tighten threshold (L2<0.06)

## 4. PureShearVelocity Experiments

- [x] 4.1 Baseline: L2(Vx)=2.03, L2(Vz)=1.98 — revealed sign bug in analytical formula
- [x] 4.2 Resolution: Nx=21,41,61 — error flat at 0.024 (inclusion-dominated)
- [x] 4.3 Without inclusion: L2=2e-8 (machine precision) — confirms bug was in formula
- [x] 4.4 Penalty: skipped — error dominated by inclusion perturbation
- [x] 4.5 Fixed sign (Vx=-ε̇·x, Vz=+ε̇·z) and array dimensions (Nx×(Nz+1), (Nx+1)×Nz)
- [x] 4.6 Summarise: sign bug caused 100% of the L2=2.0. True L2=0.024 with inclusion.
- [x] 4.7 Tightened threshold from 3.0 to 0.05 in VelocityFieldTests.cpp

## 5. ViscousDissipation Experiments

- [x] 5.1 Baseline: dT_num=0.321, dT_ana=0.173, 85% error — found TWO bugs
- [x] 5.2 Bug 1: analytical formula 2ηε̇² → 4ηε̇² (MDOODZ Wdiss = τ_II²/η)
- [x] 5.3 Bug 2: Dirichlet T BCs drain 7% heat → switched to all-Neumann
- [x] 5.4 After fix: dT_num=0.349, dT_ana=0.346, error=0.8%
- [x] 5.5 Verified stable across dt=1e11 and dt=1e12
- [x] 5.6 Summarise: two bugs (formula factor + BCs), error 85% → 0.8%
- [x] 5.7 Tightened tolerance from ±200% to ±5% in ShearHeatingTests.cpp

## 6. Documentation

- [x] 6.1 Update `TESTS/AnalyticalSolutions.md` with final measured errors for all 5 tests
- [x] 6.2 Update `skill-testing-guide/SKILL.md` with compacted findings: per-test accuracy, parameter recommendations, L1 vs L2 guidance
- [x] 6.3 Final build and run: verify all 5 tests pass with updated thresholds
