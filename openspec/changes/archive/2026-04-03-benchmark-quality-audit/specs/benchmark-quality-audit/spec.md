## ADDED Requirements

### Requirement: Experiment log for benchmark quality investigation
The system SHALL include a file `TESTS/BenchmarkExperiments.md` that documents every parameter experiment run during this audit, organized by test name.

#### Scenario: Baseline recording
- **WHEN** an experiment baseline is run for any of the 5 analytical tests
- **THEN** the log SHALL record: test name, parameters used, L2 error, L1 error (if applicable), and runtime

#### Scenario: Parameter sweep recording
- **WHEN** a parameter is varied (resolution, Nt, dt, penalty, rel_tol_KSP, physics params)
- **THEN** the log SHALL record: parameter name, values tested, resulting errors for each value, and whether an improvement was observed

#### Scenario: Finding summary per test
- **WHEN** all experiments for a test are complete
- **THEN** the log SHALL include a summary section recording: best achievable error, recommended parameters for CI, recommended threshold, and key insights

### Requirement: L1 error reporting for field-comparison tests
The tests SteadyStateGeotherm, HydrostaticPressure, and PureShearVelocity SHALL compute and print `computeL1Error()` alongside existing L2 output. No L1 assertions are required — this is informational for the experiment log.

#### Scenario: L1 printed for SteadyStateGeotherm
- **WHEN** SteadyStateGeotherm runs
- **THEN** it SHALL print both L2 and L1 error for the temperature field

#### Scenario: L1 printed for HydrostaticPressure
- **WHEN** HydrostaticPressure runs
- **THEN** it SHALL print both L2 and L1 error for the pressure field

#### Scenario: L1 printed for PureShearVelocity
- **WHEN** PureShearVelocity runs
- **THEN** it SHALL print both L2 and L1 error for Vx and Vz fields
