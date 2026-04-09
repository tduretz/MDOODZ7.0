## MODIFIED Requirements

### Requirement: ViscousDissipation threshold tightened based on experiments
The ViscousDissipation test tolerance SHALL be updated to reflect measured accuracy after parameter optimization (Nt, dt, BCs, resolution).

#### Scenario: Threshold reflects measured error
- **WHEN** the ViscousDissipation test runs with optimized parameters
- **THEN** the dT tolerance SHALL be no more than 2× the measured absolute error
