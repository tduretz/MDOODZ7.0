## MODIFIED Requirements

### Requirement: HydrostaticPressure threshold tightened based on experiments
The HydrostaticPressure test threshold SHALL be updated to reflect measured accuracy after parameter optimization (resolution, penalty, solver tolerance).

#### Scenario: Threshold reflects measured error
- **WHEN** the HydrostaticPressure test runs with optimized parameters
- **THEN** the L2(P) threshold SHALL be no more than 2× the measured L2 error
