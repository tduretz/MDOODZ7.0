## MODIFIED Requirements

### Requirement: SteadyStateGeotherm threshold tightened based on experiments
The SteadyStateGeotherm test threshold SHALL be updated to reflect measured accuracy after parameter optimization. The specific threshold value depends on experiment results.

#### Scenario: Threshold reflects measured error
- **WHEN** the SteadyStateGeotherm test runs with optimized parameters
- **THEN** the L2(T) threshold SHALL be no more than 2× the measured L2 error

### Requirement: RadiogenicHeat threshold tightened based on experiments
The RadiogenicHeat test tolerance SHALL be updated to reflect measured accuracy after parameter optimization.

#### Scenario: Threshold reflects measured error
- **WHEN** the RadiogenicHeat test runs with optimized parameters
- **THEN** the tolerance on mean temperature SHALL be no more than 2× the measured absolute error
