## MODIFIED Requirements

### Requirement: PureShearVelocity threshold tightened based on experiments
The PureShearVelocity test threshold SHALL be updated to reflect measured accuracy after parameter optimization (resolution, penalty, inclusion handling).

#### Scenario: Threshold reflects measured error
- **WHEN** the PureShearVelocity test runs with optimized parameters
- **THEN** the L2(Vx) and L2(Vz) thresholds SHALL be no more than 2× the measured L2 error
