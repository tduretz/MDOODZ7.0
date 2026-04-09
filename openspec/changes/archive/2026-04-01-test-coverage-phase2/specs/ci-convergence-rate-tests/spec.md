## ADDED Requirements

### Requirement: Newton converges faster than Picard
The test suite SHALL verify that Newton-Raphson iteration requires fewer iterations than Picard for the same non-linear problem.

#### Scenario: Newton uses fewer iterations
- **WHEN** the same power-law inclusion problem is solved with `Newton=1` and separately with `Newton=0` (Picard)
- **THEN** the Newton iteration count SHALL be less than or equal to the Picard iteration count

### Requirement: LogIsNewtonStep flag is consistent
The test suite SHALL verify that the `Iterations/LogIsNewtonStep` field correctly marks which iterations used Newton vs Picard.

#### Scenario: Newton steps are flagged
- **WHEN** `Newton=1` is set and the solver converges
- **THEN** at least one entry in `LogIsNewtonStep` SHALL be non-zero (indicating Newton steps were taken)
