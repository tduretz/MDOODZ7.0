## MODIFIED Requirements

### Requirement: Nusselt number validation
The CI test suite SHALL include a GTest case that runs the Blankenbach Case 1a simulation (Ra = 10⁴, isoviscous) to steady state and validates the top-surface Nusselt number against the published value Nu = 4.884. The steady-state test SHALL be tuned to converge in roughly 1000 steps (soft cap) by adjusting `Courant`, `dt_max`, and grid resolution. The test SHALL use `.dat` breakpoint restart if it does not converge on the first attempt.

#### Scenario: Nusselt number within tolerance
- **WHEN** the simulation reaches steady state and the Nusselt number is computed from the top-boundary temperature gradient via finite differences
- **THEN** the relative error between the computed Nu and the published value 4.884 SHALL be less than 5%

#### Scenario: Steady-state test does not crash from particle overflow
- **WHEN** the steady-state test runs to completion with particle deactivation enabled
- **THEN** the simulation SHALL NOT call `exit(190)` or `exit(1)` due to particle count overflow

## ADDED Requirements

### Requirement: Particle count stability assertion in CI smoke test
The BlankenBench CI smoke test SHALL assert that `Nb_part` at the final time step is not significantly larger than the initial particle count. The assertion SHALL use a 10% tolerance: `Nb_part_final <= 1.1 * Nb_part_initial`.

#### Scenario: Particle count bounded after N steps
- **WHEN** the CI smoke test completes its configured number of steps
- **THEN** the `Nb_part` read from the final HDF5 output SHALL be within 10% of the initial value

#### Scenario: Test fails without deactivation
- **WHEN** the deactivation logic is removed (regression) and the smoke test runs
- **THEN** the `Nb_part` stability assertion SHALL fail, catching the monotonic growth bug

### Requirement: Nb_part written to HDF5 output
The HDF5 output writer SHALL include `Nb_part` as a readable field (scalar attribute or dataset) in each output step file so that GTest cases can read and assert on it.

#### Scenario: Nb_part readable from output file
- **WHEN** a GTest reads the final HDF5 output file
- **THEN** it SHALL be able to extract `Nb_part` as an integer value

### Requirement: Steady-state test tuned for fast convergence
The `BlankenBenchSteady.txt` parameter file SHALL be tuned so thermal steady state is reached in roughly 1000 steps (soft cap). Settings to adjust include `Courant`, `dt_max`, and grid resolution. The `.dat` breakpoint restart mechanism SHALL be usable to continue the run if convergence is not reached.

#### Scenario: Steady state reached in ~1000 steps
- **WHEN** the steady-state test runs with the tuned parameters
- **THEN** the Nusselt number and Vrms SHALL converge to within 5% of the published Blankenbach values (Nu = 4.884, Vrms = 42.865) in approximately 1000 steps
