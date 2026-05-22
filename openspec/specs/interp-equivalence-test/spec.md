## ADDED Requirements

### Requirement: Field-level equivalence test between interp_mode 0 and 3
The system SHALL include a GTest case that runs the same scenario twice — once with `interp_mode = 0` and once with `interp_mode = 3` — and compares the resulting T, Vx, and Vz fields from HDF5 output.

#### Scenario: Temperature fields match
- **WHEN** both runs complete on a 51×51 grid for 5 steps
- **THEN** the max absolute difference between the T fields at the final step SHALL be less than 1e-6

#### Scenario: Velocity fields match
- **WHEN** both runs complete
- **THEN** the max absolute difference between the Vx fields and Vz fields at the final step SHALL each be less than 1e-6

#### Scenario: Test uses a thermally-coupled scenario
- **WHEN** the equivalence test runs
- **THEN** it SHALL use a scenario with `thermal = 1` and non-trivial temperature structure (e.g. SteadyStateGeotherm or GaussianDiffusion) so that interp differences in T, if any, are exercised

### Requirement: Separate parameter files for each interp mode
The system SHALL provide two `.txt` parameter files that are identical except for `interp_mode` and `writer_subfolder`.

#### Scenario: Mode 0 parameter file
- **WHEN** the interp_mode=0 run executes
- **THEN** it SHALL read a `.txt` file with `interp_mode = 0` and `writer_subfolder = InterpEquiv_Mode0`

#### Scenario: Mode 3 parameter file
- **WHEN** the interp_mode=3 run executes
- **THEN** it SHALL read a `.txt` file with `interp_mode = 3` and `writer_subfolder = InterpEquiv_Mode3`
