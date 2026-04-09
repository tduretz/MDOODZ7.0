## ADDED Requirements

### Requirement: Nusselt number validation
The CI test suite SHALL include a GTest case that runs the Blankenbach Case 1a simulation (Ra = 10⁴, isoviscous) to steady state and validates the top-surface Nusselt number against the published value Nu = 4.884.

#### Scenario: Nusselt number within tolerance
- **WHEN** the simulation reaches steady state at 41×41 resolution and the Nusselt number is computed from the top-boundary temperature gradient via finite differences
- **THEN** the relative error between the computed Nu and the published value 4.884 SHALL be less than 5%

### Requirement: RMS velocity validation
The CI test suite SHALL include a GTest case that validates the domain-averaged RMS velocity against the published Blankenbach Case 1a value Vrms = 42.865.

#### Scenario: RMS velocity within tolerance
- **WHEN** the steady-state Vx and Vz fields are read from HDF5 output and Vrms = sqrt(mean(Vx²) + mean(Vz²)) is computed
- **THEN** the relative error between the computed Vrms and the published value 42.865 SHALL be less than 5%

### Requirement: Mid-depth maximum temperature validation
The CI test suite SHALL include a GTest case that validates the maximum temperature at mid-depth against the published Blankenbach value T_max = 0.4266.

#### Scenario: Mid-depth temperature within tolerance
- **WHEN** the steady-state temperature field is read and the maximum temperature in the horizontal row closest to z = 0 (mid-depth) is extracted
- **THEN** the relative error between the computed T_max and 0.4266 SHALL be less than 5%

### Requirement: Test registered in CMakeLists
The BlankenBenchTests executable SHALL be registered in `TESTS/CMakeLists.txt` and linked against mdoodz, GTest, and HDF5.

#### Scenario: Test builds and runs
- **WHEN** cmake builds the BlankenBenchTests target
- **THEN** the executable SHALL compile without errors and be runnable from the TESTS/ directory
