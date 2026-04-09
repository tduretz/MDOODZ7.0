## ADDED Requirements

### Requirement: Numeric regression assertions in visual tests
The visual test executable SHALL compute relative L2 norm errors between simulation output and reference HDF5 fields and return a nonzero exit code if any assertion fails.

#### Scenario: L2 norm check detects regression
- **WHEN** a visual test scenario runs and the computed ε̇_II field differs from the reference by more than the configured tolerance
- **THEN** the test SHALL print a diagnostic message identifying the failed field and tolerance, and the process exit code SHALL be nonzero

#### Scenario: L2 norm check passes for correct output
- **WHEN** a visual test scenario runs and produces output matching the reference within tolerance
- **THEN** the test SHALL print a success message and the process exit code SHALL be zero

### Requirement: Scalar quantity assertions in visual tests
The visual test executable SHALL verify key scalar quantities (peak stress, strain onset, mean temperature) against expected values.

#### Scenario: Peak stress assertion
- **WHEN** the VEP visual test completes
- **THEN** the peak average stress SHALL be within ±5% of the expected value from the reference paper (approximately 2.5 MPa for VEP_Duretz18)

### Requirement: CI workflow job for visual tests
The GitHub Actions CI workflow SHALL include an optional job that builds with `VIS=ON`, runs the visual test executable, and reports pass/fail.

#### Scenario: Visual tests run on push to main
- **WHEN** a commit is pushed to the `main` branch
- **THEN** the CI SHALL trigger the visual test job and report its exit code as the job status

#### Scenario: Visual tests are opt-in for pull requests
- **WHEN** a pull request is opened or updated without the `run-visual-tests` label
- **THEN** the visual test job SHALL be skipped

#### Scenario: Visual tests run when labeled
- **WHEN** a pull request has the `run-visual-tests` label applied
- **THEN** the CI SHALL trigger the visual test job for that PR

### Requirement: Reference regeneration target
The build system SHALL provide a `make generate-refs` target that runs all visual test scenarios and copies their output HDF5 files to the reference locations.

#### Scenario: Generate references for all visual tests
- **WHEN** `make generate-refs` is executed after building with `VIS=ON`
- **THEN** the reference HDF5 files for all visual test scenarios SHALL be updated in the `VISUAL_TESTS/` directory
