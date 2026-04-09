## ADDED Requirements

### Requirement: CI test authoring guide
The skill SHALL document how to create a new GTest-based CI test including the fixture pattern, parameter file conventions, assertion types, and CMake integration.

#### Scenario: User asks how to add a CI test
- **WHEN** the user asks how to write a new CI test for MDOODZ
- **THEN** the skill SHALL provide the complete workflow: create `.cpp` with GTest fixture, create `.txt` parameter file, add to CMakeLists.txt, and show assertion patterns (iteration count, field bounds, residual checks)

### Requirement: Visual test authoring guide
The skill SHALL document how to create a new visual regression test including the C++ class pattern, reference generation, Eigen field extraction, L2 comparison, and gnuplot integration.

#### Scenario: User asks how to add a visual test
- **WHEN** the user asks how to write a new visual regression test
- **THEN** the skill SHALL provide the workflow: create test class in `.cpp`, create `.txt` config, generate reference HDF5, add gnuplot script, register in `main.cpp`, and show the `relativeL2Error` comparison pattern

### Requirement: Test execution guide
The skill SHALL document how to build and run both CI tests and visual tests locally and in CI.

#### Scenario: User asks how to run CI tests
- **WHEN** the user asks how to run the fast CI tests
- **THEN** the skill SHALL show the commands: `make build-dev TEST=ON` then `make run-tests` (which runs `ctest`)

#### Scenario: User asks how to run visual tests
- **WHEN** the user asks how to run the visual regression tests
- **THEN** the skill SHALL show the commands: `make build-dev VIS=ON` then `make run-vis`

### Requirement: Reference update guide
The skill SHALL document how to regenerate reference HDF5 files when intentional code changes alter simulation output.

#### Scenario: User needs to update references after a code change
- **WHEN** the user asks how to update visual test references
- **THEN** the skill SHALL explain: run `make generate-refs`, verify the new output is correct, then commit the updated `.h5` files

### Requirement: Test coverage map
The skill SHALL include a table mapping tested modules/features to their corresponding test files for quick navigation.

#### Scenario: User asks what is tested
- **WHEN** the user asks which features have test coverage
- **THEN** the skill SHALL provide a table with columns: Feature, CI Test File, Visual Test, and coverage status
