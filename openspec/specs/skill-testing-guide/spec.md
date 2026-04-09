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

## MODIFIED Section: SolVi Benchmark Tests

### Requirement: Document L1 norm findings
The SolVi section of the testing guide skill SHALL document that L1 norm (`mean(|num - ana|)`) was added alongside L2 to better measure pressure convergence. It SHALL record the measured L1 convergence orders at high resolution and note how they compare to L2 orders.

#### Scenario: Skills reflects current test inventory
- **WHEN** a developer reads skill-testing-guide
- **THEN** the SolVi benchmark section SHALL mention both L2 and L1 norms and list the HighResL1Convergence test alongside existing CI tests

## ADDED Section: Analytical Benchmark Authoring

### Requirement: Analytical benchmark authoring guide
The skill SHALL document how to write a CI test that compares numerical output against an analytical solution, including: reading fields with `readFieldAsArray`, reading coordinates with `readCoordArray`, computing L2 error with `computeL2Error`, setting empirical thresholds, and documenting the analytical solution in `TESTS/AnalyticalSolutions.md`.

#### Scenario: User asks how to add an analytical benchmark
- **WHEN** the user asks how to add an L2 error benchmark to an existing or new CI test
- **THEN** the skill SHALL provide the workflow: use `readFieldAsArray()` + `readCoordArray()` to get numerical data and grid coordinates, evaluate the analytical formula at each grid point, call `computeL2Error()`, assert `EXPECT_LT(L2, threshold)`, and document the analytical solution and threshold in `AnalyticalSolutions.md`

### Requirement: Grid-convergence test pattern
The skill SHALL document the pattern for writing a grid-convergence order test: running the same problem at multiple resolutions, computing L2 at each, and asserting the convergence order p = log(L2_coarse/L2_fine) / log(h_coarse/h_fine) ≥ expected_order.

#### Scenario: User asks how to verify spatial convergence order
- **WHEN** the user asks how to write a convergence-order test
- **THEN** the skill SHALL explain: create parameter files at 3+ resolutions, run each, compute L2 error at each resolution, compute the convergence rate between successive pairs, and assert the rate is ≥ expected for the method order

### Requirement: L2 threshold calibration guide
The skill SHALL explain how to set and maintain L2 error thresholds: run the test, measure the actual L2, set threshold to measured × 3, document the baseline in `AnalyticalSolutions.md`, and re-calibrate after intentional solver changes.

#### Scenario: User asks how to set an L2 threshold
- **WHEN** the user asks what threshold to use for an L2 assertion
- **THEN** the skill SHALL advise: run the test empirically, record the measured L2, set the assertion threshold to 2–5× the measured value, and document both values

### Requirement: Analytical benchmark pitfalls
The skill SHALL document common pitfalls specific to analytical benchmarks: coordinate system mismatches, staggered grid field sizes, SI unit scaling, marker-in-cell noise on convergence order, and handling zero analytical solutions (absolute vs relative L2).

#### Scenario: User encounters L2 error larger than expected
- **WHEN** the user asks why their L2 error is unexpectedly large
- **THEN** the skill SHALL list common causes: wrong sign convention for pure shear (Vx=-ε̇·x not +ε̇·x), wrong dissipation factor (Wdiss=4ηε̇² not 2ηε̇²), insufficient time stepping for steady-state tests, Dirichlet BCs causing heat leakage, and wrong staggered grid dimensions (Vx is Nx×(Nz+1), Vz is (Nx+1)×Nz)

## ADDED Section: Non-SolVi Benchmark Accuracy Findings

### Requirement: Document per-test accuracy expectations and parameter recommendations
The skill SHALL include a section documenting the measured accuracy for each analytical benchmark test, which parameters matter most, and recommended settings.

#### Scenario: Developer reads about benchmark accuracy
- **WHEN** a developer reads skill-testing-guide for benchmark information
- **THEN** the skill SHALL list for each test: measured error, recommended .txt parameters, and the key finding

### Requirement: Document L1 vs L2 findings across non-SolVi tests
The skill SHALL document whether L1 norm provides meaningfully different accuracy measurements compared to L2 for the non-SolVi benchmarks.

#### Scenario: Developer asks about norm choice
- **WHEN** a developer asks which norm to use for a new benchmark
- **THEN** the skill SHALL provide guidance summarizing L1 vs L2 findings from both SolVi and non-SolVi experiments
