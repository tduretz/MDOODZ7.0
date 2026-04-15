## ADDED Requirements

### Requirement: Analytical benchmark authoring guide
The skill SHALL document how to write a CI test that compares numerical output against an analytical solution, including: reading fields with `readFieldAsArray`, reading coordinates with `readCoordArray`, computing L2 error with `computeL2Error`, setting empirical thresholds, and documenting the analytical solution in `TESTS/AnalyticalSolutions.md`.

#### Scenario: User asks how to add an analytical benchmark
- **WHEN** the user asks how to add an L2 error benchmark to an existing or new CI test
- **THEN** the skill SHALL provide the workflow: use `readFieldAsArray()` + `readCoordArray()` to get numerical data and grid coordinates, evaluate the analytical formula at each grid point, call `computeL2Error()`, assert `EXPECT_LT(L2, threshold)`, and document the analytical solution and threshold in `AnalyticalSolutions.md`

### Requirement: Grid-convergence test pattern
The skill SHALL document the pattern for writing a grid-convergence order test: running the same problem at multiple resolutions, computing L2 at each, and asserting the convergence order p = log(L2_coarse/L2_fine) / log(h_coarse/h_fine) ≥ expected_order.

#### Scenario: User asks how to verify spatial convergence order
- **WHEN** the user asks how to write a convergence-order test
- **THEN** the skill SHALL explain: create parameter files at 3+ resolutions, run each, compute L2 error at each resolution, compute the convergence rate between successive pairs, and assert the rate is ≥ 1.8 for second-order methods

### Requirement: L2 threshold calibration guide
The skill SHALL explain how to set and maintain L2 error thresholds: run the test, measure the actual L2, set threshold to measured × 3, document the baseline in `AnalyticalSolutions.md`, and re-calibrate after intentional solver changes.

#### Scenario: User asks how to set an L2 threshold
- **WHEN** the user asks what threshold to use for an L2 assertion
- **THEN** the skill SHALL advise: run the test empirically, record the measured L2, set the assertion threshold to 2–5× the measured value, and document both values

### Requirement: Analytical benchmark pitfalls
The skill SHALL document common pitfalls specific to analytical benchmarks: coordinate system mismatches, staggered grid field sizes, SI unit scaling, marker-in-cell noise on convergence order, and handling zero analytical solutions (absolute vs relative L2).

#### Scenario: User encounters L2 error larger than expected
- **WHEN** the user asks why their L2 error is unexpectedly large
- **THEN** the skill SHALL list common causes: wrong coordinate array for the field's grid type, SI vs non-dimensional mismatch, analytical formula evaluated at wrong positions, or field includes ghost/boundary cells
