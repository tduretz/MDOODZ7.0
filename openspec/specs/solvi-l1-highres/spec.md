## ADDED Requirements

### Requirement: L1 error norm helper function
The system SHALL provide a `computeL1Error()` function in `TESTS/TestHelpers.h` that computes the absolute mean error: `(1/N) * sum(|num_i - ana_i|)` for two equal-length arrays.

#### Scenario: L1 of identical arrays
- **WHEN** numerical and analytical arrays are identical
- **THEN** `computeL1Error()` SHALL return 0.0

#### Scenario: L1 of known offset
- **WHEN** every numerical value is exactly 1.0 above the analytical value
- **THEN** `computeL1Error()` SHALL return 1.0

### Requirement: High-resolution SolVi parameter files
The system SHALL include parameter files `SolViRes101.txt`, `SolViRes151.txt`, and `SolViRes201.txt` in `TESTS/SolViBenchmark/`. Each SHALL be identical to `SolViRes41.txt` except for `Nx`, `Nz`, and `writer_subfolder`.

#### Scenario: Parameter file consistency
- **WHEN** comparing `SolViRes101.txt` to `SolViRes41.txt`
- **THEN** only `Nx`, `Nz`, and `writer_subfolder` SHALL differ

### Requirement: HighResL1Convergence CI test
The system SHALL include a `HighResL1Convergence` test in `SolViBenchmarkTests.cpp` that runs SolVi at resolutions 41, 81, 101, 151, and 201, computes both L2 and L1 norms for Vx and P at each resolution, and prints convergence-order tables.

#### Scenario: L1 pressure convergence at high resolution
- **WHEN** running SolVi at 101×101 and 201×201
- **THEN** the L1 pressure convergence order SHALL be ≥ 0.8

#### Scenario: L1 velocity convergence at high resolution
- **WHEN** running SolVi at 101×101 and 201×201
- **THEN** the L1 Vx convergence order SHALL be ≥ 0.8

#### Scenario: L2 pressure convergence (informational)
- **WHEN** running SolVi at 101×101 and 201×201
- **THEN** the L2 pressure convergence order SHALL be printed but not asserted (existing L2 assertions at 41→81 remain unchanged)

#### Scenario: Output table format
- **WHEN** the test completes
- **THEN** it SHALL print a table showing resolution, h, L2(Vx), L1(Vx), L2(P), L1(P) for each resolution, followed by convergence orders for each consecutive pair

#### Scenario: All resolutions converge
- **WHEN** the solver runs at each of the 5 resolutions
- **THEN** the solver SHALL converge (final residual < 1e-8) at every resolution

### Requirement: Existing CI test L1 reporting
The existing `GridConvergence` test in `SolViBenchmarkTests.cpp` SHALL additionally compute and print L1 norms alongside existing L2 output. No new assertions SHALL be added to the existing test.

#### Scenario: No regression in existing assertions
- **WHEN** `GridConvergence` runs
- **THEN** existing `EXPECT_GE(order_Vx, 0.7)` and `EXPECT_GE(order_P, 0.4)` assertions SHALL remain unchanged
