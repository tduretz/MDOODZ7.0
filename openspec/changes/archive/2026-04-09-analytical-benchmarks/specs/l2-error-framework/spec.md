## ADDED Requirements

### Requirement: readFieldAsArray helper
`TestHelpers.h` SHALL provide a `readFieldAsArray()` function that reads a float field from an HDF5 group/dataset and returns it as a `std::vector<double>`.

#### Scenario: Read a Centers field
- **WHEN** `readFieldAsArray("Output00000.gzip.h5", "Centers", "T")` is called
- **THEN** it SHALL return a vector of length (Nx-1)×(Nz-1) containing the temperature values as doubles

#### Scenario: Read a VxNodes field
- **WHEN** `readFieldAsArray("Output00000.gzip.h5", "VxNodes", "Vx")` is called
- **THEN** it SHALL return a vector containing the Vx velocity field values

### Requirement: readCoordArray helper
`TestHelpers.h` SHALL provide a `readCoordArray()` function that reads a 1D coordinate array from the HDF5 `Model/` group and returns it as a `std::vector<double>`.

#### Scenario: Read center x-coordinates
- **WHEN** `readCoordArray("Output00000.gzip.h5", "xc_coord")` is called
- **THEN** it SHALL return a vector of center x-coordinate values in SI units

### Requirement: computeL2Error helper
`TestHelpers.h` SHALL provide a `computeL2Error()` function that computes the relative L2 error norm between a numerical field array and analytically computed values: L2 = sqrt(Σ(num - ana)² / Σ(ana)²). When the analytical solution is identically zero, it SHALL fall back to absolute L2: sqrt(Σ(num)² / N).

#### Scenario: Relative L2 for non-zero analytical solution
- **WHEN** `computeL2Error` is called with numerical values [1.01, 2.02, 3.03] and analytical values [1.0, 2.0, 3.0]
- **THEN** it SHALL return the relative L2 norm ≈ 0.01

#### Scenario: Absolute L2 for zero analytical solution
- **WHEN** `computeL2Error` is called with numerical values [0.001, -0.001, 0.002] and analytical values [0.0, 0.0, 0.0]
- **THEN** it SHALL return the absolute L2 norm (RMS of the numerical values) rather than dividing by zero

### Requirement: Theory documentation
A file `TESTS/AnalyticalSolutions.md` SHALL document each analytical solution used in the test suite: problem statement, analytical formula with LaTeX math, derivation sketch or literature reference, how it maps to the MDOODZ test, and the measured L2 baseline with threshold.

#### Scenario: SolVi section exists
- **WHEN** a reader opens `TESTS/AnalyticalSolutions.md`
- **THEN** there SHALL be a section documenting the SolVi benchmark including the Schmid & Podladchikov (2003) reference, the complex-variable pressure/velocity formulas, and the expected L2 error at each resolution

#### Scenario: 1D analytical solutions section exists
- **WHEN** a reader opens `TESTS/AnalyticalSolutions.md`
- **THEN** there SHALL be sections for each 1D analytical solution (steady-state geotherm, hydrostatic pressure, pure shear velocity, Maxwell stress accumulation, viscous dissipation, radiogenic heating, thermal expansion) explaining the formula and how it is applied to the 2D grid

#### Scenario: Anisotropy section exists
- **WHEN** a reader opens `TESTS/AnalyticalSolutions.md`
- **THEN** there SHALL be a section documenting the Mühlhaus director evolution equation, the anisotropic constitutive tensor, and the stress-vs-angle analytical formula

### Requirement: README cross-reference
`TESTS/README.md` SHALL reference `AnalyticalSolutions.md` in a dedicated "Analytical Benchmarks" section.

#### Scenario: README links to theory document
- **WHEN** a reader opens `TESTS/README.md`
- **THEN** there SHALL be a section titled "Analytical Benchmarks" containing a link to `AnalyticalSolutions.md` and a summary of the benchmarks performed
