## ADDED Requirements

### Requirement: Peierls creep visual regression test
The visual test suite SHALL include a scenario exercising Peierls (exponential) creep with reference HDF5 comparison.

#### Scenario: Peierls creep strain rate field matches reference
- **WHEN** a Peierls creep scenario runs to completion (e.g., ShearTemplate with `expv=1`, ~5 time steps, moderate resolution)
- **THEN** the relative L2 error between the computed strain rate invariant (ε̇_II) field and the reference field SHALL be less than 1e-5

#### Scenario: Peierls creep stress field matches reference
- **WHEN** the same scenario completes
- **THEN** the relative L2 error between the computed stress invariant (τ_II) field and the reference field SHALL be less than 1e-5

### Requirement: Anisotropic VEP visual regression test
The visual test suite SHALL include a scenario combining visco-elasto-plasticity with mechanical anisotropy.

#### Scenario: Anisotropic VEP stress evolution matches reference
- **WHEN** an extension scenario runs with `elastic=1`, `anisotropy=1`, `finite_strain=1`, plasticity enabled, for ~10 time steps
- **THEN** the relative L2 error of the stress invariant field against the reference SHALL be less than 1e-4

#### Scenario: Anisotropic VEP director field matches reference
- **WHEN** the same scenario completes
- **THEN** the relative L2 error of the director components (`nx`, `nz`) against the reference SHALL be less than 1e-4

### Requirement: Visual test gnuplot output
Each new visual test scenario SHALL produce gnuplot-based comparison images alongside the numeric verification.

#### Scenario: Side-by-side images generated
- **WHEN** any visual test scenario completes
- **THEN** the test SHALL write a `.dat` file and invoke its `.gnu` script to produce a PNG comparing the run with the reference
