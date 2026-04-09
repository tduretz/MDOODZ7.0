## ADDED Requirements

### Requirement: Dilatant flow produces non-zero divergence
The test suite SHALL verify that when `compressible=1` and dilatant plasticity is active (ψ > 0), the velocity divergence `divu` is non-zero.

#### Scenario: div(v) is non-zero with dilatancy
- **WHEN** a simulation runs with `compressible=1`, `plast=1`, high strain rate, and dilation angle ψ > 0
- **THEN** `max(|divu|)` SHALL be greater than zero

### Requirement: Incompressible flow has near-zero divergence
The test suite SHALL verify that when `compressible=0`, the velocity divergence is negligible.

#### Scenario: div(v) is approximately zero without compressibility
- **WHEN** a simulation runs with `compressible=0` and `plast=1`
- **THEN** `max(|divu|)` SHALL be less than 1e-6 times the strain rate
