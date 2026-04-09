## ADDED Requirements

### Requirement: Isolated power-law dislocation creep test
The test suite SHALL verify that power-law (dislocation) creep activates in isolation when `pwlv=40` is set. The strain-rate invariant `eII_pwl` SHALL be greater than zero, and the effective viscosity SHALL be less than the constant-viscosity reference (shear-thinning for n>1).

#### Scenario: Power-law strain rate is non-zero
- **WHEN** a simulation runs with `pwlv=40` (olivine dislocation creep, n=3.5) as the only active creep mechanism
- **THEN** `max(eII_pwl)` from Centers SHALL be greater than zero

#### Scenario: Shear-thinning reduces viscosity
- **WHEN** power-law creep is active with stress exponent n>1
- **THEN** `max(eta_n)` SHALL be less than the constant-viscosity reference value (η₀ of the inclusion)
