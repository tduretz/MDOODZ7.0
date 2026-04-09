## ADDED Requirements

### Requirement: Mechanism-specific strain rates are non-zero when enabled
The test suite SHALL verify that each creep mechanism produces a non-zero strain-rate invariant (eII_pwl, eII_lin, eII_exp, eII_gbs, eII_el, eII_pl) when that mechanism is the only one active.

#### Scenario: Power-law strain rate is non-zero
- **WHEN** only `pwlv` is set (dislocation creep active)
- **THEN** `max(eII_pwl)` SHALL be greater than zero

#### Scenario: Diffusion creep strain rate is non-zero
- **WHEN** only `linv` is set (diffusion creep active)
- **THEN** `max(eII_lin)` SHALL be greater than zero

### Requirement: Composite strain rates sum correctly
The test suite SHALL verify that when multiple creep mechanisms are active simultaneously, the individual mechanism strain rates are each non-zero and their magnitudes are physically consistent (the total strain rate is partitioned among mechanisms).

#### Scenario: Both mechanisms contribute in composite flow
- **WHEN** both `pwlv=40` and `linv=40` are active simultaneously
- **THEN** both `max(eII_pwl)` > 0 and `max(eII_lin)` > 0
