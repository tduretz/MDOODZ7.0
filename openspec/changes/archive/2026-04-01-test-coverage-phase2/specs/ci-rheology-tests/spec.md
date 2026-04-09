## MODIFIED Requirements

### Requirement: Diffusion creep strain rate verification
The existing diffusion creep tests SHALL additionally verify that the mechanism-specific strain-rate invariant `eII_lin` is non-zero when diffusion creep is active, and that smaller grain size produces higher eII_lin (grain-size sensitivity check).

#### Scenario: eII_lin is non-zero for diffusion creep
- **WHEN** a simulation runs with `linv=40` (diffusion creep active)
- **THEN** `max(eII_lin)` from Centers SHALL be greater than zero

#### Scenario: Smaller grain size produces higher diffusion strain rate
- **WHEN** two DiffusionCreep runs are compared — one with gs_ref=5e-3 and one with gs_ref=1e-3
- **THEN** `max(eII_lin)` for the smaller grain size SHALL be greater than for larger grain size
