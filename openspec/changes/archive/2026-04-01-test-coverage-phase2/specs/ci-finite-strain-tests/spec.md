## ADDED Requirements

### Requirement: Deformation gradient departs from identity
The test suite SHALL verify that the deformation gradient tensor F evolves away from the identity matrix under finite deformation when `finite_strain=1`.

#### Scenario: Fxx grows under extension (pure shear)
- **WHEN** a pure-shear simulation runs for 5 steps with `finite_strain=1`
- **THEN** `max(Fxx)` SHALL be greater than 1.0 (extension in x-direction)

#### Scenario: Fzz decreases under compression (pure shear)
- **WHEN** a pure-shear simulation runs for 5 steps with `finite_strain=1`
- **THEN** `min(Fzz)` SHALL be less than 1.0 (compression in z-direction)

#### Scenario: F tensor is finite
- **WHEN** `finite_strain=1` is enabled
- **THEN** all components Fxx, Fxz, Fzx, Fzz SHALL be finite (no NaN or overflow)
