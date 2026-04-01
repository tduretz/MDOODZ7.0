## ADDED Requirements

### Requirement: Composition advection test
The test suite SHALL verify that the composition field `X` (Centers/X) is transported by the velocity field. Under deformation, the initial composition pattern SHALL evolve (not remain frozen).

#### Scenario: Composition field changes under advection
- **WHEN** a simulation runs with advection enabled and non-zero strain rate for 3 steps
- **THEN** `max(X)` or the spatial distribution of X at step 3 SHALL differ from step 0

#### Scenario: Composition remains bounded
- **WHEN** composition is advected
- **THEN** `min(X)` SHALL be ≥ 0 and `max(X)` SHALL be ≤ 1 (no overshoots from interpolation)
