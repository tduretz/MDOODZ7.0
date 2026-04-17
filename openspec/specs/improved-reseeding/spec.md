## ADDED Requirements

### Requirement: New reseed_mode=2 with CountPartCell_v2
The system SHALL provide a new reseeding mode activated by setting `reseed_mode=2` in the `.txt` parameter file. This mode SHALL call a new function `CountPartCell_v2` at all dispatch sites where modes 0 and 1 are dispatched in `Main_DOODZ.c`. The function SHALL have the same signature as `CountPartCell`. Modes 0 and 1 SHALL remain unchanged and backward compatible.

#### Scenario: Mode 2 is dispatched correctly
- **WHEN** `reseed_mode=2` is set in the `.txt` parameter file
- **THEN** `CountPartCell_v2` SHALL be called at all three dispatch sites in `Main_DOODZ.c` (initial count, post-sedimentation count, reseed+count)

#### Scenario: Modes 0 and 1 unchanged
- **WHEN** `reseed_mode=0` or `reseed_mode=1` is set
- **THEN** the existing `CountPartCell_OLD` or `CountPartCell` SHALL be called with their original (pre-deactivation-pass) behaviour

#### Scenario: Default mode unchanged
- **WHEN** `reseed_mode` is not specified in the `.txt` file
- **THEN** the default value (0) SHALL be used, preserving existing behaviour

### Requirement: Distance-based deactivation ordering in mode 2
When a fine cell contains more than `min_part_cell + 4` particles, `CountPartCell_v2` SHALL deactivate excess particles in order of decreasing distance from the fine-cell centroid. Particles farthest from the centroid SHALL be deactivated first. The distance metric SHALL be Euclidean: `(x - xc)² + (z - zc)²`.

#### Scenario: Over-populated cell deactivates farthest particles
- **WHEN** a fine cell contains 25 particles and `min_part_cell` is 16 (threshold = 20)
- **THEN** the 5 particles farthest from the fine-cell centroid SHALL be deactivated (phase set to -1), leaving 20 particles, with the 20 closest to the centroid retained

#### Scenario: Tie-breaking in distance ordering
- **WHEN** two particles have identical distance from the centroid and both are candidates for deactivation
- **THEN** either particle MAY be deactivated (tie-breaking order is not specified)

#### Scenario: Boundary cells excluded from deactivation
- **WHEN** a fine cell has `BCt_fine.type == 30` (boundary cell)
- **THEN** no deactivation SHALL occur in that cell regardless of particle count

### Requirement: Multi-particle fill-to-target reseeding in mode 2
When a fine cell contains fewer than 2 particles, `CountPartCell_v2` SHALL add enough particles in a single pass to bring the count up to 2. Each added particle SHALL be assigned properties via `AssignMarkerProperties` using the same nearest-phase search as mode 1.

#### Scenario: Empty cell receives multiple particles
- **WHEN** a fine cell contains 0 particles and the reseeding threshold is 2
- **THEN** 2 particles SHALL be added to that cell in a single pass

#### Scenario: Cell with one particle receives fill-up
- **WHEN** a fine cell contains 1 particle and the reseeding threshold is 2
- **THEN** 1 particle SHALL be added to bring the count to 2

#### Scenario: Cell at or above threshold is not reseeded
- **WHEN** a fine cell contains 2 or more particles
- **THEN** no particles SHALL be added to that cell

#### Scenario: Particle slot recycling during multi-fill
- **WHEN** multiple particles are added to a cell and dead particle slots exist in `part_reuse[]`
- **THEN** dead slots SHALL be reused before incrementing `Nb_part`, consistent with existing recycling behaviour

#### Scenario: Nb_part_max reached during multi-fill
- **WHEN** the total particle count reaches `Nb_part_max` during a multi-particle fill
- **THEN** the system SHALL exit with `exit(190)`, same as current behaviour

### Requirement: Randomized particle placement in mode 2
New particles added during reseeding in `CountPartCell_v2` SHALL be placed at uniformly random positions within the fine cell bounds, not at the cell centroid. The position SHALL be computed as `xc - dx_fine/2 + dx_fine * U` and `zc - dz_fine/2 + dz_fine * U` where `U` is a uniform random variate in [0, 1) from `rand() / RAND_MAX`, `dx_fine = mesh->dx / 2`, and `dz_fine = mesh->dz / 2`.

#### Scenario: New particle position is within cell bounds
- **WHEN** a particle is added to a fine cell with centroid (xc, zc) and half-widths (dx_fine/2, dz_fine/2)
- **THEN** the particle position (new_x, new_z) SHALL satisfy `xc - dx_fine/2 <= new_x < xc + dx_fine/2` and `zc - dz_fine/2 <= new_z < zc + dz_fine/2`

#### Scenario: Multiple particles in same cell have distinct positions
- **WHEN** 2 or more particles are added to the same cell in a single pass (via fill-to-target)
- **THEN** each particle SHALL receive an independently sampled random position (positions will differ with high probability)

### Requirement: Revert deactivation passes from modes 0 and 1
The deactivation blocks previously added to `CountPartCell` (mode 1, lines ~2895–2911) and `CountPartCell_OLD` (mode 0, lines ~3776–3851) SHALL be removed. This restores both functions to their original behaviour before the deactivation change was introduced.

#### Scenario: Mode 1 has no deactivation
- **WHEN** `reseed_mode=1` is active and a cell exceeds `min_part_cell + 4` particles
- **THEN** `CountPartCell` SHALL NOT deactivate any excess particles (original behaviour)

#### Scenario: Mode 0 has no deactivation
- **WHEN** `reseed_mode=0` is active and a cell exceeds `min_part_cell + 4` particles
- **THEN** `CountPartCell_OLD` SHALL NOT deactivate any excess particles (original behaviour)

### Requirement: Rigid body rotation advection test for reseeding quality
A GTest SHALL exist that quantitatively compares L2 composition errors across all three reseed_modes using a prescribed rigid body rotation velocity field. See `specs/rotation-advection-test/spec.md` for full requirements.
