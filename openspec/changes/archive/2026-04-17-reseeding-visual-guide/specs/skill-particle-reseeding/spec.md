## ADDED Requirements

### Requirement: DEBUG logging traces reseeding decisions per cell
`ParticleRoutines.c` SHALL include `LOG_DBG` calls in `CountPartCell`, `CountPartCell_v2`, and `CountPartCell_OLD` that trace per-cell reseeding decisions: particle count vs threshold, number of new particles injected, and slot reuse vs append.

#### Scenario: Reseeding trigger logged
- **WHEN** a cell's particle count falls below the reseeding threshold and `log_level >= 4`
- **THEN** a DEBUG log line SHALL be emitted containing the cell indices, current count, threshold, and number of new particles to inject

#### Scenario: New particle placement logged
- **WHEN** a new particle is created during reseeding and `log_level >= 4`
- **THEN** a DEBUG log line SHALL be emitted containing the particle coordinates and whether a recycled slot or a new array index was used

#### Scenario: No output at default log level
- **WHEN** `log_level < 4` (default is 2 = INFO)
- **THEN** no reseeding DEBUG log lines SHALL be emitted

### Requirement: DEBUG logging traces deactivation decisions (mode 2)
`CountPartCell_v2` SHALL include `LOG_DBG` calls that trace per-cell deactivation: excess count, number of particles deactivated, and the distance criterion used.

#### Scenario: Deactivation trigger logged
- **WHEN** a cell exceeds `min_part_cell + 4` particles in mode 2 and `log_level >= 4`
- **THEN** a DEBUG log line SHALL be emitted containing the cell indices, particle count, excess count, and number to deactivate

#### Scenario: Individual deactivation logged
- **WHEN** a particle is marked `phase = -1` during deactivation and `log_level >= 4`
- **THEN** a DEBUG log line SHALL be emitted containing the particle's array index and its distance from the cell centroid

### Requirement: DEBUG logging traces Nb_part updates
After reseeding completes, a `LOG_DBG` call SHALL emit the before/after `Nb_part` value and the breakdown of reused vs newly created slots.

#### Scenario: Nb_part summary logged
- **WHEN** reseeding finishes for a timestep and `log_level >= 4`
- **THEN** a DEBUG log line SHALL be emitted with format: `"Nb_part: <before> -> <after> (reused=<N>, created=<M>)"`
