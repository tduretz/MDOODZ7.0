## ADDED Requirements

### Requirement: Skill file exists at standard path
A Copilot skill file SHALL exist at `.github/skills/skill-particle-reseeding/SKILL.md` with valid frontmatter (name, description).

#### Scenario: Skill file discoverable
- **WHEN** a user or Copilot scans `.github/skills/` for available skills
- **THEN** `skill-particle-reseeding/SKILL.md` SHALL be present and parseable

### Requirement: Documents reseeding modes
The skill SHALL document both active reseeding modes: mode 0 (`CountPartCell_OLD`, per-thread recycling via `ipreuse[]`) and mode 1 (`CountPartCell`, global `part_reuse[]` array). It SHALL explain when each mode is selected (`reseed_mode` parameter in `.txt`).

#### Scenario: User asks about reseeding modes
- **WHEN** a user asks "how does particle reseeding work in MDOODZ"
- **THEN** the skill SHALL provide enough context to understand both modes, their data structures, and the `.txt` parameter that selects between them

### Requirement: Documents recycling and deactivation
The skill SHALL explain the recycling mechanism (dead particles with `phase == -1` are reused for new particles) and the deactivation mechanism (excess particles in over-populated cells are marked `phase = -1` when count exceeds `min_part_cell + 4`).

#### Scenario: User asks about closed-box particle stability
- **WHEN** a user asks "why does my convection simulation crash with exit(190)"
- **THEN** the skill SHALL explain the relationship between deactivation, recycling, and `Nb_part_max`, and confirm that deactivation prevents unbounded growth

### Requirement: Documents markers struct particle fields
The skill SHALL document the key `markers` struct fields for particle management: `Nb_part`, `Nb_part_max`, `Nb_part_ini`, `min_part_cell`, `phase[]`, and the 4.1× allocation multiplier.

#### Scenario: User asks about Nb_part_max
- **WHEN** a user asks "what is the particle limit"
- **THEN** the skill SHALL explain `Nb_part_max = 4.1 × Nb_part_ini` and that it is a hard allocation limit, not configurable at runtime

### Requirement: Documents common failure modes
The skill SHALL document at least: (1) `exit(190)` from `Nb_part` reaching `Nb_part_max`, (2) closed-box simulations lacking natural particle death, and (3) the relationship between `min_part_cell` and the reseeding/deactivation thresholds.

#### Scenario: Troubleshooting particle crash
- **WHEN** a user encounters `exit(190)` or "Max number of particles reached"
- **THEN** the skill SHALL guide them to check whether deactivation is active and whether the simulation is closed-box

### Requirement: Registered in copilot-instructions.md
The skill SHALL be listed in the skill table in `.github/copilot-instructions.md` with an appropriate trigger description (e.g., "Particle reseeding, recycling, deactivation, markers, Nb_part").

#### Scenario: Copilot routes particle questions to skill
- **WHEN** a user asks about particle reseeding, recycling, or particle count issues
- **THEN** the copilot-instructions.md table SHALL map that topic to `skill-particle-reseeding`

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
