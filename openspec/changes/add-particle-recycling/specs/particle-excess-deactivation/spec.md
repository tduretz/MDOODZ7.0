## ADDED Requirements

### Requirement: Excess particles deactivated in mode 1 reseeding
`CountPartCell` (mode 1) in `ParticleRoutines.c` SHALL include a deactivation pass after reseeding that scans all cells and marks particles beyond `min_part_cell + 4` as `phase = -1`. The pass SHALL use the existing `part_cell[][]` structure to identify which particles belong to each cell.

#### Scenario: Closed-box convection does not crash
- **WHEN** a BlankenBench convection simulation runs with `reseed_mode = 1` for enough steps to trigger reseeding
- **THEN** `Nb_part` SHALL remain bounded (not grow monotonically) and the simulation SHALL NOT call `exit(190)`

#### Scenario: Excess particles become available for recycling
- **WHEN** a cell accumulates more than `min_part_cell + 4` particles due to advection
- **THEN** the deactivation pass SHALL set `phase = -1` on the excess particles so they enter the `part_reuse[]` pool on the next reseeding call

### Requirement: Excess particles deactivated in mode 0 reseeding
`CountPartCell_OLD` (mode 0) in `ParticleRoutines.c` SHALL include a serial deactivation pass after the parallel reseeding loop. The pass SHALL scan all cells using `ind_per_cell[][]` and mark particles beyond `min_part_cell + 4` as `phase = -1`.

#### Scenario: Mode 0 closed-box stability
- **WHEN** a closed-box simulation runs with `reseed_mode = 0` for enough steps to trigger reseeding
- **THEN** `Nb_part` SHALL remain bounded and the simulation SHALL NOT call `exit(1)` due to particle overflow

#### Scenario: Thread safety of deactivation pass
- **WHEN** the deactivation pass runs after the OpenMP-parallel reseeding loop in `CountPartCell_OLD`
- **THEN** the pass SHALL execute serially (outside any `#pragma omp parallel` region) to avoid race conditions on `particles->phase[]`

### Requirement: Deactivation threshold uses min_part_cell + 4
The deactivation threshold SHALL be `particles->min_part_cell + 4`, matching the logic in the existing dead `CountPartCell2`. No new `.txt` parameters SHALL be added for this threshold.

#### Scenario: Default threshold with min_part_cell = 16
- **WHEN** `min_part_cell` is set to 16 in the `.txt` file
- **THEN** cells with more than 20 particles SHALL have excess particles deactivated

### Requirement: Deactivation skips boundary cells
The deactivation pass SHALL skip cells where `mesh->BCt.type[ic] == 30` (boundary condition type 30), matching the guard in the existing `CountPartCell2` logic.

#### Scenario: Boundary cells untouched
- **WHEN** a cell has `BCt.type == 30`
- **THEN** the deactivation pass SHALL NOT modify any particles in that cell, regardless of particle count

### Requirement: Dead code removed from ParticleReseeding.c
The following dead code SHALL be deleted from `ParticleReseeding.c`: the entire `CountPartCell2` function, the `inc_reuse` variable declaration, and any allocations of `ind_part_reuse` / `ind_part_reuse2` that are never consumed.

#### Scenario: CountPartCell2 removed
- **WHEN** the codebase is searched for `CountPartCell2`
- **THEN** no definition or declaration of `CountPartCell2` SHALL exist

### Requirement: Dead code removed from ParticleRoutines.c
The `AddPartCell` and `AddPartVert` functions in `ParticleRoutines.c` SHALL be deleted. These functions are never called.

#### Scenario: AddPartCell and AddPartVert removed
- **WHEN** the codebase is searched for `AddPartCell` or `AddPartVert`
- **THEN** no definition or declaration of these functions SHALL exist (note: `AddPartCell2` is actively used and SHALL NOT be removed)
