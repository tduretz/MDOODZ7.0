---
name: skill-particle-reseeding
description: MDOODZ particle reseeding and recycling — reseeding modes 0/1, particle deactivation (excess culling), recycling mechanism, markers struct fields, min_part_cell threshold, Nb_part_max hard limit, and common failure modes.
---

# MDOODZ Particle Reseeding & Recycling

## Overview

MDOODZ uses marker-in-cell advection. Particles drift with the flow and can cluster in convergent zones (downwellings, corners) while leaving divergent cells depleted. The reseeding system detects under-populated cells and injects new particles, while the deactivation system culls excess particles from over-populated cells. Dead particles (`phase = -1`) are recycled as slots for new ones.

All reseeding logic lives in `MDLIB/ParticleRoutines.c`. The file `MDLIB/ParticleReseeding.c` contained a dead function (`CountPartCell2`) that has been removed from the build.

## Reseeding Modes

Set `reseeding_mode` in the `.txt` parameter file:

| Mode | Function | Threading | Cell→particle mapping | Notes |
|------|----------|-----------|----------------------|-------|
| 0 | `CountPartCell_OLD` | Per-thread `ipreuse[]` | `ipcell[ith][kc]` (extended fine mesh) | `nthreads` hardcoded to 1 |
| 1 | `CountPartCell` | Global `part_reuse[]` | `part_cell[kc][nb]` (fine mesh) | Default; calls `exit(190)` at `Nb_part_max` |

Both modes use a **fine mesh** (2×Ncx × 2×Ncz) for cell-particle assignment, then map back to the coarse mesh for reseeding decisions.

## Key Struct Fields (`markers`)

| Field | Type | Meaning |
|-------|------|---------|
| `Nb_part` | `int` | Current active particle count |
| `Nb_part_max` | `int` | Hard allocation limit (4.1× initial count) |
| `min_part_cell` | `int` | Minimum particles per cell (from `.txt`); triggers reseeding |
| `phase[]` | `int*` | Phase index; `-1` = dead/air (available for recycling) |
| `Nx_part` | `int` | Particles per cell in x (typically 4) |
| `Nz_part` | `int` | Particles per cell in z (typically 4) |

## Reseeding Threshold

A cell is considered under-populated when its particle count drops below `min_part_cell`. The reseeding routine (`AddPartCell2`) injects new particles at sub-cell positions, assigning phase from the nearest existing neighbour.

## Deactivation (Excess Culling)

After reseeding, both modes run a **deactivation pass** that scans all cells on the fine mesh. If a cell has more than `min_part_cell + 4` particles, the excess are marked `phase = -1` (dead), making their array slots available for recycling.

**Boundary guard:** Cells with `BCt.type == 30` (Dirichlet thermal boundary) are skipped — deactivating boundary particles would corrupt boundary conditions.

### Mode 1 Deactivation (CountPartCell)

Uses the existing `part_cell[kc][nb]` mapping and `mesh->BCt_fine.type[kc]` for the boundary check. Runs inside the `if (reseed_markers==1)` block, before `part_cell` is freed.

### Mode 0 Deactivation (CountPartCell_OLD)

Builds a fresh `npc_g[]` / `ipc_g[][]` count on the fine mesh (2×Ncx × 2×Ncz), then maps fine cell indices → coarse cell indices to check `mesh->BCt.type` (since `BCt_fine` is not populated in mode 0). Runs after the parallel reseeding section.

## Recycling Mechanism

- **Mode 1:** After counting, builds `part_reuse[]` — an array of indices where `phase == -1`. New particles are assigned to these slots first. When slots run out, `Nb_part` is incremented (up to `Nb_part_max`).
- **Mode 0:** Each thread maintains `ipreuse[ith][]`. Same logic: dead slots first, then fresh indices.

## Nb_part in HDF5 Output

`Nb_part` is written as a 1-element integer dataset in the `Model` group of each HDF5 output file. Use `readScalarInt(file, "Model", "Nb_part")` from `TestHelpers.h` to read it in tests.

## Common Failure Modes

| Symptom | Cause | Fix |
|---------|-------|-----|
| `exit(190)` crash | `Nb_part` reached `Nb_part_max` | Deactivation not working or threshold too high |
| Particles cluster in corners | Convection concentrates markers in downwellings | Normal — deactivation culls excess |
| Phase smearing at boundaries | Deactivation ran on boundary cells | Check `BCt.type != 30` guard |
| `Nb_part` grows monotonically | No deactivation pass active | Verify `reseed_markers == 1` and deactivation code is present |

## .txt Parameters

```text
particles_per_cell_x = 4     # Nx_part
particles_per_cell_z = 4     # Nz_part
min_part_cell        = 16    # Reseeding threshold
reseeding_mode       = 1     # 0 or 1
```

## BlankenBench Test Coverage

The `ConvectionDevelops` test in `TESTS/BlankenBenchTests.cpp` asserts that `Nb_part` stays within 2× the initial count after 500 steps, confirming the recycling mechanism works for closed-box convection.
