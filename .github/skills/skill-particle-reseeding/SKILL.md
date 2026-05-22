---
name: skill-particle-reseeding
description: MDOODZ particle reseeding and recycling — reseeding modes 0/1, particle deactivation (excess culling), recycling mechanism, markers struct fields, min_part_cell threshold, Nb_part_max hard limit, and common failure modes.
---

# MDOODZ Particle Reseeding & Recycling

## Overview

MDOODZ uses marker-in-cell advection. Particles drift with the flow and can cluster in convergent zones (downwellings, corners) while leaving divergent cells depleted. The reseeding system detects under-populated cells and injects new particles, while the deactivation system culls excess particles from over-populated cells. Dead particles (`phase = -1`) are recycled as slots for new ones.

All reseeding logic lives in `MDLIB/ParticleRoutines.c`. The file `MDLIB/ParticleReseeding.c` contained a dead function (`CountPartCell2`) that has been removed from the build.

## Reseeding Modes

Set `reseed_mode` in the `.txt` parameter file:

| Mode | Function | Placement | Deactivation | Slot Recycling |
|------|----------|-----------|-------------|----------------|
| 0 | `CountPartCell_OLD` | Closest-neighbour (per-thread `ipreuse[]`) | None | Per-thread `ipreuse[]` |
| 1 | `CountPartCell` | Fine-cell centroid | None | Global `part_reuse[]` |
| 2 | `CountPartCell_v2` | Random within fine cell | Distance-from-centroid, farthest first | Global `part_reuse[]` |

All modes use a **fine mesh** (2×Ncx × 2×Ncz) for cell-particle assignment, then map back to the coarse mesh for reseeding decisions.

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

**Only mode 2** (`CountPartCell_v2`) actively deactivates excess particles. After reseeding, it scans all cells on the fine mesh. If a cell has more than `fine_threshold` particles (where `fine_threshold = ceil(min_part_cell / 4.0) + 4`), the excess are sorted by distance from the cell centroid and the farthest are marked `phase = -1` (dead), making their array slots available for recycling.

**Boundary guard:** Cells with `BCt.type == 30` (Dirichlet thermal boundary) are skipped — deactivating boundary particles would corrupt boundary conditions.

**Modes 0 and 1 have no deactivation.** Particles can only become dead (`phase = -1`) by exiting the domain boundaries or being above the free surface. In a closed-box simulation (e.g., Blankenbach convection), modes 0 and 1 will see monotonic growth of active particle count, eventually hitting `Nb_part_max` and crashing with `exit(190)`. Mode 2 is the only mode that guarantees bounded active particle count in closed environments.

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
| `Nb_part` grows monotonically | Using mode 0 or 1 (no deactivation) | Switch to `reseed_mode = 2` for closed-box simulations |

## .txt Parameters

```text
Nx_part          = 4     # particles_per_cell_x
Nz_part          = 4     # particles_per_cell_z
min_part_cell    = 16    # Reseeding threshold
reseed_mode      = 1     # 0, 1, or 2
```

## DEBUG Logging

With `log_level = 3` (DEBUG), all three reseeding functions emit per-cell and per-particle diagnostic messages:

- **Cell trigger**: `cell (i,j) fine: N particles (threshold=T), reseeding`
- **Particle placement**: `new particle at (x, z), reusing slot K` or `appended at index K`
- **Deactivation** (mode 2 only): `cell (i,j): N particles, excess E, deactivating E farthest` and `deactivated particle K (dist=D from centroid)`
- **Summary**: `Nb_part: <before> -> <after> (reused=R, created=C)`

Enable in `.txt`: `log_level = 3`. These messages are suppressed at the default `log_level = 2` (INFO).

See also: `VISUAL_TESTS/ReseedingGuide.md` for a visual guide with plots from real simulation output.

## BlankenBench Test Coverage

The `ConvectionDevelops` test in `TESTS/BlankenBenchTests.cpp` asserts that `Nb_part` stays within 2× the initial count after 500 steps, confirming the recycling mechanism works for closed-box convection.
