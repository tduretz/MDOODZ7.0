## Context

MDOODZ uses marker-in-cell advection with two active reseeding paths:

- **Mode 1** (`CountPartCell` in `ParticleRoutines.c`): Builds a global `part_reuse[]` array of dead particles (`phase == -1`), then reseeds under-populated cells by recycling dead slots or allocating new ones at `Nb_part`. Hits `exit(190)` if `Nb_part` reaches `Nb_part_max`.
- **Mode 0** (`CountPartCell_OLD` in `ParticleRoutines.c`): Per-thread recycling via `ipreuse[]`. Reseeds under-populated cells within OpenMP parallel regions using `AddPartCell2`.

Neither path deactivates excess particles in over-populated cells. A third function, `CountPartCell2` in `ParticleReseeding.c`, contains deactivation logic but is never called — it is dead code along with `AddPartCell`, `AddPartVert`, `ind_part_reuse`, and `inc_reuse`.

The `markers` struct tracks `Nb_part` (current count), `Nb_part_max` (4.1× initial, hard limit), and `min_part_cell` (reseeding threshold, read from `.txt`). The `.txt` files also accept `max_part_cell` but the C code never reads it.

In closed-box simulations (BlankenBench convection), no particles leave the domain, so no dead particles accumulate for recycling. Reseeding creates new particles → `Nb_part` grows monotonically → `exit(190)`.

## Goals / Non-Goals

**Goals:**
- **Priority**: Get particle recycling working in closed-box simulations — deactivate excess particles so they are reused, not accumulated
- Do not increase particle limits (`Nb_part_max`) — the fix is recycling, not raising caps
- Remove dead code that will never be called now that the active paths have deactivation
- Validate the fix with a BlankenBench particle stability test
- Document the particle reseeding system in a Copilot skill

**Non-Goals:**
- Changing the reseeding algorithm itself (how new particles are created or properties interpolated)
- Supporting dynamic reallocation of particle arrays (the fixed `Nb_part_max` allocation stays)
- Adding a third reseeding mode — both mode 0 and mode 1 get deactivation, no new mode
- Optimising particle memory layout or introducing particle pooling
- Adding new `.txt` parameters (`max_part_cell` configurability can be revisited later)

## Decisions

### D1: Use existing `min_part_cell`-based threshold for deactivation — no new parameters

The dead `CountPartCell2` uses `min_part_cell + 4` as the deactivation threshold. This is simple and already proven. The priority is getting the recycling mechanism working in the active code paths, not adding new configurable parameters.

**Decision**: Use `min_part_cell + 4` as the deactivation threshold, matching the existing `CountPartCell2` logic. Do not add `max_part_cell` to the `markers` struct or read it from `.txt`. Do not increase `Nb_part_max`. The fix is purely about recycling: deactivate excess particles (`phase = -1`) so they become available for reuse by the existing reseeding logic. This is the top priority of the change.

**Alternative considered**: Add `max_part_cell` as a user-configurable `.txt` parameter. Rejected — adds unnecessary complexity. The heuristic threshold works and can be revisited later if needed.

### D2: Insert deactivation as a separate pass after reseeding, not inline

The deactivation scan iterates over all cells and marks excess particles as `phase = -1`. This is a distinct concern from the reseeding (adding particles to under-populated cells).

**Decision**: Add the deactivation as a separate loop after the reseeding loop in both `CountPartCell` and `CountPartCell_OLD`, following the same pattern as `CountPartCell2` lines 668–685. This keeps reseeding and deactivation logic cleanly separated and easy to reason about.

**Alternative considered**: Inline deactivation into the existing per-cell loops. Rejected because it interleaves two concerns and makes the already-complex reseeding functions harder to follow.

### D3: Use `ind_per_cell[][]` for deactivation (requires building it in mode 1)

The deactivation logic needs to know which particles belong to each cell, so it can mark the extras. `CountPartCell2` uses `ind_per_cell[ic][nb]` — a per-cell list of particle indices.

- **Mode 0** (`CountPartCell_OLD`): Already builds `ind_per_cell` as part of its reseeding logic.
- **Mode 1** (`CountPartCell`): Builds `part_cell[][]` (similar structure) during reseeding. When `reseed_markers == 0` (no reseeding), neither structure is built.

**Decision**: In mode 1, reuse `part_cell[][]` for deactivation when reseeding is active. When reseeding is off but deactivation is still desired, build a lightweight cell→particle index from `mesh->nb_part_cell` and `particles->phase`. In mode 0, reuse `ind_per_cell`.

### D4: Delete dead code rather than annotating it

`CountPartCell2`, `AddPartCell`, `AddPartVert`, and `inc_reuse` are dead code. Once deactivation is ported to the active paths, `CountPartCell2` has no remaining purpose.

**Decision**: Delete the dead functions and variables outright. They are in version control if ever needed. Annotating dead code with comments doesn't prevent confusion — someone will still wonder why it exists.

### D5: CI smoke test checks `Nb_part` stability; steady-state test tuned to converge fast (~1000 steps)

The BlankenBench steady-state run currently uses 100k steps — too expensive even for local development. The simulation settings should be tuned aggressively (coarser grid, larger `Courant`, higher `dt_max`) so thermal steady state is reached in roughly 1000 steps (soft cap — may need somewhat more, but should stay in that order of magnitude).

**Decision**: The CI smoke test runs a short batch of steps and asserts that `Nb_part` at the final step is not significantly larger than the initial count (within 10%). This catches the "monotonic growth" bug.

For the steady-state test (gated by `BLANKENBACH_STEADY=1`): tune `BlankenBenchSteady.txt` to reach thermal steady state in ~1000 steps by increasing `Courant`, `dt_max`, and possibly using a coarser grid. Use `.dat` breakpoint restart (`irestart=1`, `istep=N`) so we can re-run from a checkpoint if the first attempt doesn't converge. The final steady-state values must match the published reference: Nu = 4.884, Vrms = 42.865 (within 5% tolerance).

If the CI smoke test's step count isn't enough to trigger reseeding, we can increase it or tune time stepping to accelerate evolution.

### D6: Expose `Nb_part` in HDF5 output for test assertion

The GTest tests read HDF5 output files. To assert on `Nb_part`, we need it available in the output.

**Decision**: Check if `Nb_part` is already written to HDF5 output. If not, add it as a scalar attribute (minimal change). The test reads it from the final step's output file.

**Alternative considered**: Parse stdout/log for particle count. Rejected — fragile and couples tests to log format.

## Risks / Trade-offs

**[Deactivation may remove "good" particles]** → The deactivation picks particles beyond `max_part_cell` in cell-index order, which is arbitrary. This could remove a particle with important strain history. → Mitigation: This matches the existing `CountPartCell2` approach. Particle properties are interpolated from neighbours during reseeding, so the cell retains representative coverage. A future improvement could prioritise removing the youngest particles, but this is out of scope.

**[Mode 0 thread safety]** → `CountPartCell_OLD` runs reseeding inside OpenMP parallel regions. The deactivation pass must either be serial (after the parallel section) or use thread-local counting. → Mitigation: Run deactivation as a serial pass after the parallel reseeding loop, matching the `CountPartCell2` pattern. The cost is negligible — it's a simple scan over `ncx × ncz` cells.

**[Step count sensitivity for CI]** → 500 steps may not trigger enough cell over-population to exercise deactivation. → Mitigation: Tune the BlankenBench CI `.txt` parameters (increase `Courant`, reduce grid resolution, increase `dt_max`) so reseeding fires within fewer steps. Verify empirically before finalising the test.
