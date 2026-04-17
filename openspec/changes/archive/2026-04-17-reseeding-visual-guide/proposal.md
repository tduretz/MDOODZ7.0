## Why

MDOODZ has three particle reseeding modes (0, 1, 2) with fundamentally different strategies for adding and removing markers, but no visual documentation explains how they work. Users choosing between modes — or debugging particle count issues — must read 500+ lines of C code. A visual guide with plots generated from real MDOODZ particle data would make the differences immediately obvious and serve as reference documentation.

## What Changes

- Add a Markdown guide (`VISUAL_TESTS/ReseedingGuide.md`) that embeds generated PNG plots with minimal explanatory text — "less words, more presentation"
- Add a standalone `.cpp` test runner (outside the main GTest suite) that runs short simulations with `writer_markers=1` for each reseeding mode and produces particle HDF5 snapshots
- Add gnuplot scripts (production-ready) that read the particle HDF5 output and generate PNGs. Python may be used during development, but the deliverable plotting pipeline is C++ + gnuplot
- Plots use **real particle coordinates** from `Particles<NNNNN>.gzip.h5` output (not theoretical/schematic). Each plot focuses on **one concept** — build separate plots rather than overloading a single figure:
  - Per-cell particle positions before and after reseeding (new particles appearing)
  - Per-cell particle positions before and after deactivation (particles being removed)
  - Which particles are selected for removal (index-order in mode 1 vs distance-from-centroid in mode 2)
  - Particle array slot recycling: how dead slots (`phase == -1`) get reused, shown with real array indices
  - Nb_part evolution over time for each mode (growth in modes 0/1 vs stability in mode 2)
- **Side-by-side mode comparisons** — the guide's core value is highlighting differences between modes 0, 1, and 2. Run the **same scenario** with all 3 modes and present results side-by-side:
  - **Reseeding placement**: mode 0 (regular grid within cell) vs mode 1 (regular grid) vs mode 2 (random within cell) — same cell, same particle deficit, 3 panels
  - **Deactivation strategy**: mode 0 (no deactivation) vs mode 1 (index-order removal) vs mode 2 (farthest-from-centroid removal) — same overpopulated cell, 3 panels showing which particles survive
  - **Particle count over time**: single time-series plot with 3 curves (one per mode) on the same axes — shows mode 0/1 unbounded growth vs mode 2 stability
  - **Array slot usage**: mode 0/1 (append-only, array grows) vs mode 2 (recycles dead slots) — show array index heatmaps or slot occupancy timelines
- **Edge cases** — dedicated scenarios that trigger and visualise boundary conditions:
  - Cell hits `Nb_part_max` hard limit — what happens to excess particles, how the code handles the overflow
  - Cell drops below `min_part_cell` — reseeding kicks in, show the before/after
  - Massive particle loss (e.g., advection sweeps a region empty) — bulk reseeding response
  - Array full (`Nb_part` approaches `Nb_part_max`) — recycling vs `exit(190)` behavior
  - Multi-phase cell — how reseeding picks the dominant phase for new particles
- **Test scenarios**: may reuse and tweak existing tests (e.g., rotation advection), or create entirely new purpose-built scenarios to showcase specific features and edge cases. Different plots may come from different test setups — whatever best demonstrates the concept
- **DEBUG logging in MDLIB**: add `LOG_DBG` calls inside `ParticleReseeding.c` to trace algorithm steps — e.g., "cell (i,j) has N particles, below min_part_cell → reseeding 4 new", "deactivating particle k (distance 0.73 from centroid)", "recycling slot 1042". These emit only at `log_level=4` (DEBUG) so zero cost in production, but provide narrated traces that can be included in the guide as annotated log excerpts
- The guide doubles as a **code reference**: each plot section links to the relevant `ParticleReseeding.c` variables, array names, and index semantics (e.g., `particles.x`, `particles.phase`, `Nb_part`, `Nb_part_max`, `min_part_cell`, the cell-local particle list, the free-slot recycling index)

## Capabilities

### New Capabilities
- `reseeding-visual-guide`: The standalone C++ test runner, gnuplot scripts, generated PNGs, and Markdown guide in `VISUAL_TESTS/` showing all 3 reseeding modes with real particle data and code-level references

### Modified Capabilities
- `skill-particle-reseeding`: Adding DEBUG-level logging to `ParticleReseeding.c` that traces reseeding/deactivation decisions per cell

## Impact

- **New files**: `VISUAL_TESTS/ReseedingGuide.md`, standalone `.cpp` runner, `.gp` gnuplot scripts, generated PNGs, `.txt` parameter file(s) with `writer_markers=1`
- **Existing code**: `MDLIB/ParticleReseeding.c` gets new `LOG_DBG` calls (DEBUG-level only, no behavioral change)
- **Dependencies**: gnuplot (plotting), optionally Python 3 + h5py/numpy/matplotlib during development
