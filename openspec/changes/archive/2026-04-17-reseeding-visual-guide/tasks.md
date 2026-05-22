## 1. DEBUG Logging in ParticleRoutines.c

- [x] 1.1 Add `LOG_DBG` calls in `CountPartCell_OLD` (mode 0): cell particle count vs threshold, reseeding trigger, new particle placement with slot index
- [x] 1.2 Add `LOG_DBG` calls in `CountPartCell` (mode 1): cell particle count vs threshold, reseeding trigger, new particle placement (reused slot vs appended index)
- [x] 1.3 Add `LOG_DBG` calls in `CountPartCell_v2` (mode 2): same as 1.2 plus deactivation trigger (cell indices, excess count), individual deactivation (particle index, distance from centroid)
- [x] 1.4 Add `LOG_DBG` Nb_part summary at end of each function: `"Nb_part: <before> -> <after> (reused=<N>, created=<M>)"`
- [x] 1.5 Verify: build with `cmake --build`, run existing rotation advection tests with `log_level=4`, confirm DEBUG output appears and no regressions at default `log_level=2`

## 2. Parameter Files

- [x] 2.1 Create `VISUAL_TESTS/ReseedingGuide/` directory with `CompareMode0.txt`, `CompareMode1.txt`, `CompareMode2.txt` — identical geometry (21×21, pure advection, rotation flow), differing only in `reseed_mode`, all with `writer_markers=1`, `writer=5`, `log_level=4`, 50–100 steps
- [x] 2.2 Create `EdgeLowPart.txt` — low `min_part_cell=4`, high velocity to empty cells, `reseed_mode=2`, `writer_markers=1`
- [x] 2.3 Create `EdgeHighPart.txt` — high initial particles, low `min_part_cell` to trigger immediate deactivation, `reseed_mode=2`, `writer_markers=1`
- [x] 2.4 Create `EdgeMultiphase.txt` — 3-phase setup (two circular inclusions of different phases in a matrix) where interface cells contain multiple phases, `reseed_mode=2`, `writer_markers=1`

## 3. C++ Test Runner

- [x] 3.1 Create `VISUAL_TESTS/ReseedingGuide.cpp` with `SetPhase`, `SetHorizontalVelocity`, `SetVerticalVelocity` callbacks for rigid-body rotation (same as rotation advection test)
- [x] 3.2 Implement runner function that calls `mdoodz_run()` for each `.txt` config, writing output to `visualtests-out/ReseedingGuide/Mode0/`, `Mode1/`, `Mode2/`, `EdgeLowPart/`, `EdgeHighPart/`, `EdgeMultiphase/`
- [x] 3.3 Add `ReseedingGuide.cpp` to `VISUAL_TESTS/CMakeLists.txt` `add_executable(visualtests ...)` source list
- [x] 3.4 Add `configure_file()` calls for all `.txt` and `.gnu` files to copy to `visualtests-out/`
- [x] 3.5 Build and run: verify `visualtests` compiles and the ReseedingGuide runner produces particle HDF5 files in expected subdirectories

## 4. TSV Extractor

- [x] 4.1 Implement `extractParticlesForCell()` in `ReseedingGuide.cpp` using `HDF5pp.h` — reads particle HDF5, filters by cell bounds, writes TSV with columns: `x z phase generation index`
- [x] 4.2 Implement `extractNbPartTimeseries()` — scans output directory for all `Particles*.h5`, reads `Model/Nb_part` from each, writes TSV with columns: `step Nb_part`
- [x] 4.3 Wire extractors into the runner: after each simulation completes, extract TSV files for target cells and Nb_part timeseries
- [x] 4.4 Verify: run runner, confirm TSV files are produced and contain expected data

## 5. Gnuplot Scripts

- [x] 5.1 Create `advection_primer.gnu` — 2-panel PNG showing particles before/after a few advection steps (no reseeding), illustrating cell depletion and accumulation as motivation for reseeding
- [x] 5.2 Create `reseed_placement.gnu` — 3-panel PNG, plots `generation==1` particles from modes 0/1/2 in the same depleted cell, with fine-mesh cell boundary overlay
- [x] 5.3 Create `deactivation_strategy.gnu` — 3-panel PNG, plots surviving vs deactivated particles in overpopulated cell for modes 0/1/2, color-coding `phase==-1` particles, with fine-mesh cell boundary overlay
- [x] 5.4 Create `nbpart_timeseries.gnu` — single-axes PNG, 3 curves (one per mode) showing Nb_part over time
- [x] 5.5 Create `slot_usage.gnu` — PNG showing array index usage: append-only (modes 0/1) vs recycled slots (mode 2)
- [x] 5.6 Create `edge_empty_cell.gnu` — before/after panels showing bulk reseeding response in depleted region
- [x] 5.7 Create `edge_multiphase.gnu` — zoomed cell view showing dominant-phase selection for new particles across 3 phases
- [x] 5.8 Ensure all scatter plots overlay 2×Ncx × 2×Ncz fine-mesh cell grid lines
- [x] 5.9 Verify: run all `.gnu` scripts against extracted TSV data, confirm PNGs are produced

## 6. Markdown Guide

- [x] 6.1 Create `VISUAL_TESTS/ReseedingGuide.md` with overview table: mode number, function name, placement strategy, deactivation, recycling
- [x] 6.2 Add section 0 (Why Reseeding?): embed `advection_primer.png`, brief annotation showing particle drift causes cell depletion/accumulation, referencing `particles.x`, `particles.z` advection
- [x] 6.3 Add section 1 (Reseeding Placement): embed `reseed_placement.png`, code annotation referencing `CountPartCell`/`CountPartCell_v2`, variables `particles.x`, `particles.z`, `particles.generation`
- [x] 6.4 Add section 2 (Deactivation Strategy): embed `deactivation_strategy.png`, code annotation referencing `CountPartCell_v2` deactivation block, `PartDist[]`, `particles.phase`, `min_part_cell + 4` threshold
- [x] 6.5 Add section 3 (Particle Count Over Time): embed `nbpart_timeseries.png`, code annotation referencing `Nb_part`, `Nb_part_max`, growth vs stability
- [x] 6.6 Add section 4 (Array Slot Usage): embed `slot_usage.png`, code annotation referencing `part_reuse[]`, `phase == -1` recycling, append-only vs reuse
- [x] 6.7 Add edge-case sections (Empty Cell, Multiphase 3-phase): embed corresponding PNGs, code annotations referencing `min_part_cell`, `nb_part_cell[]`, dominant phase logic
- [x] 6.8 Add annotated log excerpt section: paste selected `LOG_DBG` output from a run with `log_level=4`, showing algorithm decision flow for one cell through reseeding+deactivation
- [x] 6.9 Review: verify all images render, fine-mesh grid lines visible on scatter plots, all code references use inline code formatting, no section exceeds 3 lines of text after the image
