## Context

MDOODZ has three particle reseeding modes (`reseed_mode` 0/1/2) implemented in `MDLIB/ParticleRoutines.c` as `CountPartCell_OLD`, `CountPartCell`, and `CountPartCell_v2`. These functions share a common pattern — build a fine mesh (2×Ncx × 2×Ncz), count particles per fine cell, inject new particles where cells are depleted, optionally deactivate excess — but differ in placement strategy, deactivation logic, and array management.

There is no visual documentation of these differences. The existing rotation advection test (`TESTS/RotationAdvectionTests.cpp`) already runs all 3 modes on identical geometry and has Python plotting (`plot_rotation.py`), but focuses on L2 accuracy, not reseeding mechanics.

The `VISUAL_TESTS/` directory has an established pattern: `.cpp` scenario runners compiled into the `visualtests` executable, `.gnu` gnuplot scripts for plotting, `.txt` parameter files, and reference HDF5 data. All output goes to `visualtests-out/`.

Particle HDF5 output (`Particles<NNNNN>.gzip.h5`) writes `x`, `z`, `phase`, `generation`, `Vx`, `Vz`, `T`, `P`, `sxxd`, `sxz` under `Particles/` group and `Nb_part` under `Model/`. Dead particles (`phase == -1`) are included in the output arrays.

## Goals / Non-Goals

**Goals:**
- Visual guide that makes mode 0/1/2 differences immediately obvious through side-by-side plots from real particle data
- One concept per plot — never overload a figure
- Cover normal operation (reseeding, deactivation) AND edge cases (Nb_part_max overflow, empty cells, multi-phase)
- Code reference: each section maps plots to `ParticleRoutines.c` variables, array names, and index semantics
- DEBUG logging in MDLIB to narrate algorithm steps for annotated log excerpts in the guide
- Production-ready gnuplot pipeline; Python acceptable during development

**Non-Goals:**
- Changing reseeding algorithm behavior
- Automated CI regression testing of reseeding (that's a separate concern)
- Julia/Makie visualisation (gnuplot is the target)
- Documenting advection accuracy (rotation advection test already covers that)

## Decisions

### 1. Test runner: new .cpp file added to VISUAL_TESTS/visualtests executable

**Choice:** Add `ReseedingGuide.cpp` to the existing `VISUAL_TESTS/CMakeLists.txt` `add_executable(visualtests ...)` list, following the same pattern as `RiftingChenin.cpp`, `ShearTemplate.cpp`, etc.

**Why not a separate standalone binary?** The `visualtests` executable already links `mdoodz`, Eigen, HDF5. Adding another `.cpp` to it is zero-friction. A separate binary would need its own CMake target and link lines for no benefit.

**The runner calls `mdoodz_run()` multiple times** with different `.txt` configs (one per mode, one per edge case), each writing particle HDF5 output to a subfolder under `visualtests-out/ReseedingGuide/`.

### 2. Scenario: purpose-built simple shear / rotation, not reusing existing tests directly

**Choice:** Create new `.txt` parameter files specifically tuned for visual clarity:
- **Comparison scenario** — small grid (e.g., 21×21), pure advection (`mechanical=0, thermal=0`), rigid-body rotation or simple shear. `writer_markers=1`, `writer=5` (frequent output). `min_part_cell=16`, `particles_per_cell_x=4/z=4`. Run 50–100 steps (not 1000 — just enough to see reseeding/deactivation happen).
- **Edge-case scenarios** — tweaked configs that force specific conditions:
  - Very low `min_part_cell=4` + high velocity → cells empty out → massive reseeding
  - Very high initial particles + low `min_part_cell` → deactivation triggers immediately
  - Small `Nb_part_max` multiplier → approach array limit quickly

**Why not reuse RotationAdvection .txt directly?** Those run 1000 steps for L2 accuracy. We need short runs (few steps) with frequent particle output. Different `min_part_cell` values to trigger edge cases. Separate configs are cleaner.

**Alternatives considered:**
- Reuse rotation advection: too many steps, wrong `writer` frequency, no edge-case coverage
- Create a new SETS/ scenario: overkill, SETS/ is for full simulations

### 3. One .txt per mode per scenario, named explicitly

```
VISUAL_TESTS/ReseedingGuide/
  CompareMode0.txt    # reseed_mode=0, writer_markers=1
  CompareMode1.txt    # reseed_mode=1, writer_markers=1
  CompareMode2.txt    # reseed_mode=2, writer_markers=1
  EdgeLowPart.txt     # min_part_cell=4, high velocity
  EdgeHighPart.txt    # extra particles, low threshold
  EdgeArrayFull.txt   # small Nb_part_max multiplier
```

The `.cpp` runner iterates over these and writes output to mode-specific subdirectories.

### 4. Gnuplot scripts: one .gnu per plot concept

**Choice:** Separate `.gnu` files, each producing one PNG:

| Script | Reads | Produces | Shows |
|--------|-------|----------|-------|
| `advection_primer.gnu` | particle HDF5 (before/after advection) | 2-panel PNG | Why reseeding is needed: particle drift, cell depletion/accumulation |
| `reseed_placement.gnu` | 3 particle HDF5 (modes 0/1/2) | 3-panel PNG | New particle positions in a depleted cell |
| `deactivation_strategy.gnu` | 3 particle HDF5 | 3-panel PNG | Which particles survive in overpopulated cell |
| `nbpart_timeseries.gnu` | `perf.csv` or extracted counts | 1 plot, 3 curves | Nb_part growth over time |
| `slot_usage.gnu` | particle HDF5 (phase array) | 2-panel PNG | Array slot occupancy: append-only vs recycled |
| `edge_empty_cell.gnu` | edge-case HDF5 | before/after panels | Bulk reseeding response |
| `edge_array_full.gnu` | edge-case HDF5 + log | annotated plot | Approach to Nb_part_max |
| `edge_multiphase.gnu` | edge-case HDF5 | cell zoom | Phase selection for new particles |

**HDF5 in gnuplot:** Gnuplot 5.4+ can read HDF5 natively with `plot 'file.h5' using ...`. If the build environment has an older gnuplot, a thin C++ HDF5-to-TSV extractor (using the existing `HDF5pp.h` in VISUAL_TESTS) pre-processes data into tab-separated files that gnuplot reads.

**Why gnuplot over Python/matplotlib?** User requirement: C++ + gnuplot is production. Python is fine during development but the deliverable scripts must be gnuplot. This matches the existing VISUAL_TESTS pattern (all existing .gnu files).

### 5. Markdown guide structure: plot-first, code-reference annotations

```markdown
# Particle Reseeding Visual Guide

## Overview
(3-sentence summary of the 3 modes — table only, no prose)

## 0. Why Reseeding?
![](advection_primer.png)
**Code:** Marker-in-cell advection moves `particles.x[k]`, `particles.z[k]` each step.
Cells lose particles on one side, gain on the other → some cells empty, some overpopulate.

## 1. Reseeding Placement
![](reseed_placement.png)
**Code:** `CountPartCell` (mode 1, line N) places at fine-cell centroid.
`CountPartCell_v2` (mode 2, line M) adds random offset within cell.
**Variables:** `particles.x[k]`, `particles.z[k]`, `particles.generation[k]=1`

## 2. Deactivation Strategy
![](deactivation_strategy.png)
**Code:** Mode 2 only (`CountPartCell_v2`, line P). Builds `PartDist[]`,
sorts by distance from centroid descending, marks farthest as `phase=-1`.
**Variables:** `particles.phase[k]`, `min_part_cell`, threshold = `min_part_cell + 4`

## 3. Particle Count Over Time
![](nbpart_timeseries.png)
...

(etc.)
```

Each section: one plot image, 2–3 lines of code reference with function name + line number + variable names. No paragraphs.

### 6. DEBUG logging: targeted LOG_DBG in CountPartCell / CountPartCell_v2

**Choice:** Add ~10 `LOG_DBG()` calls at key decision points in `ParticleRoutines.c`:

| Location | Log message |
|----------|-------------|
| After counting cell particles | `"cell (%d,%d) fine: %d particles (threshold=%d)"` |
| When reseeding triggers | `"cell (%d,%d): reseeding %d new particles"` |
| New particle placement | `"  new particle at (%.4f, %.4f), reusing slot %d"` |
| New particle appended | `"  new particle at (%.4f, %.4f), appended at index %d"` |
| Deactivation trigger (mode 2) | `"cell (%d,%d): %d excess, deactivating %d farthest"` |
| Particle deactivated | `"  deactivated particle %d (dist=%.4f from centroid)"` |
| Nb_part update | `"Nb_part: %d → %d (reused=%d, created=%d)"` |
| Array full guard | `"WARNING: Nb_part %d approaching Nb_part_max %d"` |

**Why LOG_DBG, not LOG_INFO?** Per-particle logging is extremely verbose. `log_level=4` (DEBUG) ensures zero output in normal runs. The visual guide's `.txt` files set `log_level=4` to capture these traces.

**Output:** Selected log lines appear as annotated code blocks in the guide, showing the algorithm's decision flow for a specific cell.

### 7. Data extraction: C++ helper using HDF5pp.h

**Choice:** A small helper function in `ReseedingGuide.cpp` that reads particle HDF5 and writes TSV files for gnuplot:

```cpp
void extractParticlesForCell(const char* h5file, int cellI, int cellJ,
                              double dx, double dz, const char* outTsv);
void extractNbPartTimeseries(const char* outputDir, int nSteps, const char* outTsv);
```

This uses the existing `HDF5pp.h` wrapper already in `VISUAL_TESTS/`. The extractor runs after the simulation, before gnuplot.

**Why not gnuplot's native HDF5?** Not all gnuplot builds have HDF5 support. TSV is universally readable and easy to inspect.

## Risks / Trade-offs

**[Verbose particle output for short runs]** → `writer_markers=1` writes full particle arrays every output step. For small grids (21×21) and short runs (50–100 steps), file sizes are negligible (~KB each). Mitigated by keeping grid small.

**[Gnuplot HDF5 portability]** → Not all gnuplot versions support HDF5. Mitigated by the C++ TSV extractor approach — gnuplot reads plain text files.

**[LOG_DBG per-particle noise]** → Even at DEBUG level, per-particle logs can produce thousands of lines on larger grids. Mitigated by small grids in visual test configs AND by logging only in cells that trigger reseeding/deactivation (not all cells).

**[Guide maintenance]** → Line numbers in code references will drift as `ParticleRoutines.c` changes. Mitigated by referencing function names + variable names primarily, with line numbers as secondary hints.

**[Edge case: exit(190) is fatal]** → Testing the array-full edge case that triggers `exit(190)` would kill the test runner. Mitigated by configuring `Nb_part_max` just large enough that we approach the limit but add enough particles via aggressive reseeding WITHOUT actually hitting it. The guide documents the exit behavior from code inspection + log output rather than triggering it.

## Resolved Questions

- **Fine-mesh cell boundaries on scatter plots?** → **Yes.** Overlay the 2×Ncx × 2×Ncz fine-cell grid lines on all particle scatter plots. These boundaries show the actual spatial bins the algorithm uses and help the reader coordinate particle positions with reseeding/deactivation decisions.
- **Multi-phase edge case: 2-phase or 3-phase?** → **3-phase.** Shows dominant-phase selection more clearly (majority phase wins). Won't overload the plot if represented with distinct colors per phase.
- **Advection primer?** → **Yes.** Add a short introductory section that explains *why* reseeding is needed: show a single particle advecting through the grid over a few steps, cells emptying on one side and accumulating on the other. This is not the focus — it's the instrument to motivate the problem. One plot: particle trajectory + cell occupancy before/after a few advection steps.
