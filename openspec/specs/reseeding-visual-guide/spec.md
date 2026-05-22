## Purpose

Visual test suite for particle reseeding mechanisms across three operational modes (0, 1, 2). Demonstrates reseeding placement strategies, deactivation approaches, and array management via interactive gnuplot comparisons.

## ADDED Requirements

### Requirement: Test runner produces particle HDF5 for all 3 modes
The `ReseedingGuide.cpp` runner SHALL execute MDOODZ simulations for reseeding modes 0, 1, and 2 using identical geometry and flow parameters, with `writer_markers=1` enabled, writing `Particles<NNNNN>.gzip.h5` output to mode-specific subdirectories under `visualtests-out/ReseedingGuide/`.

#### Scenario: Run comparison scenarios
- **WHEN** the `visualtests` binary is executed with the ReseedingGuide runner
- **THEN** it SHALL produce particle HDF5 files in `visualtests-out/ReseedingGuide/Mode0/`, `Mode1/`, and `Mode2/` directories from the same initial configuration

#### Scenario: Particle output contains required fields
- **WHEN** a `Particles<NNNNN>.gzip.h5` file is produced
- **THEN** it SHALL contain `Particles/x`, `Particles/z`, `Particles/phase`, `Particles/generation`, and `Model/Nb_part`

### Requirement: Edge-case scenarios trigger boundary conditions
The runner SHALL include configurations that force specific reseeding/deactivation edge cases: cells emptied below `min_part_cell`, cells overpopulated beyond `min_part_cell + 4`, and `Nb_part` approaching `Nb_part_max`.

#### Scenario: Cell depletion triggers reseeding
- **WHEN** an edge-case configuration with high velocity and low `min_part_cell` is run
- **THEN** at least one cell SHALL drop below the reseeding threshold, causing new particles to appear in subsequent output

#### Scenario: Cell overpopulation triggers deactivation (mode 2)
- **WHEN** an edge-case configuration is run with `reseed_mode=2`
- **THEN** cells exceeding `min_part_cell + 4` SHALL have excess particles marked `phase = -1` in the output

#### Scenario: Multi-phase cell reseeding
- **WHEN** a cell containing 3 phases drops below the reseeding threshold
- **THEN** new particles SHALL be assigned the dominant phase, visible by filtering `generation == 1` particles in the output, with all 3 phases distinguishable by color

### Requirement: TSV extractor converts HDF5 to gnuplot-readable format
A C++ extractor (using `HDF5pp.h`) SHALL read particle HDF5 files and write TSV files containing per-cell particle coordinates, phase, generation, and array indices for gnuplot consumption.

#### Scenario: Extract particles for a specific cell
- **WHEN** the extractor is given an HDF5 file and cell indices (i, j)
- **THEN** it SHALL produce a TSV file with columns: `x z phase generation index` for all particles within that cell's bounds

#### Scenario: Extract Nb_part time series
- **WHEN** the extractor is given an output directory with multiple particle HDF5 files
- **THEN** it SHALL produce a TSV file with columns: `step Nb_part` for each output step

### Requirement: Gnuplot scripts produce one PNG per concept
Each gnuplot `.gnu` script SHALL produce exactly one PNG file focused on a single concept. Scripts SHALL read TSV data (not HDF5 directly) for portability.

#### Scenario: Reseeding placement comparison
- **WHEN** `reseed_placement.gnu` is executed with TSV data from all 3 modes
- **THEN** it SHALL produce a 3-panel PNG showing new particle positions (filtered by `generation == 1`) in the same depleted cell for modes 0, 1, and 2, with fine-mesh cell boundaries (2×Ncx × 2×Ncz grid lines) overlaid

#### Scenario: Deactivation strategy comparison
- **WHEN** `deactivation_strategy.gnu` is executed with TSV data from all 3 modes
- **THEN** it SHALL produce a 3-panel PNG showing which particles survive in an overpopulated cell — highlighting that mode 0 has no deactivation, mode 1 removes by index order, mode 2 removes farthest from centroid — with fine-mesh cell boundaries overlaid

#### Scenario: Particle count time series
- **WHEN** `nbpart_timeseries.gnu` is executed with Nb_part TSV data from all 3 modes
- **THEN** it SHALL produce a single-axes PNG with 3 curves showing Nb_part over time — modes 0/1 growing, mode 2 stable

#### Scenario: Array slot usage comparison
- **WHEN** `slot_usage.gnu` is executed
- **THEN** it SHALL produce a PNG showing append-only growth (modes 0/1) vs recycled slot reuse (mode 2)

#### Scenario: Advection primer
- **WHEN** `advection_primer.gnu` is executed with before/after advection TSV data
- **THEN** it SHALL produce a 2-panel PNG showing particle positions before and after a few advection steps, illustrating cell depletion on one side and accumulation on the other

### Requirement: Markdown guide follows plot-first structure
`VISUAL_TESTS/ReseedingGuide.md` SHALL embed each generated PNG as an image, followed by at most 3 lines of code reference text per section. The guide SHALL NOT contain prose paragraphs.

#### Scenario: Guide structure
- **WHEN** the guide is opened
- **THEN** each section SHALL consist of: a section heading, one `![](image.png)` embed, and a **Code:** annotation naming the relevant function, line range, and variable names from `ParticleRoutines.c`

#### Scenario: Mode overview table
- **WHEN** the guide's overview section is read
- **THEN** it SHALL contain a comparison table (not prose) summarising modes 0/1/2 with columns: mode number, function name, placement strategy, deactivation, recycling

#### Scenario: Edge-case sections
- **WHEN** the guide's edge-case sections are read
- **THEN** each edge case (empty cell, overpopulated cell, array-full approach, multi-phase) SHALL have its own section with a dedicated plot and code reference

### Requirement: Guide references code variables and array semantics
Each plot section in the guide SHALL reference the specific `ParticleRoutines.c` variables involved, using inline code formatting.

#### Scenario: Variable references present
- **WHEN** any plot section is read
- **THEN** it SHALL reference at least one of: `particles.x`, `particles.z`, `particles.phase`, `particles.generation`, `Nb_part`, `Nb_part_max`, `min_part_cell`, `part_reuse[]`, `nb_part_cell[]`

### Requirement: All particle scatter plots overlay fine-mesh cell boundaries
Every gnuplot script that produces a particle scatter plot SHALL overlay the 2×Ncx × 2×Ncz fine-mesh cell grid lines, showing the spatial bins the reseeding algorithm operates on.

#### Scenario: Fine-mesh grid visible
- **WHEN** any particle scatter PNG is viewed
- **THEN** the fine-mesh cell boundaries SHALL be visible as grid lines behind the particle markers

### Requirement: Build integration
`ReseedingGuide.cpp` SHALL be added to the `VISUAL_TESTS/CMakeLists.txt` `add_executable(visualtests ...)` source list. All `.txt` parameter files and `.gnu` scripts SHALL be copied to `visualtests-out/` via `configure_file()`.

#### Scenario: Build succeeds
- **WHEN** `cmake --build` is run with the VISUAL_TESTS target enabled
- **THEN** the `visualtests` binary SHALL compile without errors including the new ReseedingGuide runner

#### Scenario: Config files deployed
- **WHEN** the build completes
- **THEN** all `ReseedingGuide/*.txt` and `*.gnu` files SHALL be present in `visualtests-out/`
