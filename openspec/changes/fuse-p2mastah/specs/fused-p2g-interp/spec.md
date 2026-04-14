## ADDED Requirements

### Requirement: Field descriptor for fused interpolation
The system SHALL define a `P2MastahField` struct that describes a single field to scatter during a fused interpolation pass. The descriptor SHALL contain: a source array pointer (`double *src`), a destination grid array pointer (`double *dst`), a boundary condition type array pointer (`char *BCtype`), a flag value (`int flag`: 0=material property, 1=particle field, -1/-2/-3=anisotropy), a property index (`int prop`), and an interpolation stencil width (`int stencil`).

#### Scenario: Descriptor fully specifies one field
- **WHEN** a `P2MastahField` struct is populated
- **THEN** it SHALL contain all the information that a single `P2Mastah` call requires beyond model, particles, mesh, interp_mode, centroid, and pool — i.e. `src`, `dst`, `BCtype`, `flag`, `prop`, and `stencil`

#### Scenario: Descriptors for different flag types
- **WHEN** building a fused batch that includes material property fields (flag=0), particle fields (flag=1), and anisotropy fields (flag=-1/-2/-3)
- **THEN** each descriptor SHALL carry its own `flag` value and the fused function SHALL handle all flag types within a single particle loop

### Requirement: P2Mastah_Fused function
The system SHALL provide a `P2Mastah_Fused()` function that performs particle-to-grid interpolation for multiple fields in a single pass over the particle array per centroid type. It SHALL compute grid cell indices and bilinear weights once per particle and scatter all fields in the batch using those weights.

#### Scenario: Single particle loop for N fields
- **WHEN** `P2Mastah_Fused` is called with an array of N `P2MastahField` descriptors and a centroid type
- **THEN** the function SHALL iterate over all particles exactly once, and for each particle compute the grid indices and weights once and scatter weighted values into all N destination grids

#### Scenario: Atomic scatter for all fields
- **WHEN** `P2Mastah_Fused` executes (mode 3)
- **THEN** all accumulations into shared `WM` and `BMWM` arrays SHALL use `#pragma omp atomic` (same approach as mode 2) for each field in the batch

#### Scenario: Normalization per field
- **WHEN** the scatter phase completes for all particles
- **THEN** the function SHALL perform a normalization loop for each field: `dst[i] = BMWM_field[i] / WM[i]`, applying harmonic/geometric inversion as needed based on each field's `flag` and `prop` values

#### Scenario: WM sharing across fields within a centroid
- **WHEN** multiple fields share the same centroid type
- **THEN** the weight array `WM` SHALL be accumulated once (not per-field) since all fields at a given centroid use identical particle weights, and all fields' `BMWM` arrays SHALL be divided by this single shared `WM`

### Requirement: Fused call-site grouping
Consecutive `P2Mastah` calls in `Main_DOODZ.c` that share the same centroid type SHALL be replaceable by a single `P2Mastah_Fused` call with a field descriptor array, when `interp_mode == 3`.

#### Scenario: Centroid-1 initialization group
- **WHEN** `interp_mode == 3` and the initialization block executes (before the first timestep)
- **THEN** the centroid-1 (cell centres) calls for eta_n, noise_n, T, Cp SHALL be dispatched as a single `P2Mastah_Fused` call

#### Scenario: Centroid-0 initialization group
- **WHEN** `interp_mode == 3` and the initialization block executes
- **THEN** the centroid-0 (vertices) calls for eta_s, noise_s SHALL be dispatched as a single `P2Mastah_Fused` call

#### Scenario: Per-timestep groups
- **WHEN** `interp_mode == 3` and a timestep's interpolation block runs
- **THEN** each contiguous block of same-centroid-type P2Mastah calls SHALL be fused into a single `P2Mastah_Fused` call

#### Scenario: Mixed-centroid calls remain separate
- **WHEN** consecutive P2Mastah calls use different centroid types (e.g. Vx nodes followed by Vz nodes)
- **THEN** they SHALL be separate `P2Mastah_Fused` calls, each containing only fields for their centroid type

#### Scenario: Non-fusible calls fall back to P2Mastah
- **WHEN** `interp_mode == 3` and a P2Mastah call cannot be fused (e.g. standalone call inside a conditional, or `prop == 1` phase percentage calls)
- **THEN** it SHALL fall through to the existing `P2Mastah` mode 2 code path (atomic scatter without fusion)

### Requirement: Phase percentage calls excluded from fusion
Calls with `prop == 1` (phase percentage interpolation via `CountPartCell`) SHALL NOT be fused into `P2Mastah_Fused` batches. They SHALL continue to use the existing `P2Mastah` mode 2 code path.

#### Scenario: CountPartCell uses P2Mastah not P2Mastah_Fused
- **WHEN** `interp_mode == 3` and `CountPartCell` calls `P2Mastah` with `prop == 1`
- **THEN** the call SHALL use the mode 2 atomic scatter path of `P2Mastah`, not `P2Mastah_Fused`

### Requirement: Per-field BMWM buffer allocation in pool
When `interp_mode == 3`, the `InterpBufPool` SHALL provide per-field `BMWM` arrays for batched scatter. For each centroid type, the pool SHALL pre-allocate a configurable maximum number of `BMWM` slots (e.g. 32) to support the largest fused batch without per-call allocation.

#### Scenario: Pool provides per-field BMWM arrays
- **WHEN** `P2Mastah_Fused` is called with N fields for a given centroid type
- **THEN** the pool SHALL provide N separate `BMWM` arrays of the correct centroid grid size, plus one shared `WM` array

#### Scenario: Batch size within pool capacity
- **WHEN** a fused batch contains N fields and N exceeds the pool's pre-allocated slot count
- **THEN** the system SHALL log an error and fall back to individual `P2Mastah` mode 2 calls for the excess fields

### Requirement: Interp_mode 3 selection via .txt parameter
When `interp_mode = 3` is set in the `.txt` parameter file, the system SHALL use fused atomic scatter for all fusible P2Mastah call groups.

#### Scenario: Mode 3 selected
- **WHEN** the `.txt` file contains `interp_mode = 3`
- **THEN** the system SHALL dispatch fusible P2Mastah call groups via `P2Mastah_Fused` and non-fusible calls via `P2Mastah` mode 2

#### Scenario: Mode 3 single-thread correctness
- **WHEN** `interp_mode == 3` and `omp_nthread = 1`
- **THEN** results SHALL match `interp_mode == 0` to within floating-point machine epsilon (same as mode 2)

### Requirement: Results match mode 2 within machine epsilon
When `interp_mode == 3`, all interpolated grid values SHALL match those produced by `interp_mode == 2` within floating-point machine epsilon for the same input. The fused path SHALL NOT change the mathematical computation — only the loop structure differs.

#### Scenario: Fused vs individual atomic scatter equivalence
- **WHEN** `interp_mode == 3` produces field values and `interp_mode == 2` produces the same field values for the same simulation state
- **THEN** the values SHALL be bit-identical (both use atomic scatter in the same order within each field)

### Requirement: CI test for mode 3 via BlankenBench with all optimizations
The CI suite SHALL include a single BlankenBench test fixture that runs with `interp_mode = 3`, `thermal_solver = 1`, and all other performance optimizations enabled. This fixture validates that the full optimized code path produces correct results.

#### Scenario: BlankenBench optimized fixture
- **WHEN** the CI suite runs
- **THEN** it SHALL include a BlankenBench test with `interp_mode = 3` and `thermal_solver = 1` that completes and produces results within the BlankenBench acceptance tolerance

#### Scenario: One BlankenBench per CI run
- **WHEN** the CI suite configures BlankenBench tests
- **THEN** there SHALL be exactly one BlankenBench fixture with all optimizations enabled, not per-mode variants

### Requirement: EC2 benchmark validation
The fused interpolation mode SHALL be validated via an EC2 highres (1001×801) benchmark sweep on c5ad.4xlarge, comparing `interp_mode = 3` against `interp_mode = 2` (Experiment 6 baseline).

#### Scenario: Benchmark sweep executed
- **WHEN** the implementation is complete
- **THEN** a benchmark sweep SHALL be run with `--resolutions "highres" --threads "1 2 4 6 8 12 16" --steps 10 --sed 's/^thermal_solver.*/thermal_solver = 1/' --sed 's/^interp_mode.*/interp_mode = 3/'`

#### Scenario: Performance target met
- **WHEN** the benchmark results are collected at 16 threads
- **THEN** `avg_interp_s` SHALL be less than 1.5 seconds (reduction from ~2.56s baseline)

#### Scenario: No regression at low thread counts
- **WHEN** the benchmark results show 1-thread performance
- **THEN** `avg_wall_s` SHALL not regress by more than 10% compared to `interp_mode == 2` at 1 thread

#### Scenario: Results documented as Experiment 7
- **WHEN** the benchmark completes
- **THEN** the results SHALL be documented in `skill-benchmarking/SKILL.md` as Experiment 7 with the standard format (hypothesis, implementation, results table, conclusion)
