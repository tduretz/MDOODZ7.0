## MODIFIED Requirements

### Requirement: Interpolation mode selection via .txt parameter
The system SHALL read an `interp_mode` integer parameter from the `.txt` input file to select the particle-to-grid interpolation strategy used by `P2Mastah`.

#### Scenario: Default is legacy mode
- **WHEN** the `.txt` file does not contain `interp_mode`
- **THEN** the system SHALL default to `interp_mode = 0` and the `P2Mastah` function SHALL execute the existing code path unchanged (per-call buffer allocation, thread-local scatter, explicit reduction)

#### Scenario: User selects persistent buffers
- **WHEN** the `.txt` file contains `interp_mode = 1`
- **THEN** the system SHALL use pre-allocated persistent buffers for `P2Mastah` with thread-local scatter and explicit reduction

#### Scenario: User selects persistent buffers with atomic scatter
- **WHEN** the `.txt` file contains `interp_mode = 2`
- **THEN** the system SHALL use pre-allocated persistent buffers with atomic scatter to global arrays and no explicit reduction loop

#### Scenario: User selects fused atomic scatter
- **WHEN** the `.txt` file contains `interp_mode = 3`
- **THEN** the system SHALL use fused multi-field atomic scatter via `P2Mastah_Fused` for fusible call groups and mode 2 atomic scatter for non-fusible calls

#### Scenario: Invalid interp_mode value
- **WHEN** the `.txt` file contains `interp_mode` with a value other than 0, 1, 2, or 3
- **THEN** the system SHALL treat it as `interp_mode = 0` (legacy)

### Requirement: Persistent buffer pool lifecycle
When `interp_mode >= 1`, the system SHALL allocate an `InterpBufPool` once during simulation initialization and free it during simulation cleanup. The pool SHALL NOT be allocated when `interp_mode == 0`.

#### Scenario: Pool allocated at startup
- **WHEN** `interp_mode >= 1` and the grid has been initialized
- **THEN** the system SHALL allocate an `InterpBufPool` containing work arrays for all four centroid grid sizes (cell centres, vertices, Vx nodes, Vz nodes) before the first `P2Mastah` or `P2Mastah_Fused` call

#### Scenario: Pool freed at cleanup
- **WHEN** the simulation completes or is interrupted
- **THEN** the system SHALL free all memory held by the `InterpBufPool`, including any per-field BMWM slots allocated for mode 3

#### Scenario: No pool in legacy mode
- **WHEN** `interp_mode == 0`
- **THEN** the system SHALL pass `NULL` as the pool pointer to `P2Mastah` and no persistent buffers SHALL be allocated

### Requirement: Pre-allocated buffers for all centroid types
The `InterpBufPool` SHALL pre-allocate separate buffer sets for each of the four centroid types used by `P2Mastah`: cell centres (`centroid=1`, size (Nx-1)×(Nz-1)), vertices (`centroid=0`, size Nx×Nz), Vx nodes (`centroid=-1`, size Nx×(Nz+1)), and Vz nodes (`centroid=-2`, size (Nx+1)×Nz). When `interp_mode == 3`, the pool SHALL additionally pre-allocate per-field `BMWM` arrays for fused batches.

#### Scenario: Correct buffer size per centroid type
- **WHEN** `P2Mastah` or `P2Mastah_Fused` is called with a given `centroid` value and `interp_mode >= 1`
- **THEN** the pool SHALL provide buffers of exactly `Nx_c × Nz_c` doubles for that centroid type, where `Nx_c` and `Nz_c` are the grid dimensions corresponding to the centroid

#### Scenario: Buffer reuse across calls
- **WHEN** `P2Mastah` or `P2Mastah_Fused` is called multiple times with the same centroid type
- **THEN** the same pre-allocated buffer SHALL be reused (no new allocation) and SHALL be zeroed with `memset` before each use

#### Scenario: Per-field BMWM slots for mode 3
- **WHEN** `interp_mode == 3`
- **THEN** the pool SHALL pre-allocate a fixed number of `BMWM` slots per centroid type (sufficient for the largest fused batch, at least 32 per centroid) to avoid per-call allocation in `P2Mastah_Fused`
