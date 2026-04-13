## ADDED Requirements

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

#### Scenario: Invalid interp_mode value
- **WHEN** the `.txt` file contains `interp_mode` with a value other than 0, 1, or 2
- **THEN** the system SHALL treat it as `interp_mode = 0` (legacy)

### Requirement: Persistent buffer pool lifecycle
When `interp_mode >= 1`, the system SHALL allocate an `InterpBufPool` once during simulation initialization and free it during simulation cleanup. The pool SHALL NOT be allocated when `interp_mode == 0`.

#### Scenario: Pool allocated at startup
- **WHEN** `interp_mode >= 1` and the grid has been initialized
- **THEN** the system SHALL allocate an `InterpBufPool` containing work arrays for all four centroid grid sizes (cell centres, vertices, Vx nodes, Vz nodes) before the first `P2Mastah` call

#### Scenario: Pool freed at cleanup
- **WHEN** the simulation completes or is interrupted
- **THEN** the system SHALL free all memory held by the `InterpBufPool`

#### Scenario: No pool in legacy mode
- **WHEN** `interp_mode == 0`
- **THEN** the system SHALL pass `NULL` as the pool pointer to `P2Mastah` and no persistent buffers SHALL be allocated

### Requirement: Pre-allocated buffers for all centroid types
The `InterpBufPool` SHALL pre-allocate separate buffer sets for each of the four centroid types used by `P2Mastah`: cell centres (`centroid=1`, size (Nx-1)×(Nz-1)), vertices (`centroid=0`, size Nx×Nz), Vx nodes (`centroid=-1`, size Nx×(Nz+1)), and Vz nodes (`centroid=-2`, size (Nx+1)×Nz).

#### Scenario: Correct buffer size per centroid type
- **WHEN** `P2Mastah` is called with a given `centroid` value and `interp_mode >= 1`
- **THEN** the pool SHALL provide buffers of exactly `Nx_c × Nz_c` doubles for that centroid type, where `Nx_c` and `Nz_c` are the grid dimensions corresponding to the centroid

#### Scenario: Buffer reuse across calls
- **WHEN** `P2Mastah` is called multiple times with the same centroid type
- **THEN** the same pre-allocated buffer SHALL be reused (no new allocation) and SHALL be zeroed with `memset` before each use

### Requirement: Mode 1 — persistent buffers with thread-local scatter
When `interp_mode == 1`, `P2Mastah` SHALL use persistent thread-local buffers from the pool instead of per-call `DoodzCalloc`/`DoodzFree`, while retaining the same scatter and reduction logic as mode 0.

#### Scenario: Thread-local buffers from pool
- **WHEN** `interp_mode == 1` and `P2Mastah` is called
- **THEN** the function SHALL use `pool->Wm[centroid_idx][thread_num]` and `pool->BmWm[centroid_idx][thread_num]` instead of allocating new thread-local arrays

#### Scenario: Reduction loop retained
- **WHEN** `interp_mode == 1`
- **THEN** the function SHALL execute the same thread-local-to-global reduction loop as mode 0

#### Scenario: Bit-identical to legacy
- **WHEN** `interp_mode == 1`
- **THEN** the interpolated node values SHALL be bit-identical to those produced by `interp_mode == 0` for the same input

#### Scenario: Wm_ph allocated per-call for prop==1
- **WHEN** `interp_mode == 1` and `prop == 1` (phase percentage interpolation)
- **THEN** the `Wm_ph` thread-local arrays SHALL be allocated and freed per-call as in mode 0 (not pooled)

### Requirement: Mode 2 — atomic scatter to global arrays
When `interp_mode == 2`, `P2Mastah` SHALL accumulate particle contributions directly into global `WM[]` and `BMWM[]` arrays using `#pragma omp atomic`, eliminating thread-local buffer copies and the reduction loop.

#### Scenario: Atomic accumulation in scatter loop
- **WHEN** `interp_mode == 2` and the particle scatter loop executes
- **THEN** each `+=` to `WM[kp]` and `BMWM[kp]` SHALL use `#pragma omp atomic` to ensure thread safety

#### Scenario: No thread-local Wm/BmWm arrays
- **WHEN** `interp_mode == 2`
- **THEN** the function SHALL NOT allocate or use thread-indexed `Wm[thread_num][]` or `BmWm[thread_num][]` arrays

#### Scenario: Reduction loop eliminated
- **WHEN** `interp_mode == 2`
- **THEN** the function SHALL NOT execute the thread-local-to-global reduction loop (the global arrays already contain the accumulated result)

#### Scenario: Results match to machine epsilon
- **WHEN** `interp_mode == 2`
- **THEN** the interpolated node values SHALL match those produced by `interp_mode == 0` within floating-point machine epsilon (the order of atomic `+=` operations may differ from the sequential reduction)

### Requirement: Mode 2 phase-percentage atomic scatter
When `interp_mode == 2` and `prop == 1`, `P2Mastah` SHALL apply atomic `+=` directly to `mesh->phase_perc_s[p][]` or `mesh->phase_perc_n[p][]` instead of using thread-local `Wm_ph` arrays.

#### Scenario: Direct atomic to mesh phase arrays
- **WHEN** `interp_mode == 2` and `prop == 1` and the scatter loop processes a particle of phase `p`
- **THEN** the weight contribution SHALL be accumulated with `#pragma omp atomic` directly into `mesh->phase_perc_s[p][kp]` (for `centroid==0`) or `mesh->phase_perc_n[p][kp]` (for `centroid==1`)

#### Scenario: No Wm_ph allocation in mode 2
- **WHEN** `interp_mode == 2` and `prop == 1`
- **THEN** no `Wm_ph` thread-local arrays SHALL be allocated

### Requirement: P2Mastah function signature accepts pool
The `P2Mastah` function signature SHALL include an `InterpBufPool *pool` parameter to receive the persistent buffer pool.

#### Scenario: Pool parameter passed from callers
- **WHEN** any caller invokes `P2Mastah`
- **THEN** the caller SHALL pass the `InterpBufPool *pool` pointer (or `NULL` when `interp_mode == 0`)

#### Scenario: CountPartCell passes pool
- **WHEN** `CountPartCell` calls `P2Mastah`
- **THEN** it SHALL pass the same `InterpBufPool *pool` pointer received from `Main_DOODZ.c`

### Requirement: Backward compatibility
The existing interpolation code path SHALL remain unchanged when `interp_mode == 0`. All existing `.txt` files SHALL produce identical results without modification.

#### Scenario: Existing simulations unaffected
- **WHEN** a `.txt` file does not contain `interp_mode`
- **THEN** the simulation SHALL behave identically to the code before this change (same allocations, same scatter, same reduction, same results)

#### Scenario: Legacy code path unchanged
- **WHEN** `interp_mode == 0`
- **THEN** `P2Mastah` SHALL execute the same `DoodzCalloc`/`DoodzFree` pairs, thread-local scatter, and reduction loop as the current code with no persistent buffers used

### Requirement: OpenMP parallel correctness
All `P2Mastah` modes SHALL produce correct results regardless of the OpenMP thread count.

#### Scenario: Single-thread correctness
- **WHEN** `omp_nthread = 1` and any `interp_mode`
- **THEN** `P2Mastah` SHALL produce the same result as with multiple threads

#### Scenario: Mode 2 atomic safety
- **WHEN** `interp_mode == 2` and multiple threads execute the scatter loop
- **THEN** all `#pragma omp atomic` operations SHALL guarantee no torn or lost writes

### Requirement: CI test coverage for interp modes
The CI test suite SHALL include `interp_mode = 1` and `interp_mode = 2` variants for all existing test suites except BlankenBench.

#### Scenario: Mode 1 CI tests pass with exact match
- **WHEN** any CI test suite (except BlankenBench) runs with `interp_mode = 1`
- **THEN** results SHALL be bit-identical to the same test with `interp_mode = 0`

#### Scenario: Mode 2 CI tests pass within tolerance
- **WHEN** any CI test suite (except BlankenBench) runs with `interp_mode = 2`
- **THEN** results SHALL match the `interp_mode = 0` baseline within a relative tolerance of 1e-10

#### Scenario: BlankenBench excluded
- **WHEN** CI tests are configured
- **THEN** the BlankenBench test suite SHALL NOT have `interp_mode` variants (too long-running for CI)
