## ADDED Requirements

### Requirement: Quick verification configuration exists
A file `SETS/RiftingCombinedYield_quick.txt` SHALL exist alongside the production `.txt`. It SHALL share the same physics, phases, and domain but with coarse resolution and minimal time steps for fast smoke-testing.

#### Scenario: Quick .txt has coarse resolution and few steps
- **WHEN** `RiftingCombinedYield_quick.txt` is parsed
- **THEN** `Nx=101`, `Nz=81`, and `Nt=10` are set, while all other physics switches and phase definitions are identical to the production `.txt`

### Requirement: Medium verification configuration exists
A file `SETS/RiftingCombinedYield_medium.txt` SHALL exist alongside the production `.txt`. It SHALL use the same coarse resolution as the quick config but with more time steps to test longer-term stability.

#### Scenario: Medium .txt has coarse resolution and moderate steps
- **WHEN** `RiftingCombinedYield_medium.txt` is parsed
- **THEN** `Nx=101`, `Nz=81`, and `Nt=100` are set, while all other physics switches and phase definitions are identical to the production `.txt`

### Requirement: Quick verification runs without divergence
The `RiftingCombinedYield` executable SHALL complete all time steps when run with the quick `.txt` without producing NaN values or solver divergence.

#### Scenario: Quick smoke-test completes successfully
- **WHEN** the user runs `./RiftingCombinedYield` with the quick `.txt` copied as `RiftingCombinedYield.txt`
- **THEN** the simulation completes 10 time steps, produces HDF5 output files, and the log shows no NaN or divergence errors

### Requirement: Medium verification runs without divergence
The `RiftingCombinedYield` executable SHALL complete all time steps when run with the medium `.txt` without producing NaN values or solver divergence. This catches slow-onset instabilities (thermal runaway, melt feedback) that 10 steps may miss.

#### Scenario: Medium run completes successfully
- **WHEN** the user runs `./RiftingCombinedYield` with the medium `.txt` copied as `RiftingCombinedYield.txt`
- **THEN** the simulation completes 100 time steps, produces HDF5 output files, and the log shows no NaN or divergence errors

### Requirement: All verification configs use same executable
Both verify `.txt` files SHALL be used with the same `RiftingCombinedYield` executable (same `.c` file). No separate build targets are needed — the user copies the desired `.txt` as `RiftingCombinedYield.txt` in the run directory.

#### Scenario: Single build target for all configs
- **WHEN** the user wants to run verification
- **THEN** they copy the desired `.txt` (quick or medium) to `RiftingCombinedYield.txt` in the execution directory and run the same `RiftingCombinedYield` binary
