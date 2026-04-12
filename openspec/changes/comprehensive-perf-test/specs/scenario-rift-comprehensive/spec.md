## ADDED Requirements

### Requirement: Scenario files exist and build
The scenario SHALL consist of `SETS/RiftingComprehensive.c` and multiple `.txt` parameter files. The `.c` file SHALL be registered in `SETS/CMakeLists.txt` via `add_set(RiftingComprehensive)`. The executable SHALL compile without errors and link against the `mdoodz` library.

#### Scenario: CMake target builds successfully
- **WHEN** the user runs `cmake --build` with `-DSET=RiftingComprehensive`
- **THEN** the `RiftingComprehensive` executable is produced in `cmake-exec/RiftingComprehensive/`

### Requirement: Based on RiftingCombinedYield with all subsystems enabled
The scenario SHALL be derived from `RiftingCombinedYield` (origin/add-meltin-vp-model branch) and SHALL enable all of: elasticity (`elastic=1`), thermal solving (`thermal=1`), shear heating (`shear_heating=1`), melting (`melting=1`), combined-yield plasticity (`plast=2`), visco-plasticity (`eta_vp`), free surface (`free_surface=1`), surface processes, anisotropy, and grain-size evolution. The `.c` callback file SHALL reuse the same `SetPhase`, `SetTemperature`, `SetDensity`, `SetBCT`, and pure-shear BC functions as RiftingCombinedYield.

#### Scenario: All subsystem switches active in default .txt
- **WHEN** the default `.txt` file is parsed
- **THEN** `elastic=1`, `thermal=1`, `shear_heating=1`, `melting=1`, `free_surface=1`, `surface_processes>=1`, `aniso=1`, `aniso_fstrain=1` are all set, and at least one phase has `gsel=1`

### Requirement: Anisotropy enabled on mantle phases only
The scenario SHALL enable anisotropy (`aniso=1`, `aniso_fstrain=1`) globally. Anisotropy-specific material parameters SHALL be set on mantle phases (IDs 1, 2, 3, 5) only. Crustal phases (IDs 0, 4, 6, 7, 8) SHALL remain isotropic (no anisotropy parameters or `ani_fac_max = 1.0`).

#### Scenario: Mantle phases have anisotropy parameters
- **WHEN** the `.txt` file is parsed
- **THEN** phases 1, 2, 3, 5 have anisotropy factors (`aniso_factor > 1.0`) and phases 0, 4, 6, 7, 8 have no anisotropy activation or `aniso_factor = 1.0`

### Requirement: Grain-size evolution on dry olivine mantle
At least the dry olivine mantle phase (ID 1) SHALL have grain-size evolution enabled (`gsel=1`) with appropriate grain-size-dependent flow law parameters.

#### Scenario: GSE parameters in phase 1
- **WHEN** the `.txt` file is parsed
- **THEN** phase 1 (lithospheric mantle) has `gsel = 1` with valid grain-size parameters

### Requirement: Time-based simulation termination
All `.txt` parameter files SHALL use `t_end` (in seconds) to set the simulation end time rather than relying on `Nt` alone. `Nt` SHALL be set to a large safety cap (e.g. 99999). When `t_end > 0`, the simulation loop SHALL terminate when `model.time >= t_end / scaling.t`, regardless of step count.

The default end time SHALL be 15 Ma (`t_end = 4.73e14`).

#### Scenario: Simulation terminates at model time
- **WHEN** the simulation runs with `t_end = 4.73e14`
- **THEN** the simulation terminates when model time reaches approximately 15 Ma, regardless of how many timesteps were taken

#### Scenario: Nt safety cap prevents infinite loop
- **WHEN** `t_end` is set and `Nt = 99999`
- **THEN** if the simulation has not reached `t_end` by step 99999, it terminates anyway

### Requirement: Four resolution variants
The scenario SHALL provide four `.txt` parameter files with different grid resolutions:

| File | Nx × Nz | Purpose |
|------|---------|---------|
| `RiftingComprehensive_lowres.txt` | 100 × 80 | Quick sanity check |
| `RiftingComprehensive.txt` | 200 × 160 | Default benchmark baseline |
| `RiftingComprehensive_medres.txt` | 501 × 401 | Scaling study |
| `RiftingComprehensive_highres.txt` | 1000 × 800 | Stress test |

All four SHALL share the same `t_end`, phase definitions, and physics switches. Each MAY have individually tuned `dt` and `Courant` values for numerical stability at its resolution.

#### Scenario: All four .txt files exist and are parseable
- **WHEN** the scenario is built
- **THEN** all four `.txt` files exist in `SETS/` and the executable runs successfully with each one (at least 3 timesteps without crash)

#### Scenario: Resolution-specific dt tuning
- **WHEN** the highres variant is run
- **THEN** the initial `dt` and `Courant` number are appropriate for the 1000×800 grid (may differ from lowres values)

### Requirement: Combined-yield and melting from RiftingCombinedYield
The scenario SHALL preserve all material properties, melting parameters, and combined-yield settings from `RiftingCombinedYield`: `plast=2` on all phases, `T_st=-10e6`, `eta_vp=2.5e20`, `melt=11` for crustal phases, `melt=41` for mantle phases, `melt_weak=30.0`, `force_melt_weak=1`.

#### Scenario: Material properties preserved from base model
- **WHEN** the `.txt` file is compared to `RiftingCombinedYield.txt`
- **THEN** all non-anisotropy, non-GSE material properties match the base model

### Requirement: Model stability across all resolutions
Before the scenario is considered complete, it SHALL have been validated by running 3 timesteps at each resolution (lowres → default → medres → highres) without crash or divergence.

#### Scenario: Validation pass succeeds
- **WHEN** `benchmark.sh --validate --scenario RiftingComprehensive` is run
- **THEN** all four resolutions complete 3 timesteps with exit code 0 and no NaN values in the output
