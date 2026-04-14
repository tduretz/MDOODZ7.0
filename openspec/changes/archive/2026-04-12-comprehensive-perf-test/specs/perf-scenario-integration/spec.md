## ADDED Requirements

### Requirement: Benchmark script accepts --scenario flag
`misc/benchmark.sh` SHALL accept a `--scenario <name>` flag that specifies which MDOODZ scenario to build and run. The default SHALL remain `RiftingBasic` for backward compatibility.

#### Scenario: Custom scenario builds and runs
- **WHEN** the user runs `./misc/benchmark.sh --scenario RiftingComprehensive`
- **THEN** the script builds with `-DSET=RiftingComprehensive` and runs the `RiftingComprehensive` executable

#### Scenario: Default scenario unchanged
- **WHEN** the user runs `./misc/benchmark.sh` without `--scenario`
- **THEN** the script builds and runs `RiftingBasic` as before

### Requirement: Resolution-based benchmark sweep
`misc/benchmark.sh` SHALL accept a `--resolutions` flag listing resolution tier names (e.g. `"lowres default medres highres"`). Each name maps to `<Scenario>_<res>.txt` (with `default` mapping to `<Scenario>.txt`). When `--resolutions` is specified, the script SHALL iterate over resolution files instead of patching grid sizes from `--grids`.

#### Scenario: Resolution sweep runs all tiers
- **WHEN** the user runs `./misc/benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default medres highres" --threads "1 4 8"`
- **THEN** the script runs 4 resolutions × 3 thread counts = 12 benchmark runs, using each resolution's .txt file directly

#### Scenario: Grid size comes from .txt file
- **WHEN** `--resolutions` is used
- **THEN** the Nx and Nz values are read from each resolution's `.txt` file, not from `--grids`

### Requirement: Validation mode
`misc/benchmark.sh` SHALL accept a `--validate` flag. When set, before the full benchmark sweep, the script SHALL run each resolution for 3 timesteps with 1 thread, proceeding from lowest to highest resolution. If any resolution crashes (non-zero exit), the script SHALL stop immediately and report the failure.

By default, validation SHALL only run on `lowres` and `default` resolutions. Higher resolutions (`medres`, `highres`) are too slow for quick validation and SHALL be skipped unless explicitly included via `--validate-resolutions "lowres default medres highres"`.

#### Scenario: Validation catches a crash
- **WHEN** `--validate` is set and the medres variant crashes at step 2
- **THEN** the script reports "VALIDATION FAILED: medres crashed at step 2" and exits without running the full benchmark

#### Scenario: Validation passes then benchmark runs
- **WHEN** `--validate` is set and all resolutions pass 3 timesteps
- **THEN** the full benchmark sweep proceeds normally

### Requirement: Expanded perf.csv columns
The `perf.csv` output SHALL include per-subsystem timing columns in addition to existing ones. The new columns SHALL be: `time_ma`, `thermal_s`, `advection_s`, `free_surface_s`, `reseeding_s`, `melting_s`, `anisotropy_s`, `gse_s`, `output_s`.

The full header SHALL be: `step,wall_s,time_ma,rheology_s,assembly_s,solve_s,thermal_s,advection_s,free_surface_s,reseeding_s,melting_s,anisotropy_s,gse_s,output_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s`

Subsystems that are disabled (e.g. `thermal=0`) SHALL write `0.0000` for their timing column.

#### Scenario: All timing columns populated
- **WHEN** `RiftingComprehensive` runs with all subsystems enabled
- **THEN** `perf.csv` contains non-zero values for `thermal_s`, `advection_s`, `melting_s`, and `anisotropy_s` columns

#### Scenario: Disabled subsystems write zero
- **WHEN** `RiftingBasic` runs (no melting, no anisotropy)
- **THEN** `perf.csv` contains `0.0000` for `melting_s`, `anisotropy_s`, and `gse_s` columns

#### Scenario: Backward compatibility with report scripts
- **WHEN** `benchmark-report.sh` processes a perf.csv with the expanded columns
- **THEN** the report generates correctly, displaying the new subsystem timing columns

### Requirement: Model time in perf.csv
The `perf.csv` SHALL include a `time_ma` column containing the current model time in Ma (megayears). The value SHALL be computed as `model.time * scaling.t / (1e6 * 365.25 * 24 * 3600)`.

#### Scenario: Model time increases across steps
- **WHEN** a simulation runs for 10 timesteps
- **THEN** the `time_ma` column shows monotonically increasing values reflecting the accumulated model time

### Requirement: Model time in log metadata
When `log_metadata = 1`, the log prefix SHALL include model time alongside step and iteration: `[timestamp|S0001|T0.42Ma|N00]`. The time unit SHALL be configurable via a `time_unit` parameter in the `.txt` file: `0` = Ma (default), `1` = Ka, `2` = yr.

#### Scenario: Model time appears in log metadata
- **WHEN** `log_metadata = 1` and model time is 0.42 Ma
- **THEN** log lines include `T0.42Ma` in the metadata prefix

#### Scenario: Time unit configurable
- **WHEN** `time_unit = 1` is set in the `.txt` file
- **THEN** the metadata shows `T420.00Ka` instead of `T0.42Ma`

### Requirement: Time-based loop termination in Main_DOODZ.c
The main simulation loop SHALL support a `t_end` parameter (in seconds). When `t_end > 0`, the loop SHALL terminate when `model.time >= t_end / scaling.t`, taking priority over `Nt`. When `t_end = 0` (default), the existing `Nt`-based termination SHALL be preserved.

#### Scenario: t_end terminates loop
- **WHEN** `t_end = 4.73e14` and `Nt = 99999`
- **THEN** the simulation stops when model time reaches ~15 Ma, regardless of step count

#### Scenario: Default behavior preserved
- **WHEN** `t_end` is not set (defaults to 0)
- **THEN** the simulation runs for exactly `Nt` timesteps as before

### Requirement: summary.csv compatibility
The `summary.csv` generated by `benchmark.sh` SHALL work correctly with the expanded `perf.csv` format. The summary SHALL include average timings for all new subsystem columns.

#### Scenario: Summary includes new columns
- **WHEN** `benchmark.sh` completes a run with the expanded perf.csv
- **THEN** `summary.csv` includes average values for `thermal_s`, `advection_s`, and other new timing columns
