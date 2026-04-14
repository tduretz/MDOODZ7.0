## ADDED Requirements

### Requirement: BlankenBench .txt available in SETS/
The benchmark setup SHALL copy `TESTS/BlankenBench/BlankenBench.txt` to `SETS/BlankenBench.txt` before running `benchmark.sh --scenario BlankenBench`, since `benchmark.sh` reads `.txt` files from `SETS/`.

#### Scenario: .txt copied before benchmark
- **WHEN** preparing to run BlankenBench via `benchmark.sh`
- **THEN** `SETS/BlankenBench.txt` SHALL exist as a copy of `TESTS/BlankenBench/BlankenBench.txt`

#### Scenario: .txt cleaned up after benchmark
- **WHEN** both benchmark sweeps (baseline and optimized) have completed
- **THEN** `SETS/BlankenBench.txt` SHALL be removed to avoid polluting the SETS directory

### Requirement: Baseline benchmark sweep
The system SHALL run a baseline BlankenBench benchmark sweep with default settings (thermal_solver=0, interp_mode=0) on MacBook M1 with writer enabled and writer_step=1.

#### Scenario: Baseline sweep parameters
- **WHEN** the baseline sweep is launched
- **THEN** `benchmark.sh` SHALL be invoked with `--scenario BlankenBench --steps 1000 --writer --threads "1 2 4 6 8"` and no `--sed` overrides for thermal_solver or interp_mode

#### Scenario: Baseline uses CHOLMOD thermal solver
- **WHEN** the baseline sweep runs
- **THEN** `thermal_solver` SHALL be 0 (CHOLMOD direct solve) and `interp_mode` SHALL be 0 (legacy per-call allocation)

### Requirement: Optimized benchmark sweep
The system SHALL run an optimized BlankenBench benchmark sweep with PCG thermal solver and fused P2Mastah on MacBook M1 with writer enabled and writer_step=1.

#### Scenario: Optimized sweep parameters
- **WHEN** the optimized sweep is launched
- **THEN** `benchmark.sh` SHALL be invoked with `--scenario BlankenBench --steps 1000 --writer --threads "1 2 4 6 8" --sed 's/^thermal_solver.*/thermal_solver = 1/' --sed 's/^interp_mode.*/interp_mode = 3/'`

### Requirement: Thread counts match M1 core topology
The benchmark SHALL test thread counts 1, 2, 4, 6, 8 to capture the M1's performance-core to efficiency-core transition.

#### Scenario: Thread scaling covers P-core and E-core
- **WHEN** the thread list is configured
- **THEN** it SHALL include 1, 2, 4 (within P-cores), 6 (P+E mix), and 8 (all cores)

### Requirement: HDF5 output enabled every step
Both sweeps SHALL run with `writer = 1` and `writer_step = 1` (HDF5 output every timestep) so that `output_s` is visible in the phase breakdown.

#### Scenario: Writer active every step
- **WHEN** benchmark.sh runs with `--writer`
- **THEN** the patched .txt SHALL have `writer = 1` and `writer_step = 1`

### Requirement: Results documented as Experiment 8
The benchmark results SHALL be analyzed and documented as Experiment 8 in `skill-benchmarking/SKILL.md`, following the established format (Hypothesis, Implementation, Results table, Conclusion, Data, Code).

#### Scenario: Experiment 8 entry created
- **WHEN** both sweeps have completed and results are analyzed
- **THEN** an Experiment 8 section SHALL be appended to `skill-benchmarking/SKILL.md` with baseline vs optimized comparison tables, phase breakdowns, and M1-specific observations

#### Scenario: Cross-platform comparison noted
- **WHEN** documenting Experiment 8
- **THEN** the conclusion SHALL note key differences between M1 (small grid, many steps) and EC2 (large grid, few steps) workload profiles, without making direct timing comparisons
