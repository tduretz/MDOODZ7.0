## ADDED Requirements

### Requirement: Benchmark configuration with HDF5 output enabled
The benchmark infrastructure SHALL support running the highres (1000×800) benchmark configuration with HDF5 output enabled (`writer = 1`, `writer_step = 1`) to populate the `output_s` column in perf.csv with actual I/O timing data.

#### Scenario: Highres benchmark with output enabled
- **WHEN** the benchmark suite runs the highres configuration on EC2
- **THEN** the `.txt` parameter file SHALL have `writer = 1` and `writer_step = 1` set, and the resulting perf.csv SHALL contain non-zero values in the `output_s` column for each timestep
