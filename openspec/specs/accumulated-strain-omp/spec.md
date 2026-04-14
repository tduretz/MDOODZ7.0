## ADDED Requirements

### Requirement: AccumulatedStrainII particle loop SHALL be OpenMP-parallelized
The `AccumulatedStrainII` function SHALL execute its particle interpolation loop (grid-to-particle strain increment accumulation) using `#pragma omp parallel for` with correct `private`, `shared`, and `firstprivate` clauses.

#### Scenario: Parallel execution produces identical results to serial
- **WHEN** `AccumulatedStrainII` runs with `OMP_NUM_THREADS > 1`
- **THEN** the resulting `particles->strain*[k]` values SHALL be bit-identical to those produced with `OMP_NUM_THREADS = 1`

#### Scenario: All local variables are thread-private
- **WHEN** the particle loop executes in parallel
- **THEN** all per-iteration variables (`dE_tot`, `dE_el`, `dE_pl`, `dE_pl_vol`, `dE_pwl`, `dE_exp`, `dE_lin`, `dE_gbs`, `sumW`, `dst`, `dxm`, `dzm`, `iSW`, `iSE`, `iNW`, `iNE`, `i_part`, `j_part`, `k`, `l`) SHALL be declared `private` to avoid data races

### Requirement: No process-terminating calls inside parallel regions
The `AccumulatedStrainII` function SHALL NOT call `exit()` or `abort()` inside an OpenMP parallel region.

#### Scenario: Negative strain increment detected during parallel execution
- **WHEN** a negative `strain_inc[c1]` is detected inside the parallel grid loop
- **THEN** the function SHALL log the error using `LOG_ERR` and terminate execution outside the parallel region

#### Scenario: Thread-safe error reporting
- **WHEN** multiple threads detect errors simultaneously
- **THEN** error logging SHALL be thread-safe (no interleaved or corrupted output)
