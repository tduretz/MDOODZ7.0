## ADDED Requirements

### Requirement: PCG iterative solver for the thermal system
The system SHALL provide a Preconditioned Conjugate Gradient (PCG) solver as an alternative to CHOLMOD for solving the SPD thermal equation. The solver SHALL operate on the CSR matrix arrays produced by the existing energy matrix assembly and SHALL use Jacobi (diagonal) preconditioning.

#### Scenario: PCG solves a standard thermal system
- **WHEN** `thermal_solver = 1` in the `.txt` file and `model.polar == 0`
- **THEN** the thermal equation SHALL be solved using PCG instead of CHOLMOD, and the solution SHALL be written to `mesh->T[]` identically to the direct solver path

#### Scenario: PCG converges within tolerance
- **WHEN** the PCG solver runs and the relative residual $\|r_k\| / \|b\|$ drops below `rel_tol_thermal`
- **THEN** the solver SHALL stop, log the iteration count and final residual via `LOG_TIME`, and return success

#### Scenario: PCG warm-starts from previous temperature
- **WHEN** `mesh->T[]` contains the solution from the previous timestep
- **THEN** the PCG solver SHALL use `mesh->T[]` as the initial guess $x_0$ rather than zero

### Requirement: Solver selection via .txt parameter
The system SHALL read a `thermal_solver` integer parameter from the `.txt` input file to select between direct (CHOLMOD) and iterative (PCG) solvers.

#### Scenario: Default is direct solver
- **WHEN** the `.txt` file does not contain `thermal_solver`
- **THEN** the system SHALL default to `thermal_solver = 0` (CHOLMOD direct solver) and behaviour SHALL be identical to the current code

#### Scenario: User selects iterative solver
- **WHEN** the `.txt` file contains `thermal_solver = 1`
- **THEN** the system SHALL use the PCG solver for all thermal solves (except polar mode)

### Requirement: PCG convergence parameters
The system SHALL read `max_its_thermal` (int, default 1000) and `rel_tol_thermal` (double, default 1e-8) from the `.txt` input file to control PCG convergence.

#### Scenario: Custom iteration limit
- **WHEN** the `.txt` file contains `max_its_thermal = 500`
- **THEN** the PCG solver SHALL stop after at most 500 iterations

#### Scenario: Custom tolerance
- **WHEN** the `.txt` file contains `rel_tol_thermal = 1.0e-10`
- **THEN** the PCG solver SHALL converge to a relative residual below $10^{-10}$

### Requirement: Fallback to CHOLMOD on PCG failure
The system SHALL fall back to the CHOLMOD direct solver if the PCG solver fails to converge within `max_its_thermal` iterations.

#### Scenario: PCG exceeds iteration limit
- **WHEN** the PCG solver reaches `max_its_thermal` iterations without satisfying the tolerance
- **THEN** the system SHALL log a warning via `LOG_WARN` with the iteration count and residual, then solve the same system using the CHOLMOD direct path, and the simulation SHALL continue

#### Scenario: Fallback produces correct solution
- **WHEN** fallback to CHOLMOD is triggered
- **THEN** the solution written to `mesh->T[]` SHALL be the CHOLMOD direct solution (not the unconverged PCG iterate)

### Requirement: OpenMP parallelism in PCG
The SpMV kernel and vector operations (dot products, axpy, scaling) in the PCG solver SHALL use OpenMP parallel for directives to scale with the thread count set by `omp_nthread`.

#### Scenario: PCG scales with thread count
- **WHEN** the PCG solver runs with N OpenMP threads
- **THEN** the SpMV and vector operations SHALL execute in parallel across N threads

#### Scenario: Single-thread correctness
- **WHEN** `omp_nthread = 1`
- **THEN** the PCG solver SHALL produce the same result as with multiple threads (deterministic reduction)

### Requirement: Polar mode exclusion
When `model.polar == 1`, the system SHALL bypass PCG and use the CHOLMOD direct solver regardless of the `thermal_solver` parameter.

#### Scenario: Polar mode forces direct solver
- **WHEN** `thermal_solver = 1` and `model.polar == 1`
- **THEN** the system SHALL use CHOLMOD (the polar A·A^T formulation is not supported by PCG) and SHALL NOT log a warning (this is expected behaviour)

### Requirement: Backward compatibility
The existing CHOLMOD direct solver path SHALL remain unchanged. All existing `.txt` files SHALL produce identical results without modification.

#### Scenario: Existing simulations unaffected
- **WHEN** a `.txt` file does not contain `thermal_solver`
- **THEN** the simulation SHALL behave identically to the code before this change (same solver, same results, same performance within measurement noise)

#### Scenario: Direct solver path unchanged
- **WHEN** `thermal_solver = 0` (explicitly or by default)
- **THEN** the `FactorEnergyCHOLMOD` and `SolveEnergyCHOLMOD` functions SHALL be called with identical arguments and behaviour as before this change
