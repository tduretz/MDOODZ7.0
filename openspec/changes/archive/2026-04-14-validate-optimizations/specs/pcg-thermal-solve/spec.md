## MODIFIED Requirements

### Requirement: PCG iterative solver for the thermal system
The system SHALL provide a Preconditioned Conjugate Gradient (PCG) solver as an alternative to CHOLMOD for solving the SPD thermal equation. The solver SHALL operate on the CSR matrix arrays produced by the existing energy matrix assembly and SHALL use Jacobi (diagonal) preconditioning. When `export_pcg_residuals = 1`, the solver SHALL write per-iteration residual norms to a CSV file.

#### Scenario: PCG solves a standard thermal system
- **WHEN** `thermal_solver = 1` in the `.txt` file and `model.polar == 0`
- **THEN** the thermal equation SHALL be solved using PCG instead of CHOLMOD, and the solution SHALL be written to `mesh->T[]` identically to the direct solver path

#### Scenario: PCG converges within tolerance
- **WHEN** the PCG solver runs and the relative residual $\|r_k\| / \|b\|$ drops below `rel_tol_thermal`
- **THEN** the solver SHALL stop, log the iteration count and final residual via `LOG_TIME`, and return success

#### Scenario: PCG warm-starts from previous temperature
- **WHEN** `mesh->T[]` contains the solution from the previous timestep
- **THEN** the PCG solver SHALL use `mesh->T[]` as the initial guess $x_0$ rather than zero

#### Scenario: PCG exports per-iteration residuals when enabled
- **WHEN** `export_pcg_residuals = 1` in the `.txt` file and `thermal_solver = 1`
- **THEN** the solver SHALL write a CSV file `pcg_residuals_step<NNNNN>.csv` in the `writer_subfolder` directory with columns: `iteration,abs_residual,rel_residual` (one row per PCG iteration)

#### Scenario: Residual export disabled by default
- **WHEN** the `.txt` file does not contain `export_pcg_residuals`
- **THEN** the system SHALL default to `export_pcg_residuals = 0` and no CSV file SHALL be written

#### Scenario: Residual export has no effect on solver behaviour
- **WHEN** `export_pcg_residuals = 1`
- **THEN** the PCG solver SHALL produce the same solution and convergence behaviour as when `export_pcg_residuals = 0`
