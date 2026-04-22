# cholmod-factor-recycling Specification

## Purpose
TBD - created by archiving change add-cholmod-factor-recycling. Update Purpose after archive.
## Requirements
### Requirement: cholmod_recycle configuration parameter

The system SHALL support an integer parameter `cholmod_recycle` in the `.txt` configuration file, read via `ReadInt2(fin, "cholmod_recycle", 0)`. The value `0` (default) SHALL disable factor recycling and produce bit-identical behaviour to the pre-change codebase. The value `1` SHALL enable recycling of the CHOLMOD numeric factor across Picard iterations within a time step, subject to the staleness budget in `cholmod_recycle_max_reuse`.

#### Scenario: Parameter absent from .txt file

- **WHEN** `cholmod_recycle` is not specified in the `.txt` configuration file
- **THEN** the system SHALL default to `0` and KillerSolver SHALL call `cholmod_factorize` on every invocation, matching the pre-change behaviour exactly

#### Scenario: Parameter set to 1 with default reuse budget

- **WHEN** `cholmod_recycle = 1` is specified and `cholmod_recycle_max_reuse` is absent
- **THEN** the system SHALL enable recycling with the default reuse budget of `3`

#### Scenario: Invalid value

- **WHEN** `cholmod_recycle` is set to any value other than `0` or `1`
- **THEN** the system SHALL emit a `LOG_WARN` explaining the accepted values and SHALL treat the value as `0`

### Requirement: cholmod_recycle_max_reuse staleness budget

The system SHALL support a non-negative integer parameter `cholmod_recycle_max_reuse` in the `.txt` configuration file, read via `ReadInt2(fin, "cholmod_recycle_max_reuse", 3)`. This parameter SHALL bound the maximum number of consecutive Picard iterations that reuse the same numeric factor before a fresh `cholmod_factorize` is forced. The value `0` SHALL be equivalent to `cholmod_recycle = 0`.

#### Scenario: Default budget of 3

- **WHEN** `cholmod_recycle = 1` and `cholmod_recycle_max_reuse` is absent
- **THEN** the system SHALL reuse the factor for at most 3 consecutive Picard iterations and THEN force a refactor on the 4th iteration

#### Scenario: Budget exhausted mid-step

- **WHEN** `cholmod_recycle_max_reuse = N` and the current Picard iteration index relative to the last factorization is equal to `N`
- **THEN** KillerSolver SHALL call `cholmod_factorize` on the next invocation before entering the Powell-Hestenes loop, resetting the reuse counter to 0

#### Scenario: Budget set to zero disables recycling

- **WHEN** `cholmod_recycle_max_reuse = 0` (regardless of `cholmod_recycle`)
- **THEN** the system SHALL behave identically to `cholmod_recycle = 0`

### Requirement: Recycling triggered by Picard iteration index

The system SHALL pass an explicit Picard iteration index (or equivalent boolean `is_first_picard_of_step`) from `SolveStokesDecoupled` / `SolveStokesDefectDecoupled` into `KillerSolver`. KillerSolver SHALL call `cholmod_factorize` unconditionally when the index is `0`. For subsequent indices, factorization SHALL be skipped if and only if recycling is enabled, the reuse counter has not exceeded `cholmod_recycle_max_reuse`, and no fallback has been triggered.

#### Scenario: First Picard iteration of a time step

- **WHEN** KillerSolver is called with Picard iteration index `0`
- **THEN** `cholmod_factorize(Lcml, Lfact, &c)` SHALL be called and the reuse counter SHALL be reset to `0`

#### Scenario: Subsequent Picard iteration with recycling enabled

- **WHEN** KillerSolver is called with Picard iteration index `k > 0`, `cholmod_recycle = 1`, reuse counter `r < cholmod_recycle_max_reuse`, and the previous GCR solve converged within budget
- **THEN** KillerSolver SHALL skip `cholmod_factorize` and proceed directly to the Powell-Hestenes loop, incrementing the reuse counter to `r + 1`

#### Scenario: Operator switch (Picard ↔ Newton) within a time step

- **WHEN** the nonlinear solver switches between the Picard operator and the Newton Jacobian operator mid-step
- **THEN** the caller SHALL signal operator change by passing Picard iteration index `0` on the next KillerSolver call, forcing a fresh factorization

### Requirement: GCR-convergence-based fallback to fresh factor

The system SHALL monitor the iteration count `its_KSP` returned from the `kspgcr` call inside KillerSolver's Powell-Hestenes loop when recycling is active. If `its_KSP` equals `model.max_its_KSP` (the GCR budget) AND the residual norm did not satisfy `model.rel_tol_KSP`, the recycled factor SHALL be treated as too stale. The system SHALL then call `cholmod_factorize(Lcml, Lfact, &c)` once, reset the reuse counter to `0`, and retry the Powell-Hestenes iteration that failed.

#### Scenario: Stale factor preconditions well

- **WHEN** recycling is active and `kspgcr` converges with `its_KSP < model.max_its_KSP` and residual ≤ `model.rel_tol_KSP`
- **THEN** no fallback SHALL be triggered and the Powell-Hestenes iteration SHALL proceed normally

#### Scenario: Stale factor exhausts GCR budget

- **WHEN** recycling is active and `kspgcr` returns `its_KSP == model.max_its_KSP` with residual > `model.rel_tol_KSP`
- **THEN** the system SHALL log a `LOG_WARN` message indicating the fallback, call `cholmod_factorize` to refresh the factor, reset the reuse counter to `0`, and retry the current Powell-Hestenes iteration using the fresh factor

#### Scenario: Fallback succeeds

- **WHEN** a fallback refactorization is performed and the retried `kspgcr` call converges within budget
- **THEN** execution SHALL continue normally for the remainder of the time step, with the fresh factor available for subsequent Picard iterations up to `cholmod_recycle_max_reuse` reuses

#### Scenario: Fallback fails

- **WHEN** a fallback refactorization is performed and the retried `kspgcr` call still fails to converge within budget
- **THEN** the system SHALL emit a `LOG_ERR` and exit with a non-zero status, matching existing behaviour for unrecoverable solver failures

### Requirement: Factor lifecycle preservation

The CHOLMOD numeric factor `DirectSolver.Lfact` SHALL remain allocated across Powell-Hestenes iterations, Picard iterations, time steps, and restart operations. The factor SHALL be freed only via the existing solver-teardown path at program shutdown. No additional `cholmod_allocate_factor` or `cholmod_free_factor` calls SHALL be introduced in the recycle code path.

#### Scenario: Factor survives across Picard iterations

- **WHEN** KillerSolver completes a Powell-Hestenes iteration with recycling enabled
- **THEN** `DirectSolver.Lfact` SHALL remain allocated and non-NULL, ready for the next Picard iteration

#### Scenario: In-place refactor preserves allocation

- **WHEN** `cholmod_factorize(Lcml, Lfact, &c)` is called to refresh the factor (first Picard iteration, budget exhaustion, or fallback)
- **THEN** CHOLMOD SHALL update the factor in place without allocating a second `cholmod_factor`, and peak RSS SHALL NOT grow beyond one factor's size per time step

### Requirement: Startup logging of recycling configuration

The system SHALL log the resolved recycling configuration at startup via `LOG_INFO` to aid benchmark analysis and user diagnosis. The message SHALL report both `cholmod_recycle` and the effective `cholmod_recycle_max_reuse`.

#### Scenario: Recycling disabled

- **WHEN** `cholmod_recycle = 0` or `cholmod_recycle_max_reuse = 0` after resolution
- **THEN** the system SHALL emit `LOG_INFO("CHOLMOD factor recycling: disabled")`

#### Scenario: Recycling enabled

- **WHEN** `cholmod_recycle = 1` and `cholmod_recycle_max_reuse = N` (N ≥ 1) after resolution
- **THEN** the system SHALL emit `LOG_INFO("CHOLMOD factor recycling: enabled, max_reuse = N")`

#### Scenario: Fallback triggered during a run

- **WHEN** a fallback refactorization is performed mid-step due to GCR budget exhaustion
- **THEN** the system SHALL emit `LOG_WARN("CHOLMOD recycle: stale factor — refreshing (reuse_count = R)")` identifying how many reuses had accumulated

### Requirement: Non-interference with cholmod_threads

The `cholmod_recycle` feature SHALL be orthogonal to the existing `cholmod_threads` parameter. Both parameters SHALL remain independently configurable, and the resolved recycling behaviour SHALL NOT depend on the `cholmod_threads` value.

#### Scenario: Both parameters set

- **WHEN** `cholmod_recycle = 1` and `cholmod_threads = 4` are set in the same configuration
- **THEN** the system SHALL apply both settings: CHOLMOD SHALL use 4 internal threads AND the numeric factor SHALL be recycled across Picard iterations subject to the reuse budget

#### Scenario: Recycling without threading

- **WHEN** `cholmod_recycle = 1` and `cholmod_threads = 1` (or absent)
- **THEN** the system SHALL enable recycling and CHOLMOD SHALL operate single-threaded internally, with no interaction between the two features

### Requirement: Correctness parity with baseline

When `cholmod_recycle = 0`, the system SHALL produce bit-equivalent factorization counts, residual sequences, and final velocity/pressure fields compared to the pre-change codebase at the same git revision before the recycle change. When `cholmod_recycle = 1`, the system SHALL produce velocity and pressure fields that differ from the baseline by at most the Stokes nonlinear solver tolerance (`nonlin_abs_mom`, `nonlin_abs_div`) over any 5-step run.

#### Scenario: Default flag produces baseline behaviour

- **WHEN** `cholmod_recycle = 0` (default) is set and any existing benchmark scenario is run
- **THEN** the per-step `solve_s`, `nit`, and final state in `perf.csv` and HDF5 output SHALL match the pre-change baseline to floating-point round-off

#### Scenario: Recycling does not drift the solution

- **WHEN** a 5-step run is performed with `cholmod_recycle = 1` from `Breakpoint00020.dat` on RiftingChenin 1000×800
- **THEN** the velocity and pressure fields at the final step SHALL differ from the `cholmod_recycle = 0` run by at most `max(nonlin_abs_mom, nonlin_abs_div)` in the max norm

### Requirement: Applicability is scoped to KillerSolver

The recycling behaviour SHALL apply only to the `lin_solver = 2` / KillerSolver code path. The `DirectStokesDecoupled` (`lin_solver = 0`, silently redirected to KillerSolver), `DirectStokesDecoupledComp` (`lin_solver = -1`), and `KSPStokesDecoupled` (`lin_solver = 1`) paths SHALL NOT be modified by this change.

#### Scenario: Compressible path unaffected

- **WHEN** `cholmod_recycle = 1` and `lin_solver = -1` (compressible direct solver) are set together
- **THEN** the compressible path SHALL ignore `cholmod_recycle` and run its existing factorization-per-call behaviour, with a `LOG_INFO` noting that recycling is not supported for this solver

#### Scenario: KSP path unaffected

- **WHEN** `cholmod_recycle = 1` and `lin_solver = 1` (Krylov) are set together
- **THEN** the KSP path SHALL ignore `cholmod_recycle` and run its existing behaviour, with a `LOG_INFO` noting that recycling is not supported for this solver

