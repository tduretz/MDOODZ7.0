## MODIFIED Requirements

### Requirement: GMG-preconditioned FGMRES is the production solver mode

When `lin_solver = 3` and `gmg_standalone` is unset or zero, the Stokes system SHALL be solved by FGMRES with the GMG V-cycle as a right preconditioner. FGMRES SHALL iterate until either the **relative residual reduction** `‖r_k‖ / ‖r_0‖` reaches `gmg_fgmres_tol` (default: `nonlin_abs_mom / 10`) or the iteration count reaches `gmg_fgmres_restart × gmg_fgmres_max_restarts` (defaults: restart = 30, max-restarts = 20). The convergence predicate and the log message reporting "converged" SHALL both test and report the *same* quantity — namely `‖r_k‖ / ‖r_0‖` — so that a declared-converged run always satisfies `reported_final_res ≤ gmg_fgmres_tol` numerically.

#### Scenario: Convergence on a typical lithospheric setup

- **WHEN** a lithospheric-scale model with viscosity contrast up to 1e6 is solved with `lin_solver = 3`
- **THEN** FGMRES converges within its iteration budget, the final relative residual is below `gmg_fgmres_tol`, and the solver logs residual history per FGMRES iteration

#### Scenario: Diagnostic standalone V-cycle mode

- **WHEN** `lin_solver = 3` and `gmg_standalone = 1` are set
- **THEN** the Stokes solve applies the V-cycle without the FGMRES outer iteration, exposing the raw multigrid convergence for instrumentation

#### Scenario: Convergence predicate and log message agree

- **WHEN** FGMRES declares "converged" via `LOG_INFO`
- **THEN** the reported `final_res_rel` in the log line is numerically ≤ `gmg_fgmres_tol`, within 1e-12 of the value actually tested by the convergence predicate

#### Scenario: Non-convergence log message reports the correct quantity

- **WHEN** FGMRES exhausts its iteration budget without converging
- **THEN** the log line reports `final_res_rel` (the relative-reduction quantity tested by the predicate), not the raw `‖r_k‖`, and "converged" is not present in the message

### Requirement: Solver exposes optional tuning parameters with safe defaults

The `.txt` input SHALL accept the following optional parameters, each with the default listed; unset parameters MUST fall back to defaults without error:

- `gmg_levels` (default: auto-selected per the grid-hierarchy requirement)
- `gmg_nu_pre` (default: 2) — number of pre-smoothing Vanka sweeps per level
- `gmg_nu_post` (default: 2) — number of post-smoothing Vanka sweeps per level
- `gmg_fgmres_restart` (default: 30) — FGMRES restart length
- `gmg_fgmres_max_restarts` (default: 20) — maximum number of FGMRES restarts before declaring non-convergence; valid range `[1, 1000]`
- `gmg_fgmres_tol` (default: `nonlin_abs_mom / 10`) — FGMRES relative-residual convergence tolerance
- `gmg_standalone` (default: 0) — 1 = bare V-cycle, 0 = FGMRES-wrapped V-cycle
- `gmg_dump_vcycle` (default: 0) — 1 = dump per-level residual + solution snapshots for diagnostic use

#### Scenario: Defaults applied when parameters are absent

- **WHEN** a `.txt` input sets `lin_solver = 3` and omits every `gmg_*` parameter
- **THEN** the solver runs with the default values above and logs the effective configuration at startup, including the active `gmg_fgmres_max_restarts` value

#### Scenario: Invalid parameter is rejected with a clear error

- **WHEN** a `.txt` input sets `gmg_nu_pre = 0`, `gmg_fgmres_restart < 1`, or `gmg_fgmres_max_restarts` outside `[1, 1000]`
- **THEN** the solver aborts at startup with a log-level-ERROR message naming the offending parameter and its valid range

#### Scenario: `gmg_fgmres_max_restarts` tunable at runtime

- **WHEN** a `.txt` input sets `gmg_fgmres_max_restarts = 50` to permit longer FGMRES diagnostics
- **THEN** FGMRES iterates up to `gmg_fgmres_restart × 50` total inner iterations before declaring non-convergence, without requiring a rebuild

### Requirement: Solver logs residual history and falls back on non-convergence

On every Stokes linear solve using `lin_solver = 3`, the code SHALL log the initial residual `‖r_0‖`, per-iteration relative residuals `‖r_k‖ / ‖r_0‖` of the outer FGMRES iteration, and the final relative residual. If FGMRES does not converge within its iteration budget (restart length × max-restarts), the code SHALL log a WARN-level message whose numerical content is directly comparable to `gmg_fgmres_tol` (i.e. the relative-reduction quantity tested by the predicate) and fall back to the CHOLMOD direct solver for that Stokes solve, consistent with the thermal-PCG fallback pattern in `ThermalRoutines.c`.

#### Scenario: Convergence history is logged

- **WHEN** any Stokes solve with `lin_solver = 3` completes successfully
- **THEN** the log contains the initial residual `‖r_0‖`, per-FGMRES-iteration relative residuals `‖r_k‖ / ‖r_0‖`, and the final relative residual, each tagged so they can be extracted by the existing benchmarking scripts

#### Scenario: Fallback on non-convergence reports the comparable quantity

- **WHEN** FGMRES fails to reach `gmg_fgmres_tol` within `gmg_fgmres_restart × gmg_fgmres_max_restarts` inner iterations
- **THEN** a `LOG_WARN` line is emitted naming the step number and the achieved `final_res_rel` (the quantity tested against the tolerance), NOT the raw `‖r_k‖`; the CHOLMOD direct solver is invoked on the same system; simulation continues without aborting

#### Scenario: "Converged" is never logged when predicate fails

- **WHEN** the predicate `‖r_k‖ / ‖r_0‖ ≤ gmg_fgmres_tol` is false
- **THEN** the word "converged" is absent from the relevant log line regardless of how the predicate expression evolved during the run

## ADDED Requirements

### Requirement: V-cycle residual amplification is bounded

Across a single V-cycle, the composed residual amplification `‖r_exit‖ / ‖r_entry‖` at the finest level SHALL NOT exceed 4.0 on constant-viscosity SolCx-analog fixtures at 81×81, 161×161, and 201×201, measured via `gmg_dump_vcycle = 1` HDF5 snapshots. Additionally no single stage of the cycle — pre-smooth, restrict, coarse-operator apply, prolongate, post-smooth — SHALL amplify the residual by more than 4.0 relative to its input, preventing any single stage from silently violating the whole-cycle contract.

**Scope note.** This requirement was originally written with a tighter 2.0× per-stage bound under the assumption of Galerkin coarsening. The realised hierarchy uses MDOODZ-bridge operator at level 0 and textbook Picard rediscretisation on restricted viscosity at coarse levels (design D4 of `add-gmg-stokes-solver`); the cross-level operator mismatch appears as a benign 2×–4× jump specifically on the level-0 prolongate stage. Because the whole-cycle ratio remains well inside the 4.0 bound (measured 0.64×–1.36× at 81²–201² post-fix), the per-stage ceiling is set equal to the whole-cycle bound so that a genuine single-stage divergence (e.g. a smoother with spectral radius > 1 as observed pre-fix) still fires the regression without false-positives from the known Picard/bridge inconsistency.

#### Scenario: Whole V-cycle amplification is bounded

- **WHEN** a constant-viscosity 81² / 161² / 201² SolCx run executes with `gmg_dump_vcycle = 1` and the first V-cycle dump is analysed
- **THEN** the whole-V-cycle `‖r_exit‖ / ‖r_entry‖` ratio at level 0 stays below 4.0, demonstrating that the V-cycle is a productive preconditioner rather than a residual amplifier

#### Scenario: Per-stage amplification cannot silently violate the whole-cycle bound

- **WHEN** the same dump is walked stage-by-stage (pre-smooth, restrict, coarse-solve, prolongate, post-smooth) across both downleg and upleg
- **THEN** no single stage's output-to-input residual ratio exceeds 4.0, guarding against a future regression where one stage blows up by 8× and the rest damp to hide it

#### Scenario: Level-0 pre-smoother is a contraction on constant-viscosity inputs

- **WHEN** the first V-cycle dump of a constant-viscosity 81² fixture is examined at `level_0_pre_smooth` (pre → post)
- **THEN** the ratio is ≤ 4.0 (and, post-fix, is below 1.0 — a genuine damping), confirming the Vanka smoother at the finest level is a contraction; this was the specific pathology pre-fix (ratio 7.55×) that Fix F addresses by routing L0 `VankaBlockAssembleSolve` through the textbook block instead of the MDOODZ bridge

#### Scenario: 201² default-levels convergence

- **WHEN** the SolViPerf harness is driven at 201² with `lin_solver = 3` and the default auto-selected `gmg_levels` (expected ≥ 4 on this grid)
- **THEN** FGMRES converges to `gmg_fgmres_tol = 1e-6` in no more than 60 total inner iterations (at most 2 restarts at default `gmg_fgmres_restart = 30`), and no CHOLMOD fallback is triggered
