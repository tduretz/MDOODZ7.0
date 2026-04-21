# gmg-stokes-solver — delta spec

This change retires the geometric multigrid Stokes solver. All prior normative requirements in `openspec/specs/gmg-stokes-solver/spec.md` are removed. The main spec file is replaced at apply time with a short tombstone that points to archived changes and git history.

## REMOVED Requirements

### Requirement: GMG solver is selected via `lin_solver = 3`

**Reason**: The GMG-preconditioned FGMRES implementation is removed; `lin_solver = 3` is no longer a valid dispatch.

**Migration**: Set `lin_solver` to `0`, `-1`, `1`, or `2` per existing MDOODZ documentation. See archived change `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/` for historical behaviour.

### Requirement: Grid hierarchy is built with a bounded coarsest level

**Reason**: Multigrid hierarchy construction (`MultigridLevels`) is deleted with the GMG removal.

**Migration**: N/A — CHOLMOD direct solve does not use a grid hierarchy.

### Requirement: Restriction preserves constants and averages smooth fields

**Reason**: Restriction operators lived in the removed GMG stack.

**Migration**: N/A.

### Requirement: Prolongation is exact on fields representable on the coarse grid

**Reason**: Prolongation lived in the removed GMG stack.

**Migration**: N/A.

### Requirement: Viscosity is restricted harmonically across contrasts

**Reason**: Coarse-level viscosity restriction lived in the removed GMG stack.

**Migration**: N/A.

### Requirement: Vanka block smoother reduces the residual

**Reason**: Vanka smoothing lived in the removed GMG stack.

**Migration**: N/A.

### Requirement: Coarsest-level solve uses CHOLMOD

**Reason**: The standalone coarse-level GMG solve path is removed; CHOLMOD remains the sole Stokes linear algebra backend for production runs.

**Migration**: Use `lin_solver = 0` or `-1` for CHOLMOD-backed Stokes solves.

### Requirement: The V-cycle is grid-independent for constant viscosity

**Reason**: V-cycle implementation is deleted.

**Migration**: N/A.

### Requirement: GMG-preconditioned FGMRES is the production solver mode

**Reason**: FGMRES-over-GM preconditioning is removed; production Stokes uses existing direct / Powell–Hestenes paths.

**Migration**: Use CHOLMOD (`lin_solver = 0` or `-1`).

### Requirement: Newton Jacobian is handled without polluting the smoother

**Reason**: Newton-mode GMG integration is removed with the solver.

**Migration**: Use existing Newton + direct-solver workflows documented for MDOODZ.

### Requirement: Fine-level operator is bit-identical to MDOODZ's assembled matrix

**Reason**: GMG assembly bridge (`StokesAssemblyGMG`) is deleted.

**Migration**: N/A.

### Requirement: Newton-mode GMG is phased; interim behaviour falls back to CHOLMOD

**Reason**: There is no GMG branch; Newton Stokes uses pre-existing solver dispatch only.

**Migration**: N/A.

### Requirement: Active mask is consistent across multigrid levels

**Reason**: Multigrid levels are removed.

**Migration**: N/A.

### Requirement: GMG produces correct shear-zone localisation

**Reason**: GMG solver is removed.

**Migration**: Validate shear localisation with CHOLMOD-backed runs.

### Requirement: GMG results are consistent with CHOLMOD on existing benchmarks

**Reason**: GMG equivalence tests and solver are removed.

**Migration**: Benchmark against CHOLMOD-only baselines.

### Requirement: Memory scales near-linearly with problem size

**Reason**: This requirement described GMG memory behaviour; GMG is removed.

**Migration**: N/A — CHOLMOD memory characteristics apply instead.

### Requirement: Solver exposes optional tuning parameters with safe defaults

**Reason**: All `gmg_*` input parameters are removed from `params` and ignored at parse time (see ADDED Requirements).

**Migration**: Delete `gmg_*` keys from `.txt` files; if present, they are ignored with a warning.

### Requirement: Solver logs residual history and falls back on non-convergence

**Reason**: GMG logging and fallback-to-CHOLMOD from `lin_solver = 3` are removed; `lin_solver = 3` is now a hard error (see ADDED Requirements).

**Migration**: Do not set `lin_solver = 3`.

### Requirement: V-cycle residual amplification is bounded

**Reason**: V-cycle implementation is deleted (the upleg-fix work is historical only).

**Migration**: N/A.

## ADDED Requirements

### Requirement: GMG Stokes sources are absent from the build

The codebase SHALL NOT compile or link `MultigridStokes.c`, `MultigridLevels.c`, or `StokesAssemblyGMG.c`. The `MDLIB/CMakeLists.txt` target SHALL list none of these translation units.

#### Scenario: Clean build after removal

- **WHEN** a developer configures and builds MDOODZ from a tree containing this change
- **THEN** the link step succeeds without any object file produced from the former GMG translation units

### Requirement: `lin_solver = 3` aborts at startup with a clear error

The input reader SHALL detect `lin_solver = 3` and SHALL terminate the run before time-stepping begins, emitting a `LOG_ERR` message that states GMG was removed, lists valid `lin_solver` values (`-1`, `0`, `1`, `2`), and points readers to `openspec/specs/gmg-stokes-solver/spec.md` and the archived GMG changes.

#### Scenario: Legacy config with GMG dispatch

- **WHEN** a `.txt` input sets `lin_solver = 3`
- **THEN** the program exits with non-zero status and the error message identifies the invalid value and valid alternatives

### Requirement: `gmg_*` input parameters are ignored with a single aggregated warning

If the parser encounters any key whose name matches `gmg_*` (case-insensitive), it SHALL NOT populate removed `params` fields. After parameter parsing completes, the program SHALL emit at most one `WARN` line listing the deprecated parameter family and naming ignored keys, and the run SHALL continue when `lin_solver` is valid.

#### Scenario: Legacy config retains gmg keys with CHOLMOD

- **WHEN** a `.txt` input sets `lin_solver = 0` and includes one or more `gmg_*` keys
- **THEN** the program emits the aggregated warning and completes initialization without treating `gmg_*` keys as errors

### Requirement: Published `gmg-stokes-solver` spec is a retirement tombstone

The file `openspec/specs/gmg-stokes-solver/spec.md` SHALL be replaced with a short document that states **RETIRED**, summarises the negative-outcome rationale (reference `benchmarks/ec2/results/PERFORMANCE_REPORT.md` §§11–13 as the measurement record), lists archived changes `2026-04-19-add-gmg-stokes-solver`, `2026-04-21-add-gmg-stokes-defence`, and `2026-04-18-add-gmg-upleg-fix`, and states that CHOLMOD is the production Stokes solver. It SHALL NOT restate the removed normative requirements as active contracts.

#### Scenario: Reader opens main spec after apply

- **WHEN** a reader opens `openspec/specs/gmg-stokes-solver/spec.md` on the branch containing this change
- **THEN** the document clearly indicates retirement in the first screen and does not read as an active implementation specification
