## Why

MDOODZ's Stokes solver has one analytical CI benchmark — SolVi (viscous circular inclusion) — which covers a single regime (smooth radial viscosity contrast). The audit in [skill-testing-guide/SKILL.md](.claude/skills/skill-testing-guide/SKILL.md) showed that SolVi's P-convergence stalls at ~0.75 (L2) and most improvement attempts (marker density, `eta_average`, D-tensor harmonic averaging) were ruled out. We therefore don't know whether the solver behaves correctly under a **discontinuous viscosity jump** with spatially-varying buoyancy forcing. ASPECT ships exactly such a benchmark — SolCx (Zhong 1996 / Duretz et al. 2011) — with independent reference values from multiple codes. Adding it closes this coverage gap at low cost.

Implementing SolCx in MDOODZ revealed that the solver's default body-force path uses per-phase density from `materials.rho[p]`, not the spatially-varying `particles->rho` set by `SetDensity`. A small library addition (a `SetGridDensity` callback) is therefore required so that SolCx's `ρ(x,z) = sin(π·z)·cos(π·x)` density pattern reaches the grid.

## What Changes

- Add `TESTS/SolCxBenchmarkTests.cpp` — Zhong (1996) step-viscosity benchmark with 10⁶ jump at x=0.5, sinusoidal density; L2(Vx), L2(Vz), L1(P) assertions at 51×51 plus grid-convergence order check across 21/41/81.
- Add `TESTS/SolCx/` directory with `.txt` parameter files (one per resolution) plus the ported Velic analytical evaluator (`SolCxAnalytical.h/.cpp`, ~2900 LOC verbatim from ASPECT with GPL v2+ attribution).
- **Add a `SetGridDensity` callback to MDLIB** ([mdoodz.h](MDLIB/include/mdoodz.h) typedef + `SetParticles_ff` field) with **per-iteration grid-override semantics**: when non-null, it replaces the computed `mesh->rho_n[c]` / `mesh->rho_s[c]` on every rheology pass with the callback's return value evaluated at the cell / vertex coordinate. This gives MDOODZ an ASPECT-style "material model" hook for spatially-varying density that cannot be expressed as a per-phase constant.
- **Add startup validation** in [Main_DOODZ.c](MDLIB/Main_DOODZ.c): when `SetGridDensity != NULL`, every declared phase must have `density_model = 0`. Mixed configuration produces `LOG_ERR` + `exit` — same loud-failure pattern as existing `Slim >= C` / `coh_soft` checks.
- Wire the callback into the end of `UpdateDensity` and the `NonNewtonianViscosityGrid`/`NonNewtonianViscosityGridAniso` routines in [RheologyDensity.c](MDLIB/RheologyDensity.c) / [AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c).
- Extend [TESTS/TestHelpers.h](TESTS/TestHelpers.h) with `subtractMean` helper (added for potential future use; not required by SolCx itself).
- Wire into [TESTS/CMakeLists.txt](TESTS/CMakeLists.txt) — one new executable, added to `run-tests` target.
- Add a **log-log convergence PNG figure** for SolCx via a committed gnuplot script, embedded in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md). Reads a `solcx_convergence.dat` file that the test writes during its `GridConvergence` case.
- **Deferred**: SolKz (smooth exponential viscosity) and Donea–Huerta — initial implementation of SolKz revealed a persistent magnitude discrepancy between MDOODZ's numerical solution and the Velic analytical (3-20× on pressure, 20-30% on velocity; solver residual plateaus at ~1e-4). Root cause appears to be a subtle interaction between MDOODZ's Newton Jacobian and the override path. SolKz is filed as a dedicated follow-up change along with its required `SetViscosity` callback.

## Capabilities

### New Capabilities
- `benchmark-solcx`: Analytical Stokes benchmark with a step-function viscosity jump (10⁶ contrast at x=0.5). Validates the FD solver's handling of discontinuous η at cell faces, with a sinusoidal density pattern driving buoyant flow.
- `set-grid-density-callback`: New `SetGridDensity_f` callback in `SetParticles_ff` with per-iteration grid-override semantics — overrides `mesh->rho_n` / `mesh->rho_s` with `callback(x, z, phase)` at every iteration when set. Startup validation rejects the combination of the callback with any non-zero `density_model` to prevent silent ambiguity. Enables spatially-varying density setups that can't be expressed as per-phase constants (e.g., SolCx's sinusoidal buoyancy).

### Modified Capabilities
None.

## Impact

- **Affected code**:
  - [TESTS/](TESTS/) — new test files, new fixture `.txt` files, [TESTS/CMakeLists.txt](TESTS/CMakeLists.txt) updates.
  - [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h) — add `SetGridDensity_f` typedef + field in `SetParticles_ff`.
  - [MDLIB/mdoodz-private.h](MDLIB/mdoodz-private.h) — declarations for `RegisterSetGridDensityCallback` and `ApplySetGridDensityOverride`.
  - [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c) — state statics, register/apply functions, invocation at end of `UpdateDensity` and `NonNewtonianViscosityGrid`.
  - [MDLIB/AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) — invocation at end of `NonNewtonianViscosityGridAniso`.
  - [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c) — startup validation block + register call.
- **Backwards compatibility**: `SetGridDensity` field is `NULL` by default (zero-initialised in C); no existing scenario sets it, so behaviour is unchanged. No `.txt` schema change. No checkpoint-file change. Validated by full CI run — 23/24 existing tests pass (the one failure, `RotationAdvectionTests`, is a pre-existing flake on this branch, confirmed to fail on baseline too).
- **CI budget**: SolCx runs 4 single-step Stokes solves (21, 41, 51, 81 grids). Measured total runtime ~10 s. Under budget.
- **Documentation**:
  - Update [skill-testing-guide/SKILL.md](.claude/skills/skill-testing-guide/SKILL.md) coverage table with SolCx and a new pitfall entry documenting `SetGridDensity` semantics (per-iteration override; validation refuses mixing with `density_model ≠ 0`).
  - **Add a new section to [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md)** for SolCx, mirroring the existing `§1 SolVi` pattern: Problem Statement, Analytical Solution (with math), Implementation Details, Measured L2 Errors table (per resolution), Convergence Orders table, and Code Assertions snippets.
- **Out of scope**: SolKz (deferred — see separate follow-up), Donea–Huerta (requires a body-force callback — separate follow-up), other ASPECT benchmarks (Rayleigh-Taylor van Keken, Crameri Case 2, McKenzie–Jackson finite strain, slab detachment), and the Popov2025 ASSERT-snippet documentation fix — each should be its own change.
