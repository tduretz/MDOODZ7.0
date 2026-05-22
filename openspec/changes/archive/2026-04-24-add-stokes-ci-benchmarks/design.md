## Context

The proposal initially scoped three ASPECT Stokes benchmarks (SolCx, SolKz, Donea–Huerta) for MDOODZ's CI. A pre-design audit of [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h) and [MDLIB/mdoodz-private.h:30-42](MDLIB/mdoodz-private.h#L30-L42) revealed:

- **SolCx** needs spatially-varying density `ρ(x,z) = sin(πz)·cos(πx)`. MDOODZ's body-force path uses per-phase `materials.rho[p]` via `EvaluateDensity()` in [RheologyDensity.c](MDLIB/RheologyDensity.c) — the `SetDensity` callback writes `particles->rho` but that value is overwritten by grid-to-particle interpolation. The two-phase approach used by SolVi (distinct constant `eta0` per phase) doesn't work for SolCx's sinusoidal density. A small library addition is needed.
- **SolKz** needs exponentially-varying viscosity `η(x,z) = exp(2Bz)`. MDOODZ's viscosity lives on the grid (`mesh->eta_n/s`) recomputed per iteration from phase flow laws. Initial attempts revealed a persistent magnitude discrepancy (3-20× on pressure) between MDOODZ's numerical solution and the Velic analytical that could not be resolved in this change's debugging window. **SolKz is deferred to a follow-up change** along with its required `SetViscosity` callback.
- **Donea–Huerta** needs an arbitrary 2D body force that cannot be factored as `ρ·(gx, gz)` with constant gravity — deferred.

This change therefore scopes in **one** library addition — a `SetGridDensity` callback that overrides the per-cell grid density `mesh->rho_n[c]` / `mesh->rho_s[c]` every rheology pass — to enable SolCx.

Relevant existing surface:
- `SetPhase_f(MdoodzInput*, Coordinates) → int` — picks a material phase per marker
- `SetDensity_f(MdoodzInput*, Coordinates, int phase) → double` — writes `particles->rho` (not used by body force)
- Gravity enters momentum RHS as `-gz * mesh->rho_n` in [RheologyDensity.c:1066](MDLIB/RheologyDensity.c#L1066) via `EvaluateRHS` (callers in `StokesRoutines.c`)
- Grid density `mesh->rho_n` is set by `UpdateDensity` ([RheologyDensity.c:1879](MDLIB/RheologyDensity.c#L1879)) from per-phase flow-law-derived values

The existing SolVi benchmark ([TESTS/SolViBenchmarkTests.cpp](TESTS/SolViBenchmarkTests.cpp)) fits the callback surface because it uses two phases (matrix, inclusion) with distinct `eta0` values — the same fit SolCx would exploit for the viscosity jump. It's only the density pattern that requires new library support.

## Goals / Non-Goals

**Goals:**
- Add **SolCx** as a CI test with L2(V), L1(P) assertions at 51×51 and a 3-resolution grid-convergence-order assertion.
- Introduce a **`SetGridDensity` callback** in MDLIB with per-iteration grid-override semantics: fires inside `UpdateDensity` and `NonNewtonianViscosityGrid` to replace `mesh->rho_n/s` with the callback's return value. Matches ASPECT's "material model" pattern for spatially-varying density.
- Share one staggered-grid coordinate mapping helper across the new benchmark, using the corrected formula from the April 2026 SolVi audit.
- Keep total CI time added < 30 s (SolCx is 4 single-step Stokes solves, ~10 s total).

**Non-Goals:**
- Porting ASPECT's non-Stokes benchmarks (Rayleigh–Taylor, Crameri Case 2, slab detachment, Tosi 2015) — deferred to separate changes.
- Reaching ASPECT's absolute error levels. Their benchmarks use higher-order FE; MDOODZ is 2nd-order FD on staggered grids. We assert MDOODZ-realistic thresholds (2–5× above measured).
- Fixing the Popov2025 ASSERT-snippet doc issue (separate change).
- **SolKz** (smooth exp-viscosity Stokes) and its required `SetViscosity` callback. Debugging revealed a persistent 3-20× pressure-magnitude discrepancy that couldn't be resolved without deeper MDOODZ-internals work. SolKz + `SetViscosity` will ship as a dedicated follow-up change where the solver-interaction can be investigated on its own merits.
- **Donea–Huerta** and its required `SetBodyForce` callback — same "separate follow-up" argument.
- Per-phase opt-in for `SetGridDensity`. The override is whole-simulation scope (validation refuses mixing with any non-zero `density_model`). YAGNI for now.
- Re-evaluating `SetGridDensity` on state (strain rate, temperature, pressure). The callback signature is `(MdoodzInput*, Coordinates, int phase) → double` — pure function of position and phase.

## Decisions

### Decision 1: SolCx is in scope

SolCx specifies `η(x) = 1` for `x < 0.5` and `η(x) = 10⁶` for `x ≥ 0.5` (canonical ASPECT headline contrast), plus `ρ(x,z) = sin(πz)·cos(πx)` and gravity `(0, 1)`. This maps to:
- Two phases via `SetPhase` returning `0` for `x < 0.5` and `1` otherwise, with `eta0[0] = 1`, `eta0[1] = 1e6` in the `.txt` file.
- `SetGridDensity` returning `sin(π·z)·cos(π·x)` (see Decision 2).
- Gravity `gx = 0, gz = 1` in the `.txt` file.
- `Nt = 1`, `Nx = Nz ∈ {21, 41, 81}` for the convergence sweep, primary single-resolution case at 51×51.

**Analytical solution:** verbatim port of Velic's SolCx Chebyshev-series evaluator from ASPECT `benchmarks/solcx/solcx.h` (GPL v2+, originally from Underworld). Lives at `TESTS/SolCx/SolCxAnalytical.cpp` with full attribution — `numbers::PI → M_PI` substitution and namespace wrapping, otherwise unchanged.

**Assertions (measured values at 51×51):**
- `L2(Vx) < 2.5e-1` — measured 7.4e-2 (3.4× margin)
- `L2(Vz) < 2.5e-1` — measured 8.5e-2 (2.9× margin)
- `L1(P) < 1e-2`    — measured 2.3e-3 (4.4× margin)

**Assertions (convergence order 41→81):**
- `order(Vx) ≥ 0.7` — measured 0.98
- `order(Vz) ≥ 0.7` — measured 0.99
- `order(P_L1) ≥ 0.5` — measured 2.57 (near-cubic — consistent with Duretz et al. 2011 §4)

L1 (not L2) for pressure because the analytical pressure crosses zero inside the domain and relative L2 becomes ill-conditioned — lesson from the SolVi audit.

### Decision 2: Ship SolCx by adding a `SetGridDensity` callback

MDOODZ's body-force path reads `mesh->rho_n` / `mesh->rho_s`, which are computed in `UpdateDensity` from per-phase constants via `EvaluateDensity()` — the `SetDensity` callback's return value is stored on markers but overwritten by grid-to-particle interpolation. SolCx's sinusoidal density cannot be expressed as per-phase constants.

We add:

```c
typedef double (*SetGridDensity_f)(MdoodzInput *input, Coordinates coordinates, int phase);
```

Registered in `SetParticles_ff` alongside `SetPhase`, `SetDensity`, etc. in [mdoodz.h](MDLIB/include/mdoodz.h).

**Semantics — per-iteration grid override:**
- If `SetGridDensity == NULL` (the default for every current scenario): current behaviour preserved entirely.
- If `SetGridDensity != NULL`: every time `UpdateDensity` or `NonNewtonianViscosityGrid` finalises `mesh->rho_n[c]` and `mesh->rho_s[c]`, the computed value is replaced by `SetGridDensity(input, {x, z}, phase)` using the cell/vertex coordinates and a representative phase from that cell.
- **Startup validation**: in [Main_DOODZ.c](MDLIB/Main_DOODZ.c), if `SetGridDensity != NULL`, verify that for every declared phase, `density_model == 0` (the default constant-density branch). If any phase has `density_model ≠ 0` → `LOG_ERR` + `exit`. Same loud-failure pattern as existing `Slim >= C` / `coh_soft` checks. The two sources of density cannot sensibly coexist — the callback always wins, so mixing is silent work thrown away.

**Why grid-override (rejected alternatives):**
- **Per-marker `particles->rho` via `SetDensity`**: rejected — the marker `rho` is overwritten by grid-to-particle interpolation each step (Main_DOODZ.c:369), so the callback's value doesn't reach the body force.
- **Multi-phase staircase approximation**: rejected — bounded by the 20-phase cap in MDLIB's material arrays, and staircase artefacts would pollute the convergence measurement.
- **Temperature-EoS abuse** (set `density_model=1` with `SetTemperature` encoding the pattern): rejected — conflates temperature with density, requires thermal solver machinery, and obscures what the test is doing.

**Why this matches ASPECT conventions:** ASPECT has a material-model plugin system where `evaluate()` is called at every quadrature point to return density (and viscosity). Our `SetGridDensity` callback is functionally equivalent — a user-provided function that's consulted every time the solver needs density. The "override" framing is MDOODZ-architectural (because MDOODZ caches density on the grid); the FUNCTION is standard.

**Integration points:**
- Typedef + struct field in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h) — `SetGridDensity_f` + field in `SetParticles_ff`.
- Module-level static state + register function in [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c): `RegisterSetGridDensityCallback(input, cb)` called once from `RunMDOODZ`.
- Override applied in `ApplySetGridDensityOverride(mesh, model)` at end of `UpdateDensity` and `NonNewtonianViscosityGrid` (and the Aniso variant). Null-check guards a no-op when callback unset.
- Declarations in [MDLIB/mdoodz-private.h](MDLIB/mdoodz-private.h).
- Startup validation + register call in [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c) — placed after `ReadInputFile` returns and the `MdoodzInput` is built, before any solver work.

### Decision 3: Donea–Huerta and SolKz are deferred

DH requires a `SetBodyForce` callback (arbitrary 2D RHS); SolKz requires a `SetViscosity` callback AND a resolution of the numerical mismatch identified during debugging (Newton residual plateau, 3-20× pressure magnitude discrepancy). Both are separate follow-up changes using the `SetGridDensity` pattern established here as the template.

### Decision 4: Test harness follows SolVi precedent

- One GTest suite (`SolCxBenchmarkTests.cpp`), fixture class mirroring `SolViBenchmark`.
- Subdirectory for `.txt` fixture files: `TESTS/SolCx/SolCx21.txt`, `SolCx41.txt`, `SolCx51.txt`, `SolCx81.txt`.
- Reuse `computeL2Error`, `computeL1Error`, `readFieldAsArray`, `readCoordArray` from [TESTS/TestHelpers.h](TESTS/TestHelpers.h).
- Analytical-Dirichlet BCs on all four walls (same pattern as SolVi) — anchors the solution and gives the cleanest convergence measurement.
- Each test runs with `writer_subfolder = <TestName>`, `coh_soft = 0`, `Nt = 1`, no thermal coupling, no elasticity, no plasticity, `density_model = 0` on every phase.

### Decision 5: Use analytical BCs

Each benchmark imposes its analytical velocity as Dirichlet boundary conditions via the E/W/N/S callbacks. For SolCx this means evaluating the Velic reference solution on the boundary. No free-slip BC helpers — benchmarks need exact-analytic BCs or the interior solution is inconsistent with the BC-imposed solution and L2 no longer measures the solver's internal accuracy.

## Risks / Trade-offs

- **[SolCx pressure order may underperform]** → Same as SolVi at the jump: we expect P order ≤ 1 for relative L2, hence we use L1. Mitigation: assert `order(P_L1) ≥ 0.5`, measured 2.57 — comfortable margin.

- **[Velic SolCx analytical solution is 2900 lines of Maple-generated code]** → Non-trivial port. Mitigation: taken verbatim from ASPECT (GPL v2+) with attribution preserved; changes limited to `numbers::PI → M_PI` and namespace wrapping. Cross-validated by self-consistency test (divergence-free, no NaN/Inf, finite values).

- **[Phase attribution in the grid-override]** → `SetGridDensity` takes a `phase` argument, but `mesh->rho_n[c]` is a cell-aggregate derived from multiple markers potentially of different phases. We pass a representative phase (dominant by percentage) to the callback. For SolCx (two phases) each cell is dominated by one of them; for the callback's general contract, users writing position-dependent callbacks should not depend on the exact phase when the domain has multi-phase cells. Acceptable scope for now.

- **[CI runtime]** → SolCx: 4 resolutions × single-step solve ≈ 10 s total. Under budget.

- **[Scope expansion via library change]** → Adding `SetGridDensity` touches [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h), [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c), [MDLIB/AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c), and [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c). Diff is small (~60 LOC) — one typedef, one struct field, module-level state + register/apply functions, two call sites, one validation block. Isolated behind null-check so existing scenarios are unaffected. Confirmed backwards-compatible by full CI run.

- **[Misusing gravity direction]** → SolCx needs `gz = 1` with `ρ = +sin(πz)·cos(πx)` — empirically verified by matching the analytical (test passes with 3-20× margin on all assertions).

- **[Parameter-file stale copy]** → CMake copies fixture `.txt` files at configure time only. Documented; task list reminds to re-run `cmake -B` after editing fixtures.

## Migration Plan

`SetGridDensity` is additive and backwards-compatible. The field defaults to `NULL` because `SetParticles_ff` is typically initialised with designated initialisers (current scenarios only populate the subset they use); unset fields are zero-initialised in C, so existing scenarios continue to run unchanged. No schema migration, no `.txt` file changes, no rebuild of checkpoint files. Confirmed by running the full CI suite after the changes — 23/24 tests pass; the one failure (`RotationAdvectionTests`) is a pre-existing flake that also fails on the baseline.

## Open Questions

None — all five prior open questions have been resolved:
- **Q1 (SolCx contrast)**: 10⁶ only, primary case at 51×51, no 10³ secondary.
- **Q2 (SolCx P-convergence floor)**: 0.5 on L1(P) (permissive; measured 2.57).
- **Q3 (SolKz mode numbers)**: moot — SolKz deferred.
- **Q4 (SolKz B value)**: moot — SolKz deferred.
- **Q5 (Callback scope — thermal/elastic/plast)**: validation restricted to `density_model != 0` only. Thermal / elastic / plast flags can legitimately coexist with a prescribed-density Stokes problem in future use cases. For SolCx all three are zero anyway.
