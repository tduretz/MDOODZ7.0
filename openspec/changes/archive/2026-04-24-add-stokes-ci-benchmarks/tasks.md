## 1. MDLIB library — SetGridDensity callback

- [x] 1.1 Added `SetGridDensity_f` typedef + `SetGridDensity` field in `SetParticles_ff` in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h), placed alongside `SetDensity_f`.
- [x] 1.2 Added doc comment block next to the typedef stating per-iteration grid-override semantics, the validation rule, and the intended use.
- [x] 1.3 Audit finding recorded: `mesh->rho_n` is populated by `UpdateDensity` ([RheologyDensity.c:1879](MDLIB/RheologyDensity.c#L1879)) from per-phase constants via `EvaluateDensity()`. The `SetDensity` callback's value on `particles->rho` is overwritten by grid-to-particle interpolation and doesn't reach the body force. Callback therefore targets the grid.
- [x] 1.4 Added `RegisterSetGridDensityCallback` + `ApplySetGridDensityOverride` in [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c) — module-level static holds callback + `MdoodzInput*`; override loops over centres (`mesh->rho_n`) and vertices (`mesh->rho_s`), using `xc_coord` / `zc_coord` / `xg_coord` / `zg_coord` and `phase_perc_{n,s}` max-percentage for representative phase. Invoked at end of `UpdateDensity`, `NonNewtonianViscosityGrid`, and the Aniso variant. Null-check guards a no-op when callback unset. Declarations added to `mdoodz-private.h`.
- [x] 1.5 Added startup validation in [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c): if `SetGridDensity != NULL`, every phase must have `density_model == 0`; `LOG_ERR` + `exit(1)` on conflict. Pattern matches existing `Slim >= C` check.
- [x] 1.6 Build clean (`libmdoodz.dylib` + all test executables build). Test suite: 23/24 pre-existing tests pass. The 1 failure (`RotationAdvectionTests`, labelled `experimental`) is **pre-existing** — verified by running against the stashed baseline. Backwards compatibility preserved.

## 2. SolCx analytical solution port

- [x] 2.1 Created `TESTS/SolCx/SolCxAnalytical.h` + `TESTS/SolCx/SolCxAnalytical.cpp` — verbatim port of ASPECT `benchmarks/solcx/solcx.h` `_Velic_solCx` (GPL v2+, originally from Underworld). 2891 LOC with full attribution preserved. Transformations: `numbers::PI` → `M_PI`, dedent by 6, wrapped in `MdoodzSolCx::detail`, thin `EvalSolCx(x, z, eta_A, eta_B, xc, n) → SolCxValue` facade exposed.
- [x] 2.2 Self-consistency test inside `SolCxBenchmarkTests.cpp` — evaluates analytical at 5 sample points, asserts no NaN/Inf, divergence-free (`exx + ezz ≈ 0`), and finite values. Passes.

## 3. Test helpers

- [x] 3.1 `computeL1Error` already present at [TestHelpers.h:237-244](TESTS/TestHelpers.h#L237) — no change needed.
- [x] 3.2 Added `subtractMean(std::vector<double>&)` to [TestHelpers.h](TESTS/TestHelpers.h) — in-place spatial-mean subtraction. Useful helper (not required by SolCx itself; kept for future benchmarks).
- [x] 3.3 Staggered-grid coordinate handling already available via `readCoordArray` ([TestHelpers.h:198](TESTS/TestHelpers.h#L198)). No extraction needed.

## 4. SolCx benchmark

- [x] 4.1 Created `TESTS/SolCx/SolCx21.txt`, `SolCx41.txt`, `SolCx51.txt`, `SolCx81.txt` with `density_model=0` on every phase (required by SetGridDensity validation), all flow-law indices zero, `eta0[0]=1`, `eta0[1]=1e6`, `gz=1`, domain [0,1]², `Nt=1`, `thermal=0`, `elastic=0`, `plast=0`, `coh_soft=0`.
- [x] 4.2 Created [TESTS/SolCxBenchmarkTests.cpp](TESTS/SolCxBenchmarkTests.cpp) with `SetPhase` (0 if x<0.5 else 1), `SetGridDensity` (returns `sin(π·z)·cos(π·x)`), `SetBCVx`/`SetBCVz` Dirichlet from the Velic analytical on all four walls, and `BuildAnalytical` helper reading HDF5 coordinates to build analytical vectors matching the staggered-grid layout.
- [x] 4.3 `L2Error` GTest case at 51×51: asserts `L2(Vx) < 2.5e-1`, `L2(Vz) < 2.5e-1`, `L1(P) < 1e-2`. Measured 7.4e-2 / 8.5e-2 / 2.3e-3 (3–4× margin).
- [x] 4.4 `GridConvergence` GTest case at 21/41/81: asserts `order(Vx) ≥ 0.7`, `order(Vz) ≥ 0.7`, `order(P_L1) ≥ 0.5`. Measured 0.98 / 0.99 / 2.57 (41→81 pair) — consistent with Duretz et al. 2011 §4.
- [x] 4.5 Added `SolCxBenchmarkTests` to [TESTS/CMakeLists.txt](TESTS/CMakeLists.txt) with `file(COPY TESTS/SolCx ...)` and linked to `mdoodz` + `GTest::gtest_main`. All three GTest cases pass.

## 5. Documentation — AnalyticalSolutions.md

- [x] 5.1 Added new section "SolCx Benchmark (Step Viscosity)" as §7 in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) with Problem Statement, Analytical Solution (Velic Chebyshev reference with GPL attribution), Implementation Details, Measured L2 Errors table (21/41/51/81), Convergence Orders table, Code Assertions with inline margin comments, and embedded convergence PNG.
- [x] 5.2 Added six SolCx rows to the "Summary of L2 Assertions" table: L2(Vx) / L2(Vz) / L1(P) at 51×51 plus the three convergence-order thresholds.

## 6. Gnuplot figure

- [x] 6.1 Wrote [TESTS/SolCx/plot_solcx.gp](TESTS/SolCx/plot_solcx.gp) — reads `solcx_convergence.dat` (written by the test during `GridConvergence`) and produces a log-log convergence PNG with slope-1 and slope-2 reference lines. Scope reduced from the originally-planned 3-panel layout (numerical heatmap + error map + convergence) since the field-visualisation panels require an HDF5-to-.dat extractor step (substantial extra code). Convergence-only PNG is sufficient for CI benchmark documentation; field/error panels are a follow-up polish item.
- [x] 6.2 Added `solcx_convergence.dat` write in [SolCxBenchmarkTests.cpp](TESTS/SolCxBenchmarkTests.cpp) — `GridConvergence` appends `h L2(Vx) L2(Vz) L1(P)` per resolution at test end.
- [x] 6.3 Ran the script, generated and committed [TESTS/SolCx/solcx_convergence.png](TESTS/SolCx/solcx_convergence.png) (~42 KB).
- [x] 6.4 Embedded the PNG in `TESTS/AnalyticalSolutions.md` §7 after "Code Assertions" with regeneration instructions and a note about the L1(P) 21→41 kink (known staggered-FD behaviour at η-jumps).
- [~] 6.5 SolVi figure retrofit — **deferred** out of this change. Requires separate work on SolVi's test to write a convergence .dat file. Follow-up issue: add a matching convergence plot to the existing SolVi section for visual consistency across Stokes benchmarks.

## 7. Skill testing guide + pitfalls

- [x] 7.1 Added `SolCxBenchmarkTests` rows (3 tests) to the "Existing CI Tests" coverage table in [.claude/skills/skill-testing-guide/SKILL.md](.claude/skills/skill-testing-guide/SKILL.md).
- [x] 7.2 Added two new pitfall entries in the same file: (a) "`SetGridDensity` callback is a whole-simulation grid override" explaining `SetDensity` doesn't reach body force, `mesh->rho_n/s` is the authoritative source, callback overrides it, and validation refuses `density_model ≠ 0`. (b) "MDOODZ pressure output includes hydrostatic from initial lithostatic integration" explaining `Centers/P = total pressure`, enable `track_T_P_x_z = 1` to access `Centers/P0` for dynamic-only comparison.

## 8. Final verification

- [x] 8.1 Full build clean (`mdoodz` + all test executables).
- [x] 8.2 Full CI suite: 23/24 tests pass. Only failure is pre-existing `RotationAdvectionTests` flake (labeled `experimental`, also fails on baseline — unrelated to this change). `SolCxBenchmarkTests` passes in ~10 s with all three test cases green.
- [x] 8.3 `openspec validate add-stokes-ci-benchmarks` → clean.
- [x] 8.4 Spec-to-implementation cross-check: every `#### Scenario:` maps to a concrete artefact:
  - `benchmark-solcx/spec.md` scenarios → GTest cases in `SolCxBenchmarkTests.cpp`, `.txt` fixtures in `TESTS/SolCx/`, `AnalyticalSolutions.md` §7 section, committed PNG
  - `set-grid-density-callback/spec.md` scenarios → typedef/field in `mdoodz.h`, `Register`/`Apply` functions in `RheologyDensity.c`, validation block in `Main_DOODZ.c`, pitfall entry in skill-testing-guide
