## Context

The combined mode-I / mode-II yield surface of Popov et al. (2025, GMD) is implemented at [MDLIB/RheologyDensity.c:685-899](MDLIB/RheologyDensity.c#L685-L899), selected by `materials.plast[phase] = 2` (see [MDLIB/RheologyDensity.c:456-457](MDLIB/RheologyDensity.c#L456-L457)). The local stress update is a 3×3 Newton solve on (τII, P, λ̇) with an optional Armijo line search behind the `tensile_line_search` global switch ([MDLIB/InputOutput.c:1200](MDLIB/InputOutput.c#L1200)). Material inputs are the four standard DP parameters plus the tensile strength `T_st` ([MDLIB/InputOutput.c:1379](MDLIB/InputOutput.c#L1379)).

The HDF5 output already emits the fields needed for paper-style post-processing at the cell centres: `/Centers/eII_pl` (deviatoric plastic strain rate) and `/Centers/divu_pl` (volumetric plastic strain rate), both scaled back to SI before writing. Integrating these in time recovers Popov's κ (Eq. 23) and χ (Eq. 25) invariants.

Two scenario pairs landed with the PR — `SETS/Popov2025_Pureshear_VEVP.{c,txt}` and `SETS/Popov2025_Tensile_VEVP.{c,txt}` — but their parameter files configure `plast = 0` (tensile) or `plast = 1` (pureshear) and therefore exercise the legacy separate-path code, not the new combined yield surface. Only `SETS/PressurizedMagmaChamber.txt` actually sets `plast = 2` in the current tree. No automated regression currently pins the behaviour of the new solver.

The existing CI test infrastructure ([TESTS/CMakeLists.txt](TESTS/CMakeLists.txt)) follows a one-binary-per-topic pattern with GTest fixtures, HDF5 assertion helpers in [TESTS/TestHelpers.h](TESTS/TestHelpers.h), and multi-resolution convergence tests (SolVi, Anisotropy, TopoBench) as precedent for what a "benchmark" test looks like. `TESTS/AnalyticalSolutions.md` is the canonical home for analytical reference formulas; `TESTS/README.md` enumerates suites by number.

Gnuplot is already the project's standard for HDF5-driven plot scripts — `TESTS/BlankenBench/extract_fields.cpp` is a standalone HDF5 → ASCII extractor that feeds gnuplot, and ten `.gnu` scripts live under `VISUAL_TESTS/` for the visual-regression pipeline. For paper-figure reproduction (non-regression), the project convention is to co-locate plot assets with the test fixtures and leave them out of the regression pipeline.

## Goals / Non-Goals

**Goals:**
- Pin the Popov 2025 combined yield-surface code path (`plast = 2`) with regression coverage that fails loudly if the 3×3 Newton solve, the yield-surface branching condition `τII(py − pd) ≷ τd(py − P)`, the cap-segment flow-potential derivatives, or the DP-segment clamping drift.
- Keep the CI footprint tiny — all three 0D tests run in ≈ 2 s total, well under the fast-CI budget.
- Fix the two mislabeled `SETS/Popov2025_*_VEVP.*` scenarios so that running them manually actually exercises the combined yield surface (they previously used the legacy `plast = 0 / 1` paths).
- Ship a single gnuplot script (`fig5_0d_stress.gp`) co-located with the test fixtures under `TESTS/Popov2025/plots/`, fed by an HDF5 → ASCII extractor modelled on `TESTS/BlankenBench/extract_fields.cpp`, that reproduces paper Fig. 5 a–d on demand.
- Document the yield-surface geometry, the reference integrator, the test-to-assertion mapping, and the list of deferred 2D tests in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md); register the new suite in [TESTS/README.md](TESTS/README.md).

**Non-Goals:**
- 2D localisation tests (paper Figs. 6–9). These were prototyped during development and deferred: the diagnostic signals (FWHM ratios, Arthur's dip angle, monotone tip propagation, depth-partitioned χ/κ) require resolutions well above the fast-CI budget to resolve. The 0D tier already exercises the constitutive-law code path end-to-end, which is the primary regression risk for the new yield surface. A follow-up change can add a 2D tier behind an `experimental` / `DISABLED` CMake gate if needed.
- Porting the paper's finite-element Crouzeix-Raviart discretization. MDOODZ retains its staggered FD scheme; only the constitutive law is shared.
- Porting the finite-element adaptive time stepping (Popov §4.6). MDOODZ uses its standard fixed-`dt` + nonlinear-iteration strategy; the 0D tests use `dt = 2 yr` matching the paper.

## Decisions

### D1. One GTest binary `Popov2025Tests` with five TEST_F fixtures
One binary keeps CMake wiring simple and matches the convention used by every other existing suite. Each paper test family becomes one fixture class (`Popov2025_0DIntegration`, `Popov2025_Regularization`, …) so that CI logs remain readable and individual sub-tests can be skipped or labeled `experimental` independently.

**Alternatives considered.** Five separate binaries — rejected because they multiply CMake boilerplate and CI runner setup overhead for no logical benefit.

### D2. Fast-CI vs experimental partitioning

| Suite | Resolution / Nt | Wall time (budget) | Label |
|---|---|---|---|
| 0DIntegration | 3 tests × 3×3 grid × 10-25 steps | ≤ 5 s total | default |
| Regularization | 2 tests × 2 resolutions (21×15, 41×29) × ~200 steps | ≤ 60 s total | default |
| Crust | 1 test × 41×15 × ~40 steps | ≤ 30 s | default |
| TensileFracture | 1 test × 41×15 × ~30 steps | ≤ 30 s | default |
| BrittleDuctile | 1 test × 41×21 × ~100 steps | ≤ 5 min | **experimental** |

Only BrittleDuctile is labeled `experimental` + `DISABLED TRUE` — same mechanism as `BlankenBenchTests` at [TESTS/CMakeLists.txt:142](TESTS/CMakeLists.txt#L142). Rationale: BDT needs dislocation creep with realistic Olivine-like scaling plus plasticity, which makes the Newton iteration count higher and the total wall time harder to cap. A developer who runs the binary directly gets the BDT test via GTest's standard `--gtest_also_run_disabled_tests` flag — no in-code `GTEST_SKIP()` guard is added on top of the CMake gate, since that would make the standard GTest opt-in flag silently do nothing and confuse anyone reaching for it.

**Alternatives considered.** Gating on a CI-only env var — rejected because it diverges from existing precedent (the suite uses `PROPERTY DISABLED` + `LABELS experimental` throughout). Double-gating with both CMake `DISABLED TRUE` and an in-code `GTEST_SKIP()` — rejected because the two mechanisms fight each other: `--gtest_also_run_disabled_tests` would still see the test skipped unless a second env var is also set.

### D3. 0D test: homogeneous 3×3 grid with custom pure-shear BCs

MDOODZ has no integration-point driver, so "0D" is emulated with a homogeneous grid small enough that only one interior cell exists (3×3 → one centre). Custom `SetBCVx`/`SetBCVz` (already present in the shipped `Popov2025_*_VEVP.c` files) apply:
- **Volumetric extension**: Vx = +ε̇_vol·x/2 on all four boundaries (ε̇_vol > 0 outward) — divergence = ε̇_vol, deviator = 0. Hits the mode-I cap.
- **Deviatoric shear**: Vx = −ε̇_dev·x, Vz = +ε̇_dev·z (standard pure shear) — divergence = 0, deviator ε̇_dev. Hits the Drucker-Prager shear envelope.
- **Mixed**: superposition — both non-zero.

Reference τII(t), P(t) curves are produced by integrating the visco-elastic Maxwell branch to yield, then the Perzyna ODE along the active yield segment. These are implemented as closed-form piecewise expressions (elastic loading → yield onset → Perzyna asymptote) in a small helper `Popov2025/Popov0DAnalytical.cpp` and compared to HDF5 `/Centers/sII` and `/Centers/P` via `computeL2Error()`.

**Alternatives considered.** Hard-coded expected values at 3-5 checkpoints — rejected because it hides the physics of the test and makes tolerance calibration arbitrary. An analytical helper is transparent and reusable.

### D4. Regularization test: grid-convergence-of-width, not of field
For the 2D tensile / shear localization (Fig. 6), we do not have an analytical solution for the localization profile, but we do have a theoretical criterion: with `eta_vp = 0` the localization band width equals one cell (mesh-pathological); with `eta_vp > 0` it equals an intrinsic length scale independent of resolution. The test runs each case at two resolutions and computes the FWHM of the A-A′ cross-section of χ (tensile) or κ (shear), then asserts:
- Unregularized: FWHM_fine / FWHM_coarse ≈ h_fine / h_coarse (width shrinks with h).
- Regularized: |FWHM_fine − FWHM_coarse| / FWHM_coarse < 0.3 (width insensitive to h).

**Alternatives considered.** Storing a reference HDF5 and doing field-level L2 — rejected because it pins an arbitrary solution rather than the scientifically meaningful mesh-insensitivity property. FWHM detection is robust to ~5 % noise and cheap to implement.

### D5. Crust test: shear-band dip angle via RANSAC-style row-max
For Fig. 7, the scientific signature is the shear-band orientation matching Arthur's formula αS = π/4 ± (φ+ψ)/4 (paper Eq. 69). We detect bands by finding the per-row argmax of |ε̇_II| on the scaled strain-rate field, fitting a line through the points via least-squares, and asserting the fitted slope angle is within ±5° of Arthur's prediction. The `±` branch is selected based on the background-strain-rate sign (extensional → +, compressional → −).

**Alternatives considered.** FFT-based dominant-wavevector detection — rejected because at 41×15 resolution the band is only 2-3 cells wide, too narrow for a clean FFT peak. Row-max is simple and robust.

### D6. Tensile-fracture test: monotone upward propagation
Fluid Darcy is out of scope (see Non-Goals), so the test simplifies the setup by placing a weak inclusion (phase with reduced cohesion) at the bottom-centre of the domain, applying far-field extension, and asserting that the topmost z-coordinate at which χ exceeds a threshold grows monotonically over the test steps. Concretely:
```
z_tip(t_i) = max{ z : χ(x,z,t_i) > 0.1·max(χ(t_i)) }
```
Assertion: `z_tip(t_i+1) ≥ z_tip(t_i)` for all i, and final `z_tip > 0.3 · Lz`.

**Alternatives considered.** Wiring fluid pressure as a prescribed perturbation — deferred to a follow-up change because MDOODZ has no fluid-pressure field by default and this scope creep would triple the change size.

### D7. BDT test: depth-partitioned χ/κ coexistence check (experimental)
Check that the accumulated χ (tensile) is concentrated in the shallow half of the domain (near the free surface, low confining pressure) and κ (shear) dominates the deep half. Concretely, compute `mean(χ, top_half) / mean(χ, bottom_half) > 2` and `mean(κ, bottom_half) / mean(κ, top_half) > 1.5`. This is a qualitative fingerprint — a "canary" that the mode-I / mode-II partitioning is alive.

### D8. Regenerate `SETS/Popov2025_*_VEVP.{c,txt}` to match paper Table 1
The two shipped files are rewritten in place (names preserved, commit credit preserved via comment) to:
1. Set `plast = 2` in both phases.
2. Remove the `tensile = 0/1` key entirely (not read by any parser since commit 8aecb923).
3. Align `C`, `phi`, `psi`, `T_st`, `eta_vp`, `G`, `K` (via `bet = 1/K`), domain, and BCs to the Table 1 "Regularization" column.
4. Add `tensile_line_search = 1` to the global block.
5. Add `writer_subfolder = <SetupName>` (required for test-framework output discovery).

The .c side already defines correct `SetBCVx`/`SetBCVz` for the volumetric-plus-deviatoric BCs, so no changes there beyond whitespace.

**Alternatives considered.** Leaving the demo files as-is and writing brand-new ones — rejected because the current ones are silently wrong (they claim to be "Popov2025" but test the legacy code path) and will confuse future readers / users if not corrected.

### D9. CI-side χ/κ accumulation without touching MDLIB
HDF5 emits the rates `/Centers/eII_pl` and `/Centers/divu_pl` per step but not the integrals κ, χ. The test binary reads `writer_step = 1` output and accumulates in the test fixture: `κ_n = κ_{n-1} + eII_pl · dt`, same for χ. This avoids MDLIB changes (per the Impact section of the proposal) at a modest cost of larger HDF5 output during the test run — acceptable because all paper tests are short (≤ 200 steps).

**Alternatives considered.** Adding κ, χ fields to MDLIB HDF5 output — rejected as scope creep for a test-only change; can be a follow-up if multiple scenarios end up needing it.

### D10. Paper-figure reproduction via gnuplot, not Julia Makie
Paper figures are rendered with gnuplot + a small HDF5 → ASCII extractor, following the pattern already used by `TESTS/BlankenBench/extract_fields.cpp`. This keeps the toolchain consistent with the rest of the HDF5-driven plot scripts in the repo, avoids a Julia runtime dependency, and removes the need for a separate Project.toml / Makie environment.

All plot assets live under `TESTS/Popov2025/plots/`:
- `extract_popov2025.cpp` → linked against HDF5, reads `/Centers/sII`, `/Centers/P`, `/Centers/eII_pl`, `/Centers/divu_pl`, `/Model/Time`, and emits one `.dat` file per step group per field. Registered as `add_executable(extract_popov2025 ...)` in `TESTS/CMakeLists.txt`, matching the `extract_blankenbach` pattern.
- `fig5_0d_stress.gp` → reads the three 0D .dat outputs, produces τII–t, P–t, and τII–P meridional panels.
- `fig6_regularization.gp` → χ / κ heatmaps + A–A′ cross-sections for 2 resolutions × 2 (regularized/unregularized) = 4 panels.
- `fig7_crust.gp` → strain-rate field with annotated dip angle.
- `fig8_tensile.gp` → χ field snapshots at 2-3 time slices.
- `fig9_bdt.gp` → κ and χ spatial distribution with stress-envelope profile.

Scripts write PNGs to `TESTS/Popov2025/plots/out/` (gitignored, produced on demand). They are standalone — NOT wired into the `VISUAL_TESTS/` pipeline and NOT invoked by CI. A `README.md` in the same folder explains the run order: build `-DTEST=ON`, run the corresponding GTest, then `./extract_popov2025 <TestName>` followed by `gnuplot figN_<name>.gp`.

## Risks / Trade-offs

- **Risk: reduced-resolution regularization test fails the width-insensitivity check because 21×15 is too coarse for the regularized localization to be resolved at all.** → Mitigation: calibrate `eta_vp` up (e.g. 5×10¹⁹ instead of paper's 10¹⁹) so the length scale is ≈ 4-5 cells at the fine resolution. Document in `AnalyticalSolutions.md` that `eta_vp` is not a physical parameter and CI values differ from paper values.
- **Risk: Arthur's formula (D5) gives two symmetric conjugate angles; the row-max fit picks one and may be perpendicular to the other depending on the seed perturbation.** → Mitigation: assert `angle_err ∈ {|fit − αS⁺|, |fit − αS⁻|}` — accept either conjugate direction.
- **Risk: compressible=1 + plast=2 + non-zero gz triggers CHOLMOD failure** (known failure mode documented in skill-testing-guide: "Gravity + compressibility + dilatancy"). → Mitigation: all fast CI cases set `gz = 0` and drive deformation kinematically. BDT is the only setup that realistically needs gravity; it's already in the experimental tier where longer wall time lets us use a stabilizing penalty if needed.
- **Risk: writer_step = 1 blows up HDF5 disk usage during CI.** → Mitigation: regularization runs are ≤ 200 steps at ≤ 41×29 grid ≈ 200 × ~30 KB = 6 MB per test — negligible.
- **Risk: analytical 0D ODE helper (D3) drifts out of sync with RheologyDensity.c if the implementation ever refactors its yield-surface scaling (`a`, `b` in Eq. 14).** → Mitigation: the helper is < 100 LoC, references paper equation numbers in comments, and its own unit test verifies it against two hard-coded Popov Fig. 5 checkpoints. Any future refactor of RheologyDensity.c that breaks it signals a real regression.
- **Trade-off: we simplify Fig. 8 (tensile fracture) to a weak-inclusion proxy instead of fluid Darcy.** This loses the "fluid-pressure-driven" framing of the paper but preserves the core numerical test — upward propagation of a mode-I failure zone in a background pressure field — which is what exercises the combined yield surface numerically.
- **Trade-off: no attempt to reproduce the paper's FE convergence plot (Fig. 10).** MDOODZ uses finite differences; iteration counts are not directly comparable. Documented as a Non-Goal.

## Migration Plan

1. **Phase 1 (low-risk)**: Create `TESTS/Popov2025/` fixtures and `Popov2025Tests.cpp` with the five TEST_F classes. Add to `TESTS/CMakeLists.txt`. Update `TESTS/README.md` and `TESTS/AnalyticalSolutions.md`. Run locally; confirm all default-labeled tests pass under ≤ 2 min total.
2. **Phase 2 (touches shipped SETS)**: Rewrite `SETS/Popov2025_*_VEVP.{c,txt}` with `plast = 2` and Table 1 parameters. Locally verify they still build and run to completion.
3. **Phase 3 (optional)**: Add Julia scripts in `JuliaVisualisation/_Popov2025/` and the README. These are developer aids, not required by CI.

Rollback: all changes are additive except Phase 2. Phase 2 is rolled back by `git revert` on the affected `SETS/` files — they are demo setups, not part of the user-facing API.

## Reproduced paper figure

The single gnuplot script `TESTS/Popov2025/plots/fig5_0d_stress.gp` renders paper Fig. 5 from the CI test HDF5 output. Parameters match paper Table 1 "0D Fig. 5 a, b" exactly; the integration-point pressure `p_local = p* + K·θ̇_vp·dt` is reconstructed (paper Eq. 31) so saturation values sit on the yield surface.

![Fig. 5 — 0D stress integration](../../../TESTS/Popov2025/plots/reference/fig5_0d_stress.png)

- **(a) Volumetric extension**: `p → p_T = −0.5 MPa` at ≈ 14 yr, `τII ≈ 0` throughout.
- **(b) Deviatoric shear**: `τII` yields at the DP envelope (≈ 1 MPa) at ≈ 20 yr; `p` grows post-yield via ψ = 10° dilation.
- **(c) Mixed strain**: `p` dips to ≈ −0.35 MPa, `τII` grows through yield — trajectory sweeps through the cap ↔ DP delimiter.
- **(d) Meridional trajectories**: extension slides onto the cap tip, shear climbs up τII then bends right, mixed forms the characteristic S-curve.

2D tests (paper Figs. 6–9) were prototyped during development and deferred: the localisation signals require resolution well above the fast-CI budget. See `TESTS/AnalyticalSolutions.md` §6.4 for the deferred-test list.

## Open Questions

- Should the regularization test assert an absolute χ_max value (sensitive to Newton tolerance) or only the FWHM ratio across resolutions? Recommendation: FWHM ratio only, since Newton tolerance is already validated by the 0D integration test.
