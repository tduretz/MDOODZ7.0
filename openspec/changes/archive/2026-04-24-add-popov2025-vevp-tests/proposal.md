## Why

The combined mode-I / mode-II yield surface of Popov et al. (2025, GMD, doi:10.5194/gmd-18-7035-2025) was merged into MDLIB as `plast = 2` in commit 8aecb923 (PR #159), but the merge ships with two mislabeled `SETS/Popov2025_*_VEVP.*` scenarios that still use the legacy `plast = 0 / 1` paths, and no CI coverage at all for the combined yield surface. As a result, a non-trivial new constitutive law (circular tensile cap + Drucker-Prager shear envelope, Perzyna regularization, 3×3 Newton local solve, Armijo line search) can regress silently. The Popov paper provides a self-contained test hierarchy — 0D stress integration, 2D mesh-independent localization, crust-scale shear bands, tensile dyke propagation, and brittle-ductile transition — that is ready to be ported as analytical / regression benchmarks.

## What Changes

- **Add three CI tests** as `Popov2025_0DIntegration.{VolumetricExtension,DeviatoricShear,MixedStrain}` in `TESTS/Popov2025Tests.cpp`, reproducing paper Fig. 5 a/b/c. Each fixture runs a homogeneous 21×15 grid with paper Table 1 "0D" parameters (η^vp = 0, perfect plasticity) and asserts yield-surface membership at the final step — cap circle for volumetric extension, Drucker-Prager envelope for shear, either-branch for mixed. The tests exercise the full `plast = 2` code path: segment branching condition, 3×3 local Newton, and cap/DP flow potentials.
- **Ship a test-only reference integrator** `TESTS/Popov2025/Popov0DAnalytical.{cpp,h}` that implements the paper's §3.6 implicit stress update independently (no MDLIB link), so future regressions can be compared against an independent coding of the same algorithm. Includes a standalone `-DPOPOV0D_STANDALONE` self-check against two Fig. 5 checkpoints.
- **Fix the two shipped demo scenarios** `SETS/Popov2025_Pureshear_VEVP.{c,txt}` and `SETS/Popov2025_Tensile_VEVP.{c,txt}`: set `plast = 2` (they previously used the legacy `plast = 0 / 1` paths), add `tensile_line_search = 1`, remove the unused `tensile = 0/1` keys, and align parameters with Table 1 "Regularization" plus a top-of-file comment pointing at the paper section.
- **Register `Popov2025Tests` in `TESTS/CMakeLists.txt`** under `LABELS "popov2025"`; included in the default `make run-tests` target.
- **Document the suite** in `TESTS/README.md` (new "Suite 20 — Popov2025 (Combined Tensile Cap / DP Yield)") and the reference solutions in `TESTS/AnalyticalSolutions.md` (new "§6 Popov et al. (2025) Combined Yield Surface"), with an embedded reference PNG rendered from the CI output.
- **Reproduce paper Fig. 5** with a gnuplot script (`TESTS/Popov2025/plots/fig5_0d_stress.gp`) plus an HDF5 → ASCII extractor (`extract_popov2025.cpp`, modelled on `TESTS/BlankenBench/extract_fields.cpp`). The script reconstructs the integration-point pressure `p_local = p* + K·θ̇_vp·dt` (paper Eq. 31) from the HDF5 global `/Centers/P` and `/Centers/divu_pl` so saturation values sit on the yield surface. Developer aid only — not wired into CI.
- **Defer 2D localisation tests (paper Figs. 6–9)** to a follow-up change. They were prototyped during development but require resolution beyond the fast-CI budget to produce the diagnostic signals (FWHM ratios, Arthur's angle, monotone tip propagation, depth-partitioned χ/κ). The 0D tier gives regression coverage of the constitutive-law code path, which is the primary source of risk for the new yield surface.

## Capabilities

### New Capabilities
- `ci-popov2025-vevp-tests`: CI + regression coverage for the combined mode-I / mode-II (tensile cap + Drucker-Prager) yield surface of Popov et al. 2025 — the five test families above, their analytical / reference-based assertions, and documentation in `TESTS/README.md` and `TESTS/AnalyticalSolutions.md`.

### Modified Capabilities
None. `ci-plasticity-tests` continues to cover the legacy `plast = 1` Drucker-Prager path; the Popov 2025 combined formulation is a distinct yield surface (different constitutive law, different code path in `RheologyDensity.c`) and belongs in its own capability.

## Impact

- **Code**: no changes to `MDLIB/`. The change is test fixtures, CMake wiring, documentation, and Julia post-processing scripts. The two `SETS/Popov2025_*_VEVP.{c,txt}` files are rewritten (not deleted), which is non-breaking because they are demo scenarios, not user API.
- **CI**: adds ~1 GTest binary with 5 test groups (≈ 10-15 test cases). Fast cases target ≤ 30 s each at 21-31 grid; regularization and BDT multi-resolution runs are labeled `experimental` so the default CI pipeline stays under its current budget.
- **Dependencies**: none new. Tests use the existing `TESTS/TestHelpers.h` L2-error framework; Julia plotting uses the already-vendored Makie stack.
- **Docs**: `TESTS/README.md`, `TESTS/AnalyticalSolutions.md`, and a new `JuliaVisualisation/_Popov2025/README.md` describing the reproduction of paper figures.
- **Risk**: the five paper tests are published reference solutions and the combined yield surface is already merged, so the main risk is calibrating CI tolerances — addressed by the multi-resolution + L2 framework already used by SolVi and Anisotropy benchmarks.
