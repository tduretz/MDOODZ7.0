## 1. Base fixture file

- [x] 1.1 Create [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) — single base scenario with `Nb_phases = 1`, `anisotropy = 1`, per-phase `ani_fstrain = 1`, `elastic = 0`, `thermal = 0`, `Nx = Nz = 21`, default pure-shear settings (`shear_style = 0`), `pure_shear_ALE = 1`, `periodic_x = 1`, `reseed_markers = 0`, `ani_fac_max = 1e6`, `Nt = 5`, `dt = 0.5`. Tests override the BC mode to simple shear via MutateInput.
- [x] 1.2 Verify the bare fixture parses, executes step 1, and writes a valid `Output*.gzip.h5` with `Centers/ani_fac` populated.

## 2. MutateInput callback

- [x] 2.1 Add `AniFstrainOverride` struct + `g_aniFstrainOverride` static + `mutateAniFstrainInput(MdoodzInput *input)` callback at file scope in [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp), supporting overrides for `bkg_strain_rate`, `Nt`, `dt` (with `dt`, `dt0`, `dt_start` set together), `ani_fac_max` (per-phase, phase 0), `shear_style`, `periodic_x`, `pure_shear_ALE` (sentinel `-99`), and `writer_subfolder` (heap-safe `free`+`strdup`).
- [x] 2.2 Decided: keep file-scope (override fields specific to this test family).

## 3. Test case: simple shear, unsaturated regime

- [x] 3.1 Add `AnisotropyBenchmark.AniFstrainSimpleShear` with `bkg_strain_rate = 1.0` (γ̇ = 2), `ani_fac_max = 1e6`, `Nt = 5`, `dt = 0.5`, `shear_style = 1`, `periodic_x = 1`.
- [x] 3.2 At each output step `k ∈ {1..Nt}`, read `Centers/ani_fac` from `Output<k>.gzip.h5`, compute mean over interior cells, assert it matches `1 + γ²/2 + γ·sqrt(γ²/4+1)` evaluated at `γ = γ̇·(k-1)·dt` (output offset by one step) with relative error `< 1e-6`.
- [x] 3.3 Assert spatial L2 error `< 1e-4`.
- [x] 3.4 Assert spatial homogeneity `(max/min - 1) < 1e-4`.
- [x] 3.5 Assert `min(ani_fac) >= 0.999` (FS_AR ≥ 1, with float32 HDF5 round-trip slack) and finite.

## 4. Test case: simple shear, saturation clamp

- [x] 4.1 Add `AnisotropyBenchmark.AniFstrainSaturation` with `ani_fac_max = 4`, `bkg_strain_rate = 1.0` (γ̇ = 2), `Nt = 5`, `dt = 0.5` — γ ∈ {0, 1, 2, 3, 4}, γ_sat ≈ 1.94.
- [x] 4.2 Assert steps 1–2 (γ ∈ {0, 1}, unsaturated): mean matches `FS_AR(γ)` with relative error `< 1e-6`.
- [x] 4.3 Assert steps 3–5 (γ ∈ {2, 3, 4}, saturated): mean equals `ani_fac_max = 4` with absolute error `< 1e-9`.
- [x] 4.4 Assert spatial homogeneity at every step.

## 5. Pure-shear test (DEFERRED — MDOODZ bug)

- [ ] 5.1 ~~Add `AnisotropyBenchmark.AniFstrainPureShear` test~~ — **DROPPED FROM CI**. MDOODZ's `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` segfaults at the start of step 2 (silent SIGSEGV after Output00001 writes cleanly) under all `pure_shear_ALE` settings (0, 1, -1) and multiple strain rates. Root cause is in MDOODZ's pure-shear time-loop machinery, not in this test setup. Documented in [design.md](design.md) "Open Issues". Pure-shear analytical formula `FS_AR(t) = exp(2·ε̇·t)` is documented in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md). Investigation deferred to a follow-up change.

## 6. CI visualisation

- [x] 6.1 The saturation test writes `aniso_factor_evolution.dat` with columns `step, t, gamma, FS_AR_analytical, delta_clamped_analytical, delta_mdoodz_mean` (saturation test owns this — its data exhibits the `min()` clamp visually).
- [x] 6.2 Added `TESTS/AnisotropyBenchmark/plot_aniso_factor.gp` rendering `aniso_factor_evolution.png` with three curves: unbounded `FS_AR(γ)` (solid), clamped `δ = min(FS_AR, 4)` (dashed), MDOODZ markers (red triangles), plus reference lines at `y = ani_fac_max` and `x = γ_sat`.
- [x] 6.3 Ran gnuplot; committed `TESTS/aniso_factor_evolution.png`.

## 7. Documentation

- [x] 7.1 Added a new numbered §9 "Finite-Strain Anisotropy" to [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md): closed-form derivations for `FS_AR(t)` (pure shear) and `FS_AR(γ)` (simple shear) from polar decomposition `F = R·U`, `min(FS_AR, ani_fac_max)` saturation form with code reference, parameter values, expected L2 magnitudes, GTest assertion blocks, off-by-one output-indexing note, embedded `aniso_factor_evolution.png`, §9.7 footnote about the pure-shear segfault, §9.8 pointer to the SETS research note.
- [x] 7.2 Added 3 rows to the summary table at the bottom of `TESTS/AnalyticalSolutions.md` (simple-shear, saturation, spatial homogeneity).

## 8. Research artefact (NOT CI)

- [x] 8.1 Created directory `SETS/AnisoFstrainResearch/`.
- [x] 8.2 Added `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` — runnable simple-shear scenario with `bkg_strain_rate = 1.0`, `dt = 0.1`, `Nt = 20`, `ani_fac_max = 14` (Hansen+12 asymptote), reaches γ ≈ 4.
- [ ] 8.3 ~~Verify it executes via cmake-exec/AniFstrainSimpleShear~~ — MDOODZ does not auto-build a binary for this `.txt` (no `.c` setup file). Manual run instructions documented in `hansen2012_comparison.md` §6.
- [x] 8.4 Wrote `SETS/AnisoFstrainResearch/hansen2012_comparison.md` (~200 lines): live MDOODZ formula with code refs, Hansen+12 fit `δ(γ) = 13·tanh(0.25γ) + 1`, quantitative comparison table at γ ∈ {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 8}, form-mismatch discussion, three future-work options (erfc blend / direct tanh / keep min() and document), full citation Hansen+12 *Nature* 492, 415–418, doi:10.1038/nature11671.
- [ ] 8.5 ~~Optional PNG overlay~~ — table-only comparison sufficient; PNG can be added later if useful.

## 9. Verification

- [x] 9.1 Ran `AnisotropyBenchmarkTests` (8/8 PASS, 5.3 s wall — DirectorEvolution, DirectorDtConvergence, StressAnisotropy, StressAnisotropyL2, DirectorEvolutionInterpMode1/2, **AniFstrainSimpleShear**, **AniFstrainSaturation**) and `RheologyCreepTests` (12/12 PASS — all 11 RheologyCreep cases plus PinchSwellGSESmoke). No regressions.
- [x] 9.2 L2 thresholds calibrated based on observed FP/HDF5-float32 residuals: relative `1e-6`, absolute saturated `1e-9`, spatial L2 `1e-4`. Observed residuals are `~3e-8` to `~5e-8` per step (well inside the threshold) — float32 round-trip noise, not method error. Documented in [design.md D5](design.md).
- [x] 9.3 `aniso_factor_evolution.png` renders cleanly (clamp knee visible at γ_sat ≈ 1.94, MDOODZ markers sit on dashed clamp line, ani_fac_max=4 reference line visible). Embedded in `AnalyticalSolutions.md` §9.6.
- [x] 9.4 `SETS/AnisoFstrainResearch/` exists with `AniFstrainSimpleShear.txt` and `hansen2012_comparison.md`. NOT referenced from any `TESTS/` or CI build script — research-only by design.
