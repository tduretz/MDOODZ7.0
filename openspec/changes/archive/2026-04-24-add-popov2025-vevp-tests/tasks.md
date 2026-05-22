## 1. Test fixture skeleton & CMake wiring

- [x] 1.1 Create directory `TESTS/Popov2025/` with `plots/` and `plots/out/` (the latter gitignored)
- [x] 1.2 Create `TESTS/Popov2025Tests.cpp` with `Popov2025_0DIntegration` fixture and the three 0D `TEST_F` cases (VolumetricExtension, DeviatoricShear, MixedStrain)
- [x] 1.3 Add `Popov2025Tests` target to `TESTS/CMakeLists.txt` with `LABELS "popov2025"`; registered via `add_test` and linked against `mdoodz` + `GTest::gtest_main`
- [x] 1.4 Verify build green with `cmake -B cmake-build -DTEST=ON && cmake --build cmake-build --target Popov2025Tests`

## 2. 0D reference integrator (Popov §3.6)

- [x] 2.1 Create `TESTS/Popov2025/Popov0DAnalytical.{cpp,h}` — an implicit backward-Euler scalar integrator that mirrors the paper's local Newton solve, independent of MDLIB; exposes `tauII_trajectory(p)`, `P_trajectory(p)`, and a `geometry(p)` helper for the yield-surface constants
- [x] 2.2 Wrap with a `-DPOPOV0D_STANDALONE` self-check against two hard-coded Fig. 5 checkpoints (volumetric-extension steady state and deviatoric-shear steady state)
- [x] 2.3 Link `Popov0DAnalytical.cpp` into the `Popov2025Tests` target (not into `mdoodz`)

## 3. 0D stress-integration tests (paper Fig. 5)

- [x] 3.1 Author `TESTS/Popov2025/Popov0D_VolumetricExtension.txt` — 21×15 grid, Table 1 "0D" parameters (G=1e10, K=2e11, φ=30°, ψ=10°, C=1e6, T_st=-5e5, η^vp=0, Δt=2 yr), outward normal BCs `user2 = user3 = 1.167e-15` giving `ε̇_xx = ε̇_zz = 2.333e-15`
- [x] 3.2 Author `TESTS/Popov2025/Popov0D_DeviatoricShear.txt` — standard pure-shear BCs with `bkg_strain_rate = 7e-14` (paper value), `tensile_line_search = 0` to avoid Newton divergence on large per-step elastic overshoot
- [x] 3.3 Author `TESTS/Popov2025/Popov0D_MixedStrain.txt` — asymmetric outward normal BCs (`user2 = -3.325e-14`, `user3 = +2.573e-14`) that produce `trace(ε̇) = 7e-15` and `ε̇_II_dev = 7e-14` simultaneously, matching paper Fig. 5c
- [x] 3.4 Implement `TEST_F(Popov2025_0DIntegration, VolumetricExtension)` — reconstruct `sII` from `sxxd`/`szzd`/`sxz`, assert `divu_pl > 0` (dilatant activation), `|sII| < 15% τ_d`, and `|R̂_y − R_y|/R_y < 30%` (tensile cap membership)
- [x] 3.5 Implement `TEST_F(Popov2025_0DIntegration, DeviatoricShear)` — assert `eII_pl > 0`, `|sII − (k·P + c)| / (k·P + c) < 5%` (DP envelope membership), `|P| < 10·c_MC` (bounded dilatant growth)
- [x] 3.6 Implement `TEST_F(Popov2025_0DIntegration, MixedStrain)` — assert at least one plastic mode active, and the smaller of the cap-circle and DP-envelope residuals is < 30% of τ_d

## 4. Rewrite existing SETS/Popov2025 demo scenarios

- [x] 4.1 Rewrite `SETS/Popov2025_Tensile_VEVP.txt`: set both phases to `plast = 2`, remove obsolete `tensile = 0/1` keys, add `tensile_line_search = 1`, align (φ, ψ, C, Ce, T_st, η^vp, G, K) to paper Table 1 "Regularization" column, add top-of-file header comment referencing the paper and commit 8aecb92
- [x] 4.2 Rewrite `SETS/Popov2025_Pureshear_VEVP.txt` with the same edits for the shear case
- [x] 4.3 Add one-line header comments to `SETS/Popov2025_{Tensile,Pureshear}_VEVP.c` pointing at the paper section
- [x] 4.4 Registering the two scenarios in `SETS/CMakeLists.txt` is deferred (pre-existing gap from commit 8aecb92; out of scope for this change)

## 5. Gnuplot pipeline for paper Fig. 5

- [x] 5.1 Create `TESTS/Popov2025/plots/extract_popov2025.cpp` — HDF5 → ASCII extractor that reads `/Centers/sxxd`, `/Centers/szzd`, `/Vertices/sxz` to reconstruct `sII`, plus `/Centers/{P, eII_pl, divu_pl}`; emits one `trajectory.dat` per test output directory
- [x] 5.2 Register `extract_popov2025` as `add_executable` in `TESTS/CMakeLists.txt`, linked against `${HDF5_C_LIBRARIES}`
- [x] 5.3 Author `TESTS/Popov2025/plots/fig5_0d_stress.gp` — four-panel layout (a) P(t)/τII(t) for volumetric extension, (b) for deviatoric shear, (c) for mixed strain, (d) meridional P–τII trajectories for all three. Pressure is reconstructed as `p_local = p* + K·θ̇_vp·dt` (paper Eq. 31) so saturation values sit on the yield surface
- [x] 5.4 Author `TESTS/Popov2025/plots/README.md` describing the invocation order and overrides

## 6. Documentation

- [x] 6.1 Add "Suite 20 — Popov2025 (Combined Tensile Cap / DP Yield)" to `TESTS/README.md` with the three-test table, parameters, how-to-run, and a pointer to the deferred 2D tier
- [x] 6.2 Add "§6 Popov et al. (2025) Combined Yield Surface" to `TESTS/AnalyticalSolutions.md` covering: (a) yield-surface geometry (paper Eqs. 13–17); (b) embedded reference PNG of Fig. 5; (c) 0D reference integrator description; (d) test-case-to-assertion mapping; (e) deferred 2D tests list
- [x] 6.3 Commit a fresh `TESTS/Popov2025/plots/reference/fig5_0d_stress.png` rendered from the final calibrated run, referenced from both docs

## 7. Validation

- [x] 7.1 Run `make run-tests` / direct binary — all three 0D tests pass in ~2 s total
- [x] 7.2 Render `fig5_0d_stress.png` and visually confirm all four panels match paper Fig. 5 qualitatively (P→T_st at ~14 yr for volumetric; τII yields at ~20 yr for shear; characteristic S-curve for mixed; three trajectories in (d))
- [x] 7.3 `openspec validate add-popov2025-vevp-tests` passes

## 8. Deferred — 2D localisation tests (paper Figs. 6–9)

These were prototyped during development and deferred to a follow-up change — the CI-budget resolution doesn't resolve the diagnostic signals. The follow-up change should:

- [ ] 8.1 Produce the Fig. 6 mesh-insensitivity test (FWHM ratio across 2+ resolutions) at a higher-resolution tier behind `LABELS experimental`
- [ ] 8.2 Reproduce the Fig. 7 crustal-scale shear bands with Arthur's-angle detection
- [ ] 8.3 Reproduce the Fig. 8 tensile-fracture upward propagation (either via fluid-pressure perturbation or a weak-inclusion proxy)
- [ ] 8.4 Reproduce the Fig. 9 brittle-ductile depth partitioning with dislocation-creep rheology
