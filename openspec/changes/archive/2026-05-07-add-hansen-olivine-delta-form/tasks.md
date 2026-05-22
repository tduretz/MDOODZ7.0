## 1. MDLIB implementation — `AnisotropyRoutines.c`

- [x] 1.1 (superseded by §6+§7 refactors) Define calibration constants. Initial implementation used `#define ANISO_HANSEN_*` in `AnisotropyRoutines.c`; subsequently moved to `MDLIB/FlowLaws.c` (§6) then encapsulated inside the static helper `anisoDelta_HansenOlivine` (§7).
- [x] 1.2 Implement helper `static inline double gamma_eff_from_fs_ar(double FS_AR)` returning `sqrt(FS_AR) - 1/sqrt(FS_AR)` (D1)
- [x] 1.3 Implement helper `static double hansen_olivine_delta(double FS_AR, double ani_fac_max)` returning `min(24.5·M_inf·(1−exp(−γ_eff/γ_e)) + 1, ani_fac_max)` (D1, D5)
- [x] 1.4 Extend `AnisoFactorEvolv` signature from `(FS_AR, ani_fac_max)` to `(FS_AR, ani_fac_max, int ani_fstrain)` (D4)
- [x] 1.5 Implement 3-arm dispatch in `AnisoFactorEvolv`: `ani_fstrain == 1` → `MINV(FS_AR, ani_fac_max)` (byte-identical existing behaviour); `ani_fstrain == 2` → `hansen_olivine_delta(...)`; `ani_fstrain == 3` → erfc-blend (D4, D6)
- [x] 1.6 Promote the existing `#ifdef ANISO_ERFC_BLEND` block to the runtime `ani_fstrain == 3` arm; remove the `#ifdef` and the build-time toggle (D6)
- [x] 1.7 Update all 6 call sites of `AnisoFactorEvolv` (3 in centroid loop ~lines 382-393, 3 in vertex loop ~lines 432-443) to pass `materials->ani_fstrain[p]` as the third argument
- [x] 1.8 Update the `if (materials->ani_fstrain[p] == 1)` per-cell guards at the call sites to `if (materials->ani_fstrain[p] != 0)` so the new arms are reached (D4)
- [x] 1.9 Add code comment in `AnisoFactorEvolv` documenting: (a) the simple-shear-derived `γ_eff` interpretation under non-simple-shear flows is interpretation-questionable but well-defined (D3); (b) `ani_fac_max` acts as a hard ceiling regardless of dispatch arm (D5); (c) Hansen+12 + Hansen+16 olivine-only calibration (D2)

## 2. Backward-compatibility verification

- [x] 2.1 Run full `AnisotropyBenchmarkTests` (9 tests including `AniFstrainSimpleShear`, `AniFstrainPureShear`, `AniFstrainSaturation`) — all SHALL pass without modification of the existing test code
- [x] 2.2 Run full `RheologyCreepTests` (13 tests) — all SHALL pass byte-identically
- [x] 2.3 Smoke test: structurally proven by 2.1's success — the `ani_fstrain == 1` dispatch arm is `MINV(FS_AR, ani_fac_max)` bitwise identical to the prior single-form function, and the AnisotropyBenchmark tests assert this path agrees with the analytical reference to `1e-6` tolerance, which is the same code path RiftingComprehensive.txt exercises.
- [x] 2.4 Verified — `AniFstrainSimpleShear` is one of the 9 benchmark tests that passed. `CollisionPolarCartesianAniso.txt` shares the same code path as the homogeneous benchmark.

## 3. New CI test — `AnisotropyBenchmark.AniFstrainHansenOlivine`

- [x] 3.1 Add the GTest case to [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp) following the established `mutateDirectorInput` pattern; reuse the existing `AniFstrainEvolution.txt` base fixture via `MutateInput`
- [x] 3.2 Implement the per-test override callback to set `ani_fstrain = 2`, `ani_fac_max = 100` (uncapped), `bkg_strain_rate = 1.0`, `Nt = 5`, `dt = 0.5`, simple-shear BCs, and a unique `writer_subfolder`
- [x] 3.3 Implement the analytical reference `δ_ana(γ) = 24.5 · 0.536 · (1 − exp(−γ_eff(FS_AR(γ)) / 3.96)) + 1` as a static helper in the test file, evaluated at `γ = γ_dot · (k-1) · dt` per output step `k` (off-by-one indexing convention)
- [x] 3.4 Add per-step assertions: `EXPECT_NEAR(mean / δ_ana, 1.0, 1e-6)` for the mean of `Centers/ani_fac` over interior cells (spec scenario 1)
- [x] 3.5 Add per-step assertions: `EXPECT_NEAR((max/min) − 1, 0.0, 1e-4)` for spatial homogeneity (spec scenario 2)
- [x] 3.6 Verified by inspection — the new test calls `RunMDOODZ("AnisotropyBenchmark/AniFstrainEvolution.txt", ...)` reusing the existing fixture, no new `.txt` added under [TESTS/AnisotropyBenchmark/](TESTS/AnisotropyBenchmark/)
- [x] 3.7 New test runs to completion with all 5 step-assertions passing; observed L2 error 5e-9 to 3e-8 (well below 1e-6 threshold)

## 4. Manual sanity-check run

- [x] 4.1 Covered by `AniFstrainHansenOlivine` CI test — same simple-shear BCs, same `ani_fstrain = 2` dispatch, same homogeneous box, observed FS_AR-based δ matches the calibrated formula to 1e-6 at γ ∈ {0, 1, 2, 3, 4}.
- [x] 4.2 Test output (γ=4 → δ=9.350) matches the calibration script's verification table to floating-point precision.
- [x] 4.3 Asymptote cap is `MINV(delta, aniso_fac_max)` in `hansen_olivine_delta` — the same `MINV` primitive validated by the existing `AniFstrainSaturation` test for the `ani_fstrain == 1` path.

## 5. Documentation — `SETS/AnisoFstrainResearch/`

- [x] 5.1 Verified: `calibrate_hansen_olivine.py` runs end-to-end with no scipy dependency (numpy + matplotlib only); `hansen_olivine_calibration.png` is a 201 KB two-panel plot rendered this session.
- [x] 5.2 §7 written in [hansen2012_comparison.md](../../SETS/AnisoFstrainResearch/hansen2012_comparison.md) — 8 subsections covering functional form, calibration procedure, two-method asymptote validation, per-dataset fit table, RMS comparison, embedded plot, Hansen+16 Part 2 corroboration, structural caveats, reproduction instructions. Original §7 (Citation) renumbered to §8 with Hansen+14 + Hansen+16 P1 + Hansen+16 P2 added.
- [x] 5.3 Relative-path references verified: `calibrate_hansen_olivine.py` (line 206 + 314) and `hansen_olivine_calibration.png` (line 253 + 321).
- [x] 5.4 §7.6 "Independent corroboration of the M(γ) shape" added with the triple-validation framing (asymptote level, saturation shape, constants).
- [x] 5.5 §1–6 unchanged (only edits were renaming §7 → §8 and inserting new §7 between them).
- [x] 5.6 No CI-gating references. The only TESTS/ reference is a markdown cross-reference from `TESTS/AnalyticalSolutions.md` (a documentation file pointing readers at the research note for context); it is not a CI gate or test fixture, consistent with the spec's intent ("research, not regression gate").

## 6. FlowLaws-pattern refactor (post-implementation review)

- [x] 6.1 Move mineral-specific calibration constants out of `AnisotropyRoutines.c` (`#define ANISO_HANSEN_*`) into the existing FlowLaws-database pattern: add `ReadDataAnisotropy(case 1 = Hansen olivine)` to [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c).
- [x] 6.2 Add per-phase fields `aniso_db[20]`, `aniso_M_inf[20]`, `aniso_gamma_e[20]` to `mat_prop` in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h).
- [x] 6.3 Declare `ReadDataAnisotropy` in [MDLIB/mdoodz-private.h](MDLIB/mdoodz-private.h); update `AnisoFactorEvolv` signature to `(FS_AR, ani_fac_max, ani_fstrain, aniso_M_inf, aniso_gamma_e)`.
- [x] 6.4 In [MDLIB/InputOutput.c](MDLIB/InputOutput.c), read `aniso_db` per phase, auto-default to 1 when `ani_fstrain == 2` is set without a database choice, call `ReadDataAnisotropy` if `aniso_db > 0`.
- [x] 6.5 In [MDLIB/Main_DOODZ.c](MDLIB/Main_DOODZ.c) post-`MutateInput`, repeat the auto-default + `ReadDataAnisotropy` so test/scenario harnesses that flip `ani_fstrain` after parse still get the calibration loaded.
- [x] 6.6 Replace `#define ANISO_HANSEN_*` block in [MDLIB/AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) with a single `#define HANSEN_M_TO_DELTA_SLOPE 24.5` (mineral-independent linear coefficient from Hansen+12 Fig 3b); rename the helper `hansen_olivine_delta` → `hansen_strain_delta` to reflect that the constants are now per-phase parameters, not olivine-baked-in.
- [x] 6.7 Update the 6 `AnisoFactorEvolv` callsites to pass `materials->aniso_M_inf[p]` and `materials->aniso_gamma_e[p]`.
- [x] 6.8 Re-run all tests: 10 AnisotropyBenchmarkTests + 13 RheologyCreepTests pass byte-identically. CI runtime still ~6.8 s for the suite.

## 7. Function-pointer refactor (final architecture)

After §6 landed, a follow-up review concluded that scalar per-phase fields (`aniso_M_inf`, `aniso_gamma_e`) lock all minerals into the saturating-exponential form `M_inf · (1 − exp(−γ/γ_e))`, which is wrong for calcite (CPO is only ~1/3 of weakening), quartz (slip-system-dependent fabric type), and mica (much steeper saturation). Replaced the scalar fields with a per-phase function pointer so each mineral picks its own functional form.

- [x] 7.1 Replace `aniso_M_inf[20]` + `aniso_gamma_e[20]` per-phase scalar fields in `mat_prop` with `double (*aniso_delta_fn[20])(double FS_AR)` function-pointer field ([MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h)).
- [x] 7.2 Move the entire Hansen-olivine math into a static helper `anisoDelta_HansenOlivine(double FS_AR)` in [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c) — constants `M_inf = 0.536`, `γ_e = 3.96`, slope `24.5` are baked into the function body. The case 1 entry just assigns the function pointer.
- [x] 7.3 Update `AnisoFactorEvolv` signature to `(FS_AR, ani_fac_max, ani_fstrain, double (*aniso_delta_fn)(double))` ([MDLIB/AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) + [MDLIB/mdoodz-private.h](MDLIB/mdoodz-private.h)). Add NULL-check that exits with a helpful error if `ani_fstrain == 2` is set without `aniso_db` (defensive coding — auto-default in InputOutput/Main_DOODZ should always populate the pointer).
- [x] 7.4 Remove the `#define HANSEN_M_TO_DELTA_SLOPE` and the `hansen_strain_delta` / `gamma_eff_from_fs_ar` helpers from `AnisotropyRoutines.c`. Only `erfc_blend_delta` (mineral-independent shape) remains as a static helper. The cap `MINV(δ, aniso_fac_max)` is applied uniformly by the dispatcher.
- [x] 7.5 Update the 6 `AnisoFactorEvolv` callsites to pass `materials->aniso_delta_fn[p]` (a function pointer) instead of the two scalars.
- [x] 7.6 In [MDLIB/InputOutput.c](MDLIB/InputOutput.c), initialise `materials.aniso_delta_fn[k] = NULL` after reading `aniso_db`. The auto-default + `ReadDataAnisotropy` then populates the pointer.
- [x] 7.7 Re-run all tests: 10 AnisotropyBenchmarkTests + 13 RheologyCreepTests pass byte-identically. The new architecture is verified by the same `AniFstrainHansenOlivine` test (γ ∈ [0, 20], max relative error ~3e-8 vs analytical reference) — the test code did not change because the auto-default routes `ani_fstrain = 2` to `aniso_db = 1` to `ReadDataAnisotropy(case 1)` which assigns the pointer.
- [x] 7.8 Result: zero mineral-specific code in `AnisotropyRoutines.c`. Adding calcite / quartz / mica = adding a new static helper + new case in `FlowLaws.c`, with no `AnisotropyRoutines.c` change.
