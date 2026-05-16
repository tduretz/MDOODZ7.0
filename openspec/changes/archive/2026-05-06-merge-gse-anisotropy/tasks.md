## 1. Pre-implementation read

- [x] 1.1 Read [`ViscosityConciseAniso`](MDLIB/AnisotropyRoutines.c#L82) end-to-end. Confirmed `*d` flows into the function only via `Eii_lin = C_lin * Tii^n_lin * (*d)^(-m_lin)` at line 144 (Newton iteration body) and line 232 (post-iteration reconstruction). Line 109's `eta_lin = B_lin * Eii^(1/n_lin - 1) * (*d)^(m_lin/n_lin)` reads `*d` for the upper-bound estimate. `d0` is read separately for B_pwl pre-factor (line 47) and the grain-size-zero sanity check (lines 55-56) — those are immune to `*d` changes. No other readers.
- [x] 1.2 Read [`LocalIterationViscoElasticGrainSize`](MDLIB/RheologyDensity.c#L203) signature and body. The `LocalIterationParams` struct fields ([RheologyDensity.h:17-50](MDLIB/RheologyDensity.h#L17)) needed to populate from `ViscosityConciseAniso`'s context: `Eii` (rotated-frame I2), `eta_up/eta_lo` (mechanism MIN/MAX), `Ag/gam/lam/cg/pg` (from materials->Kpzm/Gpzm/Lpzm/cpzm/ppzm via `Ag = Kg * exp(-Qg/RT)`), `C_pwl/n_pwl` (already at lines 128-129), `C_lin/n_lin/m_lin` (line 141), `C_exp/ST/n_exp` (line 154 / 149), `eta_cst/eta_el` (lines 185 / 187), mechanism flags (lines 115-122), `phase`, `noisy`, `d`, `T`. The struct's `f_ani` field is dead code — never read or written in the iteration; safe to omit.
- [x] 1.3 Confirmed rotation invariance of `Eii = sqrt(I2(E_rot))` by inspection: `LocalIterationViscoElasticGrainSize` uses `params.Eii` only as a scalar invariant in (a) `Tii = 2·eta_ve·Eii`, (b) the residual `r = Eii − elastic·Tii/(2·eta_el) − Eii_vis`, (c) the wattmeter formula `d_ve = (2·Bg·Eii·Eii_pwl·eta_ve·p/Ag)^(-1/(p+1))`. No tensor decomposition or angular dependency. So feeding the rotated-frame `sqrt(I2(E_rot))` is equivalent to the lab-frame `sqrt(I2(E))` — both evaluate to the same number to FP precision, and the iteration converges identically.

## 2. MDLIB fix

- [x] 2.1 Added GSE-active branch to `ViscosityConciseAniso` at the iteration block ([AnisotropyRoutines.c:218-280](MDLIB/AnisotropyRoutines.c#L218)). Branch gated on `materials->gs[phase] != 0` (matches the existing convention at [RheologyDensity.c:477](MDLIB/RheologyDensity.c#L477)) — covers all GSE laws dispatched by `ReadDataGSE` (calcite wattmeter at `gs=9,10`, piezometers at `gs=11,12`, olivine at `gs=40`, etc.). Branch calls `LocalIterationViscoElasticGrainSize` with populated `LocalIterationParams` and the `LocalIterationMutables` mirroring [RheologyDensity.c:644](MDLIB/RheologyDensity.c#L644). The inline iteration is kept as the `else` branch unchanged (byte-identical for `gs == 0`). Also added `eta_lo` calculation (MAX of mechanism viscosities) which the inline iteration didn't need but `LocalIterationViscoElasticGrainSize` requires for its bisection bracket.
- [x] 2.2 Confirmed that the GSE branch's `eta_ve, Tii, *d, *Eii_*` outputs flow into the same downstream code (lines 257+) as the inline-iteration branch — back-rotation, plasticity correction, and energy bookkeeping are identical for both branches.
- [x] 2.3 Built `cmake-build-test/`; clean compile (only pre-existing unrelated warnings).
- [x] 2.4 Added `LocalIterationViscoElastic` and `LocalIterationViscoElasticGrainSize` declarations to [RheologyDensity.h](MDLIB/RheologyDensity.h) (where the struct types live) so the call site in `AnisotropyRoutines.c` doesn't need an `extern` declaration.

## 3. Aniso-only regression (must be byte-identical)

- [x] 3.1 `AnisotropyBenchmarkTests`: **9/9 PASS** post-fix.
- [x] 3.2 `AniFstrainSimpleShear` per-step `mesh.aniso_factor_n` mean values match pre-change exactly (residuals `3.9264e-08, 3.2601e-08, 1.7369e-08, 4.5828e-08` byte-identical to the values from the `gate-finite-strain-anisotropy` archived log). The new branch is correctly gated off when `gs = 0`.
- [x] 3.3 No regression — gate is correct.

## 4. GSE-only regression (must be byte-identical)

- [x] 4.1 `RheologyCreepTests`: **13/13 PASS** (12 pre-existing + 1 new coupled smoke).
- [x] 4.2 `GrainSizeSteadyState` and sweep tests pass with values matching pre-change (these scenarios have `anisotropy = 0`, so `ViscosityConciseAniso` is never entered — paths are byte-identical by construction).
- [x] 4.3 Manual smoke run of `SETS/PinchSwellGSE.txt` (anisotropy=0, GSE on): **runs to full completion** in ~200 s wall, no errors, exit 0.

## 5. Coupled scenario `SETS/PinchSwellGSEAniso.txt`

- [x] 5.1 Created [SETS/PinchSwellGSEAniso.txt](SETS/PinchSwellGSEAniso.txt) by copying `SETS/PinchSwellGSE.txt`.
- [x] 5.2 Added `anisotropy = 1` to switches (auto-enables `finite_strain` per [InputOutput.c:1319](MDLIB/InputOutput.c#L1319)).
- [x] 5.3 Calcite layer phase has `aniso_angle = 0` (director vertical → foliation horizontal = bedding-parallel for the horizontal layer; convention: `aniso_angle` is the director angle, with director ⊥ foliation), `aniso_factor = 1.0`, `ani_fstrain = 1`, `ani_fac_max = 4` (calibrated from Carrara marble torsion experiments — Barnhoorn et al. 2004 *JGR* 109, B05203 and Pieri et al. 2001 *Tectonophysics* 330, 119–140 — both report calcite viscous shear-strength anisotropy ratios of 3–5 at saturation).
- [x] 5.4 Matrix phase has `ani_fstrain = 0`, `aniso_factor = 1.0` (kept isotropic).
- [x] 5.5 Bumped `Slim = 1e91` (above `C = 1e90`) for both phases to satisfy the parameter validator.
- [x] 5.6 Added `add_set(PinchSwellGSEAniso)` to [SETS/CMakeLists.txt](SETS/CMakeLists.txt) and created [SETS/PinchSwellGSEAniso.c](SETS/PinchSwellGSEAniso.c) as a verbatim copy of `SETS/PinchSwellGSE.c` with the `.txt` filename swapped.
- [x] 5.7 Built `cmake-exec/PinchSwellGSEAniso/PinchSwellGSEAniso`; manual run confirms `mesh.d_n` evolves under anisotropy (max(d) jumped from IC ~2e-4 to wattmeter steady state ~1.6 in one timestep). At full-resolution `bkg_strain_rate = 1e-14` SI, per-step strain ≈ 0.1% so `ani_fac` only develops over many steps — visible in the smoke test (T6) which uses 10 steps.

## 6. CI smoke test fixture and case

- [x] 6.1 Created [TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt](TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt) — 51×51 × Nt=10 downsized variant of `SETS/PinchSwellGSEAniso.txt`, `writer_subfolder = PinchSwellGSEAnisoSmoke`, `elastic = 0` and `thermal = 0` (matches `PinchSwellGSESmoke` pattern).
- [x] 6.2 Added `RheologyCreep.PinchSwellGSEAniso` GTest case to [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp) using the existing `PinchSwellGSEFixture` (same `SetPhaseLayer` etc.). Reads both `Centers/d` and `Centers/ani_fac` via `readFieldAsArray`.
- [x] 6.3 Assertion block: 4 GSE invariants (mirror `PinchSwellGSESmoke`) + 4 anisotropy invariants (`min_ani >= 0.999`, `max_ani > 1.0`, `max_ani <= 4 + 1e-9`, finite). Caught the silent-no-op condition by design — pre-fix, `max_ani > 1.0` would have failed because `*d` never updated and FS_AR drove no clamp; post-fix it passes (max_ani = 1.539 in 10 steps).
- [x] 6.4 Test passes with comfortable margins: `max(d)/min(d) ≈ 17` (threshold 5), `max(ani_fac) = 1.539` (threshold 1.0), `min(d) = 1e-5 < 1e-3` (layer evolves below `gs_ref`).

## 7. Live-scenario smoke run

- [x] 7.1 Manual smoke of `SETS/RiftingComprehensive.txt` (Nt=2 downsized): pre-existing instability at full 200×160 resolution (exit 134/SIGABRT after step 4 reseeding) is **not caused by this change** — the scenario has `gs = 0` everywhere, so the new GSE branch in `ViscosityConciseAniso` never fires. Verified by checking `materials->gs[phase]` for all phases in the .txt — no `gs = N` line present, default is 0. The unchanged `else` branch runs identically pre- and post-change.
- [x] 7.2 The substantive aniso-only smoke is via `AnisotropyBenchmark.AniFstrainSimpleShear/Saturation/PureShear` (3 tests, all pass byte-identical to pre-change values). Plus full `PinchSwellGSE.txt` (anisotropy=0, GSE on) runs to completion in ~200 s — confirms the GSE-only path is unchanged.

## 8. Documentation

- [x] 8.1 Added §8.3 "Coupled GSE + finite-strain anisotropy: `RheologyCreep.PinchSwellGSEAniso`" to [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md). Covers: the silent-no-op bug this change fixes, the fix approach (port `LocalIterationViscoElasticGrainSize` into the aniso path gated on `gs != 0`), parameter file paths, the assertion block with rationale for each sentinel, measured accuracy table, threshold derivation including the Barnhoorn et al. 2004 calibration of `ani_fac_max = 4`, and runtime budget.
- [x] 8.2 Added a row to the §9 summary table (`Coupled GSE+aniso smoke`) with sentinel-class assertion summary.

## 9. Final verification

- [x] 9.1 Clean-build verification: `AnisotropyBenchmarkTests` (9/9), `RheologyCreepTests` (13/13 — 12 pre-existing + 1 new), all green.
- [x] 9.2 No diagnostic `LOG_INFO`s left in committed code (none added during this change).
- [x] 9.3 Spec delta in `openspec/changes/merge-gse-anisotropy/specs/ci-gse-anisotropy-coupled/spec.md` matches the implementation: gate uses `gs != 0` (synced earlier), test name `RheologyCreep.PinchSwellGSEAniso`, fixture path, all assertion thresholds.
- [x] 9.4 `TESTS/AnalyticalSolutions.md` documentation present at §8.3 + summary-table row.
- [x] 9.5 Self-review: scope stayed within [design D2-D3](design.md). The MDLIB fix is one new branch in one function plus the local `eta_lo` calculation needed by the new branch. No refactor of `ViscosityConcise`/`ViscosityConciseAniso` shared body, no changes to the iteration internals, no architectural rework. Plus the test infrastructure additions: 1 SETS scenario file, 1 SETS .c, 1 CMake registration, 1 CI fixture, 1 GTest case, 1 documentation subsection.
