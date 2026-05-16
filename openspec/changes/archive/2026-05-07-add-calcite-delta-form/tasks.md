## 1. MDLIB implementation — `MDLIB/FlowLaws.c` (no `AnisotropyRoutines.c` change)

- [x] 1.1 Added `static double anisoDelta_Calcite(double FS_AR)` to [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c) with `M_inf=0.50`, `gamma_e=3.50`, `slope=6.0`.
- [x] 1.2 Comment block in place documenting Pieri+01 / Bruijn+11 calibration sources, Skemer+05 J→M conversion, Tommasi+09 VPSC slope bounds, asymptote δ_∞ = 4 matching empirical PinchSwell, and HPT non-applicability.
- [x] 1.3 `case 2 :` added to `ReadDataAnisotropy` switch — assigns `aniso_delta_fn[k] = anisoDelta_Calcite`, logs the citation.
- [x] 1.4 Verified — auto-default in InputOutput.c and Main_DOODZ.c routes `ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1` (olivine). Additionally extended Main_DOODZ.c to (re-)call `ReadDataAnisotropy` whenever `aniso_db > 0` post-`MutateInput` so test/scenario harnesses that flip aniso_db after parse get the function pointer populated correctly. (This is a small fix to the post-MutateInput path; non-MutateInput scenarios are unchanged.)
- [x] 1.5 Verified — zero changes to AnisotropyRoutines.c, mdoodz.h, mdoodz-private.h. The function-pointer architecture absorbs the calcite addition entirely in FlowLaws.c (proves the §7 refactor pays off).

## 2. Backward-compatibility verification

- [x] 2.1 All 10 pre-existing AnisotropyBenchmarkTests pass byte-identically (HansenOlivine still hits 1.94e-8 max relative error at γ=20). Now 11 with the new AniFstrainCalcite test.
- [x] 2.2 All 13 RheologyCreepTests pass byte-identically.
- [x] 2.3 Verified — `ani_fstrain ∈ {0, 1, 3}` paths untouched; `aniso_db = 1` olivine helper byte-identical; auto-default still routes ani_fstrain=2+aniso_db=0 to olivine.
- [x] 2.4 PinchSwellGSEAnisoSmoke test in RheologyCreepTests passed byte-identically (`d=[1.000e-05,1.704e-04] r=17.04, ani=[1.000e+00,1.539e+00] r=1.54`). Confirms the existing calcite `ani_fstrain = 1` path is unaffected by adding case 2.

## 3. New CI test — `AnisotropyBenchmark.AniFstrainCalcite`

- [x] 3.1 GTest case `AniFstrainCalcite` added to [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp) at end-of-file. Reuses `AniFstrainEvolution.txt` base fixture via MutateInput lambda. No new fixture.
- [x] 3.2 Override callback sets `ani_fstrain[0]=2`, `aniso_db[0]=2`, `ani_fac_max=100`, `bkg_strain_rate=1.0`, `Nt=11`, `dt=0.5`, simple-shear BCs, `writer_subfolder="AniFstrainCalcite"`.
- [x] 3.3 Analytical reference `anaDeltaCalcite(gamma_dot, t)` defined with `kCalciteMInf=0.50`, `kCalciteGammaE=3.50`, `kCalciteSlope=6.0` mirroring MDLIB constants.
- [x] 3.4 Per-step `EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)` assertion in place.
- [x] 3.5 Per-step `EXPECT_LT(s.l2_rel, 1e-4)` and `(max/min)-1 < 1e-4` homogeneity assertions in place.
- [x] 3.6 Test passes at γ ∈ {0, 1, 2, ..., 10}. δ at γ=10 = 3.828 (95.7% of asymptote 4.0). Max L2 error 4.5e-8, well below 1e-6 threshold.
- [x] 3.7 Full AnisotropyBenchmarkTests suite passes — 11/11 in 6.8 s (was 10/10 in 6.4 s pre-change).

## 4. Calcite calibration artefacts — `SETS/AnisoFstrainResearch/calcite/`

- [x] 4.1 Created `SETS/AnisoFstrainResearch/calcite/`.
- [x] 4.2 [calibrate_calcite.py](SETS/AnisoFstrainResearch/calcite/calibrate_calcite.py) written — pure stdlib + numpy + matplotlib. Inline 5-point Pieri+01 + 3-point Bruijn+11 datasets. 2D grid search on (M_inf, γ_e) with slope=6.0 fixed. RMS comparison: `min(FS_AR,4)` = 1.571, Hansen olivine form = 5.66, **NEW calcite case 2 = 0.79** (50% improvement over existing default).
- [x] 4.3 [mdoodz_vs_calcite_data.png](SETS/AnisoFstrainResearch/calcite/mdoodz_vs_calcite_data.png) rendered — two-panel comparison plot showing the new calcite curve threading through Pieri+01 + Bruijn+11 lab data clouds, with min(FS_AR, 4) and Hansen olivine forms as references.
- [x] 4.4 [calcite_calibration.md](SETS/AnisoFstrainResearch/calcite/calcite_calibration.md) written — 10 sections covering functional form, calibration procedure, asymptote-uncertainty discussion (slope=6.0 calibrated to empirical PinchSwell, NOT direct lab fit), per-datapoint verification table, Pieri vs Bruijn disagreement at γ=5, γ=11 non-monotonic outlier, simple-shear caveat, HPT non-applicability, J→M Skemer+05 conversion caveat, reproduction instructions, citations.
- [x] 4.5 Verified — relative-path references `calibrate_calcite.py` (lines 39, 210) and `mdoodz_vs_calcite_data.png` (lines 123, 219) all check out.
- [x] 4.6 No CI-gating references. Single soft cross-reference from `TESTS/AnisotropyBenchmarkTests.cpp` is a code comment (`// See SETS/AnisoFstrainResearch/calcite/`) pointing at this research note for context — same pattern as the existing AnalyticalSolutions.md → hansen2012_comparison.md cross-reference, consistent with the spec's intent ("research-grade, not regression gate").
- [x] 4.7 Created [mdoodz_vs_calcite_comparison.py](SETS/AnisoFstrainResearch/calcite/mdoodz_vs_calcite_comparison.py) mirroring `mdoodz_vs_hansen_comparison.py`. Reads HDF5 outputs from `cmake-build-test/TESTS/AniFstrainCalcite/`, overlays 11 red stars on the calibrated curve. Max MDOODZ-vs-formula relative error 1.62e-7 confirms dispatch + integration + HDF5 round-trip preserves the calibration to single-precision FP.

## 5. Quartz documentation — `SETS/AnisoFstrainResearch/quartz/` (no MDLIB code)

- [x] 5.1 Created `SETS/AnisoFstrainResearch/quartz/`.
- [x] 5.2 [calibrate_quartz.py](SETS/AnisoFstrainResearch/quartz/calibrate_quartz.py) written — pure stdlib + numpy + matplotlib. Inline 17-point Pennacchioni+10 dataset (4 WDV + 5 MDV + 8 SDV samples). Saturation check: SDV regime mean J=10.81, std=0.90 — confirms Pennacchioni+10's "no further strengthening past γ ≈ 4" claim. RMS comparison: `min(FS_AR, 2.5)` = 0.516 (best), `min(FS_AR, 3.0)` = 0.616 (recommended for slight conservatism).
- [x] 5.3 [mdoodz_vs_quartz_data.png](SETS/AnisoFstrainResearch/quartz/mdoodz_vs_quartz_data.png) rendered — two-panel plot with regime-coded markers (WDV/MDV/SDV) and `min(FS_AR, ani_fac_max)` curves at 2.0, 2.5, 3.0. Visual case for "min form is fit-for-purpose" is clear.
- [x] 5.4 [quartz_calibration.md](SETS/AnisoFstrainResearch/quartz/quartz_calibration.md) written — 8 sections covering Pennacchioni+10 dataset summary, the explicit "no further strengthening past γ ≈ 4" verbatim quote, min-form fit-for-purpose RMS argument, recommendation to use `ani_fstrain = 1 + ani_fac_max = 3.0`, three independent reasons not to add quartz `case 3` (data quality + T-aware-dispatch needed + asymptote is small), T-range caveat (~500°C, prism-a slip), reproduction instructions, full Pennacchioni+10 citation.
- [x] 5.5 Verified NO new MDLIB code for quartz: `grep -i quartz MDLIB/FlowLaws.c` shows only pre-existing PowerLaw entries (cases 13, 22, 26, 27, 28 — flow laws, NOT anisotropy) plus the inline comment `// Future entries (deferred — see SETS/AnisoFstrainResearch/quartz/)` in `ReadDataAnisotropy` switch. No quartz case in switch, no new helper, no new mat_prop entry.
- [x] 5.6 Verified — `quartz_calibration.md` references `calibrate_quartz.py` (line 162) and `mdoodz_vs_quartz_data.png` (line 99) via relative paths.
- [x] 5.7 Verified — no references to `quartz_calibration` or `calibrate_quartz` from any TESTS/ file (`grep -rn ... TESTS/` returns empty). Research-grade, not regression gate.
