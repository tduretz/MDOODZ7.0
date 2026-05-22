## Why

MDOODZ has a finite-strain-driven anisotropy-strength evolution code path (`ani_fstrain = 1` per phase, [AnisotropyRoutines.c:46-54](MDLIB/AnisotropyRoutines.c#L46)) that's used live in published-research scenarios (`SETS/RiftingComprehensive.txt`, `SETS/RiftingComprehensive_highres.txt`, `SETS/CollisionPolarCartesianAniso.txt`) but has **zero CI coverage**. The form is `δ = min(FS_AR, ani_fac_max)` where FS_AR is the strain-ellipse aspect ratio computed via polar decomposition of the deformation gradient ([RheologyParticles.c:634-678](MDLIB/RheologyParticles.c#L634)). Both `FS_AR(t)` for pure shear and `FS_AR(γ)` for simple shear have closed-form analytical solutions — testable to machine precision against MDOODZ's own implementation. This change adds those gates before `merge-gse-anisotropy` (the planned follow-up that combines GSE with anisotropy state) builds on top of the δ-evolution path.

## What Changes

- Add a new GoogleTest fixture and 2 tests in `TESTS/AnisotropyBenchmark/` that exercise `ani_fstrain = 1` under prescribed simple shear (unsaturated and saturated regimes), assert that MDOODZ's measured `aniso_factor_n` matches the closed-form `min(FS_AR(γ), ani_fac_max)` to FP precision (relative error < `1e-6`, absolute error at saturation < `1e-9`).
- **Pure-shear evolution is documented but NOT in CI.** Implementation revealed an MDOODZ bug: `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` segfaults at step 2 start (silent SIGSEGV after step 1's clean output write), under all `pure_shear_ALE` settings tried (0, 1, -1) and multiple strain rates. The pure-shear analytical formula `FS_AR(t) = exp(2·ε_dot·t)` is documented in `TESTS/AnalyticalSolutions.md` for completeness. The simple-shear tests gate the same `FiniteStrainAspectRatio → AnisoFactorEvolv` code path. Investigation of the pure-shear segfault is left to a follow-up change.
- Add a CI plot `aniso_factor_evolution.png` (analytical-only, MDOODZ's own form, gnuplot-rendered from the test's `.dat` side-effect) — same pattern as `director_benchmark.png` and `grain_size_benchmark.png`.
- Add a numbered section to [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) documenting the closed-form references (FS_AR for pure/simple shear, derivation from polar decomposition of F).
- **Separate research artefact** in `SETS/AnisoFstrainResearch/` (not in CI):
  - One or two `.txt` scenarios showing `ani_fstrain = 1` in pure/simple shear
  - A markdown comparison note `hansen2012_comparison.md` documenting the form mismatch with the Hansen+12 *Nature* `δ(γ) = 13·tanh(0.25γ) + 1` empirical fit (their Eq. 7) — explicit known-limitation note motivating future work toward a smoother saturation form.
- Use `MutateInput` per-test override pattern from day one (no per-strain-rate fixture proliferation).

## Capabilities

### New Capabilities

- `ci-finite-strain-anisotropy`: gates the `ani_fstrain = 1` code path with closed-form analytical references for pure shear (`FS_AR(t) = exp(2·Eii·t)`) and simple shear (`FS_AR(γ) = 1 + γ²/2 + γ·√(γ²/4+1)`), plus the `min(·, ani_fac_max)` saturation clamp.

### Modified Capabilities

(none — this is a new capability that doesn't change existing requirements)

## Impact

- **Tests**: new `AniFstrainEvolution` base fixture (`TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt`) and 2 GTest cases (`AniFstrainSimpleShear`, `AniFstrainSaturation`) in [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp). Single shared base `.txt` (no fixture proliferation) with `MutateInput` overrides for `bkg_strain_rate`, `Nt`, `dt`, `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, `writer_subfolder`.
- **Code under test**: zero changes — gates existing live code in [AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c), [RheologyParticles.c](MDLIB/RheologyParticles.c), and [InputOutput.c](MDLIB/InputOutput.c).
- **Documentation**: new section in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) (closed-form FS_AR references for both pure shear and simple shear; `min()` clamp behaviour; visualisation embed; pure-shear segfault footnote).
- **Research artefact (NOT CI)**: new `SETS/AnisoFstrainResearch/` directory with `.txt` scenario(s) and `hansen2012_comparison.md` documenting MDOODZ's `min()` form vs Hansen+12's tanh fit — flagged as future-work motivation, not as a CI gate.
- **CI runtime**: 2 tests × 5 short timesteps × 21×21 homogeneous box → ~0.4 s total; well under the existing budget.
- **Risk**: minimal — adds tests against an analytical reference that's a textbook continuum-mechanics result, no new physics.
- **Out of scope**:
  - Changing MDOODZ's `min(FS_AR, ani_fac_max)` form to a smoother saturation (e.g., the commented-out `erfc` blend, or a tanh fit to Hansen+12). Documented as future work in the research note.
  - Coupling with GSE — that's the follow-up `merge-gse-anisotropy` change, blocked on this one.
  - Quantitative CI comparison to Hansen+12 lab data. The research note in `SETS/` makes the comparison qualitatively for visual reference.
  - **Investigation of the pure-shear (`shear_style = 0`) + `ani_fstrain = 1` step-2 segfault.** Documented as known MDOODZ bug in design.md; left to a follow-up change.
