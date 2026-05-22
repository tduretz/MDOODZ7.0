## Why

MDOODZ's pure-shear time loop crashes silently (SIGSEGV) at the start of step 2 when `shear_style = 0` is combined with `anisotropy = 1` and per-phase `ani_fstrain = 1`. Step 1 alone runs cleanly and writes `Output00001.gzip.h5`; the crash happens after that, before any "Time step 00002" log line. Reproduced under all `pure_shear_ALE` settings (0, 1, -1) and at multiple strain rates (`Eii ∈ {0.05, 0.1, 1.0}`). Discovered during the just-archived `gate-finite-strain-anisotropy` change, where it forced the pure-shear evolution test to be dropped from CI (see [archived design.md "Open Issues"](openspec/changes/archive/2026-05-06-gate-finite-strain-anisotropy/design.md)).

Three reasons to fix this now:

1. **Live research scenarios depend on this combination.** [SETS/RiftingComprehensive.txt](SETS/RiftingComprehensive.txt), [SETS/RiftingComprehensive_highres.txt](SETS/RiftingComprehensive_highres.txt), and [SETS/CollisionPolarCartesianAniso.txt](SETS/CollisionPolarCartesianAniso.txt) all set `shear_style = 0` and `ani_fstrain = 1` per phase. They run today only because their first time step has not been hit by whatever step-2 setup path triggers the crash in the test fixture — but the crash mechanism is presumably present and could be triggered by parameter perturbations.
2. **`merge-gse-anisotropy` (also in `openspec/changes/`) is gated on this.** That change combines GSE with anisotropy state including δ-evolution, and its planned `PinchSwellGSEAniso` scenario is a pure-shear setup. Building on top of a known step-2 crash means tests there inherit the same simple-shear-only limitation we just documented.
3. **The closed-form pure-shear analytical reference exists and is testable.** `FS_AR(t) = exp(2·ε̇·t)` was derived during `gate-finite-strain-anisotropy` and documented in `TESTS/AnalyticalSolutions.md` §9.2 specifically because we couldn't run it. Once the crash is fixed, the deferred `AniFstrainPureShear` GTest case (templated in [archived tasks.md §5.1](openspec/changes/archive/2026-05-06-gate-finite-strain-anisotropy/tasks.md)) drops in.

## What Changes

- **Investigate** the step-2 segfault under `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1`. Bisect MDLIB step-2 setup paths (BC re-application after step-1 mesh reshape, finite-strain particle field handling, anisotropy P2G, advection bookkeeping) until the offending code is isolated. Prefer running under valgrind / ASan / debug build to catch the actual fault rather than just the SIGSEGV.
- **Fix the bug** in MDLIB. Scope is bounded by what the bisection finds — likely a single uninitialised-state or stale-pointer issue in the pure-shear branch. Out-of-scope: redesigning the pure-shear time loop or anisotropy update order.
- **Add the deferred CI test** `AnisotropyBenchmark.AniFstrainPureShear` from `gate-finite-strain-anisotropy` §5.1: simple homogeneous pure-shear, `Nt = 5`, asserts `mean(Centers/ani_fac) ≈ exp(2·ε̇·(k-1)·dt)` to FP precision (`< 1e-6` relative). Reuses the existing `AniFstrainEvolution.txt` base fixture and `mutateAniFstrainInput` callback.
- **Update [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §9.7**: replace the "Pure shear — analytical only (NOT in CI)" footnote with active CI coverage; add the pure-shear measured-accuracy table next to the existing simple-shear one.

## Capabilities

### New Capabilities

(none)

### Modified Capabilities

- `ci-finite-strain-anisotropy`: lift the documented "pure shear is NOT exercised in CI" limitation. Add a requirement for the pure-shear analytical L2 benchmark (`FS_AR(t) = exp(2·ε̇·t)`). Update the simple-shear requirement's note about pure-shear not being tested.

## Impact

- **Code under change**: localised MDLIB fix in whatever step-2 setup path is at fault. Not yet known — bisection is the first task. Likely candidates from the conversation context: pure-shear BC re-application after `pure_shear_ALE` mesh reshape, finite-strain particle field reseeding interaction, or anisotropy P2G in the pure-shear branch.
- **Tests**: one new GTest case `AnisotropyBenchmark.AniFstrainPureShear` in [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp). Zero new fixture files — reuses `AniFstrainEvolution.txt` via the existing `mutateAniFstrainInput` callback.
- **Documentation**: §9.7 in `TESTS/AnalyticalSolutions.md` flips from "deferred" to "active CI gate"; spec delta covers the requirement change.
- **CI runtime**: +1 test × 5 short timesteps × 21×21 box → ~0.2 s. Within budget.
- **Risk**: **medium-high** until bisection completes — touching the pure-shear time loop has wider blast radius than the additive `gate-finite-strain-anisotropy` change. Mitigations: existing CI coverage (8 AnisotropyBenchmark tests + 12 RheologyCreep tests + all archived rifting/anisotropy live scenarios) catches regressions; debug-build/valgrind-first investigation reduces guesswork.
- **Out of scope**:
  - Changes to the `min(FS_AR, ani_fac_max)` saturation form. Still tracked in `SETS/AnisoFstrainResearch/hansen2012_comparison.md` as future work.
  - GSE coupling (still `merge-gse-anisotropy`).
  - Rewriting `pure_shear_ALE` semantics.
  - Cross-`shear_style` semantic harmonisation. Just fix the crash on the documented configuration.
