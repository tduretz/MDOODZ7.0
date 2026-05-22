## Context

The crash is reproducible but unexplained. Symptoms from the [`gate-finite-strain-anisotropy` debug session](openspec/changes/archive/2026-05-06-gate-finite-strain-anisotropy/design.md):

- Configuration: `shear_style = 0`, `anisotropy = 1`, per-phase `ani_fstrain = 1`. Single phase, homogeneous, `Nx = Nz = 21`. `pure_shear_ALE ∈ {0, 1, -1}` all crash. `Eii ∈ {0.05, 0.1, 1.0}` all crash. `reseed_markers = 0` does NOT prevent the crash.
- Step 1 runs cleanly — Stokes solve converges, advection updates `particles.F`, ALE resizes the mesh (when on), `Output00001.gzip.h5` is written successfully and contains valid `Centers/ani_fac` (= 1.0 everywhere, i.e. IC).
- The last log line before the crash is the per-step trailer: `Current dt = 5.00e-01 sec, Old dt = 5.00e-01 sec`. There is no "Time step 00002" log line. Process exits with code 139 (SIGSEGV) without any MDOODZ error message.
- Switching to `shear_style = 1` (simple shear) with everything else identical → all 5 steps run cleanly.

The narrow trigger (`shear_style = 0` + `ani_fstrain = 1`) means the crash is in code that is conditional on **both** of those flags, executed at the start of step 2. Candidates:

1. **BC re-application after step-1 mesh reshape (when ALE = 1)**: pure-shear BCs are recomputed from the updated `model.xmin/xmax/zmin/zmax`. Anisotropy state may carry stale pointers or assumptions.
2. **Finite-strain particle field handling at the start of step 2**: `FiniteStrainAspectRatio` runs at line ~658 of `Main_DOODZ.c` using `particles.F` from end of step 1. Pure-shear advection wrote new F values; something in the read-back path may be uninitialised on the second pass only.
3. **Anisotropy P2G after particle deformation**: `mesh.FS_AR_n` and `mesh.aniso_factor_n` are interpolated from particles each step. Pure shear means particles have moved relative to the (possibly reshaped) grid — reseeding/inflow handling may differ.
4. **`PureShearALE` boundary-velocity averaging**: under `pure_shear_ALE = 1`, the routine reads `mesh->u_in`/`v_in` after step 1's solve. Reading these on the second pass may hit a stale buffer or misaligned offset.

We do not know which of these is at fault — bisection is required.

## Goals / Non-Goals

**Goals:**

- Identify the exact MDLIB code path that crashes under `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` at start of step 2.
- Fix the bug with the minimum-scope change consistent with not regressing simple-shear, pure-shear-without-anisotropy, or pure-shear-with-anisotropy-but-static-`ani_fstrain = 0`.
- Add the deferred `AnisotropyBenchmark.AniFstrainPureShear` GTest case (templated in [archived `tasks.md` §5.1](openspec/changes/archive/2026-05-06-gate-finite-strain-anisotropy/tasks.md)) to gate the closed-form `FS_AR(t) = exp(2·ε̇·t)` reference at FP precision.
- Update `TESTS/AnalyticalSolutions.md` §9.7 to reflect active CI coverage instead of the documented deferral.

**Non-Goals:**

- **No redesign of pure-shear time-loop semantics.** If the bisection lands on something requiring an architectural rework (rather than a localised fix), this change pauses and a follow-up is opened. Scope here is "find and fix one specific crash."
- **No change to the `min(FS_AR, ani_fac_max)` form.** Still tracked in [SETS/AnisoFstrainResearch/hansen2012_comparison.md](SETS/AnisoFstrainResearch/hansen2012_comparison.md).
- **No coupling with GSE.** That's `merge-gse-anisotropy`, still pending.
- **No new `pure_shear_ALE` mode or `ParticleInflowCheck` semantic changes.**
- **No spec changes outside `ci-finite-strain-anisotropy`.**
- **No fix-via-runtime-guard.** The aim is to run pure shear cleanly, not to detect-and-exit-with-message.

## Decisions

### D1: Investigation methodology — debug build + valgrind/ASan first, gdb backtrace second

**Choice:** Reproduce the crash under tooling that produces a useful backtrace before doing any code changes. Order of attempts:

1. **Debug build with sanitisers** (`-fsanitize=address,undefined`, `-O0 -g`). ASan catches use-after-free, stack-buffer-overflow, and uninitialised reads with line numbers. UBSan catches integer overflow / null deref. This usually localises the bug to a specific source line on first run.
2. **gdb on a debug build**. If sanitisers don't fire (e.g., the crash is a segfault from a legitimate-looking but corrupt pointer), `gdb -batch -ex run -ex bt` gives the call stack at the crash point.
3. **Valgrind memcheck**, only if (1) and (2) are inconclusive. Slower (~30× wall) but catches some classes of issues sanitisers miss.

**Rationale:** the crash is **silent** (no MDOODZ assertion, no error message), so the working hypothesis is memory corruption — not a logic error MDOODZ knew to trap. Sanitisers are designed for exactly this category. Skipping straight to "read the code and guess" wastes time when we can get a labelled call stack in one CMake reconfigure.

**Alternatives considered:**

- *Reading `Main_DOODZ.c` step-2 setup line-by-line and guessing which path is broken*: rejected. Step 2 setup is hundreds of lines; without a backtrace we'd be running candidates in random order.
- *Adding `LOG_INFO` statements between every step-2 setup call*: workable but slower and noisier than ASan, and only narrows the offending function — doesn't tell us the offending line within it.

### D2: Bisection candidate ordering (when D1 doesn't immediately resolve)

If sanitisers/gdb don't pinpoint the bug (e.g., crash is in a hot loop with stripped frames), bisect by **toggling step-2 setup blocks** in order of suspected blast radius, smallest first:

1. **`FiniteStrainAspectRatio` at line ~658**: temporarily skip the second call (just for diagnosis — definitely not the fix). If skipping it makes step 2 run, the bug is in F-readback or eigenvalue logic for the post-advection F state.
2. **`InitialiseDirectorVector` / per-step director re-init**: similarly skip the second-pass anisotropy state init.
3. **`PureShearALE` boundary-velocity averaging**: temporarily set `pure_shear_ALE = 0` AND `periodic_x = 1` AND keep the box static. We already showed `pure_shear_ALE = 0` alone crashes, so this isolates the BC-reapplication path from the ALE mesh-resize path.
4. **Particle reseeding interaction**: we've already shown `reseed_markers = 0` doesn't prevent the crash, so reseeding is ruled out.

**Rationale:** narrow each diagnostic step to one suspect block. If the crash disappears with that block disabled, we have the location; reasoning about WHY then suggests the fix.

**Alternatives considered:**

- *Bisecting commits (git bisect)*: not applicable. We don't have a "last known good" version of pure-shear + ani_fstrain — it's never been tested in CI, and the live research scenarios that use this combination may have been masking the same bug latently.
- *Bisecting MDOODZ flags one at a time*: already done during `gate-finite-strain-anisotropy`. The trigger pair (`shear_style = 0` + `ani_fstrain = 1`) is fully characterised.

### D3: Fix scope — minimum change, no architecture work

**Choice:** Once the bug is localised, the fix is:

- **One-line / few-line change** if the bug is e.g. an uninitialised array element or a missing reset of a per-step buffer.
- **A single function rewrite** if the bug is e.g. a logic flaw in how `FiniteStrainAspectRatio` handles post-advection F (still localised).
- **Pause this change and open a follow-up** if the bisection lands on something architectural — e.g. requires reordering the time-loop or restructuring how anisotropy state is held across steps.

**Rationale:** scoping discipline. A "fix the segfault" change should not become a "redesign anisotropy state lifecycle" change. The pause-and-escalate clause makes the scope boundary explicit.

**Alternatives considered:**

- *Inline the fix even if it grows large*: rejected. Concentration of risk is the danger here — anisotropy + finite-strain + GSE merger is coming next; we want each piece to be reviewed in isolation.

### D4: Test for the fix — reuse the deferred `AniFstrainPureShear` template, no new fixture

**Choice:** Once the fix lands, add `AnisotropyBenchmark.AniFstrainPureShear` to [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp) using:

- Existing base fixture: [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt)
- Existing override callback: `mutateAniFstrainInput` (already supports `shear_style`, `periodic_x`, `pure_shear_ALE`, `bkg_strain_rate`, `Nt`, `dt`, `ani_fac_max`, `writer_subfolder`).
- Set `setup.SetBCs->SetBCVx = SetPureShearBCVx; SetBCVz = SetPureShearBCVz;` in the test body.
- Override: `{ Eii, Nt=5, dt=0.5, ani_fac_max=1e6, shear_style=0, periodic_x=0, pure_shear_ALE=1, "AniFstrainPureShear" }` (or whatever combo the fix validates).
- Per-step assertion: `EXPECT_NEAR(s.mean / exp(2*Eii*(k-1)*dt), 1.0, 1e-6)` plus the same homogeneity / non-finite / lower-bound checks as the simple-shear test.

**Rationale:** zero fixture proliferation. The original `gate-finite-strain-anisotropy` design already foresaw this test; reusing the existing template means adding ~30 lines of test code and one entry to the spec.

**Alternatives considered:**

- *New `AniFstrainPureShearEvolution.txt` fixture*: rejected — violates the no-fixture-proliferation rule established in the upstream change.

### D5: Strain-rate choice for the pure-shear test

**Choice:** Match the simple-shear test's `bkg_strain_rate = 1.0` non-dim, `Nt = 5`, `dt = 0.5`. This gives `t ∈ {0, 0.5, 1.0, 1.5, 2.0}` and `FS_AR ∈ {1, 2.72, 7.39, 20.09, 54.60}`. Sets `ani_fac_max = 1e6` to stay unsaturated across the whole range.

**Rationale:** symmetry with the simple-shear test simplifies side-by-side reading of `AnalyticalSolutions.md` §9. Same Nt / dt / box / particle settings makes per-step comparison meaningful.

**Caveat:** at `Eii = 1.0`, the per-step strain is 50%, which is what triggered the original `pure_shear_ALE = 1` mesh-reshape problem during `gate-finite-strain-anisotropy` debug. If the fix lets the time loop survive but still has trouble with such large per-step deformation, we drop to `Eii = 0.1` (10% per step, FS_AR_final ≈ exp(1.0) ≈ 2.72) — small range but enough to verify the formula. Decision deferable to implementation: try `Eii = 1.0` first, fall back if needed.

**Alternatives considered:**

- *Lower default strain rate (`Eii = 0.1`)*: less aggressive, but underutilises the analytical formula's range. Use only as a fallback.

### D6: Spec delta — MODIFIED `ci-finite-strain-anisotropy`

**Choice:** Add a new requirement to `ci-finite-strain-anisotropy` and modify the existing simple-shear requirement:

- **ADD** `Requirement: Analytical L2 benchmark for finite-strain anisotropy under pure shear` — closed-form `FS_AR(t) = exp(2·ε̇·t)`, gated to `< 1e-6` relative.
- **MODIFY** `Requirement: Analytical L2 benchmark for finite-strain anisotropy under simple shear` — drop the "Pure shear is NOT exercised in CI" paragraph from the requirement description; replace with a reference to the pure-shear requirement.
- **MODIFY** `Requirement: Benchmark detects NaN propagation from finite-strain integration` — already says "any of the benchmark simulations"; no edit needed but verify the wording covers both tests cleanly.
- **MODIFY** `Requirement: Single base fixture file with MutateInput per-test overrides` — update from "two test cases" to "all three test cases" (or "the test cases").
- **MODIFY** `Requirement: Benchmark generates a CI visualisation aniso_factor_evolution.png` — no edit, the saturation test still owns the `.dat` file.

**Rationale:** delta isolated to lifting the pure-shear deferral. The simple-shear test's behaviour does not change; only its sibling-test list grows.

## Risks / Trade-offs

- **[Risk] Sanitisers don't reproduce the crash.** ASan/UBSan miss certain classes (e.g. legitimate-looking but logically wrong pointer arithmetic that doesn't violate memory-model rules). → **Mitigation**: gdb backtrace as Plan B (D1 step 2). If both fail, fall back to D2 manual bisection.
- **[Risk] The fix has wider blast radius than expected.** Touching pure-shear time-loop code can affect `pure_shear_ALE` semantics, BC re-application timing, or particle field integration — paths used by [SETS/RiftingComprehensive*.txt](SETS/RiftingComprehensive.txt), [SETS/AnisoVEVP.txt](SETS/AnisoVEVP.txt), and the existing CI tests (`StressAnisotropy`, `StressAnisotropyL2`). → **Mitigation**: full `AnisotropyBenchmarkTests` + `RheologyCreepTests` + (manual) live-scenario smoke run before declaring done. Specifically run `AnisoVEVP` from `SETS/` for at least 2 timesteps. If anything regresses, investigate and resolve before merging.
- **[Risk] Bisection lands on architectural issue beyond the scope wall (D3).** → **Mitigation**: pause this change, file `redesign-pure-shear-anisotropy-state` (or similar) as a new openspec change with the full diagnostic captured. Don't try to expand scope here.
- **[Risk] At `Eii = 1.0`, the fix lets step 2 run but step 3+ still misbehaves due to large per-step deformation.** → **Mitigation**: D5 fallback to `Eii = 0.1`; document range at which the test is FP-precision-stable.
- **[Trade-off] Test parameters mirror simple-shear (Eii=1, Nt=5, dt=0.5) for symmetry**, even though for purely diagnostic-of-the-fix purposes a tiny `Eii` would suffice. The symmetry buys readability in `AnalyticalSolutions.md`; the cost is a tighter coupling to whatever the fix's stable strain-rate range turns out to be.

## Migration Plan

1. **Reproduce the crash** in a debug build with ASan enabled: `cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS=-fsanitize=address -DCMAKE_CXX_FLAGS=-fsanitize=address …`. Run `AnisotropyBenchmarkTests --gtest_filter='*AniFstrainPureShear*'` (test added back temporarily, then deleted again until the fix lands).
2. **Capture the backtrace** (or ASan report). Document in design.md "Open Questions" if not pinned.
3. **Isolate the bug** to a single MDLIB function. If D1 succeeds, this is automatic. If not, run D2 bisection.
4. **Apply the minimum-scope fix.** Verify it with the same crashing configuration in a release build.
5. **Run full regression suites**: `AnisotropyBenchmarkTests` (currently 8 tests), `RheologyCreepTests` (12), plus a 2-step smoke run of `AnisoVEVP` from `SETS/` (manual).
6. **Add `AnisotropyBenchmark.AniFstrainPureShear`** per D4, calibrate thresholds (expect ~`1e-7` to ~`1e-6` per step like simple shear).
7. **Update [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md)** §9: replace §9.7 deferral footnote with active coverage; add a measured-accuracy table next to the simple-shear one.
8. **Update spec delta** at [openspec/specs/ci-finite-strain-anisotropy/spec.md](openspec/specs/ci-finite-strain-anisotropy/spec.md) per D6 — happens at archive sync.

Rollback: revert each step in reverse order. The MDLIB fix should be a single small commit easy to revert if step 5 finds a regression.

## Open Questions (resolved during implementation)

- ~~**What is the actual root cause of the crash?**~~ **Resolved.** The crash is at [ParticleRoutines.c:483](MDLIB/ParticleRoutines.c#L483):
  ```c
  if (j_part>Nx-2) {
      LOG_ERR("Should never be here II! (Vertices2Particle)"); exit(1);
      j_part = Nx-2;  // dead code — exit(1) above prevents reaching this fallback
  }
  ```
  `Vertices2Particle` is called from the director-rotation loop in [`UpdateParticleStress`](MDLIB/RheologyParticles.c#L1453) (gated on `anisotropy == 1 && advection == 1`, NOT specifically on `ani_fstrain`). Under pure shear with ALE mesh-shrink, particles can drift outside the new mesh in x; the bound check fires.

  **Build-mode-dependent crash signature confirms the bug**:
  - Release `-O3`: silent SIGSEGV from invalid `NodeField[iNE]` access (the buffered `LOG_ERR` is lost when the OMP thread is killed before flush). Sometimes `LOG_ERR; exit(1)` is reordered/elided by the optimiser when its side effects are deemed unobservable along the executed path; the array access on the next line then segfaults.
  - Debug `-O0` (no sanitisers): explicit `LOG_ERR` fires and prints the message; `exit(1)` follows.
  - Debug `-O0` + ASan/UBSan: passes — different binary layout shifts particle positions slightly relative to the post-ALE mesh boundary; the trigger condition isn't reached. **Note**: this means the bug is silent in instrumented builds but manifests in release. Sanitisers are NOT a substitute for the fix.

  Three of the four boundary checks in the same function (cases I, IV, V at [ParticleRoutines.c:476-498](MDLIB/ParticleRoutines.c#L476)) already gracefully clamp, with their `LOG_ERR / exit` calls commented out. **Only case II (line 483) still has an active `exit(1)`** — almost certainly an oversight when the others were disabled. The dead-code clamp `j_part = Nx-2;` immediately after is the intended fallback.

  **The bug is not specific to `ani_fstrain = 1`.** Any pure-shear scenario with `anisotropy = 1` and `advection = 1` running for ≥ 2 steps with `pure_shear_ALE` mesh shrink would hit it. The reason no existing test caught it: all `AnisotropyBenchmark` tests prior to this change used simple shear (`shear_style = 1, periodic_x = 1`), where particles wrap and never leave the mesh.

- **Is the crash also present in the live research scenarios** (`SETS/RiftingComprehensive*.txt`, `SETS/CollisionPolarCartesianAniso.txt`)? **Likely yes, but masked by their specific dt × strain rates and SI scaling.** Their `bkg_strain_rate ~ 1e-15`/`-16` produces per-step displacements smaller than mesh resolution × float epsilon → particles never measurably exceed mesh bounds. Worth a 2-step smoke run after the fix anyway (Task 4.4).
- **Does the fix interact with `pure_shear_ALE = -1` + `ParticleInflowCheck`** specifically? The fix lives in `Vertices2Particle` which is called regardless of `pure_shear_ALE` mode. The clamp behaviour is mode-independent — if a particle is outside the mesh under any mode, the fix gracefully clamps. Verified during regression (Task 4).
- ~~**Can `Eii = 1.0` actually run cleanly post-fix?**~~ **Use `Eii = 0.1` per D5 fallback** — the test demonstrates the fix and the analytical formula at modest strain rate. `Eii = 1.0` with 50% per-step strain remains an open question for the user, not a CI requirement.
