## Context

`ani_fstrain = 1` ([InputOutput.c:1452](MDLIB/InputOutput.c#L1452)) is the per-phase switch that selects evolving anisotropy strength via [AnisotropyRoutines.c:46-54](MDLIB/AnisotropyRoutines.c#L46):

```c
double AnisoFactorEvolv(double FS_AR, double aniso_fac_max) {
  return MINV(FS_AR, aniso_fac_max);
}
```

`FS_AR` is the strain-ellipse aspect ratio (= ratio of largest to smallest principal stretch) computed at [RheologyParticles.c:634-678](MDLIB/RheologyParticles.c#L634) by:
1. Per-particle right Cauchy-Green tensor: `C = F^T·F`
2. Right stretch tensor: `U = √C` via closed-form 2D matrix square-root
3. Eigenvalues of U → principal stretches `(e1, e2)`
4. **`FS_AR = e1/e2`**
5. P2G interpolation onto cell-centre `mesh.FS_AR_n` and vertex `mesh.FS_AR_s`

The deformation gradient `F` is integrated per-particle by `DeformationGradient` ([Main_DOODZ.c:1275](MDLIB/Main_DOODZ.c#L1275)), which auto-enables when `model.anisotropy == 1` ([InputOutput.c:1319](MDLIB/InputOutput.c#L1319)).

This code path is used in three live research scenarios (`SETS/RiftingComprehensive*.txt`, `SETS/CollisionPolarCartesianAniso.txt`) but has no CI regression gate. The Schmalholz & Duretz 2017 GSE work (just shipped) and the Duretz et al. 2025 anisotropy work (also shipped) both validated their pieces in isolation; the next change `merge-gse-anisotropy` will combine GSE with anisotropy state including δ-evolution. Gating `ani_fstrain = 1` first lets that merge proceed against a confirmed-correct δ-evolution path.

## Goals / Non-Goals

**Goals:**

- CI tests that exercise `ani_fstrain = 1` under prescribed pure shear and simple shear and assert MDOODZ's measured `aniso_factor_n` matches the closed-form `min(FS_AR(t), ani_fac_max)` to numerical-method precision.
- Test the saturation clamp explicitly (run past `t_sat` with a low `ani_fac_max`) — this is the entire physical content of the simple `min()` form.
- Single CI plot `aniso_factor_evolution.png` showing the analytical curve and the MDOODZ measurements — same pattern as `director_benchmark.png` and `grain_size_benchmark.png`.
- Research-artefact deliverable in `SETS/AnisoFstrainResearch/`: a `.txt` scenario file plus a markdown findings record `hansen2012_comparison.md` documenting the form mismatch with Hansen et al. 2012 *Nature*'s tanh fit. The markdown is a **first-class deliverable**, not optional documentation — it records the known limitation of MDOODZ's simple `min()` form for future researchers and motivates the smoother saturation forms (`erfc` blend already commented out in `AnisoFactorEvolv`, or direct tanh).
- No fixture-file proliferation: single base `.txt` + `MutateInput` per-test override pattern from day one.

**Non-Goals:**

- Changing MDOODZ's `min(FS_AR, ani_fac_max)` saturation form. The mismatch with Hansen+12 is documented in the research note as future work; this change gates only what's already in the code.
- Coupling with grain-size evolution. That's `merge-gse-anisotropy`, blocked on this.
- Quantitative CI comparison to Hansen+12 lab data. Digitisation noise (5–10 % per point) is comparable to the MDOODZ-vs-tanh form mismatch we'd be flagging — using digitised data as a CI gate would conflate "MDOODZ has a bug" with "the model is approximate". The research note compares the two analytically (MDOODZ's `min()` curve vs Hansen+12 Eq. 7's `13 tanh(0.25γ) + 1`) without involving lab data points.
- Per-phase `ani_fstrain` heterogeneity testing. Single-phase scenarios suffice to gate the code path.
- Convergence-order tests in `dt`. The existing director benchmark already gates the time integration of state variables under simple shear; adding redundant dt-convergence here doesn't pull its weight.

## Decisions

### D1: Test fixture — single base `.txt` + `MutateInput` overrides

**Choice:** One base file `TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt` with `Nb_phases = 1`, `anisotropy = 1`, `ani_fstrain = 1`, `ani_fac_max` and other params at sensible defaults. Two GoogleTest cases share this base via `MutateInput` overrides:

- `AniFstrainSimpleShear` — overrides `bkg_strain_rate`, `Nt`, `dt`, `shear_style = 1`, `periodic_x = 1`, `writer_subfolder`. Constant simple-shear, asserts `δ(γ) ≈ 1 + γ²/2 + γ·√(γ²/4+1)` while unsaturated. Also writes `aniso_factor_evolution.dat` for the CI plot.
- `AniFstrainSaturation` — same plus `ani_fac_max = 4`, runs past `γ_sat`. Asserts `δ` clamped at `ani_fac_max` to FP precision.

The override struct supports `bkg_strain_rate`, `Nt`, `dt` (sets `dt`, `dt0`, `dt_start` together), `ani_fac_max`, `shear_style`, `periodic_x`, `pure_shear_ALE`, and `writer_subfolder`.

**Pure shear is deferred (not in CI).** Implementation revealed that MDOODZ's `shear_style = 0` + `anisotropy = 1` + `ani_fstrain = 1` configuration segfaults at the start of step 2 (after step 1 writes its output successfully) under all `pure_shear_ALE` settings tried (0, 1, -1) and at multiple strain rates. The crash is silent — no assertion message, just `exit(SIGSEGV)`. Step 1 alone runs cleanly. Investigation of this MDOODZ bug is left for a follow-up; this change covers what is testable now. The pure-shear analytical formula `FS_AR(t) = exp(2·ε_dot·t)` is documented in `TESTS/AnalyticalSolutions.md` for completeness.

**Rationale:** Matches the convention established in `gate-finite-strain-anisotropy`'s upstream changes (`fix-grain-size-evolution`, `extend-gse-ci-coverage`). One fixture file vs two is consistent with the no-fixture-proliferation rule. Each test sets `g_aniFstrainOverride = {…}` and calls `RunMDOODZ` on the shared `.txt`. Results land in `<subfolder>/Output*.gzip.h5` — no file collisions. The simple-shear formula exercises the same `FiniteStrainAspectRatio` → `mesh.FS_AR_n` → `AnisoFactorEvolv` → `mesh.aniso_factor_n` code path as pure shear would; both formulas validate the underlying `δ = min(FS_AR, ani_fac_max)` form.

**Alternatives considered:**
- *(A) One `.txt` per scenario:* would mean 2+ fixtures plus duplication for the saturation variant — same problem the GSE refactor solved.
- *(B) Reuse `DirectorEvolution.txt` as the base:* tempting since the fixture is already there, but `DirectorEvolution.txt` doesn't have `ani_fstrain = 1` set. Would require either changing its meaning (breaks `DirectorEvolution` test) or overriding too many fields per call. Cleaner to have a dedicated fixture.
- *(C) Force pure-shear to work via `pure_shear_ALE = -1` + `ParticleInflowCheck`:* tried, still segfaults at step 2.
- *(D) Force pure-shear to work via `Nt = 1`:* runs cleanly but then only verifies the IC (FS_AR = 1, δ = 1), which is already covered by step 1 of the simple-shear tests (γ = 0). Adds no coverage.

### D2: Saturation regime test parameters

**Choice:** For `AniFstrainSaturation`, set `ani_fac_max = 4`, `bkg_strain_rate = 1.0` (so `γ_dot = 2`), `Nt = 5`, `dt = 0.5`. Output indexing reads `Output<k>` to reflect `F` from end of step `k-1`, so `γ` at output step `k` is `γ_dot · (k-1) · dt ∈ {0, 1, 2, 3, 4}`. The simple-shear formula `FS_AR(γ) = 1 + γ²/2 + γ·sqrt(γ²/4+1)` gives `FS_AR ∈ {1, 2.62, 5.83, 10.91, 17.94}`. Saturation strain `γ_sat ≈ 1.94` (where `FS_AR = 4`) falls between step 2 (γ=1) and step 3 (γ=2). So output steps 1–2 are unsaturated and 3–5 are saturated. This exercises both regimes in one short run.

Assertions per output step:
- Unsaturated (`γ < γ_sat`): `δ ≈ FS_AR(γ)` to within FP cancellation (relative tolerance `1e-6`)
- Saturated (`γ > γ_sat`): `δ == ani_fac_max` to FP precision (absolute tolerance `1e-9`)
- The transition itself is automatic by construction — no explicit "transition" assertion needed.

**Rationale:** Picking a small `ani_fac_max` keeps the saturation transition inside the test's wall-clock budget. The actual numerical value (4) is arbitrary — any value > 1 works for testing the clamp.

### D3: CI plot — analytical-only, MDOODZ form

**Choice:** Single-panel plot `aniso_factor_evolution.png`:

- x: shear strain γ (or accumulated strain ε for pure shear), linear, range [0, 3]
- y: anisotropy strength δ, log-scale, range [1, ~20]
- **Solid black line**: `FS_AR(γ)` analytical (unbounded, just for reference)
- **Dashed black line**: `δ(γ) = min(FS_AR(γ), ani_fac_max)` (the MDOODZ formula, with knee)
- **Red triangles**: MDOODZ measurements (4–6 points spanning unsaturated and saturated regimes)
- **Grey dashed horizontal line** at y = `ani_fac_max`, labelled
- **Grey dashed vertical line** at x = saturation strain, labelled

**Rationale:** Same template as `grain_size_benchmark.png` — analytical curves + MDOODZ overlay. Knee shape is the visual signature of the `min()` clamp; if MDOODZ's implementation is correct, all red triangles sit on the dashed black line.

### D4: Research note — `SETS/AnisoFstrainResearch/hansen2012_comparison.md`

**Choice:** First-class markdown deliverable next to a `SETS/AnisoFstrainResearch/AniFstrainSimpleShear.txt` scenario. Contents:

1. **What MDOODZ implements**: `δ = min(FS_AR, ani_fac_max)` with the live formula and code references.
2. **What Hansen et al. 2012** *Nature* **fitted from olivine torsion experiments**: `δ(γ) = 13·tanh(0.25γ) + 1`, asymptote at 14, half-saturation at γ ≈ 2.2.
3. **Quantitative comparison** at several γ values (table or short plot, not in CI):
   - γ = 1: MDOODZ ~2.6 vs Hansen ~4.2
   - γ = 2: MDOODZ ~5.8 vs Hansen ~7.0
   - γ = 4: MDOODZ saturates (if `ani_fac_max = 14`) vs Hansen ~10.9
4. **What this means**: the `min()` form is too sharp at the knee and rises too fast at low strain compared to the empirical fit — it's a crude saturation that gets the asymptote right but the shape wrong.
5. **Future work**: replace with the commented-out `erfc` blend or with the tanh form directly, calibrated to Hansen+12 constants. **Out of scope for this change.**
6. **Citation**: full reference to Hansen, Zimmerman, Kohlstedt (2012) *Nature* 492, 415–418, doi:10.1038/nature11671.

**Rationale:** This is research-grade documentation — when someone later wonders "why is the simulated δ-vs-strain shape too sharp?", they need the trail back to Hansen+12 and the explicit acknowledgement that MDOODZ's form is approximate. Putting this in `TESTS/` would conflate "regression gate" with "scientific limitation note"; putting it in `SETS/AnisoFstrainResearch/` co-locates it with the scenario it characterises and matches the existing flat-or-grouped `SETS/` convention (e.g. `SETS/NeckinReviewSystematics/`, `SETS/ShearTemplateCI/`).

The directory contains:
- `hansen2012_comparison.md` — the findings record (this is **the** research deliverable)
- `AniFstrainSimpleShear.txt` — the scenario the markdown references; can be run by hand for visualisations
- Optional: a manually-rendered PNG showing the form mismatch overlay (not auto-generated by CI; user runs gnuplot once and commits the PNG if desired)

### D5: L2 thresholds — calibrated to observed FP-precision agreement

**Choice (after measurement, not just prediction):**
- Unsaturated mean (relative): `EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)` — observed L2 is `~3e-8` to `~5e-8` per step
- Unsaturated spatial L2: `EXPECT_LT(s.l2_rel, 1e-4)` — observed `~3e-8`, headroom for future P2G changes
- Saturated absolute: `EXPECT_NEAR(s.mean, ani_fac_max, 1e-9)` — observed exactly `4.000000e+00` (FP-equality)
- Spatial homogeneity: `EXPECT_LT(max/min - 1, 1e-4)` — homogeneous setup, P2G should preserve uniformity (small residual from float32 HDF5 write)
- Lower-bound sanity: `EXPECT_GE(s.min_v, 0.999)` — FS_AR is bounded below by 1.0; the 0.001 slack accommodates float32 round-trip through HDF5

**Rationale:** Tighter thresholds than the GSE benchmarks (`5e-3`) because there's no Newton iteration here — F evolves by simple matrix multiplication per step, and the analytical reference is from the same closed-form mathematics MDOODZ implements. The expected residual is FP cancellation noise plus single-precision round-trip via HDF5 (fields are stored as `float`). Initial `1e-10`/`1e-12` thresholds were too tight precisely because of the float32 HDF5 storage; observed errors are `~1e-8` for the relative mean, sized by `float`'s ~7-digit precision against values around 1–18.

## Risks / Trade-offs

- **[Risk] Floating-point integration of F over many timesteps drifts.** The deformation gradient is updated step-by-step via velocity gradient; long-running simulations could accumulate drift. → Mitigation: tests run `Nt ≤ 5` steps. For longer runs (research scenarios in `SETS/`), drift is a known concern but out of scope for CI.
- **[Risk] P2G interpolation introduces non-uniformity in homogeneous problems.** F is per-particle, FS_AR computed per-particle, then interpolated to grid. → Mitigation: tests use a single phase with all particles in identical strain → all particles have identical F → P2G is identity. Homogeneity check (`max/min < 1.000001`) catches any subtle interpolation noise.
- **[Risk] `MutateInput` overrides interact with `model.dt0`, `model.dt_start`** (lessons from `extend-gse-ci-coverage` debug). → Mitigation: copy the multi-field `dt`-override pattern documented there. `mutateAniFstrainInput` overrides `dt`, `dt0`, `dt_start`, `Nt`, `writer_subfolder` — all fields that get copied from `dt` at parse time.
- **[Trade-off] CI plot uses log y-axis but δ starts at 1 (=10^0).** Looks fine, just worth noting that the bottom of the y-range is `10^0 = 1` not `0`. Linear y-axis would also work but log shows the saturation knee more clearly when `ani_fac_max ~ 10`.
- **[Trade-off] Research note in `SETS/` lives outside CI.** A researcher cloning the repo to run only `make run-tests` won't see it unless they browse `SETS/`. → Acceptable: the note is for someone investigating the anisotropy code, who would naturally land in `SETS/` and find the directory by name.

## Migration Plan

1. Create `TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt` (single base fixture, simple-shear default with `pure_shear_ALE = 1`, `periodic_x = 1`, `reseed_markers = 0`).
2. Add `mutateAniFstrainInput` callback + state struct in `AnisotropyBenchmarkTests.cpp` (file-scope, not promoted to `TestHelpers.h` — the override fields are specific to this test family).
3. Add two test cases (`AniFstrainSimpleShear`, `AniFstrainSaturation`).
4. Verify both pass with calibrated thresholds (relative `1e-6` for mean, absolute `1e-9` for saturated, `1e-4` for spatial L2). Account for float32 HDF5 round-trip.
5. Add the gnuplot script (single panel, MDOODZ-form curve + measurements).
6. Render `aniso_factor_evolution.png`; commit.
7. Add new section in `TESTS/AnalyticalSolutions.md` (closed-form references for FS_AR — both pure shear and simple shear formulas — the `min()` clamp behaviour, code assertions, embedded PNG, footnote about pure-shear segfault).
8. Create `SETS/AnisoFstrainResearch/`:
   - `AniFstrainSimpleShear.txt` (scenario)
   - `hansen2012_comparison.md` (findings record — table + prose + future-work note)
9. Run full RheologyCreep + AnisotropyBenchmark suites to confirm no regressions in unrelated tests.

Rollback: revert in reverse order. Each step is a small additive change; no engine code modifications.

## Open Questions (resolved during implementation)

- ~~Should `mutateAniFstrainInput` be added to `TestHelpers.h`?~~ **Kept file-scope** in `AnisotropyBenchmarkTests.cpp` — the override fields are specific to this test family and aren't shared with other test files.
- ~~Does the existing `SetSimpleShearBCVx`/`SetSimpleShearBCVz` BC produce the expected `F(t) = [[1, γt], [0, 1]]`?~~ **Yes** — observed mean `δ` matches the analytical FS_AR(γ) to ~7 significant digits at all 5 output steps with the existing BC functions. No prescribed-velocity callback needed.
- Research note inline plot: still optional. Defaulting to a table-only comparison; the gnuplot overlay can be added later without affecting CI.

## Open Issues (out of scope for this change)

- **Pure-shear (`shear_style = 0`) + `anisotropy = 1` + `ani_fstrain = 1` segfaults at step 2 start.** Reproduced under all `pure_shear_ALE` settings (0, 1, -1) and at multiple strain rates (`Eii ∈ {0.05, 0.1, 1.0}`). Step 1 runs cleanly and writes Output00001.gzip.h5. The crash is silent (no error message, just `SIGSEGV` from MDOODZ before "Time step 00002" log appears). Bisection between simple-shear (works) and pure-shear (crashes) configurations isolates `shear_style = 0` as the trigger. Documented as a known MDOODZ bug to be addressed in a follow-up change. The simple-shear formula gates the same `FiniteStrainAspectRatio → AnisoFactorEvolv` code path as pure shear would.
