## Context

A close read of the rheology routines (during proposal-grounding) revealed a more specific bug than the proposal speculated about:

- [`NonNewtonianViscosityGrid`](MDLIB/RheologyDensity.c#L1131) (anisotropy = 0 branch) calls [`ViscosityConcise`](MDLIB/RheologyDensity.c#L385), which routes through [`LocalIterationViscoElasticGrainSize`](MDLIB/RheologyDensity.c#L203) to evolve grain size.
- [`NonNewtonianViscosityGridAniso`](MDLIB/AnisotropyRoutines.c#L485) (anisotropy = 1 branch) calls [`ViscosityConciseAniso`](MDLIB/AnisotropyRoutines.c#L82), which has its **own inline local iteration** (lines 220-250) that **never updates `d`**: `*d = d0` is set at [line 113](MDLIB/AnisotropyRoutines.c#L113), the iteration loop uses `*d` only as a read-only term in `Eii_lin = C_lin * pow(Tii, n_lin) * pow(*d, -m_lin)`, and the function returns `*d = d0` unchanged.

**Net effect: with `anisotropy = 1`, the per-phase `gs = 10` (calcite paleowattmeter) flow law is silently a no-op.** No error, no warning, no grain size evolution. This is why the live `RiftingComprehensive.txt` scenarios that have `ani_fstrain = 1` but no `gs = 10` "work" — they avoid the unsupported combination by accident, not by design.

The dispatch in [Main_DOODZ.c:392-393](MDLIB/Main_DOODZ.c#L392) is mutually exclusive (one of the two branches runs per step), so the dual-writer concern in the proposal was wrong: under anisotropy, only `NonNewtonianViscosityGridAniso` writes `mesh->d_n`, and it writes the same value back (no-op).

The merge therefore needs to **port grain-size evolution into the anisotropy-aware viscosity path**, not "resolve write conflicts" or "investigate ordering."

## Goals / Non-Goals

**Goals:**

- `ViscosityConciseAniso` SHALL update `*d` according to the active grain-size-evolution law whenever `materials->gs[phase] != 0` (matching the existing convention used by [`ViscosityConcise`](MDLIB/RheologyDensity.c#L477)), so any GSE-active phase under `anisotropy = 1` actually evolves grain size — including but not limited to `gs = 10` (calcite paleowattmeter, the case the new CI test exercises).
- A new `SETS/PinchSwellGSEAniso.txt` scenario based on `SETS/PinchSwellGSE.txt` with anisotropy + `ani_fstrain = 1` added — runnable end-to-end.
- A CI smoke test `RheologyCreep.PinchSwellGSEAniso` (sentinel-class, mirroring `PinchSwellGSESmoke`) that gates the integrated path: simulation completes, `mesh.d_n` and `mesh.aniso_factor_n` finite/positive, layer develops both grain-size heterogeneity (max/min ratio > some threshold) and anisotropy strength (max δ > 1).
- Documentation in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) (a §8.3 "Coupled GSE + finite-strain anisotropy" subsection or §10, decision deferred).

**Non-Goals:**

- **No redesign of `ViscosityConciseAniso`'s inline iteration.** The work is to evolve `*d` inside it, NOT to factor out the iteration into a shared helper. Two iterations (this one and `LocalIterationViscoElasticGrainSize`) is acceptable duplication for now; merging them is a much bigger change with broader test surface.
- **No quantitative analytical L2 gate** for the coupled regime. There's no closed-form or published reference for GSE+anisotropy combined; the smoke test gates physical-invariant sentinels.
- **No change to the `min(FS_AR, ani_fac_max)` δ-saturation form.** Stays as-is; tracked in [SETS/AnisoFstrainResearch/hansen2012_comparison.md](SETS/AnisoFstrainResearch/hansen2012_comparison.md) for future work.
- **No retroactive enabling of GSE in existing aniso scenarios.** The live scenarios (`RiftingComprehensive*`, `CollisionPolarCartesianAniso`) stay as they are — adding GSE there changes their physics and is out of scope.
- **No new mesh state.** The fix uses the existing `mesh.d0_n` (input grain size from previous step) and writes to `mesh.d_n` (output). No new fields.

## Decisions

### D1: Investigation methodology — already done

**Choice:** Source read pinned the issue precisely (no ASan or bisection needed — the inline iteration's `*d = d0` is a one-line read). This decision exists in the design only to mark that this change does NOT need the sanitiser/gdb workflow that `fix-pure-shear-fstrain-segfault` did.

### D2: Fix approach — port `LocalIterationViscoElasticGrainSize` into `ViscosityConciseAniso`

**Choice:** Replace the inline Newton iteration at [`ViscosityConciseAniso` lines 220-250](MDLIB/AnisotropyRoutines.c#L220) with a call to `LocalIterationViscoElasticGrainSize` (the same function `ViscosityConcise` uses), gated on `materials->gs[phase] != 0`. When `gs == 0` (no GSE), keep the inline iteration unchanged — this preserves the no-GSE behaviour for the existing aniso scenarios. The `!= 0` condition matches the existing GSE-active gate at [RheologyDensity.c:477](MDLIB/RheologyDensity.c#L477) (`if ( materials->gs[phase] != 0 ) gs = 1;`); both wattmeter and piezometer flow laws (`gs ∈ {9, 10, 11, 12, 40, …}`) drive the same local iteration with parameters dispatched by [`ReadDataGSE`](MDLIB/FlowLaws.c#L844).

Concretely:

```c
if (materials->gs[phase] != 0) {
    // GSE active — use the wattmeter/piezometer-aware local iteration
    LocalIterationParams params = { /* populate with rotated-frame inputs and aniso-aware terms */ };
    LocalIterationViscoElasticGrainSize(
      (LocalIterationMutables){.eta = &eta_ve, .Tii = &Tii, .Eii_cst=Eii_cst, .Eii_lin=Eii_lin,
                                .Eii_gbs=Eii_gbs, .Eii_pwl=Eii_pwl, .Eii_exp=Eii_exp, .d1 = d},
      params);
} else {
    // Existing inline iteration (no GSE)
    for (int it = 0; it < nitmax; it++) { /* unchanged */ }
}
```

**Rationale:**

- **Reuse over rewrite**: `LocalIterationViscoElasticGrainSize` is well-tested (`RheologyCreep.GrainSizeSteadyState/Sweep/SweepCoupled` × 4 strain rate decades, plus the pinch-swell smoke). Calling it from the aniso path means the existing GSE CI gates this path too.
- **Minimum-scope fix**: only one branch added, gated on `gs != 0`. Does not touch the existing iteration or any other aniso scenario.
- **Strain-rate input is rotation-invariant**: `LocalIterationViscoElasticGrainSize` uses scalar `Eii = sqrt(I2(E))`. Since I2 is invariant under proper rotation, feeding it `Eii = sqrt(I2(E_rot))` (already computed at [line 184](MDLIB/AnisotropyRoutines.c#L184)) gives the same value as the un-rotated `Eii`. No need to back-rotate.
- **Tii consistency**: the inline iteration uses `Tii = 2*eta_ve*Eii` as a scalar. `LocalIterationViscoElasticGrainSize` uses the same convention. The shear-stress component reduction `T_rot.xz = 2·eta_ve·E_rot.xz / ani_fac` is applied AFTER the iteration (lines 257-260) and doesn't affect Tii's role in the iteration.

**Alternatives considered:**

- *(A) Inline the GSE update into the existing iteration*: would duplicate the bisection-then-Newton + grain-size update logic from `LocalIterationViscoElasticGrainSize`, harder to maintain. Rejected.
- *(B) Refactor `ViscosityConcise` and `ViscosityConciseAniso` to share the iteration body*: bigger change, expands scope to non-aniso path that's already working. Rejected for this change; possible future cleanup.
- *(C) Detect-and-error if `gs > 0` and `anisotropy = 1`*: makes the silent no-op explicit but doesn't add capability. Rejected because we want capability, not just a louder error.

### D3: Strain-rate input convention

**Choice:** Pass `Eii = sqrt(I2(E_rot))` (the rotated frame's invariant — already computed) to the GSE iteration. Pass `Tii` as scalar `2·eta_ve·Eii` (consistent with the existing aniso path's convention). Pass `*d = d0` as initial; iteration updates `*d` to `dnew`.

**Rationale:** I2 is rotation-invariant, so `Eii_rot ≡ Eii_un_rot`. The aniso-induced shear-stress reduction (`T_rot.xz / ani_fac`) is applied to the OUTPUT stress tensor, not used as an input to the iteration. This keeps the GSE iteration's flow-law inputs (Tii, n_pwl, A_pwl, etc.) consistent with the non-aniso path.

### D4: Scenario `SETS/PinchSwellGSEAniso.txt`

**Choice:** Copy [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt) verbatim, then add to the calcite layer phase:

- `anisotropy = 1` (in switches)
- `aniso_angle = 0` (director vertical → foliation horizontal = bedding-parallel for a horizontal layer; the natural fresh-sediment IC. See "Open Questions" resolution below for the convention.)
- `aniso_factor = 1.0` (will be overridden by FS_AR)
- `ani_fstrain = 1` (finite-strain anisotropy active)
- `ani_fac_max = 4` (modest saturation, consistent with `AniFstrainSaturation` test)

Keep all GSE settings (`gs = 10`, calcite paleowattmeter constants, `pwlv`, `linv`) unchanged. The matrix phase keeps `aniso_factor = 1` and `ani_fstrain = 0` (no aniso for matrix, mirroring `RiftingComprehensive.txt`'s pattern).

**Rationale:** Smallest delta from the working `PinchSwellGSE.txt`. Layer-vs-matrix asymmetry (only the layer has aniso) makes the foliation development visible — the matrix flows isotropically while the layer develops a strain-aligned foliation, which physically corresponds to mica/clay alignment in a competent calcite layer surrounded by weaker pelitic matrix.

**Open question (deferred to implementation):** does adding `anisotropy = 1` to a previously-isotropic-matrix scenario need any other parameter changes (e.g., `finite_strain = 1` is auto-enabled per [InputOutput.c:1319](MDLIB/InputOutput.c#L1319), but are there others)? Defer to first run.

### D5: CI smoke test `RheologyCreep.PinchSwellGSEAniso`

**Choice:** Mirror [`RheologyCreep.PinchSwellGSESmoke`](TESTS/RheologyCreepTests.cpp) (which is the existing 51×51 × Nt=10 sentinel for GSE-only). New downsized fixture `TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt` with the same dimensions. Assertions:

```cpp
// Existing GSE invariants
EXPECT_GT(minD, 0.0);
EXPECT_TRUE(std::isfinite(minD) && std::isfinite(maxD));
EXPECT_GT(maxD / minD, 5.0);              // grain size heterogeneity from wattmeter

// NEW anisotropy invariants
EXPECT_TRUE(std::isfinite(min_ani) && std::isfinite(max_ani));
EXPECT_GE(min_ani, 0.999);                 // δ ≥ 1 (FS_AR is)
EXPECT_GT(max_ani, 1.0);                   // anisotropy actually developed somewhere
EXPECT_LT(max_ani, ani_fac_max + 1e-9);   // saturation respected
```

**Rationale:** sentinel-class assertions because the coupled physics has no closed form. The four anisotropy assertions catch the silent-no-op bug we just diagnosed: if GSE were still silently disabled under anisotropy, `max_ani` would still grow (FS_AR depends only on F, not on rheology), but `max(d) / min(d)` would be 1 (no GSE evolution → uniform grain size). And vice versa: if anisotropy weren't running, `max_ani` would equal the input `aniso_factor` everywhere, `max_ani / min_ani` ≈ 1.

The test gates BOTH subsystems running together, not just one. CI runtime estimate: ~3 s (similar to `PinchSwellGSESmoke`).

### D6: Spec scope — NEW capability `ci-gse-anisotropy-coupled`

**Choice:** New spec at `openspec/specs/ci-gse-anisotropy-coupled/spec.md` with requirements covering:

1. The smoke test exists and passes the sentinel assertions.
2. The MDLIB fix (port GSE iteration into `ViscosityConciseAniso`) is present.
3. The scenario file `SETS/PinchSwellGSEAniso.txt` exists and runs.
4. Documentation subsection in `AnalyticalSolutions.md`.

**No modifications to existing capabilities.** `ci-grain-size-evolution` and `ci-finite-strain-anisotropy` continue to gate their respective isolated regimes; the coupled regime is its own capability.

**Rationale:** The coupled regime tests something the isolated specs don't — interaction. A separate capability keeps the gates orthogonal: regressions in GSE-only or aniso-only fail their existing tests; regressions in the coupling fail the new test. No requirement is duplicated.

## Risks / Trade-offs

- **[Risk] The GSE iteration converges to a different `eta_ve` under aniso vs non-aniso for the same input.** Even with rotation-invariant Eii, the existing inline aniso iteration may have different convergence behaviour than `LocalIterationViscoElasticGrainSize` (different bisection bounds, different starting guess, different tolerance). → **Mitigation**: when porting, run `AnisotropyBenchmark.AniFstrainSimpleShear/Saturation/PureShear` after the change to confirm aniso-only behaviour is unchanged. If `eta_vep` outputs differ, that's a regression even if both paths converge.
- **[Risk] Code duplication grows.** With this change, GSE iteration logic lives twice (one in `LocalIterationViscoElasticGrainSize`, one inline in `ViscosityConciseAniso` for non-GSE). → **Mitigation**: explicitly out of scope per D2; future cleanup change can refactor.
- **[Risk] Some assumption I haven't found.** `ViscosityConciseAniso` does many things (plasticity, peierls, gbs). The GSE port must integrate cleanly with all of them, not just dislocation+diffusion. → **Mitigation**: keep the GSE branch isolated to `gs != 0`; existing branches untouched. Also: regression suites (`RheologyCreepTests` 12 / `AnisotropyBenchmarkTests` 9) run both isolated regimes' full test surface.
- **[Risk] The new scenario doesn't actually develop visible coupling.** Layer phase might saturate `δ` immediately (if `ani_fac_max` is too small) or stay near-isotropic (if not enough strain accumulates by `Nt = 10`). → **Mitigation**: tune `ani_fac_max` and `Nt` empirically during D5 calibration. The proposal's `ani_fac_max = 4` and `Nt = 10` are starting points.
- **[Trade-off] Sentinel tests don't quantify accuracy.** The smoke test catches obvious regressions (NaN, no evolution) but won't catch subtle ones (e.g., grain size evolves at the wrong rate due to the rotated-frame Eii feed). Acceptable because the coupled regime has no analytical reference; future quantitative tests could come later if needed.

## Migration Plan

1. **Read `ViscosityConciseAniso` end-to-end** to confirm my reading. Specifically: trace what `*d` flows into in the rest of the function (besides `Eii_lin` at line 144) and confirm changing it doesn't break anything else.
2. **Add the GSE branch** to `ViscosityConciseAniso` at the iteration block (line 220 area). Gate on `materials->gs[phase] != 0`. Keep all existing inputs/outputs the same.
3. **Verify aniso-only regression**: run `AnisotropyBenchmark.AniFstrain*` (3 tests) — they have `gs = 0`, should be byte-identical post-change. If `eta_vep` outputs differ, abort and investigate.
4. **Verify GSE-only regression**: run `RheologyCreep.GrainSize*` and `PinchSwellGSESmoke` — they have `anisotropy = 0`, should be byte-identical.
5. **Run a manual smoke on `PinchSwellGSE.txt`** for ≥ 2 timesteps at standard parameters. Should be byte-identical.
6. **Create `SETS/PinchSwellGSEAniso.txt`**, hand-run it for ≥ 2 timesteps, inspect output. Confirm `mesh.d_n` and `mesh.aniso_factor_n` both evolve.
7. **Create `TESTS/RheologyCreep/PinchSwellGSEAnisoSmoke.txt`** (51×51 × Nt=10 downsized).
8. **Add `RheologyCreep.PinchSwellGSEAniso`** GTest case with the assertion block from D5.
9. **Calibrate** `ani_fac_max` and `Nt` for stable observable coupling (max_ani > 1, max/min(d) > 5).
10. **Run full regression**: `RheologyCreepTests` (13/13 expected post-add) + `AnisotropyBenchmarkTests` (9/9).
11. **Manual smoke**: `RiftingComprehensive.txt` ≥ 2 timesteps (uses aniso without GSE — must be unaffected).
12. **Document** in `TESTS/AnalyticalSolutions.md` §8.3 or §10.
13. **Sync delta spec** at archive time.

Rollback: revert the `ViscosityConciseAniso` GSE branch and the new test/scenario. Each commit small and reversible.

## Open Questions

- ~~**Does `LocalIterationViscoElasticGrainSize`'s strain-rate input take `Eii` directly, or some derived quantity?**~~ **Resolved by inspection.** Takes scalar `params.Eii` directly. The struct stores it as a single `double` ([RheologyDensity.h:20](MDLIB/RheologyDensity.h#L20)) and the iteration uses it only as an I2-invariant scalar in (i) the viscoelastic relation `Tii = 2·eta_ve·Eii`, (ii) the strain-rate residual `r = Eii − elastic·Tii/(2·eta_el) − Eii_vis`, and (iii) the wattmeter formula `d_ve = (2·Bg·Eii·Eii_pwl·eta_ve·p/Ag)^(-1/(p+1))`. **No tensor decomposition, no frame rotation, no angular dependency** anywhere in the iteration body. So D3's plan to feed `params.Eii = sqrt(I2(E_rot))` is mathematically equivalent to the lab-frame `sqrt(I2(E))` — already computed at [AnisotropyRoutines.c:184](MDLIB/AnisotropyRoutines.c#L184) — because I2 is rotation-invariant.

  One unrelated note: the struct has an `f_ani` field ([RheologyDensity.h:21](MDLIB/RheologyDensity.h#L21)) that is **declared but never read or written anywhere in the iteration**. Vestigial — possibly intended for a future aniso-aware iteration that wasn't completed. Safe to omit or set to any default when populating from `ViscosityConciseAniso`.
- ~~**What does `dnew` get used for downstream of `ViscosityConciseAniso`?**~~ **Resolved by inspection.** `dnew` is purely an output parameter required by the `ViscosityConcise(Aniso)` signature. Three call sites discard it (`NonNewtonianViscosityGrid(Aniso)` vertices and `FD_Jacobian.PhaseRheologyLoop_v1`); only the centroid passes (one in each `NonNewtonianViscosityGrid*`) actually use it, accumulating harmonically into `mesh->d_n[c0]` (`mesh->d_n += phase_perc * 1.0/dnew`, then `mesh->d_n = 1/mesh->d_n` after the per-phase loop). Once the fix lets `ViscosityConciseAniso` evolve `*d`, the closure is: `particles.d` → `mesh->d0_n` (P2G) → `dnew` per phase (rheology pass with the new GSE branch) → `mesh->d_n` (harmonic average) → `particles.d` (`UpdateParticleGrainSize` at [Main_DOODZ.c:1104](MDLIB/Main_DOODZ.c#L1104)). The vertex `dnew` and `FD_Jacobian.dnew` remain harmless discards.
- **Does the inline iteration's exit-on-failure (`exit(0)` at line 248) match `LocalIterationViscoElasticGrainSize`'s behaviour?** The latter has a bisection-then-Newton fallback added in `fix-grain-size-evolution`; the former does pure Newton with a 20-iter cap and exits on non-convergence. Under coupled regime, the bisection-then-Newton should help — but verify it doesn't fire unexpectedly during regression.
- **Is `aniso_factor_max = 4` a good default for the new scenario?** Smaller → faster saturation, more visible clamp, but may saturate before any layer features develop. Larger → slower to develop. Defer to D5 calibration.
- ~~**Does the foliation orientation (`aniso_angle = 90`) make physical sense for pinch-and-swell?**~~ **Resolved (per user clarification).** Convention: `aniso_angle` is the angle of the **director** vector, which is **perpendicular to the layering/foliation plane**. `aniso_angle = 90°` puts the director horizontal — for a horizontal layer this means the **foliation is vertical** (perpendicular to the director). For a horizontal pinch-and-swell setup with the goal of representing pre-existing layer-parallel fabric, the convention to use is `aniso_angle = 0°` (director vertical → foliation horizontal = bedding-parallel). **Update [SETS/PinchSwellGSEAniso.txt](SETS/PinchSwellGSEAniso.txt) to `aniso_angle = 0` for the layer phase** during implementation (task §5.3). The `ani_fstrain = 1` machinery then evolves the strain-ellipse anisotropy on top of this layer-parallel IC.
