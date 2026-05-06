## Context

The just-archived `fix-grain-size-evolution` change repaired the GSE local Newton iteration in [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c) and added two CI tests:

- `RheologyCreep.GrainSizeSteadyState` — single-point regression gate at one strain rate.
- `RheologyCreep.GrainSizeSweep` — 5-point sweep over 4 decades of `bkg_strain_rate`.

Both run the calcite paleowattmeter scenario with `pwlv = 15` (Renner et al. 2002 dislocation creep) only — `linv = 0` (diffusion off). They validate the analytical wattmeter steady state $d_{ss} = (B_g\,\dot\varepsilon\,\tau_{II}\,p/A_g)^{-1/(p+1)}$ to within 0.07% across the sweep, with L2 ≈ 7e-4.

This left two pieces of MDOODZ's GSE behavior unverified by CI:

1. The **coupled `Ėii_pwl + Ėii_lin(d_ve)` Jacobian path** in `LocalIterationViscoElasticGrainSize`. When diffusion creep is active, `Eii_lin = C_lin · τII^n · d^(-m)` re-introduces `d_ve` into the residual, and the Newton step picks up the term `dfdeta += params.m_lin * Eii_lin * dddeta / d_ve` ([RheologyDensity.c, the Newton branch in the unified loop]). The original PinchSwellGSE crash happened *exactly* in this regime — `pwlv = 15` AND `linv = 15` AND `gs = 10` together. Today's L2 sweep with `linv = 0` doesn't probe this path.
2. The **integrated 2D code path**: marker → grid interpolation (P2G), spatial advection, harmonic averaging into `mesh.d_n`, regridding back to particles. A regression in any of these would not move the homogeneous-box L2 test.

Schmalholz & Duretz (2017) studied calcite necking specifically because the combined-mechanism regime is where wattmeter physics matters most: dislocation creep dominates at swells (large `d`), diffusion creep takes over at pinches (small `d`). Their setup is what `PinchSwellGSE.txt` reproduces.

## Goals / Non-Goals

**Goals:**

- Add a CI test that exercises `LocalIterationViscoElasticGrainSize` with both `pwlv` and `linv` active, asserting the same analytical wattmeter steady state. Catches Jacobian-coupling regressions that today's sweep doesn't.
- Add a CI test that runs a downsized `PinchSwellGSE` scenario in 2D, asserting a small set of "is the integrated path healthy?" invariants. Catches regressions in advection / interpolation / harmonic averaging that today's tests miss.
- Document both extensions in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8 as a "Coverage extensions" subsection, not a new section — they belong with the headline benchmark conceptually.

**Non-Goals:**

- Reproducing Schmalholz & Duretz 2017 figures quantitatively (Fig. 5 / Fig. 8 / etc.). That's research-grade validation, runs minutes, doesn't belong in fast CI.
- Adding paleopiezometer alternatives (`gs = 11` Rutter, `gs = 12` Schmid). Those are different flow-law dispatch entries; the iteration kernel they invoke is the same as `gs = 10`. Marginal coverage benefit.
- Multiple temperatures or multi-phase configurations. The headline scenario already pins those parameters.
- Changing the iteration kernel itself. This change only adds CI coverage of paths the kernel already supports.

## Decisions

### D1: Coupled-mechanism sweep with analytical fixed-point reference

**Choice:** A single fixture file `GrainSizeSweepCoupledBase.txt` (`linv = 15`, Calcite Herwegh 2003 diffusion creep) serves as the base for all 5 strain rates. Per-iteration overrides of `bkg_strain_rate` and `writer_subfolder` are injected via MDLIB's `MutateInput` callback hook ([mdoodz.h:338](MDLIB/include/mdoodz.h#L338); runs at [Main_DOODZ.c:107](MDLIB/Main_DOODZ.c#L107) after parse, before sim init). The dislocation-only `GrainSizeSweep` test similarly shares `GrainSizeSteadyState.txt` as its base. The new test `RheologyCreep.GrainSizeSweepCoupled` solves the **coupled fixed-point** for the analytical reference at each strain rate (since `Eii_pwl ≠ Eii_total` once diffusion is on) and L2-compares MDOODZ's `mesh.d_n` to that solution.

**The coupled equilibrium.** With both dislocation and diffusion creep active, MDOODZ imposes `Eii_total = bkg_strain_rate` at the boundary, and the iteration partitions it between mechanisms at the same stress:

$$\dot\varepsilon_{II}^{tot} = \dot\varepsilon_{II}^{pwl}(\tau_{II}) + \dot\varepsilon_{II}^{lin}(\tau_{II}, d_{ss})$$
$$\dot\varepsilon_{II}^{pwl} = \left(\tau_{II}/(2 B_{pwl})\right)^{n_{pwl}}, \quad \dot\varepsilon_{II}^{lin} = \left(\tau_{II}/(2 B_{lin})\right)^{n_{lin}} d_{ss}^{-m_{lin}}$$
$$d_{ss} = \left(B_g\,\dot\varepsilon_{II}^{pwl}\,\tau_{II}\,p / A_g\right)^{-1/(p+1)}$$

This is a 1D root-find on `τ_II`. The test runs ~30 lines of fixed-point/secant iteration to solve it independently of MDOODZ. **Engine code is untouched** — the test re-derives the same equations using a *different* iteration scheme (simple secant on τ_II vs MDOODZ's bisection-then-Newton on η_ve), so a self-consistent-but-wrong MDOODZ converged state would still be exposed by L2 disagreement.

**Why this matters physically.** At T=350°C, calcite Eii=1e-14, the strain rate splits roughly 62% dislocation / 38% diffusion. So `Eii_pwl ≈ 0.62 · Eii_total` and naive use of the dislocation-only formula `d_ss = (B_g · Eii_total · ...)^(-1/(p+1))` would give a wrong answer by ~10% — well above the 5e-3 L2 threshold. The dislocation-only sweep `GrainSizeSweep` is a degenerate special case useful for isolating the wattmeter formula. The coupled-sweep version reproduces the regime calcite *actually inhabits* at 350°C — which is exactly the regime Schmalholz & Duretz 2017 study.

**Threshold:** same `L2 < 5e-3` as `GrainSizeSweep`. The systematic ~0.07% MDOODZ-vs-strict-formula offset (from F_pwl / scaling round-trip prefactors) should carry over since the same prefactor chain is in play. Verify during implementation.

**Alternatives considered:**
- *(A) Solve coupled fixed-point in test (CHOSEN):* tight CI gate, paper-faithful, ~30 lines of test code. Independent solver guards against tautology.
- *(B) Drop L2, use physics invariants only* (`d > d_ss_disloc_only`, finite/positive, range): catches NaN/gross regressions but a 10% bug in the coupled Jacobian would pass. Rejected — too loose for the regime that was actually broken pre-fix.
- *(C) Cross-check coupled converged d to dislocation-only converged d at same Eii_total:* runs MDOODZ twice, asserts ratio in a sensible range. No analytical reference; physical sanity only. Rejected for the same tightness reason.
- *(D) Rederive d_ss formula assuming Eii_pwl = Eii_total (the wrong way):* would silently produce a wrong analytical reference; CI would fail with the right behavior. Was the design's original assumption — corrected here.

### D1b: Why the dislocation-only sweep stays as a separate test

Splitting `GrainSizeSweep` (dislocation-only) and `GrainSizeSweepCoupled` (coupled) into two tests rather than parameterising one makes failures attributable: if only the coupled test fails, the regression is in the diffusion-coupling Jacobian (specifically the `dfdeta += params.m_lin * Eii_lin * dddeta / d_ve` term and the implicit `d_ve = f(Tii, Eii_pwl)` substitution). If the dislocation-only also fails, the regression is in the wattmeter formula itself.

Two **base** fixture files are also required (not just one) because `linv` is a flow-law dispatch index resolved at parse time via `ReadDataLinear`: changing `linv` from 0 → 15 in `MutateInput` would not load the corresponding flow-law constants. Keeping `linv = 0` and `linv = 15` as separate fixtures is mandatory.

### D1c: `MutateInput` ownership note

`writer_subfolder` is heap-allocated by `ReadChar` ([InputOutput.c:1096](MDLIB/InputOutput.c#L1096)) and freed at the end of `RunMDOODZ` ([Main_DOODZ.c:1537](MDLIB/Main_DOODZ.c#L1537)). The `MutateInput` callback must therefore `free()` the existing pointer before replacing with a `strdup`'d override — naively assigning a string literal would leak the original AND crash `free()` on shutdown. The test's callback handles this correctly.

### D2: 2D smoke test runs a downsized PinchSwellGSE with loose invariants

**Choice:** New `PinchSwellGSESmoke.txt` is a copy of `PinchSwellGSE.txt` with two parameters reduced: `Nx = Nz = 51` (was 101), `Nt = 10` (was 100). Everything else — phases, flow laws, `bkg_strain_rate`, T, the cosine layer perturbation in `SetPhase` — stays identical. New test `RheologyCreep.PinchSwellGSESmoke` runs it, reads `Output00010.gzip.h5`, asserts:

- The simulation completed (`getStepsCount` returns ≥ 0, file readable).
- All `d_n` values are finite and strictly positive (catches NaN propagation in the integrated code path).
- `max(d) / min(d) > 5` (catches "scenario silently produced uniform field" regressions — the wattmeter MUST drive heterogeneity given the pinch-swell stress field).
- `min(d) < gs_ref` for the layer phase, i.e. grain reduction occurred where stress is highest. Looser version: `min(d) < 5e-4 m` (≈ half the layer's `gs_ref = 1e-3`).

**Rationale:** The point of this test is *integration smoke* — run the actual paper scenario through the entire MDOODZ time-stepping pipeline, not validate any specific quantity. The chosen invariants are loose enough that floating-point drift / minor non-deterministic interpolation differences won't trip them, but tight enough to catch genuine regressions:

- "Healthy d field" (positive + finite) catches NaN in the integrated path.
- "Heterogeneity emerged" catches the scenario silently producing a flat d field (would happen if, e.g., `mesh.d_n` accidentally got reset to `gs_ref` everywhere by a regression in `UpdateParticleGrainSize`).
- "Grain reduction in layer" catches the wattmeter being silently bypassed in the layer phase (would happen if, e.g., the `gs[phase] != 0` switch got broken again, or `gs_ref` became the converged answer for some reason).

51×51 × 10 steps runs in ~5 s on a typical CI machine. 101×101 × 100 takes ~67 s — too slow for fast CI.

**Alternatives considered:**
- *(A) Full PinchSwellGSE (101×101, 100 steps):* expensive, ~67 s. Doesn't add proportional value over 51×51 × 10 — the bug classes are caught at low resolution too.
- *(B) Compare against a reference HDF5 (regression-on-equality):* fragile to FP non-determinism, requires committing a binary blob. Loose-invariant approach is more robust for CI.
- *(C) Assert specific `d` values at specific cells (e.g. center cell):* tighter signal but breaks under any iteration-tolerance drift. Loose invariants survive ulp-level perturbations.

### D3: Documentation lives in the existing §8, not a new section

**Choice:** Append a new "## Coverage extensions" subsection to §8 in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md), describing both new tests, their analytical/invariant references, and their measured residuals. Two rows added to the summary table.

**Rationale:** These are extensions of the same analytical reference (the paleowattmeter), not a separate benchmark. A new top-level §9 would suggest they validate something different. They don't — they validate the same physics on a different code path and a different scenario.

## Risks / Trade-offs

- **[Risk]** The coupled sweep produces a different L2 than the dislocation-only one because the iteration's converged value drifts. → Mitigation: measure during implementation. If the offset moves materially (say from 0.07% to >0.5%), document the discrepancy in `AnalyticalSolutions.md` rather than masking it. The Newton iteration is deterministic so the new measurement is repeatable.
- **[Risk]** The 2D smoke test's "heterogeneity > 5" threshold is too loose and catches nothing. → Mitigation: measure the actual ratio at 51×51 × 10 steps during implementation. Tighten to ~75% of measured if observed value is much larger (e.g. measured 50× → threshold of ~37×). Keep at "> 5" if measured is in the 5-10× range.
- **[Risk]** Build / runtime increase pushes CI past acceptable limits. → Mitigation: ~5 s for the smoke test + ~0.3 s for the coupled sweep. Both are well within the budget — the existing `RheologyCreepTests` binary runs in ~1 s today. Total post-change: ~7 s.
- **[Trade-off]** Smoke test asserts on heterogeneity but not on *correctness* of the spatial distribution. A regression that produces wrong-but-heterogeneous `d` values would pass. → Acceptable. Spatial correctness is what the paper figures verify; that's research-grade, out of scope. The smoke test is a sentinel, not a validator.

## Migration Plan

1. Land all six new `.txt` files first — pure data, no test changes.
2. Add `RheologyCreep.GrainSizeSweepCoupled` to [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp) — verify all 5 strain rates pass with the same 5e-3 threshold, calibrate if needed.
3. Add `RheologyCreep.PinchSwellGSESmoke` — measure the actual `max/min` ratio and `min(d)` value, set thresholds at 75% of measured (guarding against headroom under variation).
4. Document both in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8 with measured numbers.

Rollback: revert in reverse order. Each test is self-contained.

## Open Questions

- Do `Renner+02 dislocation` and `Herwegh+03 diffusion` agree on which `tpwl` correction applies? Renner uses `tpwl = 1` (axial-compression correction), and we should check the Herwegh constants in `ReadDataLinear` — the diffusion-creep `B_lin` formula may have its own correction factor. Inspect `MDLIB/FlowLaws.c case 15` in `ReadDataLinear` during implementation; if there's a meaningfully different prefactor, the analytical formula in the test needs to match exactly what MDOODZ computes for `B_lin`.
- Is `linv = 15` actually what `PinchSwellGSE.txt` uses today? Sanity-check matches what the design assumes. (We confirmed earlier that PinchSwellGSE has `linv = 15` for both phases, so this should hold.)
