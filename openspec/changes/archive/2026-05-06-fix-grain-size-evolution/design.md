## Context

The grain-size code path is exercised only by `PinchSwellGSE` (calcite paleowattmeter, `gs = 10`) and the disabled `StrainLocalization_SH_GSE` (olivine, `gs = 40`). PinchSwellGSE crashes at step 0 inside the initial viscosity solve: `LocalIterationViscoElasticGrainSize` ([RheologyDensity.c:203-424](MDLIB/RheologyDensity.c#L203-L424)) Newton-diverges → `eta_ve` becomes NaN → the NaN guard at [line 412-414](MDLIB/RheologyDensity.c#L412-L414) only emits a `LOG_INFO` and lets the value through → harmonic averaging at [line 1294,1318-1322](MDLIB/RheologyDensity.c#L1294) inverts NaN → `mesh.d_n` is NaN → `LOG_ERR("Cell went empty!!! Exiting...")` and `exit(345)`. There is no CI gate — the only existing GSE-adjacent test (`RheologyCreep.DiffusionCreepGrainSize`) runs with `gs = 0` and never exercises the iteration.

The paleowattmeter formulation is **steady-state per cell** — there is no `dd/dt` integration, the local solver returns the equilibrium grain size for the current stress state directly. So the analytical reference is the steady-state formula itself, applied to a homogeneous box.

## Goals / Non-Goals

**Goals:**

- Make `LocalIterationViscoElasticGrainSize` converge for at least the homogeneous pure-shear paleowattmeter case (calcite, `gs = 10`).
- Make NaN cause a loud, diagnosable failure (not silent propagation into `mesh.d_n` followed by an unrelated "Cell went empty" exit).
- Add a CI test that L2-compares grid `d_n` to the analytical $d_{ss}$ in a homogeneous setup, gated to ~iteration tolerance.
- Remove the dead `gsel` switch from the input parser so users do not believe they have armed GSE when they have not.
- Restore PinchSwellGSE to a non-crashing baseline (success criterion: it runs to completion at the configured `Nt`; physical correctness validated by the new analytical test, not by Pinch&Swell itself).

**Non-Goals:**

- Reformulating GSE as a true ODE (`dd/dt` integration). The current steady-state formulation is what the paper specifies and what the rest of the code assumes; adding time-evolution is a research task, not a bug fix.
- Re-enabling `StrainLocalization_SH_GSE` (olivine flow laws 40/41/42). Out of scope for this change; the analytical test on `gs = 10` validates the iteration kernel, and olivine constants can be tested in a follow-up if needed.
- Fixing `PlotGSE()` in [VISUAL_TESTS/GSE.cpp](VISUAL_TESTS/GSE.cpp#L53) — its directory pattern (`./GSE/GSE00*`) doesn't match `Output00xxx` output. Cosmetic visual-test issue, not blocking.
- Cross-code validation against ASPECT's `grain_size_pinned_state` benchmark. ASPECT uses Mulyukova & Bercovici (2018) pinned-state damage, not Austin & Evans paleowattmeter — formulations diverge enough that "agreement" would require constant-tuning rather than physics validation.

## Decisions

### D1: Fix the Newton iteration with bisection-then-Newton in log-space

**Choice:** Bracket `eta_ve` with `[eta_up, eta_lo]` (already computed at [line 630-637](MDLIB/RheologyDensity.c#L630-L637)), do 20 **log-space** (geometric-mean) bisection steps to land near the root, then 30 Newton iterations to finish. `nitmax` raised from 20 to 50 to accommodate both phases.

**Rationale:** The bracket `[eta_up, eta_lo]` for stiff scenarios like PinchSwellGSE spans 20-30 orders of magnitude (eta_up = MIN of mechanism viscosities, eta_lo = MAX — confusing naming, but eta_up is the *small* end). Linear bisection halves the bracket *width* per step, so a 30-dex bracket needs 70+ linear steps to narrow to the root scale. Log-space bisection (`eta_ve = sqrt(eta_up * eta_lo)`) halves the bracket *width in dex* per step, so 6 steps narrow 30 dex to 0.5 dex. We use 20 to leave plenty of margin; the 6-step number quoted in earlier design drafts assumed an idealized 30-dex bracket and didn't account for FP noise during fall-off.

The bracket-update logic uses monotonicity of `r(eta_ve)` directly (r decreases with eta_ve because Tii ∝ eta_ve and Eii_vis grows superlinearly with Tii). The earlier disabled code stored `r_eta_lo` once before the loop and compared signs against that stale value — that was the bug that prevented bisection convergence in the first place. The new logic uses only the sign of `r_eta_ve` against zero.

Newton handles the quadratic-convergence tail. Without bisection, Newton starts from `eta_ve = eta_up` (smallest mechanism viscosity) and overshoots into negative/NaN territory for stiff brackets — exactly the bug observed in PinchSwellGSE.

**Alternatives considered:**
- *(A) Damped Newton (line search):* would also work but adds a second tunable. Bisection-then-Newton is the established pattern in the existing commented code, and it's already partially written — minimum reinvention.
- *(B) Full bisection only:* slow (~10× more iterations to hit `tol = 1e-11`); rejected on perf grounds.
- *(C) Improve the initial guess only:* fragile — works for some scenarios but won't generalise across all the scenarios `gs ≠ 0` enables.

### D2: Make the NaN guard fatal, not silent

**Choice:** Replace the `LOG_INFO` at [line 412-414](MDLIB/RheologyDensity.c#L412-L414) with `LOG_ERR(...) ; exit(...)` carrying the cell index, phase, `Eii`, `T`, and last `eta_ve`. Do this *before* the NaN reaches the harmonic average, so the diagnostic points at the rheology not at downstream symptoms.

**Rationale:** The current behaviour ("log a generic message, let NaN through, exit much later with an unrelated error message in a different file") is the worst of both worlds — it's not silent enough to be invisible, but it's misleading enough that the actual root cause is hard to find. Making it fatal at the source point trades one failure mode for a clearer one.

**Alternatives considered:**
- *(A) Recover with a safe fallback (e.g. clamp `eta_ve = eta_up` and continue):* hides the bug, masks regressions, and lets downstream physics drift. Rejected — this is the kind of "silent fallback" that creates the next "why is this scenario unstable" investigation.
- *(B) Keep as warning, fix the iteration so it doesn't trigger:* the warning is still hit in some edge cases. Better to keep the assertion-style failure even after D1.

### D3: Remove `gsel` rather than wire it up

**Choice:** Delete the local variable, `ReadMatProps` call, and log line at [InputOutput.c:1077,1419,1487](MDLIB/InputOutput.c#L1487). Existing `.txt` files keep `gsel = N` as a benign no-op (the parser warn-logs unknown keys). Document the removal in the change's specs.

**Rationale:** The real switch is `materials.gs[phase] != 0` (the per-phase grain-size *flow-law index*). Wiring `gsel` up as a separate boolean would (a) duplicate state, (b) introduce a "what if `gsel = 0` but `gs = 10`?" ambiguity, and (c) require a deprecation cycle. Removing the dead code is the simpler, cheaper option.

**Alternatives considered:**
- *(A) Wire `gsel` to gate `gs`:* introduces redundant control plane; one of the two becomes vestigial. Rejected.
- *(B) Document `gsel` as deprecated, leave parser intact:* deprecation churn for a switch nobody depends on (because it doesn't do anything). Rejected.

### D4: Benchmark scenario — homogeneous pure shear, viscous-only, single phase

**Choice:** New `.txt` with no inclusion, `Nb_phases = 1`, `elastic = 0`, `thermal = 0`, `Nt = 1`, dislocation flow law (`pwlv`) compatible with calcite paleowattmeter (`gs = 10`). Strain rate, temperature, and grain-growth constants chosen so $d_{ss}$ falls in the physical range (~10⁻⁵–10⁻³ m).

**Rationale:** Homogeneous → uniform analytical scalar → spatial L2 collapses to iteration-tolerance scatter. Same shape as §3.4 stress-anisotropy L2 (4.4e-16 measured). Viscous-only avoids the Maxwell transient (where $\tau_{II}(t)$ would not be closed-form for power-law), and the steady-state $d_{ss}$ is the right reference. Single phase (no inclusion) means no marker-boundary noise, so the threshold can sit at machine precision rather than the 1e-2 levels seen in inclusion benchmarks (SolVi, etc.).

**Concrete analytical reference** (derivable from Austin & Evans 2002 Eq. 6 in one line):

$$\tau_{II} = 2 B_{pwl}\,\dot\varepsilon_{II}^{1/n_{pwl}}, \qquad d_{ss} = \left(\frac{B_g\,\dot\varepsilon_{II}\,\tau_{II}\,p}{A_g}\right)^{-1/(p+1)}$$

with $B_g = \lambda/(c_g\,\gamma)$, $A_g = K_g\,\exp(-Q_g/R_g T)$. Constants from [FlowLaws.c:866-875](MDLIB/FlowLaws.c#L866).

**Alternatives considered:**
- *(A) Pinch&Swell once fixed:* heterogeneous, no closed-form analytical, only "looks similar to a published figure". Useful as a visual sanity check, not as a CI gate.
- *(B) Smooth temperature variation $T(x) = T_0 + \Delta T\sin(2\pi x/L_x)$ for a grid-convergence study:* deferred to a follow-up §8.2; too much surface area for the first CI test.
- *(C) Maxwell transient + GSE:* $\tau_{II}(t)$ for power-law + elasticity isn't closed form. Rejected.

### D5b: Gnuplot visualisation — analytical curve with overlaid MDOODZ point

**Choice:** The test writes a single-row `.dat` file (`grain_size_benchmark.dat`) with `(Eii, T, tau_II, d_ss_analytical, d_n_mean_mdoodz)`. A checked-in gnuplot script `TESTS/RheologyCreep/plot_grain_size.gp` computes the analytical $d_{ss}(\tau_{II})$ curve over a stress range (e.g. $10^4$–$10^9$ Pa) using the same constants the test uses, and overlays the MDOODZ point on log-log axes. Output PNG `grain_size_benchmark.png` is embedded in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8.

**Rationale:** Matches the established pattern (§3.2 director, §7 SolCx) — test writes `.dat`, gnuplot script renders PNG, PNG checked in and referenced from the doc. The plot shape is the canonical paleowattmeter (Austin & Evans 2002 Fig. 2-style) — stress on x, grain size on y, log-log. A single MDOODZ point sitting on the analytical curve is a strong visual sanity check that complements the L2 number; reviewers can see "yes the test point lands on the line" at a glance.

**Why a single point and not a sweep:** The L2 test uses one parameter file. A multi-point sweep would need multiple `.txt` files and multiple `RunMDOODZ` invocations — too much surface area for this change. The analytical curve gives the parameter-space context for free; the MDOODZ point validates that we're on it. A multi-point sweep can be added later if richer validation is wanted.

**Alternatives considered:**
- *(A) Heatmap of `mesh.d_n`:* the field is homogeneous, so the heatmap is one colour. Not informative.
- *(B) Convergence-rate plot like §7 SolCx:* requires a grid-refinement study; deferred to follow-up §8.2.
- *(C) Skip the plot entirely:* breaks parity with every other §3+ benchmark; reviewers have to read raw numbers from the assertion. Rejected on consistency grounds.

### D5: L2 threshold = 1e-4 with 1% mean sanity bound

**Choice:**
```cpp
EXPECT_LT(L2_d, 1e-4);                      // ~iteration noise + interp
EXPECT_NEAR(meanD, d_ss_SI, d_ss_SI*0.01);  // 1% absolute
EXPECT_GT(minD, 0.0);                       // no NaN propagation
EXPECT_LT(maxD/minD - 1.0, 1e-3);           // homogeneity
```

**Rationale:** The Newton tolerance is `tol = 1e-11` ([line 206](MDLIB/RheologyDensity.c#L206)). After P2G/G2P interpolation in the homogeneous case the residual on `mesh.d_n` is dominated by FP cancellation, ~1e-9 to 1e-7 in practice. The 1e-4 threshold gives 3-5 orders of margin (matching the calibration methodology in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md#L468) — "2-5× above measured"). The 1% absolute bound catches scaling/units bugs that the relative L2 wouldn't (e.g. constants off by a factor of 10 give same relative scatter but wrong mean).

**Calibration:** measure once on a clean machine, set thresholds to 5× the observed value if observed values turn out larger than this initial estimate; 1e-4/1% are the *upper bounds* below which we expect to land.

## Risks / Trade-offs

- **[Risk] D1 fix works for `gs = 10` but not olivine `gs = 40`/`41`/`42`** → Mitigation: validation harness uses `gs = 10` only; olivine code path remains exercised solely by the disabled `StrainLocalization_SH_GSE` and is not regressed by this change. If olivine breaks under D1, that's a separate scenario-level issue and gets its own follow-up.
- **[Risk] D2's fatal-on-NaN may trip in scenarios that today silently produce NaN and continue** → Mitigation: those scenarios are already broken (NaN → "Cell went empty" → exit). Making the failure earlier and clearer is strictly better. If a real production scenario that "mostly works" surfaces, we can add a configurable `--allow-nan-gse` CLI flag in a follow-up; not needed for this change.
- **[Risk] Removing `gsel` (D3) silently changes input file semantics** → Mitigation: `gsel` was always a no-op; "silently changing" a no-op into a different no-op (warn-log instead of info-log) is benign. Called out as **BREAKING** in the proposal for visibility.
- **[Risk] D4 benchmark passes despite a real bug** because homogeneous problems don't exercise spatial coupling → Mitigation: homogeneous catches scaling, NaN, scope/disable bugs (the four most common GSE failure modes in this codebase). A heterogeneous follow-up (D4 alternative B) tests spatial coupling and is queued for §8.2 once §8.1 is green.
- **[Trade-off] D5 thresholds are empirical** → Standard practice in this codebase ([TESTS/AnalyticalSolutions.md §"Threshold Calibration Methodology"](TESTS/AnalyticalSolutions.md#L468)). Documented at the test site so a future tightening is straightforward.

## Migration Plan

1. Land D3 (`gsel` removal) first — it's a pure cleanup, no test dependency.
2. Land D1 (Newton fix) + D2 (fatal NaN guard) together — they share the same iteration code and tests.
3. Land D4 + D5 (CI test) last — the test depends on (1)+(2) for green.
4. Spot-check PinchSwellGSE runs to completion at the configured `Nt` (success criterion only — no assertion added; analytical test owns the correctness gate).

Rollback: revert in reverse order. The CI test going red after a partial revert is the desired safety net.

## Open Questions

- Which `pwlv` index pairs cleanly with `gs = 10` (calcite paleowattmeter)? `ReadDataPowerLaw` in [FlowLaws.c](MDLIB/FlowLaws.c) lists numbered flow laws; calcite candidates need to be inspected during implementation. If no calcite power-law is in the database, use a generic `cstv = 1` (constant viscosity) phase — the wattmeter formula is well-defined as long as `Eii_pwl > 0`, which we get by computing $Eii_{pwl} = \tau_{II}/(2\eta_{cst})$ from the constant viscosity. Decide during task implementation; either path satisfies the spec.
- Should the test use SI inputs or non-dim? The existing rheology tests in [TESTS/RheologyCreep/](TESTS/RheologyCreep/) use `eta=1e22, L=1e4` non-dim scaling. Match that for consistency.
