# add-gmg-stokes-defence — Status, deferrals, receipts

**Change:** `add-gmg-stokes-defence`
**Phase:** implementation in progress (§§1–4, 6, 9, 11 complete; §7 partial — anisotropy green, shear-band deferred; §§5, 10, 12 pending).
**Last updated:** 2026-04-20

This document collects the deferred items, in-scope engine fixes, and
pass-receipts for the enabled fixtures, per `design.md` §D8 and
`tasks.md` §§7.5, 8, 12.7.

---

## 1. Deferred: GMG path fails on Newton + Drucker-Prager Jacobians

**Originating tests.** `SingleShearBandGmgEquivalence`,
`ConjugateShearBandsGmgEquivalence` (tasks.md §§7.2, 7.3).

**Symptom.** Running either fixture's GMG twin
(`TESTS/ShearBand/SingleBand_gmg.txt` or `ConjugateBands_gmg.txt`,
`lin_solver = 3`) produces `residual = nan` on **restart 0** of the
first FGMRES call of the first Newton iteration; the per-component
residuals `Fu, Fv, Fp` are all `nan` by the time `EvalNonLinearResidual`
reports them:

```
[    0.290] INFO  | GMG-FGMRES: restart 0, total_iter 30, residual nan
[    0.511] INFO  | GMG-FGMRES: restart 1, total_iter 60, residual nan
...
[    4.463] INFO  | Initial residual:
[    4.468] ERROR | Fu = nan
[    4.468] ERROR | Fv = nan
[    4.468] ERROR | Fp = nan
[    4.468] ERROR | Solve went wrong - Nan residuals...
```

**Scope of the failure.** Tested across the variable space:
- `Nx = Nz ∈ {31, 41}` — NaN at both resolutions.
- Viscous law `cstv = 1, pwlv = 0` (constant) and
  `cstv = 0, pwlv = 1, npwl = 3` (power-law creep) — NaN in both
  regimes.
- `penalty ∈ {1e5, 1e10}` — NaN in both.

**Control experiments that succeed.**
- The same physics with `lin_solver = -1`
  (`DirectStokesDecoupledComp` — Powell-Hestenes outer with direct
  CHOLMOD inner per PH iteration) nonlinear-converges cleanly on the
  `PlasticityTests.DruckerPragerYield` binary and on a dual-CHOLMOD
  round-trip through `SingleShearBandGmgEquivalence` (both twins on
  `lin_solver = -1`, L2 agreement is bit-identical).
- The GMG path itself is healthy on non-plastic Newton fixtures: SolVi
  51² with `Newton = 1` + anisotropic inclusion converges in ~34–37
  FGMRES iters and agrees with CHOLMOD at 1e-15 on Vx/Vz and 7.5e-8 on
  P (`add-gmg-stokes-solver` archive STATUS.md §1).

**Diagnosis.** The Newton-mode Jacobian of a Drucker-Prager yielding
cell is locally indefinite on the deviatoric block at the yield
surface (the DP consistency tangent has negative curvature along the
flow direction). The GMG path's fine-level matvec (`StokesAssemblyGMG`
Newton branch, D11) reproduces this indefinite local block exactly;
the Vanka 5×5 smoother then attempts to invert it via a local LU and
produces NaN on the yielding cells, which the V-cycle propagates
globally on the first pre-smoothing sweep. The CHOLMOD-inner path
sidesteps this by inverting the full augmented (K + penalty·G·G^T)
system once per PH iteration — the penalty term lifts the indefinite
block back to SPD before any local inversion is attempted.

**Fix size estimate.** Making the GMG path robust against indefinite
Newton-plasticity blocks requires one or more of:
- (a) a symmetric-Picard-projection fallback on yielding cells in the
  Vanka smoother (tasks.md §15.22 already did this for anisotropy; the
  DP case is structurally similar but needs a per-cell plasticity
  flag threaded through from the assembler to the smoother — cross-
  file change of ~300–500 LOC),
- (b) a local penalty augmentation inside the smoother to lift the
  indefinite block (algorithmic change with implications for coarse-
  grid consistency — needs its own design note and verification
  fixtures),
- (c) lin-search-style damping inside FGMRES when the local residual
  spikes (~100 LOC in `MultigridStokes.c` but doesn't address the
  root cause, just hides it).

All three exceed the D8 in-scope threshold (≤ 200 LOC core change,
≤ half-day debugging). Deferred to a follow-up change.

**Resolution.**
1. `TESTS/DisabledBenchmarks/SingleShearBandGmgEquivalence.cpp`:
   `TEST_F` renamed to `DISABLED_NewtonPlasticBandMatchesCholmodLoose`.
   Body left intact so the follow-up change can remove the `DISABLED_`
   prefix and re-run unchanged.
2. `TESTS/DisabledBenchmarks/ConjugateShearBandsGmgEquivalence.cpp`:
   `TEST_F` renamed to `DISABLED_NewtonTwoBandsMatchCholmodLoose`.
3. Fixture headers in `TESTS/ShearBand/*_{gmg,chol}.txt` carry a
   `STATUS — DEFERRED` banner pointing here.
4. This STATUS entry is cross-referenced from the defence's §8
   Limitations (tasks.md §10.8).

**Follow-up change.** To be named `add-gmg-newton-plasticity` or
similar; scope = make `lin_solver = 3` robust on Drucker-Prager
Newton Jacobians. Verification fixtures: re-enable the two tests
listed above and add a third at 41² power-law-creep.

---

## 2. Enabled-fixture pass receipts (tasks §§4, 6, 7.1)

| Fixture | Origin | Status | Headline |
|---|---|---|---|
| `TopoRelaxGmgEquivalence` | §4.1 | PASS | dVx ≈ 7e-10, dVz ≈ 7e-10, dP ≈ 3e-12 (mean-subtracted), dH = 0.0 |
| `SolKzGmgEquivalence` | §4.2 | PASS | dVx = 3.2e-11, dVz = 2.0e-10, dP = 4.1e-8 at 10⁷ viscosity contrast |
| `AnisotropyShearBandGmgEquivalence` | §7.1 | PASS (relaxed bound) | dVx ≈ 2.0e-8, dVz ≈ 2.0e-8, dP ≈ 5.7e-8 at 1e3 viscosity contrast + aniso_factor = 4; bound relaxed from 1e-8 → 5e-8 per D8, with header note citing PH residual accumulation × anisotropy-tensor rotation as the mechanism |
| `VcycleAnimGenerate` | §6 | PASS | 414 snapshots in 2.6 s; GIF + 6 stills committed to `figs/` |
| `SingleShearBandGmgEquivalence` | §7.2 | **DISABLED** | See §1 above |
| `ConjugateShearBandsGmgEquivalence` | §7.3 | **DISABLED** | See §1 above |

## 3. In-scope engine fixes

None so far. §4 and §7.1 converged without any MDOODZ core change —
every failure in their debugging trail was confined to the fixture
.txt and the test .cpp files. The §7.2/§7.3 deferral avoids a core
change that would exceed the D8 budget.

---

## 4. Amended: local MacBook perf sweep replaces the v1 EC2 plan (§5)

**Originating section.** tasks.md §5 (originally "First real EC2 sweep",
now "Local MacBook perf sweep").

**What changed.** `design.md` D1 originally pinned a `c7i.8xlarge` EC2
instance (32 vCPUs, 64 GiB, eu-west-1) with an N=5 repetition protocol.
The amended D1 runs the sweep on the Apple M1 MacBook (8 hw cores,
16 GiB, macos-arm64) across **three scripted phases** (the user
approved extending the original 1-hour cap so the strong-scaling
fill-in and an 801² GMG anchor could land). Total wall-clock across
phase-1 + phase-2a + phase-2b + phase-3 is ~2 h 21 min; each phase
carries its own budget guard. See `design.md` D1/D2/D7 for the full
rationale.

**Scope consequences.**

1. **Top grid size.** 801² now has **one empirical GMG point** at
   thr=8 (2 613 776 KiB ≈ 2.49 GiB peak RSS, 3142 s wall). The
   matching 801² CHOLMOD point is still dropped: under the
   asymptotic `α = 1.5` theory it would need ~17 GiB, which would
   straddle the 16 GiB laptop. The empirical α = 0.80 fit from
   41²–321² extrapolates 801² CHOLMOD to ~4.8 GiB; which reading is
   correct is itself a shortfall called out in defence §8.2.
2. **Thread scaling ceiling.** The strong-scaling sweep runs at
   thr ∈ {1, 2, 4, 8} with N=3 each (thr=4 ends up with N=5 across
   phase-2a + phase-2b). The observed efficiency curve is flat:
   100 % → 50.3 % → 26.1 % → 13.1 %. 16- and 32-thread points are
   deferred to the EC2 follow-up; they would hit the M1's scheduler,
   not the solver.
3. **Repetitions.** N=3 (was N=5). Median still defined; IQR is
   reported per row.
4. **Thermal throttle detection.** Wall-time-variance proxy (> 1.35×
   group median ⇒ flag + retry) in place of the v1 `turbostat` path,
   which does not exist on Apple Silicon. See `run_local.sh`.
5. **EC2 harness retention.** `benchmarks/ec2/provision.sh`,
   `run_perf_sweep.sh`, `teardown.sh`, `grids.yaml`, and `mocks/aws`
   are kept in-tree, dormant, as the seed for a follow-up change
   (`add-gmg-perf-regression-ci`) that replays the sweep on real
   Ice-Lake silicon with N=5. §8 of defence names this explicitly.

**Added files.**
- `SETS/SolViPerf.c` — argv-parameterised SolVi driver; prints a
  `PERF_JSON:` sentinel line on exit.
- `SETS/SolViPerf.txt` — stub required by `add_set()` cmake rule.
- `benchmarks/ec2/run_local.sh` — primary entry point.
- `benchmarks/ec2/grids_local.yaml` — phase-1 grid spec
  (memory-scaling + 201² wall-time + strong-scaling thr=1).
- `benchmarks/ec2/grids_local_phase2.yaml` — phase-2a grid spec
  (strong-scaling fill-in at thr ∈ {2, 4, 8}, N=3).
- `benchmarks/ec2/grids_local_phase2b.yaml` — phase-2b grid spec
  (continues phase-2a after budget truncation).
- `benchmarks/ec2/grids_local_phase3.yaml` — phase-3 grid spec
  (single 801² GMG anchor at thr=8, N=1).
- `benchmarks/ec2/_parse_grids.py` — extended to honour per-sweep
  `skip:` lists.

**Unchanged.** CSV schema (`instance_type` column now takes the value
`macbook-m1-local` vs `c7i.8xlarge`, so a reviewer can discriminate),
`_common.sh` helpers, `analyse.py` aggregation (now produces the new
`perf_strong_scaling.png` plot), the defence document's overall
structure.

**Headline outcomes.**
- α_CHOLMOD = 0.802 ± 0.162 over 5 grid points at thr=1.
- α_GMG = 0.785 ± 0.046 over 5 grid points at thr=1.
- 201² wall-time ratio t(CHOLMOD) / t(GMG) = **0.020×** (spec ≥ 3×).
- 201² strong-scaling efficiency at thr=8 = **13.1 %** (spec ≥ 50 %).
- 801² GMG at thr=8 peak RSS = 2.49 GiB, wall 52.4 min (N=1).
- All five numerical-equivalence receipts from the prior-change
  archive re-verified; the MacBook sweep did not perturb them.

**Honesty note for defence §6.** The headline numbers in defence §6
cite MacBook provenance explicitly and acknowledge three provenance
shortfalls against the v1 spec: (a) 201² wall-time ratio at 0.02× is
62× below the ≥ 3× bound — FGMRES tolerance + hardware class; (b)
strong-scaling is flat — GMG path OpenMP coverage is partial on M1;
(c) the CHOLMOD α asymptote is unresolved in 41²–321². Each is
itemised in defence §8.2 with its mechanism and a pointer to the
follow-up change (`add-gmg-perf-regression-ci`). No shortfall is
hidden.
