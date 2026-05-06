## Why

The `fix-grain-size-evolution` change just landed a CI gate (`RheologyCreep.GrainSizeSteadyState` + `GrainSizeSweep`) that catches NaN-propagation regressions in MDOODZ's paleowattmeter local iteration. But that gate covers only a *subset* of the regime studied by the paper that motivated the scenario — Schmalholz & Duretz (2017), J. Struct. Geol. 103. Two specific gaps exist:

1. The L2 sweep uses **dislocation creep only** (`linv = 0`), so the coupled `Ėii_pwl + Ėii_lin(d_ve)` Jacobian path — which is exactly the regime that originally crashed `PinchSwellGSE` — is not regression-gated.
2. The 2D code path (advection + interpolation + harmonic averaging on top of the iteration) has **no CI coverage at all**. A change that breaks the spatial coupling but leaves the local iteration intact would slip past today's tests.

Closing these gaps now is cheap (~10 s of additional CI), and it's the right time — the iteration kernel is fresh and the verification scaffolding is already in place.

## What Changes

- Add `RheologyCreep.GrainSizeSweepCoupled`: same 5-strain-rate sweep as `GrainSizeSweep`, but with `linv = 15` (Calcite Herwegh 2003 diffusion creep) **active alongside** dislocation. Asserts the coupled iteration converges to the same analytical wattmeter steady state, exercising the Jacobian terms that involve `dEii_lin/dd_ve · dd_ve/deta_ve` (the term that was disabled in the broken-Newton-only loop).
- Add `RheologyCreep.PinchSwellGSESmoke`: runs a downsized PinchSwellGSE (`Nx = Nz = 51`, `Nt = 10`) and asserts the integrated 2D path stays healthy — exit code 0, `d_n` finite/positive everywhere, heterogeneity emerged (`max(d)/min(d) > 5`), grain reduction occurred (`min(d) < gs_ref` in the layer phase).
- Add 5 new `.txt` files for the coupled sweep (`GrainSizeSweepCoupled_E{12,13,14,15,16}.txt`), 1 new `.txt` for the smoke test (`PinchSwellGSESmoke.txt`).
- Update [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8 with a "Coverage extensions" subsection documenting the two new tests, what they catch beyond the headline sweep, and their measured residuals.

## Capabilities

### New Capabilities

(none — this extends an existing capability)

### Modified Capabilities

- `ci-grain-size-evolution`: adds two new requirements (coupled-mechanism sweep + 2D smoke) on top of the existing dislocation-only sweep. No existing requirements change behavior; the new ones are strict additions.

## Impact

- **Tests**: new fixtures `RheologyCreep.GrainSizeSweepCoupled` (5 RunMDOODZ × ~80 ms each) and `RheologyCreep.PinchSwellGSESmoke` (1 RunMDOODZ × ~5 s on 51×51 × 10 steps). Combined CI runtime increase: ~6–7 s.
- **Test fixtures**: 6 new `.txt` files in [TESTS/RheologyCreep/](TESTS/RheologyCreep/). All small derivatives of existing files (one parameter change each).
- **Code under test**: no new application code. These tests gate the existing `LocalIterationViscoElasticGrainSize` (coupled-Jacobian branch) and the integrated 2D code path (advection, P2G, harmonic average).
- **Docs**: [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) §8 gains a "Coverage extensions" subsection and two new rows in the summary table.
- **Risk**: low. The coupled sweep should produce the same analytical d_ss (the wattmeter equilibrium is independent of *which* mechanisms contribute, only of total Eii). The smoke test's thresholds are loose enough that minor floating-point drift won't trip them.
- **Out of scope**: paleopiezometer alternatives (`gs = 11`/`gs = 12`); reproducing Schmalholz & Duretz Fig. 5 / Fig. 8 in CI (those are research-grade validations, too heavy for fast CI).
