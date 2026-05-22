## 1. Pre-flight: confirm Herwegh+03 diffusion-creep prefactor

- [x] 1.1 Read `case 15` in `ReadDataLinear` in [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c). Confirm the constants (n, m, Q, A, prefactors) and check whether `tlin` (analogous to `tpwl`) introduces a correction factor different from Renner's `tpwl = 1` axial-compression form
- [x] 1.2 If the diffusion prefactor handling matters for the analytical L2 (e.g. introduces an extra factor in `B_lin`), document this in the test inline comments and adjust the analytical formula to match. Otherwise note that no adjustment is needed.
- [x] 1.3 Verify [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt) sets `linv = 15` for the layer phase (sanity check the design assumption)

## 2. Coupled sweep: single base `.txt` fixture (refactored)

Originally 5 per-strain-rate `.txt` files were created (and the dislocation-only sweep had 4 of its own — 9 total redundant fixtures). After implementation, refactored to use MDLIB's `MutateInput` hook for per-iteration overrides, keeping only one base file per regime. Net 9 fewer committed fixtures.

- [x] 2.1 Create [TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt](TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt) — single coupled-regime base (`pwlv=15`, `linv=15`, `gs=10`); `bkg_strain_rate` and `writer_subfolder` are injected per-iteration via `MutateInput`
- [x] 2.2 Use existing [TESTS/RheologyCreep/GrainSizeSteadyState.txt](TESTS/RheologyCreep/GrainSizeSteadyState.txt) as the dislocation-only base for `GrainSizeSweep` — no new file needed
- [x] 2.3 Verify both base files exist; verify no per-strain-rate `GrainSizeSweep_E*.txt` or `GrainSizeSweepCoupled_E*.txt` remain in source

## 3. Coupled sweep: GTest case

- [x] 3.1 Add `TEST_F(RheologyCreep, GrainSizeSweepCoupled)` in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp). Iterate over 5 strain rates, set per-iteration override (Eii + writer_subfolder) on a static struct, run MDOODZ on the shared base `.txt`. **Replaced original "Option B" L2 reference (which would have used the dislocation-only formula) with "Option A" — coupled fixed-point solved by inline secant iteration on τ_II inside the test (~30 lines). The test's solver uses a different scheme than MDOODZ's iteration kernel, so a self-consistent-but-wrong MDOODZ result would be exposed by L2 disagreement, not masked.**
- [x] 3.2 Add `MutateInput` callback `mutateSweepInput` and static `g_sweepOverride` struct holding (Eii_SI, subfolder). Callback overrides `model.bkg_strain_rate` (non-dim'd by `scaling.E`) and `model.writer_subfolder` (with `free()` of original + `strdup()` of override — `writer_subfolder` is heap-owned by `ReadChar`, freed at sim end).
- [x] 3.3 Refactored both `GrainSizeSweep` (dislocation-only) and `GrainSizeSweepCoupled` to use the same MutateInput pattern.
- [x] 3.4 Measured L2 across the coupled sweep: 5.93–6.75e-4 (slightly tighter than dislocation-only's 7.05e-4 — same prefactor offset). Threshold of `L2 < 5e-3` carries over without change.
- [x] 3.5 All 12 RheologyCreep + PinchSwellGSE tests pass post-refactor (sweep tests run in 274 ms / 250 ms each, total suite 2.5 s).
- [x] 3.6 Deleted 9 redundant fixture files: `GrainSizeSweep_E{12,13,15,16}.txt` and `GrainSizeSweepCoupled_E{12,13,14,15,16}.txt`. From 11 fixtures down to 3 (`GrainSizeSteadyState.txt`, `GrainSizeSweepCoupledBase.txt`, `PinchSwellGSESmoke.txt`).

## 4. 2D smoke test: fixture file

- [x] 4.1 Create [TESTS/RheologyCreep/PinchSwellGSESmoke.txt](TESTS/RheologyCreep/PinchSwellGSESmoke.txt) by copying [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt) and changing: `Nx = 51` (was 101), `Nz = 51` (was 101), `Nt = 10` (was 100), `writer_step = 5` (so we get a couple of intermediate files plus the final), `writer_subfolder = PinchSwellGSESmoke`
- [x] 4.2 Verify the file is well-formed by running it manually: copy to `cmake-build-test/TESTS/RheologyCreep/`, build, and execute via the existing fixture pattern. Confirm exit 0 and `PinchSwellGSESmoke/Output00010.gzip.h5` exists

## 5. 2D smoke test: GTest case

- [x] 5.1 Add `TEST_F(RheologyCreep, PinchSwellGSESmoke)` in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp). The fixture's `SetPhase` callback in `RheologyCreep` returns 0 outside an inclusion of radius `user1`; for PinchSwellGSE we need a layered phase distribution instead. Two paths:
   - (a) Add a *new* fixture class `PinchSwellGSEFixture` in the same .cpp with `SetPhase` matching the cosine-layer pattern from [SETS/PinchSwellGSE.c](SETS/PinchSwellGSE.c)
   - (b) Use the existing `RheologyCreep` fixture and accept that `SetPhase` returns 0 everywhere → single-phase smoke test (no inclusion). This is weaker because it doesn't exercise multi-phase P2G.
   - **Decision: go with (a)**. Multi-phase P2G is part of the integrated path we're trying to gate.
- [x] 5.2 Implement the test: run MDOODZ, read `mesh.d_n` from `PinchSwellGSESmoke/Output00010.gzip.h5`, compute min, max, finite-check
- [x] 5.3 Assertions per spec scenarios:
   - `EXPECT_GT(minD, 0.0)` and `std::isfinite(minD) && std::isfinite(maxD)`
   - `EXPECT_GT(maxD / minD, 5.0)` (heterogeneity)
   - `EXPECT_LT(minD, 1e-3)` (grain reduction below layer's `gs_ref = 1e-3`)
- [x] 5.4 Run the test, measure actual `max/min` ratio and `min(d)`. Tighten thresholds to ~75% of measured if observed values are much larger (say measured ratio of 50 → tighten threshold to `> 35`, etc.)
- [x] 5.5 Verify total runtime: `time ./RheologyCreepTests --gtest_filter='*Smoke'` should be under 30 s
- [x] 5.6 Confirm all RheologyCreep tests pass: `./RheologyCreepTests` reports `[ PASSED ]` for the full suite (10 prior + 2 new = 12)

## 6. Documentation: §8 "Coverage extensions" subsection

- [x] 6.1 Append a `### Coverage extensions` subsection at the end of §8 in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) (after the visualisation block, before the next section divider). Document:
   - **Coupled sweep** (`GrainSizeSweepCoupled`): same analytical reference as headline sweep, but with `linv=15` active. State the measured L2 across the 5 strain rates and note that the test's value is in catching coupled-Jacobian regressions specifically.
   - **2D smoke** (`PinchSwellGSESmoke`): explain it's a sentinel test (not a quantitative validator), list the 4 invariants asserted, state the measured `max/min` and `min(d)` values, and explain what each invariant catches.
- [x] 6.2 Append two rows to the summary table:
   - `Grain size coupled L2 | RheologyCreepTests.cpp | (B_g·Eii·τ·p/A_g)^(-1/(p+1)) (§8) | EXPECT_LT(L2_d, 5e-3) | 5e-3`
   - `Grain size 2D smoke | RheologyCreepTests.cpp | finite/positive + heterogeneity invariants (§8) | EXPECT_GT(maxD/minD, 5) | (loose) |`

## 7. CI verification

- [x] 7.1 Run full `make run-tests` from `cmake-build-test/` (or rebuild dir matching `make run-tests`) and confirm the new tests are discovered and pass
- [x] 7.2 Local ctest dry-run: RheologyCreepTests passes (3.24 s, 12 tests including the 2 new ones). Pre-existing failure: `RotationAdvectionTests` (labelled "experimental") — verified to fail on the clean `add-performance-metrics-clean` branch without our changes, so unrelated to GSE. User to push branch and observe CI.
