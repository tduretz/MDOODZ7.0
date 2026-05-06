## Why

MDOODZ's grain-size evolution (GSE) is currently broken: the only scenario that arms it (`PinchSwellGSE`) crashes at step 0 with NaN propagating from a diverging local Newton iteration into `mesh.d_n`, and there is no CI gate that would catch the regression. On top of that, the per-phase `gsel` switch in `.txt` files is read but never wired into the simulation — it only emits a log line — so users believe they are activating GSE when they are not.

## What Changes

- **Fix the local Newton iteration in `LocalIterationViscoElasticGrainSize`** so it converges (or fails loudly) for the homogeneous pure-shear paleowattmeter case. The current `isnan(eta_ve)` guard at [RheologyDensity.c:412-414](MDLIB/RheologyDensity.c#L412-L414) only logs and lets NaN propagate; it must either recover or `LOG_ERR`/`exit` with a diagnostic instead of letting "Cell went empty" fire downstream.
- **Remove the dead `gsel` switch** from [InputOutput.c:1077,1419,1487](MDLIB/InputOutput.c#L1487). The real switch is the per-phase flow-law index `gs` ([RheologyDensity.c:522](MDLIB/RheologyDensity.c#L522)); `gsel` is read into a local variable and used only for logging. **BREAKING**: existing `.txt` files that set `gsel = N` will warn-log "parameter not found" instead — physically a no-op, but visible in logs.
- **Add a CI analytical benchmark** for steady-state paleowattmeter grain size (`gs = 10` calcite, homogeneous pure shear, single phase, viscous-only). Closed-form reference $d_{ss} = (B_g \dot\varepsilon \tau_{II} p / A_g)^{-1/(p+1)}$ derivable independently from Austin & Evans (2002). L2-compared against `mesh.d_n` like the §3.4 stress-anisotropy benchmark in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md).
- **Restore PinchSwellGSE as a green visual test** once the iteration is fixed (no functional change to the scenario itself, just unblocking it).

## Capabilities

### New Capabilities

- `ci-grain-size-evolution`: Analytical L2 benchmark gate for the paleowattmeter steady-state grain size. Defines the homogeneous pure-shear test case, the closed-form reference, the L2 metric, and the threshold. Lives alongside `ci-rheology-tests` and `ci-anisotropy-stress-l2` in the GoogleTest CI binary.

### Modified Capabilities

(none — the iteration fix and `gsel` removal are implementation details below the spec layer)

## Impact

- **Code**: [MDLIB/RheologyDensity.c](MDLIB/RheologyDensity.c) (`LocalIterationViscoElasticGrainSize`, NaN guard); [MDLIB/InputOutput.c](MDLIB/InputOutput.c) (drop `gsel` local variable, remove its `ReadMatProps` call and the diagnostic log line).
- **Tests**: new `RheologyCreep.GrainSizeSteadyState` (or new fixture) in [TESTS/RheologyCreepTests.cpp](TESTS/RheologyCreepTests.cpp); new `TESTS/RheologyCreep/GrainSizeSteadyState.txt`; new section in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) and entry in the summary table.
- **Visual tests**: [VISUAL_TESTS/main.cpp](VISUAL_TESTS/main.cpp) `class GSE` becomes runnable again; the broken `PlotGSE()` directory-pattern issue ([VISUAL_TESTS/GSE.cpp:53-57](VISUAL_TESTS/GSE.cpp#L53)) is **out of scope** for this change.
- **Scenarios**: [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt), [cmake-exec/PinchSwellGSE/PinchSwellGSE.txt](cmake-exec/PinchSwellGSE/PinchSwellGSE.txt), [VISUAL_TESTS/PinchSwellGSE.txt](VISUAL_TESTS/PinchSwellGSE.txt) — these stop crashing, no edits required.
- **Docs**: [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) gains a new "Section 8: Grain Size Evolution" between Popov2025 and SolCx, plus a row in the summary table.
- **Dependencies**: none.
- **Risk**: the iteration fix changes numerical results for any scenario with `gs ≠ 0`. Today only the broken PinchSwellGSE and the disabled `StrainLocalization_SH_GSE` exercise this code path, so impact is bounded.
