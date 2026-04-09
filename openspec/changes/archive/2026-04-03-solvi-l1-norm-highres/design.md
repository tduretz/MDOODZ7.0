## Context

The SolVi benchmark (Schmid & Podladchikov 2003) tests MDOODZ's Stokes solver against an analytical solution for a viscous inclusion under pure shear. Current convergence analysis uses relative L2 norm and resolutions up to 81×81, yielding Vx order ~0.97 and P order ~0.75.

Two improvements are expected to recover the theoretical first-order convergence: (1) L1 norm `mean(|P_num - P_ana|)` gives a cleaner convergence signal for pressure; (2) resolutions below ~80 cells don't resolve the inclusion well enough to see the expected slope of ~1.0.

The existing infrastructure includes `computeL2Error()` in `TestHelpers.h` (relative L2), `SolViBenchmarkTests.cpp` (CI convergence test at 21/41/81), and `SolViMarkerComparison.cpp` (experimental tests for marker density and eta_average).

## Goals / Non-Goals

**Goals:**
- Add `computeL1Error()` to `TestHelpers.h` — absolute mean error `mean(|num - ana|)`
- Create SolVi parameter files for 101, 151, 201 resolutions
- Add `HighResL1Convergence` test measuring both L1 and L2 at 41/81/101/151/201
- Determine whether P convergence order reaches ~1.0 with L1 norm at high resolution
- Add L1 reporting to the existing CI `GridConvergence` test (informational)

**Non-Goals:**
- Changing existing L2 assertions or thresholds
- Modifying the Stokes solver or stencil
- Adding new error norms beyond L1 (e.g., L∞)
- Changing marker density or `eta_average` settings

## Decisions

### 1. L1 norm definition: absolute mean error

**Decision**: `computeL1Error(num, ana)` computes `(1/N) * sum(|num_i - ana_i|)` — the absolute mean error, not a relative norm.

**Rationale**: The recommended norm is `mean(|P_num - P_ana|)` — the absolute mean error. Using absolute mean error rather than relative L1 avoids division-by-zero issues when the analytical pressure is zero (which it is inside the inclusion). For convergence-order computation `p = log(E1/E2) / log(h1/h2)`, any consistent norm works — the ratio cancels normalisation.

**Alternative considered**: Relative L1 (`sum|diff| / sum|ana|`). Rejected because the analytical pressure field has zero-crossings and the inside-of-inclusion region is exactly zero, making relative norms ill-conditioned for pressure.

### 2. Place `computeL1Error()` alongside `computeL2Error()` in TestHelpers.h

**Decision**: Add the function as a `static` inline in `TestHelpers.h`, right after `computeL2Error()`.

**Rationale**: Same pattern as the existing L2 helper — header-only, no additional compilation units. Every test file that includes `TestHelpers.h` gets both norms.

### 3. Resolution set: 41, 81, 101, 151, 201

**Decision**: Use 5 resolutions spanning a 5× range in h. Compute convergence orders for pairs: 41→81, 81→101, 101→151, 151→201, and also 101→201 as the "well-resolved" pair.

**Rationale**: Resolutions >80 cells are needed to properly resolve the inclusion. The 101→201 pair (both >80) should give the cleanest convergence order. Including 41 and 81 as reference points shows how under-resolution affects the measured order. Odd grid sizes ensure the inclusion boundary doesn't align with cell faces.

### 4. Test placement: new test in SolViBenchmarkTests.cpp (CI suite)

**Decision**: Add `HighResL1Convergence` as a new test case in `SolViBenchmarkTests.cpp` rather than `SolViMarkerComparison.cpp`. Runtime of ~15–30s for the 201×201 case is acceptable for CI given SolVi's importance as the primary validation benchmark.

**Rationale**: The high-resolution convergence test validates a key claim (first-order P convergence) that should be verified on every commit. The existing CI suite already runs 21+41+81 (SolVi) taking ~1s; adding 101+151+201 adds ~30–60s total, which is reasonable.

**Alternative considered**: Experimental-only in `SolViMarkerComparison.cpp`. Rejected because this is a fundamental accuracy property, not an experimental investigation.

### 5. Assertions: L1 P order ≥ 0.9 on 101→201

**Decision**: Assert `order_P_L1 >= 0.9` for the 101→201 pair. Also assert `order_Vx_L1 >= 0.8`. Keep existing L2 assertions unchanged.

**Rationale**: Theory predicts ~1.0 for P at well-resolved grids. Setting threshold at 0.9 provides 10% headroom for marker noise. The Vx L1 order should be ≥1.0 theoretically; threshold at 0.8 accounts for the inclusion boundary effect.

### 6. Parameter files: minimal copies with only Nx/Nz/subfolder changed

**Decision**: Create `SolViRes101.txt`, `SolViRes151.txt`, `SolViRes201.txt` as copies of `SolViRes41.txt` with only `Nx`, `Nz`, and `writer_subfolder` changed.

**Rationale**: Same pattern as existing resolution variants. All other parameters (penalty, tolerances, markers/cell) remain identical for a clean convergence study.

## Risks / Trade-offs

**[Risk] 201×201 may take longer than estimated on slow CI machines**
→ Mitigation: If CI timeout becomes an issue, the 201 run can be gated behind a label or moved to `SolViMarkerComparison.cpp`. Monitor first CI runs.

**[Risk] L1 norm may NOT give ~1.0 P convergence order even at high resolution**
→ Mitigation: This is an experiment — if the measured order is still ~0.75, that's a valid finding worth documenting. Set the assertion conservatively (≥ 0.9) and adjust after measuring.

**[Risk] Marker-in-cell noise at the inclusion boundary may cause non-monotonic convergence at some resolution pairs**
→ Mitigation: Report all resolution pairs in the output table. Use the 101→201 pair for the assertion (both well-resolved). Individual pair fluctuations are expected.
