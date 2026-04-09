## Why

Previous SolVi convergence analysis (41→81) measured P order ~0.75 with relative L2 norm, which seemed low. Two improvements are expected to recover the theoretical first-order convergence: (1) use the L1 norm `mean(|P_num - P_ana|)` instead of relative L2, which gives a cleaner convergence signal for pressure; (2) go to higher resolutions (>80 cells) because the inclusion is under-resolved at coarser grids, suppressing the expected slope of ~1.0.

The current test suite only goes up to 81×81 and uses relative L2 exclusively. Adding L1 norm computation and higher-resolution runs (101, 151, 201) will determine whether MDOODZ actually achieves first-order P convergence at sufficient resolution — a key claim for publication credibility.

## What Changes

- Add an `computeL1Error()` helper function in `TESTS/TestHelpers.h` that computes `mean(|f_num - f_ana|)` (absolute L1 norm)
- Create SolVi parameter files for resolutions 101, 151, and 201
- Add a new experimental test `HighResL1Convergence` in `TESTS/SolViMarkerComparison.cpp` that:
  - Runs SolVi at 41, 81, 101, 151, 201 resolutions
  - Computes both L2 and L1 norms for Vx and P at each resolution
  - Prints convergence-order tables for both norms across multiple resolution pairs
  - Asserts P convergence order ≥ 0.9 with L1 norm on the 101→201 pair (expected ~1.0)
- Update `SolViBenchmarkTests.cpp` grid-convergence test to also report L1 norms (informational, no new assertions)
- Document findings in skills and `AnalyticalSolutions.md`

## Capabilities

### New Capabilities
- `solvi-l1-highres`: L1 norm helper, high-resolution SolVi parameter files, and the `HighResL1Convergence` experimental test

### Modified Capabilities
- `skill-testing-guide`: Add L1 norm findings and high-resolution convergence table
- `skill-solvers`: Update convergence-order section with L1 results and resolution dependence

## Impact

- **Test helpers**: `TESTS/TestHelpers.h` — add `computeL1Error()` alongside existing `computeL2Error()`
- **Test parameter files**: `TESTS/SolViBenchmark/` — 3 new .txt files (101, 151, 201)
- **Experimental test**: `TESTS/SolViMarkerComparison.cpp` — new `HighResL1Convergence` test case
- **CI test**: `TESTS/SolViBenchmarkTests.cpp` — add informational L1 output (no assertion changes)
- **Documentation**: `TESTS/AnalyticalSolutions.md`, skills
- **Runtime**: 201×201 SolVi takes ~15–30s on a single core; acceptable for CI given SolVi's importance as the primary accuracy benchmark
- **No breaking changes**: existing assertions and tests unchanged
