## Why

The pure shear velocity L2 test (`VelocityField.PureShearVelocity`) currently uses a threshold of 5e-2 against measured values of ~0.024, giving only a 2.1× margin. More importantly, the test runs on a single 21×21 grid with no convergence study — it cannot distinguish between a code regression and expected coarse-grid error. Adding a multi-resolution convergence test (similar to the SolVi and anisotropy benchmarks) will tighten the assertions and verify that the Stokes solver produces the correct convergence rate for free-slip pure shear with a viscosity inclusion.

## What Changes

- Add a grid-convergence test for the pure shear velocity field at 2–3 resolutions (e.g., 21×21, 41×41, 81×81)
- Tighten the existing L2 threshold based on measured convergence behavior
- Measure and assert the convergence order (should be ~1–2 depending on the inclusion)
- Optionally add a homogeneous pure shear test (no inclusion) with machine-precision threshold as a solver sanity check

## Capabilities

### New Capabilities
- `velocity-convergence-benchmark`: Grid-convergence study for pure shear velocity L2 error across multiple resolutions, with convergence order assertion

### Modified Capabilities
- `ci-velocity-field-tests`: Tighten L2 thresholds for Vx and Vz from 5e-2 to a value justified by the convergence study; add reference to convergence test

## Impact

- `TESTS/VelocityFieldTests.cpp`: New convergence test, tightened thresholds on existing test
- `TESTS/VelocityField/`: New parameter files for additional resolutions (41×41, 81×81)
- `TESTS/AnalyticalSolutions.md`: Updated §2.5 with convergence table and measured orders
- CI time: ~1–3 seconds additional (two extra Stokes solves on small grids)
