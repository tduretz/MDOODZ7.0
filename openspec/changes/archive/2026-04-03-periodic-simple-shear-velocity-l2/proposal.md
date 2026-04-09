## Why

The CI test suite verifies velocity fields only for pure shear (Dirichlet BCs on all sides). There is no test for periodic simple shear — a fundamentally different BC mode that uses periodic boundaries in x and Dirichlet (top/bottom) shear in z. A bug in the periodic BC implementation or simple shear velocity setup would go undetected.

For homogeneous linear viscosity, the analytical simple shear velocity profile Vx(z) = γ̇·z is linear, so the finite-difference scheme should reproduce it to near machine precision. This makes it an ideal regression test: any deviation from ~1e-8 L2 immediately flags a broken periodic or shear BC implementation.

## What Changes

- Add a new GTest `SimpleShearVelocity` in `VelocityFieldTests.cpp` reusing the existing `VelocityField` test fixture — the same `SetPureOrSimpleShearBCVx/Vz` callbacks already dispatch to simple shear BCs when `shear_style=1` in the parameter file.
- Add a parameter file `TESTS/VelocityField/SimpleShearVelocity.txt` based on `PureShearVelocity.txt` with `shear_style=1`, `periodic_x=1`, `pure_shear_ALE=0`, homogeneous viscosity (both phases η=1, no inclusion), 21×21 grid.
- The new test computes L2 error of Vx against the analytical linear profile Vx(z) = γ̇·z and verifies Vz ≈ 0. With homogeneous viscosity the FD scheme reproduces the linear profile to machine precision, so threshold can be very tight (~1e-6).
- Document the analytical solution in `TESTS/AnalyticalSolutions.md`.

## Capabilities

### New Capabilities
(none)

### Modified Capabilities
- `ci-velocity-field-tests`: Add simple shear velocity L2 test alongside the existing pure shear test, with a new parameter file and analytical solution.

## Impact

- `TESTS/VelocityFieldTests.cpp` — new test case added
- `TESTS/VelocityField/` — new parameter file
- `TESTS/AnalyticalSolutions.md` — new section documenting the simple shear analytical solution
- CI time: negligible increase (~100ms for a single 21×21 step)
