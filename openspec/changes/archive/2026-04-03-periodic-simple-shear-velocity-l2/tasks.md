## 1. Parameter File

- [x] 1.1 Create `TESTS/VelocityField/SimpleShearVelocity.txt` by copying `PureShearVelocity.txt` and changing: `shear_style=1`, `periodic_x=1`, `pure_shear_ALE=0`, `user1=0` (no inclusion), Phase 1 `eta0=1` (homogeneous), `writer_subfolder=SimpleShearVelocity`

## 2. Test Implementation

- [x] 2.1 Add `TEST_F(VelocityField, SimpleShearVelocity)` in `VelocityFieldTests.cpp`: run simulation, read Vx/Vz from HDF5, compute L2(Vx) against Vx_ana(z) = γ̇·z, assert L2 < 1e-6, assert max(|Vz|) < 1e-6, check boundary Vx values ≈ ±0.5

## 3. Build and Validate

- [x] 3.1 Build and run: verify test passes and record measured L2 values — L2(Vx) = 2.97e-8, L2(Vz) = 0, both tests pass (215ms total)

## 4. Documentation

- [x] 4.1 Add §2.8 to `TESTS/AnalyticalSolutions.md` documenting the simple shear analytical solution and measured accuracy
- [x] 4.2 Update `skill-testing-guide` with the new test entry
