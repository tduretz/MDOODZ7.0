## 1. Parameter Files

- [x] 1.1 ~~Create parameter files for 41×41 and 81×81~~ DROPPED: convergence test removed (see §2)

## 2. Convergence Test

- [x] 2.1 ~~Add PureShearConvergence test~~ DROPPED: L2 error is dominated by the viscosity inclusion perturbation (~0.024) which does not converge with grid refinement when compared against the homogeneous analytical solution. The SolVi benchmark already covers inclusion convergence with the proper analytical solution.

## 3. Tighten Existing Test

- [x] 3.1 Change L2 threshold in `PureShearVelocity` from 5e-2 to 3e-2 for both Vx and Vz

## 4. Build and Validate

- [x] 4.1 Rebuild and run velocity test — passes with L2 = 2.43e-2 < 3e-2 (1.23× margin)
- [x] 4.2 Measured: L2(Vx) = 2.429e-2, L2(Vz) = 2.429e-2, L1(Vx) = 4.494e-3

## 5. Documentation

- [x] 5.1 Update `TESTS/AnalyticalSolutions.md` §2.5 with tightened threshold and rationale
- [x] 5.2 Update `skill-testing-guide` with updated accuracy numbers
