## Why

The existing `StressAnisotropy` test validates anisotropic stress using the **mean** of each stress component, which can mask spatial errors that cancel out. Upgrading to a proper spatial L2 norm against the closed-form analytical solution (rotate ε̇ into anisotropy frame, apply δ-scaled viscosity, rotate back) gives a rigorous pointwise accuracy metric. This is the 5th analytical benchmark, extending anisotropy test coverage from director evolution and mean-stress to full-field stress verification.

## What Changes

- Add a new `StressAnisotropyL2` test to the existing `AnisotropyBenchmarkTests.cpp` that computes spatial L2 error norms for sxxd, szzd, and sxz against the closed-form analytical values
- Reuse the existing `StressAngle.txt` parameter file (11×11, η=1, θ=30°, δ=6, pure shear)
- Reuse the existing `analyticalStressSxxd`, `analyticalStressSxz`, `analyticalStressTII` helper functions
- Add the test to the existing `AnisotropyBenchmark` fixture — no new executable needed
- Document the L2 stress benchmark in `TESTS/AnalyticalSolutions.md` §3.4

## Capabilities

### New Capabilities
- `ci-anisotropy-stress-l2`: Spatial L2 error norm for anisotropic stress components (sxxd, szzd, sxz) under homogeneous pure shear with a fixed director angle, validated against the closed-form rotation formula

### Modified Capabilities
- `benchmark-anisotropy`: Add the stress L2 test case to the anisotropy benchmark capability

## Impact

- **Code:** `TESTS/AnisotropyBenchmarkTests.cpp` — one new TEST_F added
- **Docs:** `TESTS/AnalyticalSolutions.md` — new §3.4 section
- **Skill:** `skill-testing-guide` — add entry to test suite table
- **Build:** No CMakeLists changes needed (same executable)
- **CI time:** Negligible — reuses existing 11×11 single-step simulation (already runs in <100ms)
