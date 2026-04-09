## Why

The test suite has no benchmark verifying that free-surface topography evolves quantitatively correctly over time. The existing `SinkingBlock` test only checks that the surface deflects (non-zero), not that the magnitude follows the expected analytical decay. The TopoBench exponential relaxation benchmark — a sinusoidal perturbation h₀·cos(kx) decaying as exp(−t/τᵣ) under gravity — is the standard verification test for free-surface codes (Crameri et al., 2012) and provides a transient, multi-step L2 check against a known closed-form solution.

## What Changes

- Add a new GTest `TopoBenchRelaxation` in a new `TopoBenchTests.cpp` that runs the TopoBench Case 1 scenario (sinusoidal 7 km perturbation, λ=2800 km, η=1e21 Pa·s, ρ=3300 kg/m³) for multiple time steps and compares the amplitude decay against the analytical formula h(t) = h₀·exp(−t/τᵣ) where τᵣ = 4πη/(ρgλ)·coth(2πH/λ) (finite-depth correction).
- Add a grid-convergence test `TopoBenchConvergence` at 2–3 resolutions (e.g. 31×31, 51×51, 101×101) verifying that the relaxation error decreases monotonically and convergence order ≥ 0.5.
- Add parameter files `TESTS/TopoBench/TopoBenchRelaxation.txt` (and resolution variants) derived from `SETS/TopoBenchCase1.txt`.
- Add the `.c` setup (SetSurfaceZCoord + SetPhase callbacks) directly in the test file.
- Document the analytical solution in `TESTS/AnalyticalSolutions.md`.

## Capabilities

### New Capabilities
- `ci-topobench-relaxation`: Verification of free-surface topographic relaxation against the exponential decay analytical solution, covering amplitude L2 accuracy over multiple time steps.

### Modified Capabilities
(none — the existing `ci-free-surface-tests` SinkingBlock test is unchanged)

## Impact

- `TESTS/TopoBenchTests.cpp` — new test file
- `TESTS/TopoBench/` — new parameter file directory (base + resolution variants)
- `TESTS/CMakeLists.txt` — register new test executable
- `TESTS/AnalyticalSolutions.md` — new section documenting the relaxation formula
- CI time: ~15–30s additional (multi-step simulations at 2–3 resolutions)
