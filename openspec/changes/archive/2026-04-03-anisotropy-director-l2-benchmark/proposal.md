## Why

The two anisotropy tests in `AnisotropyBenchmarkTests.cpp` are qualitative smoke tests — `DirectorEvolution` checks `meanAngle < 44°` and `StressAnisotropy` checks `EXPECT_GT(meanTauII, 0.1)`. Neither compares against an analytical solution. The director evolution under simple shear has a known closed-form solution (θ(t) = arctan(tan(θ₀) − γ̇·t)), and the existing `benchmark-anisotropy` spec already requires L2 error < 1e-2 — but that requirement is not yet implemented. This change upgrades the tests to quantitative L2 benchmarks.

## What Changes

- **Upgrade `DirectorEvolution` test** to compute the analytical director angle θ(t) = arctan(tan(θ₀) − γ̇·t) at the final time step and measure L2 error against the numerical director field. Assert L2(θ) < threshold.
- **Add dt-convergence order test** for the director evolution — run at 2–3 different time step sizes and verify the integration converges at first order (forward Euler scheme in `RheologyParticles.c`).
- **Upgrade `StressAnisotropy` test** to compare the computed stress invariant against the analytical value from the transverse-isotropy constitutive law at known angle θ, factor δ, and strain rate ε̇. Assert relative error < 5%.
- **Tighten or remove** the existing qualitative assertions (`meanAngle < 44°`, `EXPECT_GT(meanTauII, 0.1)`) once the quantitative L2/error checks are in place.

## Capabilities

### New Capabilities
- `anisotropy-director-l2`: Quantitative L2 benchmark for director evolution under simple shear — analytical ODE integration, L2 error measurement, dt-convergence order verification, and analytical stress comparison for anisotropic constitutive law

### Modified Capabilities
- `benchmark-anisotropy`: The existing spec already requires L2 error < 1e-2 but is not implemented. This change fulfills those requirements and adds dt-convergence order verification and stress-vs-angle accuracy not yet specified.

## Impact

- **`TESTS/AnisotropyBenchmarkTests.cpp`**: Major rewrite of both test cases — add analytical solution computation, L2 error measurement, convergence order check
- **`TESTS/AnisotropyBenchmark/`**: New parameter files for different dt values (convergence test) and different angles (stress test). Grid sizes may be increased (e.g. 21×21 or 31×31) for better marker statistics.
- **`TESTS/AnalyticalSolutions.md`**: Add anisotropy section documenting the director ODE, closed-form solution, stress formula, and measured L2 errors
- **`.github/skills/skill-testing-guide/SKILL.md`**: Add anisotropy benchmark to the test coverage map and analytical benchmark sections
- **`.github/skills/skill-anisotropy/SKILL.md`**: Add reference to the quantitative benchmark tests and measured accuracy
- **CI runtime**: Grid size increase is acceptable if needed for accuracy
