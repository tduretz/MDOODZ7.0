## Why

Currently only 1 GTest unit test (8 cases) runs in CI, all exercising variations of the same ShearTemplate scenario. Major subsystems — diffusion creep, Peierls creep, grain-boundary sliding, isolated plasticity, pure thermal diffusion, phase transitions, and periodic boundary conditions — have zero test coverage. The 7 visual regression tests exist but are not integrated into CI and lack pass/fail assertions. Expanding the test suite is essential for long-term code quality.

## What Changes

- **New CI unit tests** (TESTS/): Add 8–12 new GTest cases covering untested rheology mechanisms, plasticity, thermal solver, boundary conditions, and solver modes. Each test runs a short simulation (1–5 steps) and asserts on HDF5 output values (iteration counts, field magnitudes, invariants).
- **New visual regression tests** (VISUAL_TESTS/): Add 3–5 new integration test scenarios for Peierls creep, anisotropic VEP, pure thermal diffusion, and phase transitions — each with a reference HDF5 and gnuplot comparison.
- **Test documentation skill**: Create a Copilot skill documenting how to write and run both test types, including the parameter file conventions, assertion patterns, and reference file generation workflow.
- **CI integration improvements**: Add numeric regression assertions to visual tests (L2 norm of field differences against reference HDF5, scalar quantity checks) and wire them into CI as an optional "long" test job with concrete pass/fail exit codes. Existing gnuplot image generation is kept for human review but is not the verification mechanism.

## Capabilities

### New Capabilities

- `ci-rheology-tests`: GTest unit tests isolating each creep mechanism (diffusion/linv, Peierls/expv, GBS/gbsv, composite pwl+lin) and verifying convergence behaviour
- `ci-plasticity-tests`: GTest unit tests for Drucker-Prager yield, strain softening evolution, and stress limiter enforcement
- `ci-thermal-tests`: GTest unit tests for pure thermal diffusion (analytical comparison), radiogenic heat production, and thermal boundary condition types
- `ci-boundary-condition-tests`: GTest unit tests for periodic BCs, free-slip vs no-slip verification, and outflow conditions
- `ci-solver-mode-tests`: GTest unit tests for Picard-only iteration and adaptive time stepping behaviour
- `visual-advanced-rheology`: Visual regression tests for Peierls creep and anisotropic VEP scenarios with reference HDF5 comparison
- `visual-thermal-phase`: Visual regression tests for pure thermal diffusion and quartz-coesite phase transition
- `ci-visual-integration`: Add numeric regression checks to visual tests (relative L2 norm `‖field - field_ref‖₂ / ‖field_ref‖₂ < tol` and scalar quantity assertions like peak stress/strain values) using existing Eigen matrices, convert visual test executable to return nonzero on failure, and add a CI workflow job that runs them as an optional matrix step
- `skill-testing-guide`: Copilot skill documenting how to write, run, and maintain both CI and visual tests

### Modified Capabilities

(none — existing specs are documentation skills, not test-related)

## Impact

- **TESTS/**: New `.cpp` test files + `.txt` parameter files for each CI test scenario. Additions to `TESTS/CMakeLists.txt`.
- **VISUAL_TESTS/**: New `.cpp` test classes, `.txt` configs, `.gnu` plot scripts, and reference `.h5` files. Additions to `VISUAL_TESTS/CMakeLists.txt`.
- **.github/workflows/ci.yml**: Extended with optional visual test job.
- **.github/skills/skill-testing-guide/**: New skill directory with `SKILL.md`.
- **Dependencies**: No new external dependencies. Existing GTest (fetched via CMake) and Eigen3/HDF5 (already in visual tests).
- **Build time**: CI will remain fast (unit tests <2 min). Visual test job adds ~10–15 min as optional step.
