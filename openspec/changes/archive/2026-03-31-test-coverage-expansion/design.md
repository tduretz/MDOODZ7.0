## Context

MDOODZ has two test tiers: fast GTest-based CI tests in `TESTS/` and longer integration tests in `VISUAL_TESTS/`. The CI tests currently have a single test file (`NewtonIterationConvergence.cpp`) with 8 cases all using the same ShearTemplate scenario, asserting only on Newton iteration counts. The visual tests run 7 scenarios but produce only side-by-side images with no automated pass/fail. This design covers how to expand both tiers and wire them together.

### Current test architecture

- **CI tests**: GTest fixture class → `RunMDOODZ(txtFile, &setup)` → read HDF5 output → assert on scalar values. Each test case gets its own `.txt` parameter file under `TESTS/ShearTemplate/`. CTest discovers the executable. Tests run in ~30s each.
- **Visual tests**: C++17 executable → `RunMDOODZ()` for each scenario → extract fields with Eigen → write `.dat` → gnuplot → PNG. References are pre-computed `.h5` files. No assertions — human compares images. Requires Eigen3 dependency.

## Goals / Non-Goals

**Goals:**
- Cover each major rheology mechanism (linv, expv, gbsv, composite) with at least one fast CI test
- Add CI tests for plasticity, thermal solver, boundary condition types, and solver modes
- Add 3–5 new visual regression scenarios (Peierls, anisotropic VEP, thermal diffusion, phase transition)
- Add numeric pass/fail assertions to visual tests so they can run in CI
- Create a skill documenting how to write both test types
- Keep CI fast (<3 min for unit tests)

**Non-Goals:**
- Full parametric coverage of material property space (a few representative values per mechanism suffice)
- Replacing the visual test image generation (gnuplot plots stay for human review)
- Testing the Julia visualisation scripts
- Performance benchmarking or profiling
- Achieving code-path coverage metrics

## Decisions

### 1. One GTest file per capability, not one monolithic file

**Decision**: Create separate `.cpp` files for each test domain (e.g., `RheologyCreepTests.cpp`, `PlasticityTests.cpp`, `ThermalTests.cpp`, `BoundaryConditionTests.cpp`, `SolverModeTests.cpp`).

**Why**: The existing pattern puts all cases in one file with one fixture. Separate files allow independent compilation, clearer failure attribution, and easier navigation. Each file gets its own GTest fixture with appropriate setup callbacks.

**Alternative considered**: Adding all new tests to the existing `NewtonIterationConvergence.cpp`. Rejected — it would grow unwieldy and mix unrelated concerns.

### 2. Assertion strategy for CI tests: iteration count + field magnitude bounds

**Decision**: CI tests assert on:
1. **Iteration count** — same as existing pattern (`getStepsCount()` from HDF5)
2. **Field magnitude bounds** — read a scalar from HDF5 (e.g., max τ_II, max T, mean P) and assert `EXPECT_GT(value, lower)` / `EXPECT_LT(value, upper)` with physically motivated bounds
3. **Convergence residuals** — read `/Iterations/rp_abs` and `/Iterations/rx_abs`, assert final residual < tolerance

**Why**: Iteration count alone doesn't verify the solution is correct — a bug could produce the wrong answer in the right number of iterations. Field bounds catch gross errors without being brittle to floating-point changes. Residual checks verify the solver actually converged.

**Alternative considered**: Comparing full fields against reference data in CI tests. Rejected — too slow, too brittle to platform differences, and belongs in visual tests.

### 3. Each CI test uses a minimal parameter file (1–3 time steps, small grid)

**Decision**: CI test `.txt` files use the smallest grid that exercises the mechanism meaningfully (typically 21×21 to 51×51) with 1–3 time steps. `advection = 0` where possible to reduce runtime.

**Why**: CI tests must complete in seconds. The goal is to exercise the code path, not produce scientifically accurate results. Small grids + few steps keep individual tests under 5 seconds.

### 4. Visual test verification: L2 norm comparison via Eigen

**Decision**: Add a comparison function to the visual test infrastructure:

```cpp
double relativeL2Error(const MatrixXf& computed, const MatrixXf& reference) {
    return (computed - reference).norm() / reference.norm();
}
```

Each visual test scenario calls this after computing fields and asserts `relativeL2Error(eII, eII_ref) < tolerance`. Additionally, key scalar quantities are checked (e.g., peak stress, strain onset threshold). The `main()` function returns nonzero if any assertion fails.

**Why**: The Eigen matrices are already computed for both run and reference — adding the norm comparison is a few lines per test. L2 relative error is robust to minor floating-point differences across platforms. Tolerance per test: ~1e-5 for deterministic runs, relaxed to ~1e-3 if platform sensitivity is observed.

**Alternative considered**: Image-based pixel comparison. Rejected — brittle to rendering differences, font changes, gnuplot versions, and doesn't actually verify numeric correctness.

### 5. Visual tests run as optional CI job, not blocking

**Decision**: Add a separate CI workflow job `visual-tests` that runs after the build, gated by `if: github.event_name == 'push' || github.event.label == 'run-visual-tests'`. It installs Eigen3 + gnuplot, builds with `VIS=ON`, runs `./visualtests`, and checks the exit code.

**Why**: Visual tests take 10–15 minutes. Making them blocking would slow every PR. Running them on push to main + opt-in for PRs (via label) gives regression coverage without friction.

### 6. Parameter file naming convention

**Decision**: CI test parameter files follow the pattern `TESTS/<TestSuite>/<TestCaseName>.txt`, matching the existing `ShearTemplate/` convention. The test name derives the path automatically via `testing::UnitTest::GetInstance()->current_test_info()->name()`.

**Why**: Existing pattern already works. Each test case maps 1:1 to a `.txt` file, making it trivial to adjust parameters per test.

### 7. Reference HDF5 generation for new visual tests

**Decision**: New visual test references are generated by running the scenario once on a known-good build, then committing the output `.h5` as a reference file. A `make generate-refs` target automates this.

**Why**: This matches the existing pattern for ShearTemplateReference.h5, RiftingCheninReference.h5, etc. Committing binary references is acceptable given their small size (~1–5 MB each).

## Risks / Trade-offs

**[Platform-dependent floating point]** → Relaxed L2 tolerances (1e-3 to 1e-5) per test; if a test is platform-sensitive, widen tolerance rather than skip.

**[Reference file staleness]** → When code changes intentionally alter results, references must be regenerated. The `make generate-refs` target and skill documentation make this straightforward. CI failure message will suggest re-running reference generation.

**[Visual test build dependency (Eigen3)]** → Already required for existing visual tests. CI job installs it via `apt-get install libeigen3-dev`. Not needed for the fast CI tests.

**[Test maintenance burden]** → Each new test adds a `.cpp` + `.txt` pair. Mitigated by consistent patterns and the testing skill documenting the workflow.

**[Binary reference files in git]** → HDF5 files are ~1–5 MB each. Acceptable for a handful of files. If the count grows significantly, consider git-lfs.
