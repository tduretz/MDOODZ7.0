## Context

Experiments 1–7 introduced two optimised code paths: a Jacobi-preconditioned PCG thermal solver (`thermal_solver=1`) and a fused particle-to-grid interpolation loop (`interp_mode=3`). Both are off by default and backward-compatible. The existing CI tests (`ThermalTests.cpp`) already run CHOLMOD and PCG variants independently, but never compare them head-to-head against an analytic solution, and there is no interp_mode field-level equivalence check. The PCG solver currently returns only an iteration count — no per-iteration residual data is available for convergence analysis.

VISUAL_TESTS already provides a proven pattern: C++ programs read HDF5 output into Eigen matrices, write `.dat` files, and invoke gnuplot via `popen()` to produce PNG figures. We follow this pattern for all visualisation output.

## Goals / Non-Goals

**Goals:**
- Demonstrate that PCG produces the same L2 error as CHOLMOD against analytic thermal solutions (steady-state geotherm + transient Gaussian diffusion).
- Export PCG per-iteration residual norms so convergence rate can be plotted (cold start step 0 vs warm start step N).
- Demonstrate that `interp_mode=3` produces field-identical results to `interp_mode=0` within floating-point tolerance.
- Produce publication-ready gnuplot figures for all three validations.

**Non-Goals:**
- Performance benchmarking (covered by Experiments 1–8 in skill-benchmarking).
- Testing other interp modes (1, 2) — only 0 vs 3 matters for the fused P2Mastah validation.
- Modifying the PCG algorithm itself (Jacobi preconditioner is fixed).
- Blankenbach steady-state validation (would require hours of runtime; existing BlankenBenchTests cover this separately).

## Decisions

### D1: PCG residual export via file, not struct modification

**Decision**: `SolveThermalPCG` writes per-iteration residual norms to a text file `pcg_residuals.csv` (iteration, absolute_norm, relative_norm) when an environment variable or parameter flag is set. The function signature does not change.

**Alternatives considered**:
- *Add `double *residual_history` parameter* — changes the call signature in `EnergyDirectSolve`, requires plumbing through multiple layers, risks breaking the API.
- *Write to HDF5 Iterations group* — adds coupling between thermal solver and I/O layer.

**Rationale**: A simple CSV file is decoupled from the solver API, trivially parsed by gnuplot, and controlled by a `.txt` parameter (`export_pcg_residuals = 1`). Zero impact on callers who don't enable it.

### D2: Gaussian diffusion analytic L2 in existing ThermalTests fixture

**Decision**: Add `GaussianDiffusionL2` and `GaussianDiffusionL2PCG` test cases to the existing `Thermal` fixture in `ThermalTests.cpp`, reusing the same `SetTemperature` callback (Gaussian blob). The analytic solution is:

$$T(r,t) = T_{bg} + \Delta T \cdot \frac{R^2}{R^2 + 4\kappa t} \cdot \exp\!\left(-\frac{r^2}{R^2 + 4\kappa t}\right)$$

where $\kappa = k/(\rho C_p)$, computed from the `.txt` material properties.

**Rationale**: Reuses existing fixture, callbacks, and parameter file patterns. The analytic solution is well-known and tests time-dependent accuracy (complementing the steady-state geotherm L2).

### D3: Interp equivalence as a new test file `InterpEquivalenceTests.cpp`

**Decision**: Create a new GTest file that runs the same scenario twice (once with `interp_mode=0`, once with `interp_mode=3`) using separate `.txt` parameter files and `writer_subfolder` values. After both runs, read T, Vx, Vz from the final HDF5 output and assert max absolute difference < 1e-10.

**Alternatives considered**:
- *Single test that patches interp_mode at runtime* — would require exposing `model.interp_mode` for runtime modification, complicating the test.
- *Add to existing ThermalTests* — interp equivalence is not thermal-specific; it validates the P2G/G2P interpolation loop which affects all fields.

**Rationale**: Two separate runs with two `.txt` files is the simplest approach and matches the existing pattern (e.g., `FiniteStrainTests` with InterpMode1/InterpMode2 variants). A dedicated file makes the test discoverable.

### D4: Visual test programs in VISUAL_TESTS/ following existing pattern

**Decision**: Three new C++ source files added to the existing `visualtests` target in `VISUAL_TESTS/CMakeLists.txt`:
1. `ThermalAccuracy.cpp` — reads CHOLMOD and PCG HDF5 outputs, computes L2 against analytic, writes `thermal_accuracy.dat`, invokes gnuplot.
2. `PCGConvergence.cpp` — reads `pcg_residuals.csv` files from step 0 and step N, writes `pcg_convergence.dat`, invokes gnuplot.
3. `InterpEquivalence.cpp` — reads two HDF5 outputs (mode 0 and mode 3), computes per-cell difference, writes `interp_diff.dat`, invokes gnuplot heatmap.

Each has a companion `.gnu` script.

**Rationale**: Follows the established pattern (ShearTemplate.cpp, ShearHeatingDuretz14.cpp). Shared build target, Eigen + HDF5pp for data extraction, gnuplot for plotting. No new dependencies.

### D5: Test scenario choice — SteadyStateGeotherm grid (51×51)

**Decision**: Use the existing 51×51 grid from `SteadyStateGeotherm.txt` for the thermal accuracy comparison and interp equivalence tests. For Gaussian diffusion L2, use the existing `GaussianDiffusion.txt` grid (also 51×51). For PCG convergence export, use `GaussianDiffusion.txt` with 5 steps (step 0 = cold start, step 4 = warm start).

**Rationale**: 51×51 is large enough for meaningful L2 errors but small enough for fast CI execution (~2-5s per test case).

## Risks / Trade-offs

- **[Risk] PCG residual CSV file left on disk after tests** → Mitigation: Write to `writer_subfolder` directory so it's cleaned up with other test output. Or write to `/tmp/` with a predictable name.
- **[Risk] Interp equivalence fails at > 1e-10 due to floating-point non-associativity in OpenMP reductions** → Mitigation: Start with a loose tolerance (1e-6), tighten after measuring actual differences. The fused loop changes summation order, so bit-exact equality is not expected.
- **[Risk] Gaussian analytic solution drifts from numerical due to marker-grid interpolation diffusion** → Mitigation: Use a generous tolerance (L2 < 1e-3) and verify both CHOLMOD and PCG get the same L2 to within 1e-6 of each other (relative comparison is what matters).
- **[Trade-off] CSV file export adds ~10 lines to ThermalSolver.c** → Acceptable: the code is guarded by a parameter check (`model.export_pcg_residuals`) and has zero overhead when disabled.

### D6: User-facing documentation as `OPTIMIZATIONS.md` at repo root

**Decision**: Create `OPTIMIZATIONS.md` at the repository root containing:
1. **What** — overview of the two optimisations (PCG thermal solver, fused P2Mastah interpolation).
2. **Performance** — summary table from Experiments 1–8 (wall-time reduction, phase breakdown, thread scaling).
3. **Correctness** — embedded validation figures (L2 error comparison, PCG convergence plot, interp difference heatmap) with brief interpretation.
4. **How to enable** — `.txt` parameter settings (`thermal_solver = 1`, `interp_mode = 3`), when each kicks in, and how they combine.
5. **When NOT to use** — small grids (< ~100×100 where overhead dominates), scenarios where CHOLMOD factorisation cache (D4 in Exp 5) already makes thermal solve near-free, single-thread runs where PCG offers no advantage.

**Alternatives considered**:
- *Add to README.md* — README is already long and project-overview focused; a dedicated file is more discoverable.
- *Put in JuliaVisualisation/ or docs/* — neither directory exists as a natural home; repo root is most visible.

**Rationale**: A single, self-contained document serves both as a PR description artefact and as ongoing user documentation. Figures are referenced from `VISUAL_TESTS/img/` (generated by the validation visual tests).
