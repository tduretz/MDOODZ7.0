## Why

Experiments 1–7 delivered a 59% wall-time reduction on large grids, but we have no artefact that demonstrates **correctness** of the optimised code paths (`thermal_solver=1` PCG, `interp_mode=3` fused P2Mastah). Before merging to main we need easy-to-interpret evidence that the optimisations produce the same numerical results as the baseline. This also serves as documentation for users evaluating whether to enable the new options.

## What Changes

- Add a **thermal accuracy comparison** test that runs the same problem with CHOLMOD (thermal_solver=0) and PCG (thermal_solver=1) and prints L2 errors against an analytic solution side-by-side.
- Add a **GaussianDiffusion analytic L2** test using the transient exact solution $T(r,t) = T_{bg} + \Delta T \cdot R^2/(R^2+4\kappa t) \cdot \exp(-r^2/(R^2+4\kappa t))$ to verify time-dependent accuracy for both solvers.
- Add **PCG convergence logging** — expose per-iteration residual norms so convergence can be plotted (step 0 cold start vs step N warm start).
- Add an **interp_mode equivalence** test that runs a short simulation with `interp_mode=0` and `interp_mode=3`, then compares T, Vx, Vz fields to confirm they are identical within floating-point tolerance.
- Provide **C++ data-extraction programs and gnuplot scripts** that produce publication-ready figures: PCG residual convergence curve, L2 error bar chart (analytic vs CHOLMOD vs PCG), and interp_mode difference heatmap.
- Produce a **user-facing documentation file** (`OPTIMIZATIONS.md` at repo root) summarising the new optimisations, benchmark results from Experiments 1–8, validation evidence, embedded visualisations, usage instructions (`thermal_solver`, `interp_mode`), and guidance on when **not** to enable them (e.g. small grids < 100×100).

## Capabilities

### New Capabilities
- `thermal-accuracy-comparison`: Side-by-side L2 error test (CHOLMOD vs PCG vs analytic) for steady-state geotherm and transient Gaussian diffusion.
- `pcg-convergence-export`: Per-iteration residual norm export from the PCG thermal solver, enabling convergence-rate visualisation.
- `interp-equivalence-test`: Field-level comparison (T, Vx, Vz) between interp_mode=0 and interp_mode=3, asserting differences are at machine-epsilon level.
- `validation-visualisation`: C++ data-extraction programs and gnuplot scripts producing PCG convergence plots, L2 error bar charts, and interp difference heatmaps from HDF5 test output.
- `optimization-docs`: User-facing markdown document (`OPTIMIZATIONS.md`) covering what the optimisations are, benchmark results, validation evidence, how to enable them, and when not to use them.

### Modified Capabilities
- `pcg-thermal-solve`: Add per-iteration residual norm logging to `SolveThermalPCG` (library modification accepted — write residuals to array/file for external plotting).

## Impact

- **MDLIB/ThermalSolver.c**: Extend `SolveThermalPCG` to optionally write per-iteration residual norms.
- **MDLIB/ThermalRoutines.c**: Plumb residual array through `EnergyDirectSolve` and write to HDF5 output (or separate log).
- **TESTS/ThermalTests.cpp**: New test cases for L2 comparison and Gaussian analytic solution.
- **TESTS/Thermal/**: New `.txt` parameter files for the comparison tests.
- **TESTS/InterpTests.cpp** (new or extended): Field comparison test for interp_mode equivalence.
- **VISUAL_TESTS/**: New C++ data-extraction programs and gnuplot scripts for validation figures.
- **OPTIMIZATIONS.md** (new, repo root): User-facing documentation with embedded figures.
- No API changes. No breaking changes. All new functionality is opt-in. MDLIB modification for PCG iteration logging is minimal and backward-compatible.
