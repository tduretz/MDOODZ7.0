## Context

The existing `VelocityField.PureShearVelocity` test runs a single 21×21 pure shear simulation with a 10:1 viscosity inclusion (radius 0.05) and asserts Vx/Vz L2 error < 5e-2 against the analytical solution Vx = −ε̇·x, Vz = +ε̇·z. The measured error is ~0.024 (2.1× margin). There is no convergence study — the test cannot distinguish genuine regressions from expected coarse-grid error.

The SolVi benchmark already demonstrates the pattern: multi-resolution runs (21→41→81→101→151→201) with pair-wise convergence order assertions. The velocity field test should follow the same approach.

## Goals / Non-Goals

**Goals:**
- Add a multi-resolution convergence test for pure shear velocity (2–3 resolutions)
- Measure and assert convergence order for Vx and Vz
- Tighten the single-resolution L2 threshold based on measured values
- Document measured L2 values and orders in AnalyticalSolutions.md

**Non-Goals:**
- Changing the Stokes solver or boundary condition implementation
- Adding simple shear or other deformation modes (separate change)
- Testing without the viscosity inclusion (the inclusion is the whole point — homogeneous gives machine precision)

## Decisions

### Resolution sequence: 21×21, 41×41, 81×81

Three resolutions give two convergence pairs, sufficient to confirm the order. The grids are all odd-sized (matching existing convention) and double roughly in resolution. 81×81 is the finest — small enough to run in <1 second.

**Alternatives considered:**
- 4+ resolutions (like SolVi's 6): Overkill for a simple inclusion problem; adds CI time without much benefit.
- Only 2 resolutions: One pair gives no cross-check; three is the minimum for confidence.

### Convergence order threshold: ≥ 0.8

The viscosity inclusion creates a jump that limits convergence to ~first order at the interface. Expecting O(h¹)–O(h²) depending on how much the inclusion affects the global L2. A threshold of 0.8 is conservative (same as SolVi and anisotropy benchmarks).

### Tighten existing threshold from 5e-2 to 3e-2

The measured L2 at 21×21 is ~0.024. Setting 3e-2 gives 1.25× margin — tight but safe since this is a deterministic problem (no marker noise affects velocity on a single step). If measurements show more headroom at finer grids, we can adjust.

**Alternative considered:** Keep 5e-2 and only add convergence test. But the current 2.1× margin is too loose to catch subtle solver regressions.

### Test structure: separate convergence test, keep existing test

Add `VelocityField.PureShearConvergence` alongside the existing `PureShearVelocity`. The existing test stays as a fast sanity check (21×21 only); the convergence test runs all three resolutions.

### Parameter files: one per resolution

New files `PureShearVelocity41.txt` and `PureShearVelocity81.txt` alongside the existing `PureShearVelocity.txt`. Only `Nx`, `Nz`, and `writer_subfolder` differ.

## Risks / Trade-offs

- **[CI time increase]** → Two additional Stokes solves (41×41, 81×81). Expected ~0.5–2s total. Acceptable given the test suite runs in ~10s overall.
- **[Inclusion size vs grid]** → The r=0.05 inclusion on a 21×21 grid is poorly resolved (1–2 cells). This limits convergence order. If measured order is too low, we may need to increase the inclusion radius. → Mitigation: measure first, adjust if needed.
- **[Threshold too tight]** → 3e-2 with 1.25× margin could be flaky on different platforms. → Mitigation: run the test, check actual margin, adjust if <1.5×.
