## Context

The existing `StressAnisotropy` test (Test 3 in `AnisotropyBenchmarkTests.cpp`) validates anisotropic stress under pure shear by comparing the **mean** of the numerical stress field against the closed-form analytical solution. The analytical formula rotates the strain rate tensor into the anisotropy frame, applies the δ-scaled viscosity (normal stress at η, shear at η/δ), and rotates back.

The setup is homogeneous (11×11, η=1, θ=30°, δ=6, 1 step of pure shear with Exx=−1). Because the problem is spatially uniform, the analytical solution is a single constant at every grid point.

## Goals / Non-Goals

**Goals:**
- Add a `StressAnisotropyL2` test that computes the spatial L2 error norm for sxxd, szzd, and sxz against the analytical prediction at every grid point
- Detect spatially-varying errors that mean-based checks miss (e.g., boundary artifacts, interpolation bias between centers and vertices)
- Reuse the existing simulation run, analytical helpers, and parameter file

**Non-Goals:**
- No new parameter files or resolutions (reuse `StressAngle.txt`)
- No grid convergence test (the homogeneous case is already near machine precision)
- No changes to the existing `StressAnisotropy` mean-value test

## Decisions

### 1. Compute L2 for three separate stress components

Compute `computeL2Error(numerical, analytical)` independently for sxxd, szzd, and sxz. This gives more diagnostic granularity than a single combined τ_II L2 norm.

**Alternative considered:** Single L2 on τ_II. Rejected because τ_II mixes components and can hide sign errors where sxxd and sxz cancel.

### 2. Construct analytical vectors at full grid size

Since the analytical solution is spatially constant, build `std::vector<double>(N, value)` for each component. Centers fields (sxxd, szzd) have size (Nx-1)×(Nz-1); vertex field (sxz) has size Nx×Nz.

**Alternative considered:** Masked L2 skipping boundary-tagged cells. Unnecessary here — the homogeneous setup has no air cells, and the existing `computeL2Error` handles the full array.

### 3. Reuse the same simulation output

Run the simulation once via `RunMDOODZ` with `StressAngle.txt` (identical to existing Test 3). Read stress fields from the same HDF5 output. This avoids duplicate simulation time.

**Alternative considered:** Separate parameter file. Rejected — identical physics, no benefit.

### 4. Threshold: L2 < 1e-6 for all three components

The homogeneous case with constant viscosity and a single time step should produce L2 errors near machine precision. A threshold of 1e-6 is conservative while catching any real discretization or interpolation bug.

**Alternative considered:** 1e-4 (matching existing relErr threshold). Too loose for a homogeneous exact-solution test where we expect ~1e-10.

## Risks / Trade-offs

- **[Risk] sxz lives on vertices while sxxd/szzd live on centers** → The analytical value is the same everywhere, so grid-location differences don't matter. The L2 simply uses different array sizes.
- **[Risk] Boundary nodes may have Dirichlet-imposed values** → In the 11×11 pure-shear setup, boundary nodes are consistent with the analytical solution. If boundary artifacts appear, the L2 will catch them (which is the point).
- **[Risk] Duplicate simulation run if both `StressAnisotropy` and `StressAnisotropyL2` run** → GTest runs them sequentially in the same process. The simulation output from the first test persists on disk and the second can read it. Alternatively, the L2 test can call `RunMDOODZ` itself — the 11×11 run takes <100ms so duplication is negligible.
