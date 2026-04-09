## Context

The CI velocity test suite (`VelocityFieldTests.cpp`) currently verifies only pure shear (`shear_style=0`). The test fixture already uses `SetPureOrSimpleShearBCVx/Vz` callbacks that dispatch to simple shear BCs when `shear_style=1` in the parameter file, so no new BC code is needed.

For homogeneous viscosity under simple shear, the analytical velocity profile is exactly linear: Vx(z) = γ̇·z on domain [−0.5, 0.5], Vz = 0. The staggered FD scheme reproduces a linear function to machine precision, giving L2 ~ O(1e-8). This makes the test a strict regression guard on the periodic BC + simple shear pathway.

## Goals / Non-Goals

**Goals:**
- Add `SimpleShearVelocity` test in the existing `VelocityField` fixture
- Verify Vx against analytical Vx(z) = γ̇·z with a tight L2 threshold (~1e-6)
- Verify Vz ≈ 0 (L2 of Vz_numerical vs zero)
- Document the analytical solution in AnalyticalSolutions.md

**Non-Goals:**
- Grid-convergence study (linear solution is exact at any resolution)
- Non-linear or anisotropic simple shear tests
- Testing with a viscosity inclusion (pure shear test already covers that)

## Decisions

### 1. Reuse existing test fixture
The `VelocityField` class already registers `SetPureOrSimpleShearBCVx/Vz`. Setting `shear_style=1` in the .txt file is sufficient to activate simple shear BCs. No new callbacks or fixtures needed.

**Alternative considered:** Separate test fixture for simple shear — rejected, adds unnecessary code duplication for identical setup logic.

### 2. Homogeneous viscosity, no inclusion
Both phases use η₀ = 1 and `user1 = 0` (zero inclusion radius). This ensures the analytical linear profile is exact, so the L2 error measures only discretisation + solver noise, not a physical perturbation.

**Alternative considered:** Include a viscosity inclusion like the pure shear test — rejected, because we learned from the pure-shear-with-inclusion experience that the L2 error then measures the perturbation, not discretisation. A homogeneous setup gives a much tighter and more useful threshold.

### 3. Analytical solution: Vx(z) = bkg_strain_rate × z
On a domain [−0.5, 0.5] with `bkg_strain_rate=1`, the simple shear BCs impose Vx = −0.5 at the bottom and Vx = +0.5 at the top. The interior solution is the linear interpolation Vx(z) = z. Vz = 0 everywhere.

The Vx staggered grid has Nx × (Nz+1) entries. For each node (ix, iz): x_coord = xmin + ix·dx, z_coord = zmin + iz·dz − dz/2 (cell-face offset). The analytical Vx at that node is simply `bkg_strain_rate × z_coord`.

For Vz: the (Nx+1) × Nz grid should be uniformly zero.

### 4. L2 threshold: 1e-6
The FD scheme should reproduce a linear function exactly (truncation error is zero for linear polynomials). The residual L2 will come from solver tolerance and floating-point arithmetic. A threshold of 1e-6 provides ample margin above machine epsilon while catching any real regression.

### 5. Parameter file based on PureShearVelocity.txt
Copy `PureShearVelocity.txt` and change:
- `shear_style = 1`
- `periodic_x = 1`
- `pure_shear_ALE = 0`
- `user1 = 0` (no inclusion)
- Phase 1: `eta0 = 1` (same as matrix, homogeneous)
- `writer_subfolder = SimpleShearVelocity`

## Risks / Trade-offs

- **[Risk] Vx staggered-grid z-coordinate offset** → Must verify the exact z positions of Vx nodes. If off by half a cell, L2 will be ~O(dz) ≈ 0.05 instead of ~1e-8. Mitigation: cross-check with HDF5 coordinate arrays or test at high resolution first.
- **[Risk] Periodic x BCs may introduce extra ghost-node behaviour** → Mitigation: the same BC pathway is used by `SimpleShearAnisoHomo` scenario, which runs routinely. The test itself would catch any regression.
- **[Tradeoff] Homogeneous test doesn't exercise the solver as hard as an inclusion would** → Acceptable because the pure shear test already covers that, and the value here is testing the periodic + simple shear BC pathway at high precision.
