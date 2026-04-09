## Context

The anisotropy benchmark tests in `TESTS/AnisotropyBenchmarkTests.cpp` currently have two test cases:
- **DirectorEvolution**: Runs 10 steps of simple shear with initial director at 45°. Asserts `meanAngle < 44°` — a qualitative check that the director rotated.
- **StressAnisotropy**: Runs 1 step of pure shear with δ=6, θ=30°. Asserts `meanTauII > 0.1` — just checks the solver produced non-zero stress.

The director evolution in `RheologyParticles.c` uses a forward Euler update:
```
nx += dt*(-(dudx - dvdz)*nx*nz - dvdx*nz² + dudz*nx²)*nz
nz += dt*( (dudx - dvdz)*nx*nz + dvdx*nz² - dudz*nx²)*nx
normalize(nx, nz)
```

For simple shear (dudz = γ̇, all other gradients = 0), this reduces to:
```
nx += dt * γ̇ * nx² * nz²
nz += dt * (-γ̇ * nx³ * nz)    [then normalize]
```

The analytical ODE for the director angle θ = atan2(nz, nx) under constant simple shear is:
```
dθ/dt = -γ̇ * cos²(θ)
```
with closed-form solution: **θ(t) = arctan(tan(θ₀) − γ̇·t)**.

The stress computation in `ViscosityConciseAniso` for constant viscosity with anisotropy:
1. Rotate strain rate to anisotropy frame using director angle
2. Apply anisotropic viscosity: σ_nn = 2η·ε_nn, σ_ns = 2η·ε_ns/δ  
3. Rotate stress back to lab frame

## Goals / Non-Goals

**Goals:**
- Replace qualitative smoke tests with quantitative L2 error benchmarks
- Verify director evolution against analytical θ(t) = arctan(tan(θ₀) − γ̇·t)
- Verify forward Euler dt-convergence order (~1.0) for the director integration
- Verify anisotropic stress against the closed-form transverse isotropy formula
- Keep existing test infrastructure (same fixture class, same parameter file layout)

**Non-Goals:**
- Modifying the director evolution algorithm itself (it's forward Euler; that's the scheme we're benchmarking)
- Testing anisotropy in combination with plasticity, elasticity, or multi-phase setups
- Testing spatial convergence (director evolution is performed on particles, not the grid)

## Decisions

### 1. Director analytical solution

**Decision:** Use θ(t) = arctan(tan(θ₀) − γ̇·t) for homogeneous simple shear.

Under simple shear with shear rate γ̇ (velocity gradient dudz = γ̇), the director evolves by dθ/dt = −γ̇·cos²(θ). Integrating: tan(θ) = tan(θ₀) − γ̇·t.

The current test uses `bkg_strain_rate=0.5` with `shear_style=1`. In MDOODZ, simple shear applies dudz = 2·bkg_strain_rate = 1.0. So γ̇ = 1.0.

With θ₀ = 45° (π/4), dt = 0.05, Nt = 10: θ_final = arctan(tan(45°) − 1.0·0.5) = arctan(0.5) ≈ 26.57°.

**Alternative considered:** Numerically integrating the full ODE with RK4 in the test code. Rejected — the ODE has a known closed form, and the code uses forward Euler, so matching against the exact solution measures the actual time integration error.

### 2. L2 error metric for director angle

**Decision:** Compute L2 error on the angle field: L2 = sqrt(mean((θ_num − θ_ana)²)) / |θ_ana|.

Read `nx`, `nz` from HDF5 Centers, compute θ = atan2(nz, nx) for each cell, compare against the single analytical value θ_ana (field is spatially homogeneous).

The forward Euler error scales as O(dt), so at dt = 0.05 with 10 steps, expect errors of order ~dt² · γ̇ ≈ 0.0025 radians ≈ 0.14°. Set threshold at L2(θ) < 0.01 (in radians), with ~10× margin.

**Alternative considered:** Comparing nx, nz components separately. Rejected — angular L2 is more physically meaningful and avoids issues with the ±n director ambiguity.

### 3. Dt-convergence order test

**Decision:** Run the same problem at 3 different dt values (0.1, 0.05, 0.025) for the same total time (t = 0.5), giving Nt = 5, 10, 20 steps respectively. Compute L2(θ) at each, then convergence order = log(L2₁/L2₂)/log(dt₁/dt₂).

The forward Euler scheme is first-order, so expect order ≈ 1.0. Assert order ≥ 0.8 (conservative).

This requires two additional parameter files: `DirectorEvolution_dt01.txt` (dt=0.1, Nt=5) and `DirectorEvolution_dt025.txt` (dt=0.025, Nt=20). The existing file (`dt=0.05, Nt=10`) serves as the middle resolution.

### 4. Analytical stress for StressAnisotropy

**Decision:** Compute the analytical stress using the same rotation+scaling formula as the code.

For pure shear (Exx = -ε̇, Ezz = +ε̇, Exz = 0) with constant viscosity η and anisotropy factor δ at director angle θ:

1. `tet = θ - π/2` (the code stores angle of the *normal* to foliation)
2. `lx2 = cos²(tet)`, `lxlz = cos(tet)*sin(tet)`, `nz2 = sin²(tet)` = 1 - lx2
3. Rotate: `E_rot.xx = lx2*Exx + nz2*Ezz + 2*lxlz*Exz`
4. Stress: `T_rot.xx = 2*η*E_rot.xx`, `T_rot.xz = 2*η*E_rot.xz/δ`
5. Back-rotate and compute τ_II

This is fully deterministic for given (θ, δ, η, ε̇). Assert relative error < 2% (tighter than the proposed 5% since the setup is homogeneous constant-viscosity with no iteration uncertainty).

**Alternative considered:** Deriving a simplified τ_II(θ, δ) formula. This would be error-prone and duplicate the rotation logic. Using the same rotation sequence as the code (verified against rotation matrix formulas per the comment at line 175) is cleaner and directly validates the code's output.

### 5. Parameter files and naming

**Decision:** Add parameter files in `TESTS/AnisotropyBenchmark/`:
- `DirectorEvolution_dt01.txt` — dt=0.1, Nt=5
- `DirectorEvolution_dt025.txt` — dt=0.025, Nt=20

Keep existing `DirectorEvolution.txt` (dt=0.05, Nt=10) and `StressAngle.txt` unchanged.

## Risks / Trade-offs

**[Director angle wrap-around]** → The analytical solution θ(t) = arctan(tan(θ₀) − γ̇·t) is well-defined for the parameter range used (θ₀ = 45°, γ̇·t = 0.5, giving θ_final ≈ 26.6° — no wrap-around). For larger strains, atan2 handling would be needed. **Mitigation:** Keep total shear strain ≤ 1.0 in the test.

**[Marker-to-grid interpolation noise]** → The HDF5 director field is interpolated from markers to cell centers. Individual cells may have slightly different angles due to statistical noise from the marker average. **Mitigation:** Use L2 over all cells (averaging out noise) and set threshold with margin. Increase grid to 21×21 or 31×31 if 11×11 is too noisy — CI time increase is acceptable.

**[Forward Euler order < 1.0 at large dt]** → With dt = 0.1, the forward Euler error may not be cleanly first-order because θ changes significantly per step. **Mitigation:** Use the finest pair (dt=0.05 → dt=0.025) for the order assertion, where both are well in the asymptotic regime.

**[Simple shear velocity gradient interpretation]** → In MDOODZ, `shear_style=1` with `bkg_strain_rate=0.5` gives dudz = 1.0 (full shear rate = 2 × bkg_strain_rate). This needs to be verified during implementation by checking the numerical velocity gradient against the expected value. **Mitigation:** Print the actual velocity gradient from the HDF5 output in the test for debugging.

## Open Questions

1. **What is the exact relationship between `bkg_strain_rate` and dudz for `shear_style=1`?** The design assumes dudz = 2·bkg_strain_rate = 1.0, but this should be empirically verified from the Vx field in the HDF5 output during implementation. If wrong, the analytical θ_final will need adjustment.

2. **Should the StressAnisotropy test also check multiple angles?** The proposal mentions running at θ=0° and θ=45° to verify angle dependence. This could be a separate test case or combined into the existing one with a loop. Decide during implementation based on complexity.
