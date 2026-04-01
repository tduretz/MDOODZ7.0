---
name: skill-code-glossary
description: Code-to-physics glossary for MDOODZ — variable names, struct fields, function names mapped to physical meaning, equations, SI units, grid naming conventions, creep mechanism abbreviations, and solver process stages.
---

# Code-to-Physics Glossary

## Stress and Strain Rate Tensors

| Code variable | Physical quantity | Symbol | Units | Notes |
|---------------|-------------------|--------|-------|-------|
| `sxxd` | Deviatoric stress xx | τ'_xx | Pa | "d" = deviatoric |
| `szzd` | Deviatoric stress zz | τ'_zz | Pa | |
| `sxz` | Shear stress | τ_xz | Pa | |
| `exxd` | Deviatoric strain rate xx | ε̇'_xx | s⁻¹ | |
| `ezzd` | Deviatoric strain rate zz | ε̇'_zz | s⁻¹ | |
| `exz` | Shear strain rate | ε̇_xz | s⁻¹ | |
| `wxz` | Vorticity | ω_xz | s⁻¹ | ½(∂vx/∂z - ∂vz/∂x) |
| `Tii` | Stress 2nd invariant | τ_II | Pa | √(½(τ'_xx² + τ'_zz²) + τ_xz²) |
| `Eii` | Strain rate 2nd invariant | ε̇_II | s⁻¹ | √(½(ε̇'_xx² + ε̇'_zz²) + ε̇_xz²) |
| `sxxd0` / `szzd0` / `sxz0` | Previous-step stress | τ⁰ | Pa | Used in viscoelastic update |
| `div_u` | Velocity divergence | ∇·v | s⁻¹ | Mass conservation residual |

### Tensor2D Struct

| Field | Meaning |
|-------|---------|
| `.xx` | xx component |
| `.zz` | zz component |
| `.yy` | yy (out-of-plane) component |
| `.xz` | xz component |
| `.zx` | zx component |
| `.ii` | Second invariant |
| `.ii2` | Squared second invariant |

## Material Properties

| Code variable | Physical quantity | Symbol | Units | Typical range |
|---------------|-------------------|--------|-------|---------------|
| `rho` | Density | ρ | kg/m³ | 2700–3400 |
| `eta` / `eta_n` / `eta_s` | Effective viscosity | η | Pa·s | 1e18–1e25 |
| `G` | Shear modulus | G | Pa | ~1e10 |
| `Cp` | Heat capacity | Cp | J/kg/K | 1000–1050 |
| `k` / `kx` / `kz` | Thermal conductivity | k | W/m/K | 2–4 |
| `Qr` | Radiogenic heat production | Qr | W/m³ | 0–2e-6 |
| `alp` | Thermal expansivity | α | K⁻¹ | 3e-5 |
| `bet` | Compressibility | β | Pa⁻¹ | 1e-11 |
| `C` / `C_n` / `C_s` | Cohesion | C | Pa | 1e6–5e7 |
| `phi` / `fric_n` / `fric_s` | Friction angle | φ | deg (.txt) / rad (internal) | 15–35 |
| `psi` / `dil_n` / `dil_s` | Dilation angle | ψ | deg (.txt) / rad (internal) | 0–15 |
| `d` / `d_n` / `d_s` | Grain size | d | m | 1e-4–1e-2 |
| `phi_n` (grid) | Porosity | φ | - | 0–0.3 |
| `X` | Chemical composition / melt fraction | X | - | 0–1 |

## Grid Naming Conventions

| Suffix | Meaning | Grid location | Example |
|--------|---------|---------------|---------|
| `_n` | Cell centres (nodes) | Pressure points | `eta_n`, `rho_n`, `T` |
| `_s` | Cell vertices (staggered) | Shear stress points | `eta_s`, `rho_s`, `sxz` |
| `_0` | Previous time step value | Same location | `sxxd0`, `T0_n`, `p0_n` |
| `u` | Horizontal velocity (Vx) | Vertical cell edges | `mesh.u` |
| `v` | Vertical velocity (Vz) | Horizontal cell edges | `mesh.v` |
| `p` | Pressure | Cell centres | `mesh.p` |
| `ru` / `rv` / `rp` | Residuals | Same as u/v/p | Force balance error |

### Boundary Condition Prefixes

| Prefix | Meaning |
|--------|---------|
| `BCu` | BC for Vx |
| `BCv` | BC for Vz |
| `BCp` | BC for pressure |
| `BCt` | BC for temperature |
| `BCg` | BC for vertices |

## Creep Mechanism Abbreviations

| Abbreviation | Full name | Mechanism | Equation type |
|-------------|-----------|-----------|---------------|
| `pwl` | Power law | Dislocation creep | η ∝ A^(-1/n) · exp(Q/nRT) · ε̇^(1/n-1) |
| `lin` | Linear | Diffusion creep | η ∝ d^(m/n) · exp(Q/nRT) |
| `exp` | Exponential | Peierls creep | Low-T, high-stress mechanism |
| `gbs` | Grain boundary sliding | GBS creep | Intermediate mechanism |
| `cst` / `cstv` | Constant | Constant viscosity | η = η₀ (no T/P dependence) |

### mat_prop Field Naming Pattern

For each mechanism `X` = {pwl, lin, exp, gbs}:

| Prefix | Meaning | Example |
|--------|---------|---------|
| `tX` | Type flag | `tpwl[phase]` |
| `QX` | Activation energy [J/mol] | `Qpwl[phase]` |
| `VX` | Activation volume [m³/mol] | `Vpwl[phase]` |
| `nX` | Stress exponent | `npwl[phase]` |
| `mX` | Grain size exponent | `mpwl[phase]` |
| `AX` | Pre-exponential factor | `Apwl[phase]` |
| `aX` | Water fugacity exponent | `apwl[phase]` |
| `fX` | Fugacity value | `fpwl[phase]` |
| `rX` | Fugacity power | `rpwl[phase]` |
| `FX` | Correction factor | `Fpwl[phase]` |

## Particle (Marker) Fields

| Field | Physical quantity | Units |
|-------|-------------------|-------|
| `x`, `z` | Position | m (scaled) |
| `Vx`, `Vz` | Velocity | m/s (scaled) |
| `P` | Pressure | Pa (scaled) |
| `T` | Temperature | K (scaled) |
| `phase` | Material phase ID | integer |
| `d` | Grain size | m (scaled) |
| `phi` | Porosity | - |
| `X` | Composition / melt fraction | - |
| `progress` | Reaction progress | - |
| `strain` | Total accumulated strain | - |
| `strain_el` | Elastic strain | - |
| `strain_pl` | Plastic strain | - |
| `strain_pwl` | Dislocation creep strain | - |
| `strain_exp` | Peierls strain | - |
| `strain_lin` | Diffusion creep strain | - |
| `strain_gbs` | GBS strain | - |
| `Fxx`, `Fxz`, `Fzx`, `Fzz` | Deformation gradient F | - |
| `nx`, `nz` | Director (fabric normal) | - |
| `aniso_angle` | Fabric angle θ | radians |
| `T0`, `P0`, `x0`, `z0` | Initial values | (for P-T path tracking) |
| `Tmax`, `Pmax` | Peak T and P | K, Pa |
| `generation` | Particle generation | integer (reseeding ID) |
| `noise` | Random perturbation | - |

## Solver Process Stages

| Function | Physical process |
|----------|-----------------|
| `SetBCs` | Apply boundary conditions (plate velocities, thermal BCs) |
| `P2Mastah` | Particle → grid interpolation (homogenise heterogeneous medium) |
| `RheologyDensity` / `NonNewtonianViscosityGrid` | Evaluate constitutive law: compute η from T, P, ε̇, composition |
| `BuildStokesOperator` | Assemble momentum + continuity equations into sparse matrix |
| `SolveStokes` | Solve linear system for velocity and pressure |
| `UpdateNonLinearSolution` | Update solution with line search (Newton step) |
| `ThermalSolver` | Solve energy equation: heat conduction + sources |
| `AdvectFreeSurface` | Move free surface marker chain (topography evolution) |
| `RogerGunther` / particle advection | Advect material particles with RK integration |
| `Interp_Grid2P` | Grid → particle interpolation (update particle fields) |
| `CountPartCell` / `ParticleReseeding` | Maintain particle distribution (numerical housekeeping) |
| `WriteOutputHDF5` | Save snapshot to disk (HDF5 format) |
| `MarkerChainPolyFit` | Smooth free surface with polynomial fit |
| `CellFlagging` | Flag cells as internal, boundary, or air |

## Scaling and Units

| Scale struct field | Quantity | How to derive |
|-------------------|----------|---------------|
| `scaling.eta` | Viscosity [Pa·s] | Input |
| `scaling.L` | Length [m] | Input |
| `scaling.V` | Velocity [m/s] | Input |
| `scaling.T` | Temperature [K] | Input |
| `scaling.t` | Time [s] | L/V |
| `scaling.S` | Stress [Pa] | η·V/L |
| `scaling.rho` | Density [kg/m³] | S/(g·L) |
| `scaling.m` | Mass [kg] | ρ·L³ |
| `scaling.F` | Force [N] | S·L² |
| `scaling.E` | Energy [J] | S·L³ |
| `scaling.Cp` | Heat capacity [J/kg/K] | S·L/(ρ·T·V) |
| `scaling.k` | Conductivity [W/m/K] | Derived |

**Converting**: `physical_value = internal_value × scaling.X`

## Constants

| Code name | Value | Physical meaning |
|-----------|-------|------------------|
| `zeroC` | 273.15 | 0°C in Kelvin |
| `Rg` | 8.314510 | Gas constant R [J/mol/K] |
| `PI` / `M_PI` | 3.14159... | π |
| `Rad_Earth` | 6370000 | Earth radius [m] |

## BC Type Codes

| Code | Meaning |
|------|---------|
| 0 | Dirichlet (physical boundary) |
| 11 | Dirichlet (non-physical boundary) |
| 2 | Neumann (non-physical) |
| 13 | Neumann (physical) |
| -2 | Periodic |
| -1 | Interior (not a BC) |
| 30 | Air (inactive) |

## Averaging Modes

| Code | `eta_average` value | Type |
|------|---------------------|------|
| `ARITHMETIC` | 0 | Arithmetic mean |
| `HARMONIC` | 1 | Harmonic mean |
| `GEOMETRIC` | 2 | Geometric mean |
