---
name: skill-physics-theory
description: MDOODZ governing equations and numerical methods — conservation of momentum (Stokes), energy, and mass, constitutive relations, non-dimensionalisation and scaling, marker-in-cell advection, free surface algorithm, pure shear ALE mode, and theoretical background.
---

# Physics and Theory

## Governing Equations

### Conservation of Momentum (Stokes)

∂τ_ij/∂x_j - ∂P/∂x_i + ρg_i = 0

Inertia is neglected (infinite Prandtl number assumption, appropriate for geological time scales). This gives the Stokes equations for slow viscous flow.

**In component form:**
- x: ∂τ_xx/∂x + ∂τ_xz/∂z - ∂P/∂x + ρg_x = 0
- z: ∂τ_xz/∂x + ∂τ_zz/∂z - ∂P/∂z + ρg_z = 0

### Conservation of Mass (Continuity)

**Incompressible**: ∂v_i/∂x_i = 0

**Weakly compressible** (when β > 0): ∂v_i/∂x_i = -1/K · DP/Dt

where K = 1/β is the bulk modulus, controlled by the `bet` parameter per phase.

### Conservation of Energy

ρCp · DT/Dt = ∂/∂x_i(k · ∂T/∂x_i) + H_shear + H_adiab + Qr

| Source term | Expression | Switch |
|-------------|-----------|--------|
| Shear heating | H_shear = τ_ij · ε̇_ij | `shear_heating = 1` |
| Adiabatic heating | H_adiab = α·T·v_z·ρ·g | `adiab_heating = 1` |
| Radiogenic production | Qr (per phase, W/m³) | Always active if Qr > 0 |

## Constitutive Relations

### Viscous

τ_ij = 2η · ε̇_ij

where ε̇_ij is the deviatoric strain rate tensor.

### Visco-Elastic (Maxwell Model)

τ_ij + (η/G) · Dτ_ij/Dt = 2η · ε̇_ij

Discretised using a Jaumann objective stress rate:

τ^new_ij = 2η_ve · (ε̇_ij + τ⁰_ij,rotated / (2·G·Δt))

where:
- η_ve = η·G·Δt / (η + G·Δt) is the effective visco-elastic viscosity
- τ⁰_ij,rotated is the old deviatoric stress rotated by the Jaumann rate
- G is the shear modulus, Δt is the time step

### Visco-Elasto-Plastic

Above visco-elastic model with Drucker-Prager yield criterion:

τ_II ≤ C·cos(φ) + P·sin(φ)

When yield is exceeded, viscosity is reduced:
η_pl = τ_yield / (2·ε̇_II)

### Anisotropic

With anisotropy enabled, the viscosity depends on the angle between the director (foliation normal) and the stress tensor. The stress invariant is modified:

Y2 = ½(τ_xx² + τ_zz² + τ_yy²) + (δ·τ_xz)²

where δ is the anisotropy factor. See the skill-anisotropy documentation for details.

## Non-Dimensionalisation

MDOODZ uses a non-dimensional formulation. Four independent scales are chosen in the `SCALES` section of the `.txt` file:

| Scale | Parameter | Symbol | Units |
|-------|-----------|--------|-------|
| Viscosity | `eta` | η₀ | Pa·s |
| Length | `L` | L₀ | m |
| Velocity | `V` | V₀ | m/s |
| Temperature | `T` | T₀ | K |

### Derived Scales

All other scales are derived from these four:

| Quantity | Expression | Code field |
|----------|-----------|------------|
| Time | t₀ = L₀/V₀ | `scaling.t` |
| Stress | S₀ = η₀·V₀/L₀ | `scaling.S` |
| Strain rate | ε̇₀ = V₀/L₀ | (= 1/t₀) |
| Density | ρ₀ = S₀/(g·L₀) | `scaling.rho` |
| Mass | m₀ = ρ₀·L₀³ | `scaling.m` |
| Force | F₀ = S₀·L₀² | `scaling.F` |
| Energy | E₀ = S₀·L₀³ | `scaling.E` |
| Heat capacity | Cp₀ | `scaling.Cp` |
| Conductivity | k₀ | `scaling.k` |

### Converting Between Physical and Scaled

- **Physical → Scaled**: value_scaled = value_physical / scale_factor
- **Scaled → Physical**: value_physical = value_scaled × scale_factor

The code handles this automatically at input/output. Internal computations use scaled values (order ~1).

## Marker-in-Cell (PIC) Advection

### Method

Material properties are carried by Lagrangian particles (markers) and the governing equations are solved on a fixed Eulerian grid.

### Workflow per Time Step

1. **Particle → Grid** (`P2Mastah`): Interpolate particle properties (phase, η, ρ, Cp, k, T, etc.) to grid nodes using weighted averaging
2. **Solve on grid**: Stokes + thermal equations
3. **Grid → Particle** (`Interp_Grid2P`): Interpolate updated velocity, stress, temperature back to particles
4. **Advect particles**: Move particles with the velocity field using Runge-Kutta integration

### Runge-Kutta Order

Controlled by the `RK` parameter in `.txt`:

| RK | Order | Accuracy | Cost |
|----|-------|----------|------|
| 1 | First-order Euler | O(Δt) | Cheapest |
| 2 | Second-order midpoint | O(Δt²) | Moderate |
| 4 | Fourth-order RK4 | O(Δt⁴) | Best (default) |

### Interpolation

The `interp_stencil` parameter controls the interpolation kernel width (1 = nearest cell, larger = smoother).

### Particle Reseeding

Cells can become depleted or overpopulated as particles advect. Reseeding maintains a healthy particle distribution:

| Parameter | Description |
|-----------|-------------|
| `reseed_markers` | Enable reseeding (0/1) |
| `reseed_mode` | 0 = old counting, 1 = new counting |
| `min_part_cell` | Minimum particles per cell |
| `Nx_part`, `Nz_part` | Initial particles per cell (h × v) |

## Free Surface Algorithm

### Marker Chain

The free surface is represented by a marker chain (polyline) at the top of the domain. Each chain marker has an (x, z) position.

### Evolution

1. The chain is advected by the velocity field at each time step
2. Polynomial fit smooths the chain (`MarkerChainPolyFit`)
3. Grid cells above the surface are flagged as "air" (BC type 30)

### Stabilisation

`free_surface_stab` parameter damps oscillations that can arise from the free surface coupling.

### Surface Processes

| `surface_processes` | Mode |
|---------------------|------|
| 0 | None |
| 1 | Diffusion only (coefficient `surf_diff`) |
| 2 | Instantaneous fill to base level |
| 3 | Diffusion + source term |

| Parameter | Description |
|-----------|-------------|
| `surf_diff` | Topography diffusion coefficient [m²/s] |
| `surf_sedirate` | Sedimentation rate [m/s] |
| `surf_baselev` | Base level for sedimentation [m] |
| `surf_ised1`, `surf_ised2` | Sediment phase indices |

## Pure Shear ALE Mode

When `pure_shear_ALE = 1`:

- Domain boundaries move at the background strain rate (`bkg_strain_rate`)
- Grid is remeshed to the new domain
- Particles stay inside; domain stretches/compresses with them
- Negative `bkg_strain_rate` = extension; positive = compression
- Area is conserved (horizontal extension ↔ vertical thinning)

This is the standard mode for modelling lithospheric rifting at constant extension rate.

## Key Source Files

- Main solver loop: `MDLIB/Main_DOODZ.c`
- Particle advection: `MDLIB/AdvectionRoutines.c`
- Particle routines: `MDLIB/ParticleRoutines.c`
- Particle reseeding: `MDLIB/ParticleReseeding.c`
- Free surface: `MDLIB/FreeSurface.c`
- Grid routines: `MDLIB/GridRoutines.c`
- Thermal routines: `MDLIB/ThermalRoutines.c`
