## ADDED Requirements

### Requirement: Skill documents conservation equations
The skill SHALL present the governing conservation equations solved by MDOODZ:
- **Conservation of momentum** (Stokes): ∂τᵢⱼ/∂xⱼ - ∂P/∂xᵢ + ρgᵢ = 0 (inertia neglected)
- **Conservation of mass** (continuity): ∂vᵢ/∂xᵢ = 0 (incompressible) or ∂vᵢ/∂xᵢ = -1/K · DP/Dt (compressible, controlled by `bet` parameter)
- **Conservation of energy**: ρCp · DT/Dt = ∂/∂xᵢ(k · ∂T/∂xᵢ) + H_shear + H_adiab + Qr
where H_shear = τᵢⱼε̇ᵢⱼ (shear heating), H_adiab = αTvzρg (adiabatic heating), Qr = radiogenic production

#### Scenario: User asks about the governing equations
- **WHEN** the user asks what equations MDOODZ solves
- **THEN** the skill SHALL present the momentum (Stokes), mass (continuity), and energy equations with notation matching the codebase variables

### Requirement: Skill documents constitutive relations
The skill SHALL explain how deviatoric stress relates to strain rate through the constitutive law:
- **Viscous**: τᵢⱼ = 2η · ε̇ᵢⱼ
- **Visco-elastic** (Maxwell): τᵢⱼ + η/G · Dτᵢⱼ/Dt = 2η · ε̇ᵢⱼ, discretised as τᵢⱼ = 2η_ve · (ε̇ᵢⱼ + τ⁰ᵢⱼ/(2G·Δt))
- **Visco-elasto-plastic**: above with Drucker-Prager yield: τ_II ≤ C·cos(φ) + P·sin(φ)
- **Anisotropic**: viscosity depends on angle between director and stress, modifying the effective viscosity in the director frame

#### Scenario: User asks about viscoelastic stress update
- **WHEN** the user asks how stress is updated with elasticity
- **THEN** the skill SHALL explain the objective Jaumann stress rate discretisation: τ^new = 2η_ve(ε̇ + τ⁰/(2GΔt)) where τ⁰ is the rotated old stress, and η_ve = η·G·Δt/(η + G·Δt)

### Requirement: Skill documents non-dimensionalisation
The skill SHALL explain the scaling system defined by the `scale` struct and the `SCALES` section of the `.txt` file. Four independent scales are chosen: viscosity (η), length (L), velocity (V), temperature (T). All other scales are derived:
- Time: t = L/V
- Stress: S = η·V/L
- Density: ρ = S/(g·L)
- Mass: m = ρ·L³
- Force: F = S·L²
- Energy: E = S·L³
- Heat capacity: Cp = E/(m·T) = S·L/(ρ·T·V)
- Thermal conductivity: k = S·L²·V/(L·T) ... etc.

#### Scenario: User asks why input values look wrong in output
- **WHEN** the user is confused by non-dimensional values in the code
- **THEN** the skill SHALL explain scaling: physical values are divided by the corresponding scale factor on input and multiplied back on output. Show that e.g., `eta = 1e22` in `.txt` means the viscosity scale, and internal computations use η/η_scale ≈ 1

### Requirement: Skill documents marker-in-cell advection
The skill SHALL explain the particle-in-cell (PIC / marker-in-cell) method used for advection:
- Particles carry material properties (phase, temperature, stress, strain, grain size, anisotropy angle)
- Particle positions are advected using Runge-Kutta integration (order controlled by `RK` parameter: 1, 2, or 4)
- Properties are interpolated particle→grid (`P2Mastah`) before the mechanical solve and grid→particle after
- Particle reseeding (`reseed_markers`, `reseed_mode`) prevents cells from becoming empty or overpopulated
- The `interp_stencil` parameter controls interpolation kernel width

#### Scenario: User asks about material advection
- **WHEN** the user asks how materials are transported
- **THEN** the skill SHALL explain the marker-in-cell approach: particles are Lagrangian tracers carrying material properties, advected by velocity field, with RK4 integration as default (`RK = 4`)

### Requirement: Skill documents the free surface algorithm
The skill SHALL explain the free surface implementation:
- Represented by a marker chain (polyline) at the top of the domain
- Advected by the velocity field at each time step
- Stabilised by `free_surface_stab` parameter to damp oscillations
- Cells above the surface are flagged as "air" (type 30 BC)
- Surface processes: diffusion (`surface_processes = 1`, coefficient `surf_diff`), sedimentation (`surf_sedirate`, `surf_baselev`)
- Topography update controlled by `topo_update`

#### Scenario: User asks about erosion and sedimentation
- **WHEN** the user asks about surface processes
- **THEN** the skill SHALL explain: set `free_surface = 1`, `surface_processes = 1` (diffusion) or `2` (instant fill) or `3` (diffusion + source), configure `surf_diff` for diffusion coefficient, `surf_sedirate` for deposition rate, `surf_baselev` for base level, and assign sediment phase IDs via `surf_ised1`/`surf_ised2`

### Requirement: Skill documents the pure shear ALE mode
The skill SHALL explain the Arbitrary Lagrangian-Eulerian (ALE) mode for pure shear: when `pure_shear_ALE = 1`, the box boundaries move at the background strain rate (`bkg_strain_rate`) while the grid is remeshed, simulating bulk extension or compression without particles leaving the domain.

#### Scenario: User wants to model lithospheric extension
- **WHEN** the user asks about modelling rifting with a constant extension rate
- **THEN** the skill SHALL explain: set `pure_shear_ALE = 1` and `bkg_strain_rate = -1e-15` (negative for extension), which stretches the domain horizontally and compresses vertically to conserve area
