## ADDED Requirements

### Requirement: Skill provides variable-to-physics mapping for tensor quantities
The skill SHALL provide a lookup table mapping code variable names for stress and strain rate tensors to their physical meaning, equation notation, and SI units:
- `sxxd` / `szzd` / `sxz` → deviatoric stress components τ'_xx, τ'_zz, τ_xz [Pa]
- `exxd` / `ezzd` / `exz` → deviatoric strain rate components ε̇'_xx, ε̇'_zz, ε̇_xz [s⁻¹]
- `Tii` / `Eii` → second invariant of deviatoric stress τ_II = √(½(τ'_xx² + τ'_zz²) + τ_xz²) [Pa] / strain rate ε̇_II [s⁻¹]
- `sxxd0` / `szzd0` / `sxz0` → old (previous time step) stress, used in viscoelastic update
- `Tensor2D` struct: `.xx`, `.zz`, `.yy`, `.xz`, `.zx`, `.ii` (invariant), `.ii2` (squared invariant)
- `wxz` → vorticity component ω_xz = ½(∂vx/∂z - ∂vz/∂x) [s⁻¹]

#### Scenario: User encounters sxxd in the code
- **WHEN** the user asks what `sxxd` means
- **THEN** the skill SHALL explain it is the xx-component of deviatoric stress (τ'_xx), the "d" suffix stands for "deviatoric" (trace removed), stored in Pa (or scaled units internally)

### Requirement: Skill provides variable-to-physics mapping for material properties
The skill SHALL map material property variable names to physical quantities:
- `rho` → density ρ [kg/m³]
- `eta` / `eta_n` / `eta_s` → effective viscosity η [Pa·s], at cell centres (_n) or vertices (_s)
- `G` → elastic shear modulus G [Pa]
- `Cp` → specific heat capacity [J/kg/K]
- `k` / `kx` / `kz` → thermal conductivity [W/m/K], directional components
- `Qr` → radiogenic heat production [W/m³]
- `alp` → thermal expansivity α [K⁻¹]
- `bet` → compressibility β [Pa⁻¹]
- `C` / `C_n` / `C_s` → cohesion [Pa]
- `phi` / `fric_n` / `fric_s` → friction angle φ [degrees in .txt, radians internally]
- `psi` / `dil_n` / `dil_s` → dilation angle ψ [degrees/.txt, radians internally]
- `d` / `d_n` / `d_s` → grain size [m]
- `phi_n` (on grid, not friction) → porosity φ (context-dependent: distinguish from friction angle)

#### Scenario: User asks what alp and bet mean
- **WHEN** the user encounters `alp` and `bet` in a phase definition
- **THEN** the skill SHALL explain: `alp` = thermal expansivity α (typically 3.2e-5 K⁻¹), `bet` = compressibility β (typically 1e-11 Pa⁻¹), used in the equation of state ρ = ρ₀(1 - α(T-T₀))(1 + β(P-P₀))

### Requirement: Skill provides grid array naming conventions
The skill SHALL explain the systematic naming convention for grid arrays:
- `_n` suffix → cell centres (nodes), e.g., `eta_n`, `T`, `P`, `rho_n`
- `_s` suffix → cell vertices (staggered), e.g., `eta_s`, `rho_s`, `sxz`
- `_0` suffix → value from previous time step, e.g., `sxxd0`, `T0_n`, `p0_n`
- `u` / `v` → horizontal (Vx) and vertical (Vz) velocity, stored on staggered edges
- `p` → pressure, stored at cell centres
- `BCu` / `BCv` / `BCp` / `BCt` → boundary condition structs for Vx, Vz, P, T
- `phase_perc_n` / `phase_perc_s` → phase volume fractions per cell (for averaging)
- `eII_pwl`, `eII_exp`, `eII_lin`, `eII_gbs`, `eII_cst`, `eII_el`, `eII_pl` → strain rate partitioned by mechanism

#### Scenario: User asks what _n and _s mean
- **WHEN** the user sees `eta_n` and `eta_s` and asks the difference
- **THEN** the skill SHALL explain: `_n` = cell-centroid value (used for normal stress), `_s` = vertex value (used for shear stress), reflecting the staggered grid layout where different physical quantities live at different grid locations

### Requirement: Skill provides creep mechanism abbreviation glossary
The skill SHALL decode all flow law abbreviations used in the code:
- `pwl` → power law (dislocation creep): η_pwl = f(A, n, Q, V, T, P, ε̇)
- `lin` → linear (diffusion creep): η_lin = f(A, n, m, Q, V, T, P, d)
- `exp` → exponential (Peierls creep): low-T, high-stress mechanism
- `gbs` → grain boundary sliding: η_gbs = f(A, n, m, Q, V, T, P, d)
- `cst` / `cstv` → constant viscosity (η = η₀, no T/P/ε̇ dependence)
- `eta_ve` → visco-elastic effective viscosity: η_ve = η·G·Δt/(η + G·Δt)
- `eta_vep` → visco-elasto-plastic viscosity (after yield correction)
- `B_pwl`, `C_pwl` → pre-exponential factor intermediates in viscosity calculation
- `Fpwl`, `Flin`, `Fgbs` → correction factors for flow laws
- Subscript meaning in `mat_prop`: `t` prefix (e.g., `tpwl`) = type flag, `Q` = activation energy, `V` = activation volume, `n` = stress exponent, `m` = grain size exponent, `A` = pre-exponential, `a`/`f`/`r` = fugacity-related parameters

#### Scenario: User asks what pwlv=40 and linv=40 mean together
- **WHEN** the user sees both `pwlv = 40` and `linv = 40` for a mantle phase
- **THEN** the skill SHALL explain: both dislocation (`pwl`) and diffusion (`lin`) creep are active using Hirth & Kohlstedt 2003 Dry Olivine parameters; the code computes both viscosities and combines them harmonically, so the weaker mechanism dominates at each point

### Requirement: Skill provides particle/marker field glossary
The skill SHALL map marker (particle) field names to physical quantities:
- `x`, `z` → particle position [m scaled]
- `Vx`, `Vz` → interpolated velocity [m/s scaled]
- `P` → pressure [Pa scaled]
- `T` → temperature [K scaled]
- `phase` → material phase index (integer)
- `d` → grain size [m scaled]
- `phi` → porosity
- `X` → chemical composition / melt fraction
- `progress` → reaction progress variable
- `strain`, `strain_el`, `strain_pl` → total, elastic, plastic accumulated strain
- `strain_pwl`, `strain_exp`, `strain_lin`, `strain_gbs` → strain partitioned by creep mechanism
- `Fxx`, `Fxz`, `Fzx`, `Fzz` → finite deformation gradient tensor F components
- `nx`, `nz` → director (fabric normal) components
- `aniso_angle` → fabric angle θ [radians]
- `generation` → particle generation (for reseeding tracking)
- `noise` → random perturbation field
- `T0`, `P0`, `x0`, `z0` → initial values (for P-T path tracking)
- `Tmax`, `Pmax` → max T and P experienced (for metamorphic path)

#### Scenario: User asks about finite strain tracking on particles
- **WHEN** the user asks what `Fxx`, `Fxz`, `Fzx`, `Fzz` are
- **THEN** the skill SHALL explain these are components of the deformation gradient tensor F, tracking accumulated finite strain on each particle, used to compute finite strain ellipse aspect ratio (`FS_AR`) and to evolve anisotropy when `ani_fstrain = 1`

### Requirement: Skill provides solver process stage glossary
The skill SHALL describe what is physically happening at each major stage of the time loop, mapping function names to physical processes:
- `SetBCs` → apply mechanical boundary conditions (plate velocities, free slip, periodic)
- `P2Mastah` → interpolate particle properties to grid (particles → grid: homogenisation of heterogeneous medium)
- `RheologyDensity` → evaluate constitutive law: compute effective viscosity from T, P, ε̇, composition at each grid point
- `BuildStokesOperator` / `StokesAssemblyDecoupled` → assemble discretised momentum + continuity equations into sparse matrix
- `SolveStokes` → solve the linear system for velocity and pressure increments
- `UpdateNonLinearSolution` → update velocity/pressure with line search (Newton step)
- `EvaluateRheology` → re-evaluate viscosity with updated strain rates (nonlinear iteration)
- `ThermalSolver` → solve the energy equation: heat conduction + sources → temperature update
- `AdvectFreeSurface` → move the free surface marker chain with the velocity field (topography evolution)
- `RogerGunther` / particle advection → move material particles with RK integration (material transport)
- `Interp_Grid2P` → interpolate updated fields back to particles (grid → particles)
- `CountPartCell` / `ParticleReseeding` → maintain particle distribution (numerical housekeeping)
- `WriteOutputHDF5` → save snapshot to disk

#### Scenario: User asks what P2Mastah does physically
- **WHEN** the user encounters `P2Mastah` in the code
- **THEN** the skill SHALL explain: this is the particle-to-grid interpolation step — it computes grid-level material properties (viscosity, density, conductivity, etc.) by averaging over particles within each cell, effectively homogenising the heterogeneous marker cloud into continuous fields the finite difference solver can use

### Requirement: Skill provides scaling and unit convention reference
The skill SHALL document the unit conventions:
- All internal computations use non-dimensional (scaled) values
- Physical value = internal value × scale factor
- The `scale` struct fields: `eta` (viscosity), `L` (length), `V` (velocity), `T` (temperature), `t` (time = L/V), `a` (acceleration), `E` (energy), `S` (stress = η·V/L), `m` (mass), `rho` (density), `F` (force), `J` (Joule), `W` (Watt), `Cp` (heat capacity), `rhoE` (energy density), `k` (conductivity)
- Common conversion patterns: `MinMaxArray(array, scaling.X, N, "label")` prints min/max in physical units

#### Scenario: User wants to convert internal values to SI
- **WHEN** the user reads an internal value and wants the physical equivalent
- **THEN** the skill SHALL explain: multiply by the corresponding scale factor, e.g., stress_Pa = stress_internal × scaling.S, where S = η_scale · V_scale / L_scale

### Requirement: Skill provides constants and special values reference
The skill SHALL document the physical constants defined in the code:
- `zeroC = 273.15` → 0°C in Kelvin
- `Rg = 8.314510` → universal gas constant R [J/mol/K]
- `PI = 3.14159...` → π
- `Rad_Earth = 6370000` → Earth radius [m]
- BC type codes: `0` = Dirichlet (physical boundary), `11` = Dirichlet (non-physical), `2` = Neumann (non-physical), `13` = Neumann (physical), `-2` = periodic, `-1` = interior, `30` = air/inactive
- Averaging modes: `ARITHMETIC = 0`, `HARMONIC = 1`, `GEOMETRIC = 2`

#### Scenario: User asks why 273.15 appears everywhere
- **WHEN** the user sees `zeroC` or `273.15` in temperature calculations
- **THEN** the skill SHALL explain this is the Celsius-to-Kelvin offset: T_Kelvin = T_Celsius + 273.15, and the code uses Kelvin internally
