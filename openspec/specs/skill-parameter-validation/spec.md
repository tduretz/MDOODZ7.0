## ADDED Requirements

### Requirement: Skill validates scaling parameters
The skill SHALL check that the four scaling parameters in the `SCALES` section fall within physically reasonable bounds:
- `eta` (viscosity scale): SHALL be between 1e15 and 1e28 Pa·s. Typical Earth values: 1e19–1e24 Pa·s.
- `L` (length scale): SHALL be between 1e2 and 1e7 m. Typical: 1e3–1e6 m for crustal/lithospheric models.
- `V` (velocity scale): SHALL be between 1e-15 and 1e-5 m/s. Typical: 1e-10–1e-8 m/s (mm/yr to cm/yr scale).
- `T` (temperature scale): SHALL be between 1 and 1e4 K. Typical: 1e2–1e3 K.

#### Scenario: User has an unreasonable viscosity scale
- **WHEN** `eta` in `SCALES` is outside 1e15–1e28
- **THEN** the skill SHALL warn that this is outside the physically meaningful range for geological materials and suggest typical values (1e20–1e23 for crust, 1e18–1e21 for asthenosphere)

### Requirement: Skill validates domain and resolution parameters
The skill SHALL check `SPACE-TIME` parameters:
- `Nx`, `Nz` SHALL be positive integers, typically 50–2000. Warn if < 20 (too coarse) or > 4000 (very expensive).
- Domain bounds (`xmin`, `xmax`, `zmin`, `zmax`) SHALL have xmax > xmin and zmax > zmin.
- Domain extent SHALL be physically reasonable: horizontal 1e3–1e7 m, vertical 1e3–1e6 m for typical geodynamic models.
- Cell aspect ratio `dx/dz` SHOULD be between 0.1 and 10 to avoid numerical issues.
- `dt` SHALL be positive. Warn if dt > 1e14 s (~3 Myr) which may cause large Courant numbers.
- `Courant` SHALL be between 0 and 1, typically 0.3–0.5.
- `Nt` SHALL be positive.

#### Scenario: User sets Nx=10 for a 300 km wide model
- **WHEN** `Nx = 10` with domain width 300 km
- **THEN** the skill SHALL warn that cell size is 30 km, which is too coarse to resolve crustal structures, and suggest at minimum Nx=100–200 for 1.5–3 km resolution

### Requirement: Skill validates material properties per phase
The skill SHALL check per-phase properties in `PHASE PROPERTIES` against physically reasonable ranges:
- `rho` (density): 1000–4500 kg/m³ (water=1000, sediments=2000–2500, crust=2700–2900, mantle=3200–3400, eclogite=3500+)
- `G` (shear modulus): 1e9–1e11 Pa (typical ~1e10 Pa for rocks)
- `Cp` (heat capacity): 500–1500 J/kg/K (typical 1000–1050 for silicates)
- `k` (thermal conductivity): 0.5–10 W/m/K (typical 2–3 for crust, 3–4 for mantle)
- `Qr` (radiogenic heating): 0–5e-6 W/m³ (typical 0.5–1.5e-6 for crust, ~0 for mantle)
- `C` (cohesion): 1e5–1e9 Pa (typical 1e6–5e7 for rocks)
- `phi` (friction angle): 0–45° (typical 15–35° for rocks)
- `alp` (thermal expansivity): 1e-6–1e-4 K⁻¹ (typical 2–4e-5)
- `bet` (compressibility): 1e-13–1e-9 Pa⁻¹ (typical 1e-11)
- `aniso_factor`: 1.0–20.0 (1.0 = isotropic, >1 = anisotropic)
- `aniso_angle`: 0–180° (initial foliation orientation)

#### Scenario: User sets density of 330 instead of 3300
- **WHEN** a phase has `rho = 330`
- **THEN** the skill SHALL warn this is far below any geological material density and suggest the user likely meant 3300 kg/m³

### Requirement: Skill validates flow law index consistency
The skill SHALL check that flow law indices (`pwlv`, `linv`, `gbsv`, `expv`) reference valid entries in the `FlowLaws.c` database. The skill SHALL warn if:
- All creep mechanisms are disabled (`cstv=0, pwlv=0, linv=0, gbsv=0, expv=0`) — this means no viscosity is defined
- `expv` (Peierls) is enabled without `pwlv` (dislocation) — unusual and likely unintended
- `linv` (diffusion) is enabled but no grain size initialisation is provided
- Phase uses `cstv > 0` (constant viscosity) together with creep laws — constant viscosity overrides creep

#### Scenario: User has a phase with no active creep mechanism
- **WHEN** all of `cstv`, `pwlv`, `linv`, `gbsv`, `expv` are 0 for a phase
- **THEN** the skill SHALL warn that no viscosity mechanism is defined for this phase, which will cause undefined behaviour, and suggest setting at least `cstv = 1` for constant viscosity

### Requirement: Skill validates switch combinations
The skill SHALL check for incompatible or problematic switch combinations:
- `thermal = 0` with `shear_heating = 1`: warn that shear heating has no effect without thermal solver
- `elastic = 1` without `G` defined for all phases: error — shear modulus required
- `free_surface = 1` with `periodic_x = 1`: warn about potential incompatibility
- `anisotropy = 1` requires at least one phase with `aniso_factor > 1`: warn if all phases have factor = 1 (effectively isotropic)
- `melting = 1` requires phases with `melt` flag: warn if no phase has melting enabled
- `pure_shear_ALE = 1` requires `bkg_strain_rate != 0`: warn if strain rate is zero

#### Scenario: User enables shear heating without thermal solver
- **WHEN** `shear_heating = 1` but `thermal = 0`
- **THEN** the skill SHALL warn that shear heating requires `thermal = 1` to have any effect

### Requirement: Skill validates boundary condition consistency
The skill SHALL check:
- `bkg_strain_rate` sign convention: negative = extension, positive = compression. Warn if sign seems inconsistent with intended deformation.
- Gravity: `gz` should typically be -9.81 (negative = downward). Warn if positive or zero for models expecting gravity.
- `gx` is typically 0 unless modelling tilted problems.
- `penalty` parameter: typical range 1e0–1e5. Warn if < 1 (weak coupling) or > 1e7 (potential ill-conditioning).
- `min_eta` / `max_eta` bounds: typical 1e17–1e25. Warn if range is narrower than 4 orders of magnitude.

#### Scenario: User sets gravity as positive
- **WHEN** `gz = 9.81` (positive)
- **THEN** the skill SHALL warn that positive gz means upward gravity, which is non-physical for most models, and suggest `gz = -9.81`

### Requirement: Skill validates time step and Courant number
The skill SHALL check for potential Courant number violations:
- Compute effective CFL: Courant_eff = V_max · dt / min(dx, dz)
- Warn if `dt` combined with typical velocities (from `V` scale and `bkg_strain_rate`) would exceed `Courant` limit
- Recommend `constant_dt = 0` (adaptive time stepping) for most simulations
- Warn if `dt` > 1e14 s with `constant_dt = 1` (fixed large steps risk instability)

#### Scenario: User uses a very large fixed time step
- **WHEN** `constant_dt = 1` and `dt` > L/(10·V) where L is domain size and V is velocity scale
- **THEN** the skill SHALL warn about potential Courant violations and recommend adaptive time stepping (`constant_dt = 0`)

### Requirement: Skill provides a checklist-style validation workflow
The skill SHALL present a structured validation checklist that users can follow when preparing or reviewing a `.txt` file:
1. Check scaling parameters are in reasonable ranges
2. Verify domain bounds and resolution are consistent
3. Confirm `Nb_phases` matches the number of `ID` blocks
4. Validate each phase's material properties
5. Check flow law indices reference valid database entries
6. Verify switch combinations are compatible
7. Confirm boundary conditions are consistent with the intended setup
8. Estimate Courant number for the chosen time step

#### Scenario: User asks to validate their configuration
- **WHEN** the user provides a `.txt` file for review
- **THEN** the skill SHALL systematically walk through the checklist, flagging any values outside recommended ranges and suggesting corrections
