## ADDED Requirements

### Requirement: Skill documents the viscosity computation pipeline
The skill SHALL explain the full viscosity evaluation chain: for each material phase, the effective viscosity is computed by combining contributions from multiple creep mechanisms (dislocation, diffusion, Peierls/exponential, grain boundary sliding, constant) with elastic and plastic branches. The skill SHALL describe that these are evaluated in `ViscosityConciseAniso` (in `AnisotropyRoutines.c` / `RheologyDensity.c`) and combined via harmonic averaging.

#### Scenario: User asks how viscosity is computed
- **WHEN** the user asks about the viscosity pipeline
- **THEN** the skill SHALL present the chain: identify active mechanisms per phase (from `.txt` switches `cstv`, `pwlv`, `linv`, `gbsv`, `expv`) → compute each mechanism's viscosity → harmonic average → apply elastic correction (Maxwell model: η_ve = η·G·dt/(η + G·dt)) → apply plastic correction (if yield stress exceeded) → clamp to [min_eta, max_eta]

### Requirement: Skill documents power law (dislocation) creep
The skill SHALL explain the power law (dislocation creep) flow law: η_pwl = A^(-1/n) · exp(Q + P·V)/(n·R·T) · ε̇_II^(1/n - 1), where `pwlv` in the `.txt` file selects from a database of published flow laws in `FlowLaws.c`. The skill SHALL list the key database entries including at minimum:
- Case 1: User-defined power law (eta0, n, Q)
- Case 10: Wet Quartzite (Ranalli 1995)
- Case 11: Maryland Diabase (Ranalli 1995, crust)
- Case 12: Dry Granulite (Ranalli 1995)
- Case 19: Westerly Granite (Hansen & Carter 1983)
- Case 25: Plagioclase An75 (Ranalli 1995)
- Case 40: Dry Olivine (Hirth & Kohlstedt 2003, mantle)
- Case 41: Wet Olivine (Hirth & Kohlstedt 2003)

#### Scenario: User asks what pwlv=40 means
- **WHEN** the user encounters `pwlv = 40` in a `.txt` file
- **THEN** the skill SHALL explain this selects Dry Olivine dislocation creep from Hirth & Kohlstedt (2003), with specific A, n, Q, V values, suitable for lithospheric mantle rheology

### Requirement: Skill documents linear (diffusion) creep
The skill SHALL explain diffusion creep: η_lin = A^(-1/n) · d^(m/n) · exp(Q + P·V)/(n·R·T) · ε̇_II^(1/n - 1), where grain size `d` matters (m > 0). The `linv` switch selects from the database. Key entries include Case 40 (Dry Olivine diffusion, Hirth & Kohlstedt 2003).

#### Scenario: User asks about grain size sensitive creep
- **WHEN** the user asks about diffusion creep or grain size dependence
- **THEN** the skill SHALL explain that `linv` activates diffusion creep, the grain size exponent `m` controls sensitivity, and initial grain size is set via `SetGrainSize` callback

### Requirement: Skill documents Peierls (exponential) creep
The skill SHALL explain Peierls mechanism (low-temperature, high-stress plasticity): activated by `expv` switch, applicable mainly to olivine below ~1200°C. The skill SHALL note the exponential stress dependence and that it is capped at `TmaxPeierls = 1200 + 273.15 K`.

#### Scenario: User wants to include Peierls creep for mantle
- **WHEN** the user asks about low-temperature creep
- **THEN** the skill SHALL explain `expv = 40` activates Peierls mechanism for dry olivine, effective only below ~1200°C

### Requirement: Skill documents elastic and plastic branches
The skill SHALL explain:
- **Elasticity**: Maxwell visco-elastic model where η_ve = η·G·Δt/(η + G·Δt), activated by `elastic = 1` in `.txt`, requiring shear modulus `G` per phase
- **Plasticity**: Drucker-Prager yield criterion with cohesion `C`, friction angle `phi`, dilation angle `psi`, and optional strain softening (`phi_soft`, `coh_soft`, `pls_start`, `pls_end`), tensile strength `sig_tens`
- **Viscoplastic regularisation**: `eta_vp` parameter for regularised plasticity

#### Scenario: User asks how to enable viscoelasticity
- **WHEN** the user asks about elastic behaviour
- **THEN** the skill SHALL explain: set `elastic = 1` globally and provide `G` (shear modulus, Pa) for each phase in the `PHASE PROPERTIES` section

#### Scenario: User asks about strain softening
- **WHEN** the user asks about weakening during deformation
- **THEN** the skill SHALL explain the softening parameters: `phi_end` (final friction angle), `C_end` (final cohesion), `pls_start`/`pls_end` (accumulated plastic strain range over which softening occurs)

### Requirement: Skill documents the mat_prop structure
The skill SHALL explain the `mat_prop` struct that holds per-phase material properties (up to 20 phases), including density (`rho`), thermal properties (`Cp`, `k`, `Qr`, `alp`, `bet`), elastic modulus (`G`), plasticity (`C`, `phi`, `psi`, `Slim`), and flow law indices (`cstv`, `pwlv`, `linv`, `gbsv`, `expv`) with their associated creep parameters.

#### Scenario: User needs to understand material property arrays
- **WHEN** the user asks about data structures for materials
- **THEN** the skill SHALL explain that `mat_prop` contains arrays of size 20 (max phases) indexed by phase ID, populated from the `PHASE PROPERTIES` section of the `.txt` file

### Requirement: Skill documents density models
The skill SHALL explain how density is computed: base density `rho` with thermal expansion (`alp`: α) and compressibility (`bet`: β), following ρ = ρ₀(1 - α(T-T₀))(1 + β(P-P₀)). The `density_model` switch allows alternative formulations including phase-diagram-based density.

#### Scenario: User asks how density varies with temperature
- **WHEN** the user asks about thermal buoyancy
- **THEN** the skill SHALL explain the equation of state with `alp` (thermal expansivity, typically 3.2e-5 K⁻¹) and `bet` (compressibility, typically 1e-11 Pa⁻¹)
