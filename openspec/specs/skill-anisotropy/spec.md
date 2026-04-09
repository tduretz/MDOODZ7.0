## ADDED Requirements

### Requirement: Skill documents the director-based anisotropy formulation
The skill SHALL explain that MDOODZ tracks mechanical anisotropy via a director field (n₁, n₂) representing the orientation of the weak plane (foliation) on each marker particle. The director evolves with the flow (rotation by the vorticity tensor) and optionally via finite strain tracking (`ani_fstrain`). The skill SHALL reference the anisotropic rifting paper's formulation of viscous anisotropy using the director approach.

#### Scenario: User asks what the anisotropy model is
- **WHEN** the user asks how anisotropy is represented
- **THEN** the skill SHALL explain: each particle carries an anisotropy angle (`aniso_angle` on markers, initialised via `SetAnisoAngle` callback or `aniso_angle` in `.txt`), which defines the director orientation. The director (n₁, n₂) = (cos θ, sin θ) tracks foliation orientation and rotates with material flow.

### Requirement: Skill documents anisotropy factor definitions
The skill SHALL explain the anisotropy factors δ_n and δ_s that quantify the degree of viscous anisotropy:
- `aniso_factor` in `.txt`: the ratio of normal-to-shear viscosity contrast
- `ani_fac_max`: maximum anisotropy factor (saturation limit)
- When `ani_fstrain = 1`: the anisotropy factor evolves with accumulated finite strain via `AnisoFactorEvolv(FS_AR, ani_fac_max)`, clamped at `ani_fac_max`
- The function `Y2` computes the anisotropic second invariant of deviatoric stress as 0.5(τ_xx² + τ_zz² + τ_yy²) + (δ·τ_xz)²

#### Scenario: User asks what aniso_factor means
- **WHEN** the user encounters `aniso_factor = 2.0` in a `.txt` file
- **THEN** the skill SHALL explain this sets the initial ratio of directional viscosity contrast — a value of 1.0 means isotropic, values > 1 indicate increasing anisotropy with the weak direction aligned to the director

### Requirement: Skill documents coupling with the Stokes solver
The skill SHALL explain that when `anisotropy = 1` is set in the `.txt` file, the code automatically activates the Newton solver context (`input.model.Newton = 1`) because the anisotropic constitutive law requires consistent linearisation. The anisotropic viscosity is evaluated in `ViscosityConciseAniso` which rotates stress/strain-rate tensors into the director frame before applying the viscosity law.

#### Scenario: User enables anisotropy and sees Newton activated
- **WHEN** the user sets `anisotropy = 1` and notices Newton iteration is forced on
- **THEN** the skill SHALL explain that anisotropic rheology requires the Newton linearisation for convergence, so the code forces `model.Newton = 1` regardless of the user's Newton setting

### Requirement: Skill documents grid-level anisotropy fields
The skill SHALL explain the grid arrays related to anisotropy: `aniso_factor_n`, `aniso_factor_s` (anisotropy factor at cell centres and vertices), `d1_n`, `d2_n`, `d1_s`, `d2_s` (director components), `angle_n`, `angle_s` (director angle), and `FS_AR_n`, `FS_AR_s` (finite strain aspect ratio). These are interpolated from particles to the grid at each time step.

#### Scenario: User wants to visualise anisotropy
- **WHEN** the user asks what anisotropy fields are available for output
- **THEN** the skill SHALL list the grid-level fields: anisotropy factor, director angle, director components (d1, d2), and finite strain aspect ratio, at both cell centres (_n) and vertices (_s)

### Requirement: Skill documents relevant anisotropy scenarios
The skill SHALL catalogue the anisotropy-related scenarios in `SETS/`, including at minimum:
- `RiftingAnisotropy` / `RiftingCheninAniso`: anisotropic continental rifting
- `SimpleShearAnisoHomo`: homogeneous simple shear with anisotropy
- `AnisoVEVP` / `AnisoHomoVEVP`: visco-elasto-visco-plastic anisotropy tests
- `AnisoViscTest_*`: viscous anisotropy test variants (MWE, HomoRandomAngle, evolving, ellipses, debug)
- `AnisotropyDabrowski` / `AnisotropyLayer`: layered anisotropy benchmarks
- `CollisionAnisotropy` / `CollisionPolarCartesianAniso`: collision with anisotropy
- `CrustalShearBandingAnisotropy`: crustal shear band formation
- `ShearTemplateAniso` / `ShearInclusionAniso`: shear zone anisotropy tests
- `WedgeAnisotropy`: accretionary wedge with anisotropy
- `TestAnisotropy`: basic anisotropy verification

#### Scenario: User wants an example of anisotropic rifting
- **WHEN** the user wants to run an anisotropic rifting model
- **THEN** the skill SHALL recommend `RiftingAnisotropy` or `RiftingCheninAniso` as starting points, noting these implement the setup described in the anisotropic rifting paper

### Requirement: Skill references the anisotropic rifting paper methodology
The skill SHALL summarise the key methodological points from the attached anisotropic rifting paper: the use of director-based transverse isotropy to model inherited crustal fabric, the role of anisotropy factor in controlling rift symmetry/asymmetry, and how initial fabric orientation influences strain localisation patterns during extension.

#### Scenario: User asks about the scientific background
- **WHEN** the user asks for theoretical context on the anisotropy implementation
- **THEN** the skill SHALL reference the anisotropic rifting paper's formulation and explain how the code implements transverse isotropy via the director field, with viscosity depending on the angle between the director and the principal stress axes
