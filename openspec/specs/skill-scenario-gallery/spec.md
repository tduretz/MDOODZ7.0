## ADDED Requirements

### Requirement: Skill organises scenarios by geodynamic application
The skill SHALL catalogue all `.c` scenario files in `SETS/` grouped by geodynamic application category. Each category SHALL include a brief description of the physical process being modelled and list all scenarios in that category.

#### Scenario: User asks what scenarios are available
- **WHEN** the user asks what models they can run
- **THEN** the skill SHALL present the full categorised list with brief descriptions of each category's purpose

### Requirement: Skill documents rifting scenarios
The skill SHALL list and describe rifting scenarios:
- `RiftingBasic`: minimal continental rifting setup
- `RiftingChenin`: multi-layer lithosphere rifting with Chenin-style rheology
- `RiftingCheninAniso`: rifting with anisotropic crust (anisotropic rifting paper setup)
- `RiftingAnisotropy`: anisotropic rifting variant
- `RiftingMelting`: rifting with partial melting
- `RiftingRoman`: rifting configuration by Roman
- `NeckingReview` / `NeckingRoman`: necking instability studies

#### Scenario: User wants to model continental rifting
- **WHEN** the user asks for a rifting starting point
- **THEN** the skill SHALL recommend `RiftingChenin` as a well-tested baseline with realistic multi-layer rheology (upper crust, lower crust, lithospheric mantle, asthenosphere), and suggest `RiftingAnisotropy` if anisotropic fabric effects are needed

### Requirement: Skill documents subduction scenarios
The skill SHALL list and describe subduction scenarios:
- `SubductionBenchmarkCase2`: benchmark subduction problem
- `SubductionMelting`: subduction with melting
- `ThanushikaSubduction`: specific subduction configuration
- `AnneloreSubduction`: subduction variant
- `DoubleSubduction`: two-sided subduction
- `ThermalStructureSubduction`: thermal structure analysis

#### Scenario: User wants to model subduction
- **WHEN** the user asks about subduction setups
- **THEN** the skill SHALL recommend `SubductionBenchmarkCase2` as a validated starting point

### Requirement: Skill documents collision and orogenesis scenarios
The skill SHALL list and describe collision scenarios:
- `CollisionZone`: generic collision zone
- `CollisionIra`: collision variant
- `CollisionAnisotropy`: collision with anisotropic rheology
- `CollisionPolarCartesian` / `CollisionPolarCartesianAniso`: polar-Cartesian collision formulations

#### Scenario: User wants to model mountain building
- **WHEN** the user asks about collision or orogeny setups
- **THEN** the skill SHALL recommend `CollisionZone` as a general starting point and note `CollisionAnisotropy` for fabric-sensitive orogenesis

### Requirement: Skill documents shear deformation scenarios
The skill SHALL list and describe shear scenarios:
- `ShearTemplate`: basic simple shear template
- `ShearTemplateAniso`: simple shear with anisotropy
- `ShearInclusionAniso`: shear around anisotropic inclusion
- `ShearHeatingDuretz14`: shear heating benchmark (Duretz 2014)
- `ShearBoxLaeti`: shear box variant
- `SimpleShearAnisoHomo`: homogeneous simple shear with anisotropy
- `StressBC_ShearTemplate` / `StressBC_FreeSurf`: stress-driven shear

#### Scenario: User wants a simple test case
- **WHEN** the user wants a minimal shear test
- **THEN** the skill SHALL recommend `ShearTemplate` as the simplest starting point for testing rheological behaviour

### Requirement: Skill documents anisotropy test scenarios
The skill SHALL list and describe anisotropy-specific test scenarios:
- `AnisoVEVP` / `AnisoHomoVEVP`: visco-elasto-visco-plastic anisotropy
- `AnisoViscTest_MWE`: minimal working example for viscous anisotropy
- `AnisoViscTest_HomoRandomAngle`: random initial fabric orientations
- `AnisoViscTest_evolv`: evolving anisotropy factor
- `AnisoViscTest_Ellipses`: elliptical inclusions with anisotropy
- `AnisoViscTest_Debug`: debug configuration
- `AnisotropyDabrowski` / `AnisotropyLayer`: layered anisotropy benchmarks
- `TestAnisotropy`: basic verification
- `CrustalShearBandingAnisotropy`: crustal shear banding
- `WedgeAnisotropy`: accretionary wedge

#### Scenario: User wants to verify anisotropy implementation
- **WHEN** the user wants to test anisotropy features
- **THEN** the skill SHALL recommend `AnisoViscTest_MWE` as minimal working example and `AnisotropyDabrowski` as a published benchmark

### Requirement: Skill documents thermal and chemical scenarios
The skill SHALL list and describe thermal/chemical scenarios:
- `ThermalDiffusion`: pure thermal diffusion test
- `OceanicCooling`: oceanic lithosphere cooling
- `ThermoElastic`: thermo-elastic coupling
- `ChemicalDiffusion`: chemical diffusion test
- `ChemicalProduction`: chemical production test
- `CrystallisationStuewe1995`: crystallisation benchmark

#### Scenario: User wants to test thermal features
- **WHEN** the user asks about thermal-only tests
- **THEN** the skill SHALL recommend `ThermalDiffusion` for verifying thermal solver correctness

### Requirement: Skill documents benchmark and validation scenarios
The skill SHALL list and describe validation scenarios:
- `ViscousInclusion` / `ViscousInclusionTestSuite`: viscous inclusion benchmark
- `ViscousSinker`: falling block benchmark
- `TopoBenchCase1`: topography benchmark
- `FreeSurfaceBenchmarkFEM`: free surface validation
- `ElasticFlexureLineLoad`: elastic flexure under load
- `VEP_Duretz18`: visco-elasto-plastic benchmark (Duretz 2018)
- `Popov2025_Tensile_VEVP` / `Popov2025_Pureshear_VEVP`: Popov 2025 benchmarks
- `QuartzCoesite` / `QuartzCoesite_Simple` / `QuartzCoesite_paper`: phase transformation benchmarks

#### Scenario: User wants to validate the code
- **WHEN** the user asks about benchmarks
- **THEN** the skill SHALL recommend `ViscousInclusion` (simple analytical solution), `ViscousSinker` (Stokes flow), and `VEP_Duretz18` (full rheology) as a validation suite

### Requirement: Skill documents magmatic and melting scenarios
The skill SHALL list and describe melting/magmatic scenarios:
- `MeltingOverpressure`: melt overpressure effects
- `CoolingChamber`: magma chamber cooling
- `PressurizedMagmaChamber`: pressurised chamber
- `TCMagmaticSystem`: thermo-chemical magmatic system

#### Scenario: User wants to model magmatic processes
- **WHEN** the user asks about melting or magma chambers
- **THEN** the skill SHALL list available magmatic scenarios and note that `melting = 1` must be enabled in the `.txt` file

### Requirement: Skill documents miscellaneous scenarios
The skill SHALL list remaining scenarios:
- `Multilayers` / `MLPS_Ellipses`: multilayer folding
- `PinchSwellGSE`: pinch-and-swell with grain size evolution
- `PlasticityLayers`: plasticity in layered media
- `StrainLocalization_SH_GSE`: strain localisation with shear heating and grain size
- `FaultSmearing`: fault zone processes
- `ChristmasTree`: decorative/fun test
- `RiverTom` / `RiverTomLowerCrust`: river incision models
- `Shrinking*` / `Anais*`: shrinking/loading experiments

#### Scenario: User asks about a scenario not in the main categories
- **WHEN** the user asks about a specific miscellaneous scenario
- **THEN** the skill SHALL provide its category and brief description based on the filename and setup pattern
