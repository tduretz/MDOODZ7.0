## ADDED Requirements

### Requirement: Skill explains the paired-file setup pattern
The skill SHALL document that every MDOODZ simulation is defined by two files: a C callback file (`.c`) in `SETS/` and a plain-text parameter file (`.txt`) in `SETS/`. The skill SHALL explain that the `.c` file contains callback functions invoked by the solver to set initial conditions, while the `.txt` file contains all runtime parameters (scaling, resolution, switches, material properties).

#### Scenario: User asks how to create a new simulation
- **WHEN** the user asks how to set up a new model
- **THEN** the skill SHALL explain the `.c`/`.txt` pair pattern and provide a step-by-step workflow: copy an existing pair, rename, modify callbacks in `.c`, tune parameters in `.txt`, build with `make build SET=<Name>`, run with `make run SET=<Name>`

### Requirement: Skill documents all callback function signatures
The skill SHALL list and explain every callback function signature that a setup `.c` file can implement, including at minimum:
- `SetPhase(MdoodzInput *input, Coordinates coordinates)` → returns `int` phase ID
- `SetTemperature(MdoodzInput *input, Coordinates coordinates)` → returns `double` temperature
- `SetGrainSize(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double` grain size
- `SetSurfaceZCoord(MdoodzInput *input, double x_coord)` → returns `double` topography elevation
- `SetHorizontalVelocity(MdoodzInput *input, Coordinates coordinates)` → returns `double` Vx
- `SetVerticalVelocity(MdoodzInput *input, Coordinates coordinates)` → returns `double` Vz
- `SetBCPType(MdoodzInput *input, POSITION position)` → returns `char` pressure BC type
- `SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates)` → returns `SetBC`
- `SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates)` → returns `SetBC`
- `SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double gridTemperature)` → returns `SetBC`
- `SetDualPhase(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `int`
- `SetPorosity(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double`
- `SetDensity(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double`
- `SetXComponent(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double`
- `SetPressure(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double`
- `SetNoise(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `double`
- `SetAnisoAngle(MdoodzInput *input, Coordinates coordinates, int phase, double predefined)` → returns `double`
- `SetDefGrad(MdoodzInput *input, Coordinates coordinates, int phase)` → returns `Tensor2D`
- `MutateInput(MdoodzInput *input)` → void, for programmatic parameter overrides

#### Scenario: User wants to set initial phase geometry
- **WHEN** the user asks how to define material phases spatially
- **THEN** the skill SHALL explain `SetPhase` with examples showing coordinate-based conditionals (e.g., layered crust/mantle, inclusion geometries using `IsEllipseCoordinates` / `IsRectangleCoordinates`)

#### Scenario: User wants to set boundary conditions
- **WHEN** the user asks about boundary conditions
- **THEN** the skill SHALL explain the `POSITION` enum (`NE`, `NW`, `SE`, `SW`, `N`, `S`, `W`, `E`, `INTERNAL`, `free_surface`), the `SetBC` struct (value + type), and the BC type codes (0 = Dirichlet physical, 11 = Dirichlet non-physical, 2/13 = Neumann, -2 = periodic, -1 = interior, 30 = air)

### Requirement: Skill documents the Coordinates struct
The skill SHALL explain the `Coordinates` struct fields: `x` (horizontal position), `z` (vertical position), `l` (distance parameter), `k` (index parameter), and clarify that coordinates are in scaled (non-dimensional) units within callbacks.

#### Scenario: User is confused about coordinate units in callbacks
- **WHEN** the user encounters unexpected coordinate values in a callback
- **THEN** the skill SHALL explain that coordinates are non-dimensionalised by the scaling length `L`, and show how to convert to physical units using `input->scaling.L`

### Requirement: Skill documents the MdoodzSetup struct
The skill SHALL explain the `MdoodzSetup` struct that groups all callback function pointers (`BuildInitialTopography`, `SetParticles`, `SetBCs`, `MutateInput`) and how the `main()` function in each `.c` file populates this struct and calls `RunMDOODZ`.

#### Scenario: User wants to understand the entry point
- **WHEN** the user asks how a simulation starts
- **THEN** the skill SHALL show the canonical `main()` pattern: allocate `MdoodzSetup`, assign callback pointers, call `RunMDOODZ(inputFileName, &setup)`

### Requirement: Skill documents the .txt configuration file sections
The skill SHALL document all configuration sections in the `.txt` parameter file, including at minimum: `RESTART`, `OUTPUT FILES`, `INPUT FILE FOR PARTICLES`, `SCALES`, `SPACE-TIME`, `SWITCHES`, `SETUP DEPENDANT`, `GRAVITY`, and `PHASE PROPERTIES`. The skill SHALL list key parameters in each section with their meaning and units.

#### Scenario: User wants to know what parameters are available
- **WHEN** the user asks about `.txt` file options
- **THEN** the skill SHALL provide a section-by-section reference with parameter names, types, units, and brief descriptions

#### Scenario: User wants to configure material phases
- **WHEN** the user asks how to set rheological properties for a phase
- **THEN** the skill SHALL explain the `PHASE PROPERTIES` section format: `Nb_phases`, followed by per-phase blocks keyed by `ID`, with properties like `rho`, `G`, `Cp`, `k`, `C`, `phi`, `cstv`, `pwlv`, `linv`, `gbsv`, `expv`, `aniso_factor`, `aniso_angle`

### Requirement: Skill provides helper geometry utilities
The skill SHALL document the built-in geometry helpers `IsEllipseCoordinates(coordinates, ellipse, scalingL)` and `IsRectangleCoordinates(coordinates, rectangle, scalingL)` with their parameter structs (`Ellipse`: centreX, centreZ, radiusX, radiusZ, angle; `Rectangle`: centreX, centreZ, sizeX, sizeZ, angle).

#### Scenario: User wants to place an elliptical inclusion
- **WHEN** the user needs to define an elliptical phase region
- **THEN** the skill SHALL show how to define an `Ellipse` struct and use `IsEllipseCoordinates` inside `SetPhase` to assign a phase ID within that region

### Requirement: Skill provides built-in BC templates
The skill SHALL document the pre-built boundary condition templates: `SetPureShearBCVx`, `SetPureShearBCVz`, `SetSimpleShearBCVx`, `SetSimpleShearBCVz`, `SetPureOrSimpleShearBCVx`, `SetPureOrSimpleShearBCVz`.

#### Scenario: User wants standard pure shear boundary conditions
- **WHEN** the user needs pure shear BCs
- **THEN** the skill SHALL show they can assign `.SetBCVx = SetPureShearBCVx` and `.SetBCVz = SetPureShearBCVz` directly without writing custom callbacks
