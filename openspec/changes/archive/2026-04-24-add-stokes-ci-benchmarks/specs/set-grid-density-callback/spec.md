## ADDED Requirements

### Requirement: SetGridDensity callback typedef and struct field
MDLIB SHALL expose a `SetGridDensity_f` function-pointer type in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h) and a corresponding field in the `SetParticles_ff` struct, following the convention of the existing `SetDensity_f`, `SetTemperature_f`, `SetPressure_f` typedefs.

#### Scenario: Typedef has the documented signature
- **WHEN** `mdoodz.h` is included by a scenario source file
- **THEN** the public API SHALL contain `typedef double (*SetGridDensity_f)(MdoodzInput *input, Coordinates coordinates, int phase);`

#### Scenario: Struct field exists and defaults to NULL
- **WHEN** `SetParticles_ff` is constructed via a C designated initializer that does not set `.SetGridDensity`
- **THEN** the `SetGridDensity` field SHALL be `NULL` (via standard C zero-initialization of unset struct fields) and behaviour SHALL be indistinguishable from MDOODZ before this change

### Requirement: Per-iteration grid-override semantics
When `SetGridDensity != NULL`, the grid-level density arrays `mesh->rho_n` (cell centres) and `mesh->rho_s` (vertices) SHALL have their final values (after the existing `UpdateDensity` and `NonNewtonianViscosityGrid` pipelines complete) replaced by the callback's return value at each cell / vertex coordinate. This replacement SHALL happen on every rheology pass where the mesh-density computation runs, not once at initialization. The replacement loop SHALL be guarded by a null check so scenarios without the callback are unaffected.

#### Scenario: Grid density is overridden at cell centres
- **WHEN** a simulation runs with `SetGridDensity != NULL` and at least one rheology pass completes
- **THEN** for every cell-centre index `c` with coordinates `(x_c, z_c)` and representative phase `p`, `mesh->rho_n[c]` SHALL equal the result of `SetGridDensity(input, {x_c, z_c}, p)`

#### Scenario: Grid density is overridden at vertices
- **WHEN** a simulation runs with `SetGridDensity != NULL` and at least one rheology pass completes
- **THEN** for every vertex index `v` with coordinates `(x_v, z_v)` and representative phase `p`, `mesh->rho_s[v]` SHALL equal the result of `SetGridDensity(input, {x_v, z_v}, p)`

#### Scenario: Callback is never invoked for null pointer
- **WHEN** a scenario uses the default `SetGridDensity = NULL`
- **THEN** `mesh->rho_n` and `mesh->rho_s` SHALL contain exactly the values computed by the pre-existing rheology density pipeline, and the callback pointer SHALL NOT be dereferenced

### Requirement: Startup validation rejects conflict with active density models
When `SetGridDensity != NULL`, MDOODZ SHALL validate at startup (after the per-phase material-property parsing loop, in `Main_DOODZ.c::RunMDOODZ`) that for every declared phase, `density_model = 0` (the default constant-density branch). Any non-zero `density_model` on any phase while the callback is active SHALL cause MDOODZ to log a clear error and exit before any solver work begins.

#### Scenario: Clean configuration is accepted
- **WHEN** a scenario sets `SetGridDensity != NULL` and every phase in the `.txt` file has `density_model = 0`
- **THEN** startup validation SHALL pass and the simulation SHALL proceed normally

#### Scenario: Mixed configuration is rejected with a clear error
- **WHEN** a scenario sets `SetGridDensity != NULL` and any phase has a non-zero `density_model`
- **THEN** MDOODZ SHALL emit an `LOG_ERR` message naming the offending phase and the value, and call `exit` with a non-zero status before running any solver work

#### Scenario: Null callback does not trigger validation
- **WHEN** a scenario uses the default `SetGridDensity = NULL` (including every pre-existing scenario)
- **THEN** the new validation block SHALL do nothing — no error, no exit, behaviour is unchanged

### Requirement: Backwards compatibility
This capability SHALL NOT alter the behaviour of any existing scenario that does not set `SetGridDensity`. All pre-existing CI tests, visual tests, and `SETS/*.c` scenarios SHALL continue to run without source or parameter-file modification.

#### Scenario: Existing test suite passes without modification
- **WHEN** the MDOODZ CI suite is run after this capability is implemented, with no edits to pre-existing test parameter files or scenario source
- **THEN** every pre-existing GTest case SHALL pass with outputs equivalent to the pre-change baseline

### Requirement: SetGridDensity documentation
The `SetGridDensity` callback SHALL be documented in two locations: adjacent to the typedef in [MDLIB/include/mdoodz.h](MDLIB/include/mdoodz.h), and as a new pitfall entry in [skill-testing-guide/SKILL.md](.claude/skills/skill-testing-guide/SKILL.md).

#### Scenario: Header comment explains override semantics and the validation rule
- **WHEN** a developer reads the `SetGridDensity_f` typedef in `mdoodz.h`
- **THEN** the adjacent comment block SHALL state that (a) the callback overrides `mesh->rho_n` / `mesh->rho_s` on every rheology pass, (b) startup validation refuses the combination of the callback with any non-zero `density_model`, (c) the callback is intended for benchmarks and prescribed-density studies

#### Scenario: Skill testing guide has a dedicated pitfall entry
- **WHEN** [skill-testing-guide/SKILL.md](.claude/skills/skill-testing-guide/SKILL.md) is read
- **THEN** it SHALL contain a pitfall entry explaining that `SetGridDensity` is a whole-simulation override, that the startup validator rejects mixing with any non-zero `density_model`, and that benchmarks using the callback (e.g., SolCx) therefore require `density_model = 0` on every phase
