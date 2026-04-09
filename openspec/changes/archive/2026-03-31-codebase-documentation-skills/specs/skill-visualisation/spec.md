## ADDED Requirements

### Requirement: Skill documents HDF5 output reading
The skill SHALL explain that MDOODZ writes simulation output in HDF5 format (via `HDF5Output.c`), with files named `Output<step>.gzip.h5` in the run directory. Output frequency is controlled by `writer_step` in `.txt`. The skill SHALL explain the HDF5 file structure: grid arrays (velocity, pressure, temperature, stress, strain rate, viscosity, phase percentages), particle data (if `writer_markers = 1`), and metadata.

#### Scenario: User asks how to access simulation results
- **WHEN** the user asks where output data is stored
- **THEN** the skill SHALL explain the HDF5 output naming convention, the `writer_step` parameter, and that Julia's `read_bin.jl` or any HDF5 reader can open the files

### Requirement: Skill documents the Julia/Makie visualisation framework
The skill SHALL explain that the primary visualisation toolkit is in `JuliaVisualisation/` using Julia + CairoMakie/GLMakie. The main script is `Main_Visualisation_Makie_MD7.jl`. The skill SHALL explain the Julia project setup: `Project.toml` defines dependencies, activate with `] activate JuliaVisualisation`.

#### Scenario: User wants to get started with visualisation
- **WHEN** the user asks how to visualise results
- **THEN** the skill SHALL provide: (1) install Julia, (2) activate the JuliaVisualisation project, (3) instantiate dependencies, (4) edit the script's path/step variables, (5) run the script

### Requirement: Skill documents available visualisation field types
The skill SHALL list all fields available for plotting, including at minimum:
- **Phases**: material phase distribution (compositional field)
- **Temperature**: thermal field in °C or K
- **Velocity**: Vx, Vz components and magnitude
- **Pressure**: total pressure, dynamic pressure, lithostatic pressure
- **Stress**: deviatoric stress components (τ_xx, τ_zz, τ_xz) and second invariant τ_II
- **Strain rate**: deviatoric strain rate components (ε̇_xx, ε̇_zz, ε̇_xz) and second invariant ε̇_II
- **Viscosity**: effective viscosity (η_s at vertices, η_n at centres)
- **Density**: ρ at centres and vertices
- **Grain size**: d field
- **Fabric/Anisotropy**: anisotropy angle, anisotropy factor, director components, finite strain aspect ratio
- **Topography**: free surface elevation

#### Scenario: User wants to plot strain rate
- **WHEN** the user asks to visualise strain localisation
- **THEN** the skill SHALL show how to select the strain rate second invariant field, recommend log-scale colour mapping, and suggest overlaying phase contours for context

### Requirement: Skill documents overlay and multi-panel options
The skill SHALL explain the overlay capabilities: phase contours over any field, temperature contours, velocity arrows, topography profile, and fabric tick marks for anisotropy. The skill SHALL also explain multi-panel layouts available in `Main_Visualisation_Makie_MD7_3x3.jl`.

#### Scenario: User wants a composite figure
- **WHEN** the user asks to create a multi-panel figure with different fields
- **THEN** the skill SHALL explain the 3x3 layout script and how to configure which fields appear in each panel

### Requirement: Skill documents marker trajectory analysis
The skill SHALL explain the marker history tracking tools: `Read_Marker_History_MD7.jl` reads particle trajectory data, `Select_Markers_MD7.jl` filters markers by region/phase, and `T_Select_Markers_MD7.jl` provides temperature-path analysis. These require `track_T_P_x_z = 1` in the `.txt` file.

#### Scenario: User wants P-T paths for specific particles
- **WHEN** the user asks about pressure-temperature paths
- **THEN** the skill SHALL explain: enable `track_T_P_x_z = 1`, run the simulation, then use `Select_Markers_MD7.jl` to filter markers by initial position/phase and `T_Select_Markers_MD7.jl` to extract and plot their P-T trajectories

### Requirement: Skill documents batch export workflows
The skill SHALL explain how to export visualisations as PNG files in batch mode (iterating over time steps), and how to create animations from exported frames. The skill SHALL also mention the `_VisualTests/` directory for visual regression testing.

#### Scenario: User wants to create an animation
- **WHEN** the user asks about making a movie of the simulation
- **THEN** the skill SHALL explain: loop over time steps in the Julia script, save each frame as PNG with consistent colour scale, then assemble frames using ffmpeg or similar tool
