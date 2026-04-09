## ADDED Requirements

### Requirement: Temperature field colour map
A gnuplot script SHALL produce a 2D temperature-field plot from the Blankenbach simulation HDF5 output using a heat-palette colormap (blue-to-red or similar thermal palette) with `pm3d` or `image` mode.

#### Scenario: Temperature field rendered as PNG
- **WHEN** the gnuplot script is run on a completed simulation output
- **THEN** it SHALL produce a PNG file showing the 2D temperature distribution in the unit square with a colour bar labelled "Temperature"

### Requirement: Velocity vector overlay
The gnuplot visualization SHALL overlay velocity vectors on the temperature field to show the convective flow pattern.

#### Scenario: Vectors show convection cell
- **WHEN** the plot is generated from steady-state output
- **THEN** velocity arrows SHALL be visible showing the characteristic single-cell convection pattern (upwelling on one side, downwelling on the other)

### Requirement: Standalone script
The gnuplot script SHALL be runnable independently (not part of GTest CI) and SHALL include comments documenting how to extract data from HDF5 and run the script.

#### Scenario: Manual execution
- **WHEN** a user runs `gnuplot TESTS/BlankenBench/plot_convection.gp` from the TESTS/ build directory after a simulation
- **THEN** the script SHALL produce the visualization without errors, reading data from the simulation output directory
