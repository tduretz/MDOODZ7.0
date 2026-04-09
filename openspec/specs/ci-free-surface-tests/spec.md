## ADDED Requirements

### Requirement: Topography deflection under gravity
The test suite SHALL verify that the free surface deflects when a dense body sinks under gravity, reading the `Topo/z_grid` field from HDF5 output.

#### Scenario: Surface deflects above a sinking block
- **WHEN** a simulation runs with `free_surface=1`, a denser inclusion, and gravity for 3 steps on a 41×41 grid
- **THEN** the maximum of `Topo/z_grid` at step 3 SHALL differ from the initial flat surface (topography evolved)

#### Scenario: Topography field is well-defined
- **WHEN** `free_surface=1` is enabled
- **THEN** all values in `Topo/z_grid` SHALL be finite (no NaN or Inf)
