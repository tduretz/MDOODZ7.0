## ADDED Requirements

### Requirement: Thermo-elastic coupling verification
The CI test suite SHALL include a GTest case that runs a thermo-elastic simulation with linearly increasing temperature and verifies the resulting pressure response.

#### Scenario: CHOLMOD solver produces correct thermal pressure
- **WHEN** a simulation runs with `elastic=1`, `compressible=1`, `fix_temperature=1`, `thermal_solver=0` (CHOLMOD), a FixTemperature callback that returns T = 273.15/scaling.T + 0.1·time (linearly increasing), for Nt time steps
- **THEN** the mean pressure change in the domain SHALL be consistent with the analytical estimate ΔP ≈ −K·α·ΔT within 5% relative error, where K is the bulk modulus, α is the thermal expansion coefficient, and ΔT is the imposed temperature change

#### Scenario: PCG solver produces correct thermal pressure
- **WHEN** the same simulation runs with `thermal_solver=1` (PCG) and all other parameters identical
- **THEN** the mean pressure change SHALL satisfy the same analytical tolerance (5% relative error)

#### Scenario: CHOLMOD and PCG produce equivalent thermo-elastic results
- **WHEN** both CHOLMOD and PCG simulations complete on the same thermo-elastic setup
- **THEN** the relative L2 difference of the full 2D temperature fields SHALL be less than 1e-6
- **AND** the relative difference of mean pressure between the two runs SHALL be less than 1e-4

### Requirement: Thermo-elastic parameter files
The test suite SHALL provide `.txt` parameter files for the thermo-elastic scenario in `TESTS/Thermal/`.

#### Scenario: CHOLMOD parameter file exists
- **WHEN** the test harness looks for `Thermal/ThermoElastic.txt`
- **THEN** the file SHALL exist and configure `elastic=1`, `compressible=1`, `fix_temperature=1`, `thermal_solver=0` (or absent), with `writer_subfolder=ThermoElastic`

#### Scenario: PCG parameter file exists
- **WHEN** the test harness looks for `Thermal/ThermoElasticPCG.txt`
- **THEN** the file SHALL exist and be identical to the CHOLMOD variant except for `thermal_solver=1` and `writer_subfolder=ThermoElasticPCG`

### Requirement: Thermo-elastic gnuplot visualisation
The repository SHALL include a gnuplot script for visual inspection of thermo-elastic test results.

#### Scenario: Gnuplot script produces a PNG
- **WHEN** a user runs the gnuplot script with paths to the CHOLMOD and PCG CSV outputs
- **THEN** the script SHALL produce a PNG with two panels: (1) temperature vs time showing the linear ramp, and (2) mean pressure vs time showing CHOLMOD, PCG, and analytical curves
