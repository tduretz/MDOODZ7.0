## ADDED Requirements

### Requirement: Oceanic cooling analytical verification
The CI test suite SHALL include a GTest case that runs a half-space cooling simulation and compares the final vertical temperature profile against the analytical erf solution.

#### Scenario: CHOLMOD solver matches analytical erf profile
- **WHEN** a simulation runs with `mechanical=0`, `thermal=1`, `thermal_solver=0` (CHOLMOD), uniform initial temperature T_m, Dirichlet T=0 at the top surface, Dirichlet T=T_m at the bottom, Neumann (zero flux) on sides, for Nt time steps with constant dt
- **THEN** the relative L2 error of the mid-column vertical temperature profile against the analytical solution T(z,t) = T_m · erf(−z / √(4tκ)) with κ = k/(ρ·Cp) SHALL be less than 5e-2

#### Scenario: PCG solver matches analytical erf profile
- **WHEN** the same simulation runs with `thermal_solver=1` (PCG) and all other parameters identical
- **THEN** the relative L2 error against the analytical erf solution SHALL be less than 5e-2

#### Scenario: CHOLMOD and PCG produce equivalent results
- **WHEN** both CHOLMOD and PCG simulations complete on the same oceanic cooling setup
- **THEN** the absolute difference of their L2 errors against the analytical solution SHALL be less than 1e-6

### Requirement: Oceanic cooling parameter files
The test suite SHALL provide `.txt` parameter files for the oceanic cooling scenario in `TESTS/Thermal/`.

#### Scenario: CHOLMOD parameter file exists
- **WHEN** the test harness looks for `Thermal/OceanicCooling.txt`
- **THEN** the file SHALL exist and configure `mechanical=0`, `thermal=1`, `thermal_solver=0` (or absent, defaulting to 0), with `writer_subfolder=OceanicCooling`

#### Scenario: PCG parameter file exists
- **WHEN** the test harness looks for `Thermal/OceanicCoolingPCG.txt`
- **THEN** the file SHALL exist and be identical to the CHOLMOD variant except for `thermal_solver=1` and `writer_subfolder=OceanicCoolingPCG`

### Requirement: Oceanic cooling gnuplot visualisation
The repository SHALL include a gnuplot script that plots the numerical vs analytical temperature profiles for visual inspection.

#### Scenario: Gnuplot script produces a PNG
- **WHEN** a user runs the gnuplot script with paths to the CHOLMOD and PCG CSV outputs
- **THEN** the script SHALL produce a PNG showing depth on the y-axis, temperature on the x-axis, with three curves: analytical (dashed), CHOLMOD (solid), and PCG (dotted), and a legend identifying each
