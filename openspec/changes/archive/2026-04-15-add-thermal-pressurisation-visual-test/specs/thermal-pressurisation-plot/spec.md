## ADDED Requirements

### Requirement: ThermoElastic gnuplot shows P-vs-T with error panel

The `ThermoElastic.gnu` gnuplot script SHALL render a two-panel layout:
- **Top panel**: Pressure (Pa) vs Temperature (K) showing CHOLMOD points, PCG points, and the analytical line $P_\text{ana} = -K \alpha (T - T_0)$.
- **Bottom panel**: Relative pressure error (%) vs time step for both CHOLMOD and PCG solvers.

The output SHALL be written to `thermoelastic.png`.

#### Scenario: Visual test generates P-vs-T plot
- **WHEN** `PlotThermoElastic()` runs via the visual test executable
- **THEN** `thermoelastic.png` contains a top panel with P vs T scatter+line and a bottom panel with error bars

#### Scenario: Error values are visible in data file
- **WHEN** `thermoelastic.dat` is written by `PlotThermoElastic()`
- **THEN** the file contains columns for meanT, meanP (both solvers), analytical P, and relative error (both solvers)

### Requirement: PlotThermoElastic writes error columns

The `PlotThermoElastic()` function in `Validation.cpp` SHALL compute the relative pressure error as $|P_\text{num} - P_\text{ana}| / |P_\text{ana}| \times 100$ for each step and each solver, and write these as additional columns in `thermoelastic.dat`.

#### Scenario: Error column written for each step
- **WHEN** `PlotThermoElastic()` processes step $i$
- **THEN** the output line includes `err_cholmod` and `err_pcg` columns as percentage values

### Requirement: Thermal plots embedded in VISUAL_TESTS README

The `VISUAL_TESTS/readme.md` SHALL include sections for:
1. **ThermoElastic** — description of thermal pressurisation benchmark with embedded `thermoelastic.png`
2. **OceanicCooling** — description of half-space cooling benchmark with embedded `oceanic_cooling.png`
3. **ThermalAccuracy** — description of solver accuracy comparison with embedded `thermal_accuracy.png`

#### Scenario: README displays thermal pressurisation plot
- **WHEN** a reader opens `VISUAL_TESTS/readme.md`
- **THEN** there is a section titled "ThermoElastic" with an embedded image reference to `img/thermoelastic.png`

#### Scenario: README displays oceanic cooling plot
- **WHEN** a reader opens `VISUAL_TESTS/readme.md`
- **THEN** there is a section titled "OceanicCooling" with an embedded image reference to `img/oceanic_cooling.png`

#### Scenario: README displays thermal accuracy plot
- **WHEN** a reader opens `VISUAL_TESTS/readme.md`
- **THEN** there is a section titled "ThermalAccuracy" with an embedded image reference to `img/thermal_accuracy.png`
