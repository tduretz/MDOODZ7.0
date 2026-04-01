## ADDED Requirements

### Requirement: Pure thermal diffusion test
The CI test suite SHALL include a GTest case that runs a thermal-only simulation (no mechanical deformation) with a known initial temperature distribution and verifies diffusion behaviour.

#### Scenario: Gaussian perturbation decays
- **WHEN** a simulation runs with `mechanical=0`, `thermal=1`, a homogeneous medium with conductivity k and heat capacity Cp, and an initial Gaussian temperature perturbation centred in the domain, for multiple time steps
- **THEN** the peak temperature SHALL decrease monotonically between steps and the final temperature field SHALL be smoother than the initial field (max T at final step < max T at initial step)

#### Scenario: Steady state with fixed boundary temperatures
- **WHEN** a simulation runs with `mechanical=0`, `thermal=1`, Dirichlet temperature BCs (T_top=0°C, T_bottom=1300°C), for enough steps to reach steady state
- **THEN** the temperature profile SHALL be approximately linear with depth (max deviation from linear < 1% of ΔT)

### Requirement: Radiogenic heat production test
The CI test suite SHALL include a GTest case verifying that nonzero `Qr` increases temperature over time.

#### Scenario: Temperature increases with radiogenic heating
- **WHEN** a simulation runs with `thermal=1`, `Qr=1e-6` W/m³ for all phases, insulating (Neumann) side BCs and fixed top/bottom BCs, for multiple time steps
- **THEN** the mean temperature in the domain SHALL increase between the first and last time step

### Requirement: Thermal boundary condition types test
The CI test suite SHALL include a GTest case exercising Dirichlet thermal BCs on top and bottom boundaries.

#### Scenario: Fixed temperature BCs are respected
- **WHEN** a simulation runs with `thermal=1` and Dirichlet BCs setting T=273.15 K at the top and T=1600 K at the bottom
- **THEN** the temperature values at the top and bottom cell rows SHALL match the prescribed BC values within 0.1%
