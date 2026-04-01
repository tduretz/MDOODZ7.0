## ADDED Requirements

### Requirement: Pure thermal diffusion visual regression test
The visual test suite SHALL include a thermal-only scenario with a known analytical or reference solution.

#### Scenario: Thermal diffusion temperature field matches reference
- **WHEN** a pure thermal diffusion scenario runs with `mechanical=0`, `thermal=1`, fixed top/bottom T BCs, homogeneous medium, for ~10 time steps at moderate resolution
- **THEN** the relative L2 error between the computed temperature field and the reference field SHALL be less than 1e-5

#### Scenario: Thermal diffusion convergence to steady state
- **WHEN** the scenario runs for enough steps to approach steady state
- **THEN** the temperature profile SHALL converge toward the analytical linear geotherm (max pointwise error < 1% of ΔT)

### Requirement: Quartz-coesite phase transition visual regression test
The visual test suite SHALL include a phase transition scenario based on the existing `QuartzCoesite.c` setup with reference comparison.

#### Scenario: Phase transition density field matches reference
- **WHEN** a quartz-coesite phase transition scenario runs to completion
- **THEN** the relative L2 error of the density field (`/Centers/rho_n`) against the reference SHALL be less than 1e-4

#### Scenario: Phase boundary location matches reference
- **WHEN** the scenario completes
- **THEN** the phase composition field (`/VizGrid/compo`) SHALL match the reference (phase boundaries within ±2 cells of reference position)
