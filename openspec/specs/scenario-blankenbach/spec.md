## ADDED Requirements

### Requirement: Blankenbach Case 1a scenario .c file
The SETS/ directory SHALL contain a `BlankenBench.c` file implementing the Blankenbach Case 1a setup with callbacks for SetPhase (single phase), SetTemperature (conductive profile plus cosine perturbation to seed convection), SetBCVx/SetBCVz (free-slip on all walls), and SetBCT (Dirichlet T=0 top, T=1 bottom, Neumann on sides).

#### Scenario: Single-phase homogeneous domain
- **WHEN** SetPhase is called for any coordinate in the domain
- **THEN** it SHALL return phase 0 (single isoviscous material)

#### Scenario: Initial temperature with perturbation
- **WHEN** SetTemperature is called
- **THEN** it SHALL return a linear conductive profile T(z) = 1 - (z - z_bot)/H plus a small cosine perturbation to break symmetry and seed a single convection cell

#### Scenario: Free-slip velocity BCs
- **WHEN** SetBCVx or SetBCVz is called for boundary positions
- **THEN** free-slip conditions SHALL be applied on all four walls (zero normal velocity, zero tangential stress)

#### Scenario: Dirichlet thermal BCs on top/bottom
- **WHEN** SetBCT is called for north (top) positions
- **THEN** it SHALL return type=1 (Dirichlet) with T = 0 (scaled)
- **WHEN** SetBCT is called for south (bottom) positions
- **THEN** it SHALL return type=1 (Dirichlet) with T = 1 (scaled)
- **WHEN** SetBCT is called for east/west positions
- **THEN** it SHALL return type=0 (Neumann, zero heat flux)

### Requirement: Blankenbach parameter file
A `.txt` parameter file SHALL configure the Blankenbach Case 1a benchmark with: non-dimensional unit scaling (η=1, L=1, V=1, T=1), Ra = 10⁴ via gravity parameter, constant viscosity, thermal solver enabled, no elasticity/plasticity/free-surface, and sufficient time steps to reach steady state.

#### Scenario: Rayleigh number equals 10⁴
- **WHEN** the scaling parameters and material properties are combined
- **THEN** the effective Rayleigh number Ra = αρgΔTL³/(κη) SHALL equal 10⁴

#### Scenario: Sufficient steps for steady state
- **WHEN** the simulation runs for the configured number of time steps
- **THEN** the thermal convection SHALL have reached a steady single-cell pattern
