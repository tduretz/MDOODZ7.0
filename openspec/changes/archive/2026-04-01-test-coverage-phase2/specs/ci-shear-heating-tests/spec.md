## ADDED Requirements

### Requirement: Viscous dissipation heating test
The test suite SHALL verify that shear heating (`shear_heating=1`) raises the temperature when mechanical work is done, following H_s = 2η·ε̇²_II.

#### Scenario: Temperature increases under shear heating
- **WHEN** a simulation runs with `shear_heating=1`, constant viscosity, high strain rate, and `thermal=1`
- **THEN** `max(T)` at the final step SHALL be greater than `max(T)` at the initial step

#### Scenario: Temperature increase is order-of-magnitude correct
- **WHEN** shear heating is active with known η, ε̇, Δt, ρ, Cp
- **THEN** the temperature increase ΔT SHALL be within a factor of 2 of the analytical estimate 2η·ε̇²·Δt/(ρCp)
