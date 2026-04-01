## MODIFIED Requirements

### Requirement: Steady-state geotherm linearity check
The existing SteadyStateGeotherm test SHALL additionally verify that the temperature gradient is approximately linear (consistent with no internal heat sources and fixed-T boundaries).

#### Scenario: Temperature gradient is approximately linear
- **WHEN** a SteadyStateGeotherm simulation reaches near steady state after 20 steps
- **THEN** `mean(T)` SHALL be approximately the average of T_top and T_bot (within 20%)

### Requirement: Radiogenic heat quantitative check
The existing RadiogenicHeat test SHALL additionally verify the temperature increase is order-of-magnitude consistent with ΔT ≈ Qr·Δt/(ρCp).

#### Scenario: Temperature increase matches analytical estimate
- **WHEN** RadiogenicHeat runs with known Qr, Δt, ρ, Cp
- **THEN** the mean temperature increase SHALL be within a factor of 5 of the analytical Qr·Nt·Δt/(ρCp)
