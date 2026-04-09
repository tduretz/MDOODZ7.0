## ADDED Requirements

### Requirement: Drucker-Prager yield test
The CI test suite SHALL include a GTest case that verifies the Drucker-Prager yield criterion activates when stress exceeds the yield envelope.

#### Scenario: Plastic yielding bounds stress
- **WHEN** a simulation runs with plasticity enabled (`plast=1`), cohesion C=20e6 Pa, friction angle φ=30°, on a small grid with a background strain rate high enough to exceed yield
- **THEN** the stress second invariant (τ_II) at every cell centre SHALL NOT exceed the Drucker-Prager yield stress `C·cos(φ) + P·sin(φ)` by more than 1%

#### Scenario: No yielding below yield envelope
- **WHEN** a simulation runs with the same plasticity parameters but a low background strain rate (1e-18 s⁻¹) such that elastic/viscous stress stays below yield
- **THEN** the plastic strain rate invariant (`/Centers/eII_pl`) SHALL be zero (or effectively zero, <1e-30) at all cell centres

### Requirement: Strain softening evolution test
The CI test suite SHALL include a GTest case that verifies cohesion and friction angle evolve with accumulated plastic strain when softening is enabled.

#### Scenario: Cohesion decreases with plastic strain
- **WHEN** a simulation runs for 3 time steps with `coh_soft=1`, `C=20e6`, `C_end=5e6`, `pls_start=0`, `pls_end=1`
- **THEN** the cohesion field (`/Centers/cohesion`) SHALL contain values less than the initial C=20e6 in regions where plastic strain has accumulated

### Requirement: Stress limiter enforcement test
The CI test suite SHALL include a GTest case verifying that the stress limiter (`Slim`) caps deviatoric stress.

#### Scenario: Stress does not exceed Slim
- **WHEN** a simulation runs with `Slim=100e6` Pa for a phase
- **THEN** the stress second invariant (τ_II) SHALL NOT exceed `Slim` at any cell centre for that phase
