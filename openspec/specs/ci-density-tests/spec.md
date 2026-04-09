## ADDED Requirements

### Requirement: Thermal expansion density test
The test suite SHALL verify that density decreases with increasing temperature when thermal expansion coefficient α > 0, following ρ(T) = ρ₀(1 - α(T - T_ref)).

#### Scenario: Hot phase has lower density
- **WHEN** two phases have identical reference density but different temperatures (T_hot=1600K vs T_cold=600K) with α=3e-5 K⁻¹
- **THEN** `min(rho_n)` (corresponding to the hot region) SHALL be less than `max(rho_n)` (cold region)

### Requirement: Hydrostatic pressure test
The test suite SHALL verify that pressure increases with depth following P ≈ ρgz for a static column with no imposed strain rate.

#### Scenario: Pressure increases with depth
- **WHEN** a simulation runs with `bkg_strain_rate=0`, `pure_shear_ALE=0`, `gz=-10`, uniform density, and no-slip BCs
- **THEN** `max(P)` SHALL be greater than zero and within 50% of the expected ρgH value
