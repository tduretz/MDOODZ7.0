## MODIFIED Requirements

### Requirement: Maxwell stress accumulation test
The test suite SHALL verify that deviatoric stress accumulates over multiple time steps when `elastic=1` (Maxwell visco-elastic model). The elastic strain rate field `eII_el` SHALL be non-zero, and the maximum deviatoric stress `sxxd` SHALL increase between step 1 and step 5.

#### Scenario: Stress grows under sustained loading
- **WHEN** a pure-shear simulation runs for 5 steps with `elastic=1`, constant viscosity, and `G=1e10` Pa
- **THEN** `max(sxxd)` at step 5 SHALL be greater than `max(sxxd)` at step 1

#### Scenario: Elastic strain rate is non-zero
- **WHEN** `elastic=1` is enabled in a deforming simulation
- **THEN** `max(eII_el)` SHALL be greater than zero

#### Scenario: Maxwell stress accumulation L2 error
- **WHEN** a pure-shear simulation runs for Nt steps with `elastic=1`, constant viscosity η, shear modulus G, and strain rate ε̇
- **THEN** the relative L2 error of the sxxd field at the final step against the analytical Maxwell solution σ(t) = 2η·ε̇·(1 − exp(−G·t/η)) SHALL be less than 5e-2
