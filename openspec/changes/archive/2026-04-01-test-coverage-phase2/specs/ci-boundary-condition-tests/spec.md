## MODIFIED Requirements

### Requirement: Velocity field verification for existing BC tests
The existing boundary condition tests SHALL additionally read Vx/Vz from HDF5 and verify that velocity fields are consistent with the imposed boundary conditions.

#### Scenario: Pure shear Vx matches imposed strain rate
- **WHEN** a FreeSlipPureShear simulation runs with ε̇=1 on domain [-0.5, 0.5]
- **THEN** `max(Vx)` from VxNodes SHALL be approximately 0.5 (within 20%)

#### Scenario: No-slip velocity is near zero at boundaries
- **WHEN** a NoSlip simulation runs
- **THEN** `max(P)` SHALL remain bounded (existing) AND `max(Vx)` SHALL be less than a reasonable bound given the gravity-driven flow
