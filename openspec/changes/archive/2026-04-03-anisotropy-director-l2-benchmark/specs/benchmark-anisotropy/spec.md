## MODIFIED Requirements

### Requirement: Director evolution smoke test upgraded to L2
_Modifies: "Director evolution under simple shear matches analytical ODE" from main spec_

The existing `DirectorEvolution` test SHALL be upgraded from a qualitative check (`meanAngle < 44°`) to a quantitative L2 error check against the closed-form solution θ(t) = arctan(tan(θ₀) − γ̇·t). The test parameters (θ₀=45°, δ=2, γ̇=1.0, 10 steps, dt=0.05) remain unchanged.

#### Scenario: L2 replaces qualitative threshold
- **WHEN** the `DirectorEvolution` test runs
- **THEN** it SHALL assert `L2_director < 0.01` (radians) instead of `meanAngle < 44°`
- **AND** the existing `EXPECT_LT(meanAngle, 44.0)` assertion SHALL be removed

### Requirement: Stress invariant upgraded to analytical comparison
_Modifies: "Stress invariant at known director angle" from main spec_

The existing `StressAnisotropy` test SHALL be upgraded from a qualitative check (`meanTauII > 0.1`) to a quantitative comparison against the analytical stress tensor computed from the transverse isotropy rotation formula. The current parameters (θ=30°, δ=6, pure shear) remain unchanged.

#### Scenario: Relative error replaces magnitude check
- **WHEN** the `StressAnisotropy` test runs
- **THEN** it SHALL assert that τ_II relative error < 2% instead of `EXPECT_GT(meanTauII, 0.1)`

### Requirement: Parameter spec alignment
_Modifies: "Parameter files exist and are loadable" from main spec_

The main spec scenario for director evolution references δ=10 and 50 steps; the actual test uses δ=2 and 10 steps. The spec SHALL be updated in the main spec to reflect the actual parameters used: θ₀=45°, δ=2, γ̇=1.0, dt=0.05, Nt=10.

#### Scenario: Main spec matches implementation
- **WHEN** the main spec `benchmark-anisotropy/spec.md` is synced
- **THEN** the director evolution scenario SHALL reference δ=2, Nt=10, dt=0.05, total time t=0.5
