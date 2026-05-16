## ADDED Requirements

### Requirement: `ani_fstrain = 3` enables temperature-dependent δ-relaxation

The `ani_fstrain` enumeration SHALL define value `3` as: the `ani_fstrain == 2` δ-dispatch (mineral-calibrated `aniso_delta_fn(FS_AR)` selected by `aniso_db`) **plus** a temperature-dependent kinetic relaxation of the anisotropy factor δ toward the isotropic limit. The relaxation physics, state, integration scheme, inputs, and CI coverage are specified by the `aniso-delta-relaxation` capability. The previous meaning of `ani_fstrain == 3` (an upstream `T < ani_T_threshold` deformation-gradient freeze) and the associated `ani_T_threshold` per-phase parameter SHALL no longer exist.

#### Scenario: ani_fstrain enumeration documents value 3

- **WHEN** the `ani_fstrain` enumeration is documented (code comments and `TESTS/AnalyticalSolutions.md`)
- **THEN** value `3` SHALL be described as "`ani_fstrain == 2` dispatch + temperature-dependent δ-relaxation"
- **AND** `ani_T_threshold` SHALL NOT appear as a per-phase parameter

#### Scenario: ani_fstrain values 0, 1, 2 are unaffected

- **WHEN** this change is applied
- **THEN** the dispatch arms and CI tests for `ani_fstrain` in {0, 1, 2} SHALL be byte-identical to their prior behaviour

#### Scenario: T-threshold-freeze tests are removed

- **WHEN** the anisotropy CI suite is built after this change
- **THEN** the `AniFstrainT_Threshold_Freeze_Sentinel` and `AniFstrainT_Threshold_Freeze_Active` GoogleTest cases SHALL no longer exist
- **AND** the suite SHALL instead contain the `aniso-delta-relaxation` analytical-unit-test cases
