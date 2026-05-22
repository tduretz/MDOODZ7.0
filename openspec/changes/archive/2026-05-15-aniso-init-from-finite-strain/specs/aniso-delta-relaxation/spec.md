## MODIFIED Requirements

### Requirement: Cold-limit behaviour recovers the `ani_fstrain == 2` form

When the relaxation term is inactive (τ_eff much longer than the run time, e.g. because the marker is cold, or because the strain-rate gate `f_gate → 0` suppresses DiSRX, or both), the δ produced under `ani_fstrain == 3` SHALL track the δ produced under `ani_fstrain == 2` for the same inputs, **with no permanent additive offset, independent of the initial `aniso_factor[phase]` value**. This invariant SHALL hold because the per-marker `F`-tensor and `aniso_delta_fs_prev` are initialised consistent with `aniso_factor` (see capability `aniso-init-from-finite-strain`), so the cold-limit telescope `aniso_delta_N → aniso_delta_fn(FS_AR_N) + (aniso_delta_0 − aniso_delta_fs_prev_0)` evaluates to `aniso_delta_N → aniso_delta_fn(FS_AR_N)` because the bracketed initial-condition difference is zero by construction.

#### Scenario: Cold marker under ani_fstrain == 3 tracks the ani_fstrain == 2 form

- **WHEN** a marker is run under `ani_fstrain == 3` at a temperature cold enough that τ_eff is much greater than the run time
- **THEN** its δ SHALL track `aniso_delta_fn(FS_AR)` identically to an `ani_fstrain == 2` marker with the same inputs and the same initial `F` tensor

#### Scenario: Inherited δ does not exceed the Hansen ceiling in cold + actively-deformed cells

- **WHEN** a marker on a phase with `aniso_factor = 4`, `aniso_db = 1` (Hansen 38-point) is in a strain-rate-gated cell (gate closes, τ_eff → ∞) and accumulates strain γ_total
- **THEN** at all times `aniso_delta[k] ≤ aniso_delta_fn(FS_AR(γ_total)) ≤ 14.13`
- **AND** the `+(aniso_factor − 1) = +3` overshoot above the Hansen ceiling that prior to this change appeared in v2 rifting runs (observed `δ_max = 16.04`) SHALL NOT occur

#### Scenario: aniso_factor = 1 phase is unchanged by the init-consistency requirement

- **WHEN** a phase has `aniso_factor = 1.0`, `ani_fstrain == 3`
- **THEN** the per-marker `F` SHALL remain the identity tensor at simulation start
- **AND** the cold-limit telescope SHALL evaluate to `aniso_delta_N → aniso_delta_fn(FS_AR_N)` with `aniso_delta_fs_prev_0 = 1`
- **AND** behaviour SHALL be byte-identical to the pre-change `aniso_factor = 1` initial state

### Requirement: δ is a per-marker integrated state variable

The relaxing anisotropy factor SHALL be carried as a per-marker state field (δ, or the equivalent CPO-strength index `M`), allocated for the lifetime of the run, advected with the marker, copied on reseeding and inflow, and persisted across restart. **At simulation start, this state SHALL be initialised to `aniso_delta_fn(FS_AR(F_init))` where `F_init` is set by the `aniso-init-from-finite-strain` pathway, so that `aniso_delta_0 = aniso_factor[phase]` and `aniso_delta_fs_prev_0 = aniso_factor[phase]` are mutually consistent.** It SHALL NOT be a memoryless function of the instantaneous `FS_AR`.

#### Scenario: History-dependence

- **WHEN** two markers reach the same current `FS_AR` but one spent time hot (relaxing) and the other stayed cold
- **THEN** their stored δ values SHALL differ

#### Scenario: Restart round-trip

- **WHEN** a run is checkpointed and restarted
- **THEN** the per-marker δ state SHALL be bit-identical across the restart

#### Scenario: Initialisation consistency invariant

- **WHEN** the model is at step 0 for any phase with `ani_fstrain == 3` and `aniso_factor[phase] > 1`
- **THEN** `aniso_delta_fs_prev[k]` SHALL equal `aniso_delta_fn(FS_AR(F_init[k]))`
- **AND** `aniso_delta[k]` SHALL also equal `aniso_delta_fn(FS_AR(F_init[k]))`
- **AND** both SHALL equal `aniso_factor[phase]` to within floating-point precision
