# aniso-init-from-finite-strain Specification

## Purpose
Establish consistency between inherited fabric magnitude (`aniso_factor`),
orientation (`aniso_angle`), and the per-marker deformation-gradient tensor F
by initialising F to encode the finite-strain history that would produce the
configured anisotropy strength under the chosen analytic δ-form. Eliminates
the "+offset above Hansen ceiling" failure mode in the operator-split update
when an inherited cratonic fabric (`aniso_factor > 1`) is prescribed.
## Requirements
### Requirement: Initial fabric magnitude encoded as past finite strain via the per-marker F-tensor

For any phase with `ani_fstrain >= 1`, an analytic `aniso_delta_fn` populated by the configured `aniso_db`, and `aniso_factor[phase] > 1`, the per-marker deformation gradient `F` (fields `Fxx, Fxz, Fzx, Fzz`) SHALL be initialised at simulation start such that the FS_AR computed from `F` by the closed-form 2x2 SVD (see `aniso-delta-relaxation` requirement "Cold-limit behaviour", `MDLIB/RheologyParticles.c:712-732`) satisfies `aniso_delta_fn(FS_AR_init) = aniso_factor[phase]` to within floating-point precision. The principal extension direction of `F_init` SHALL be aligned with the easy-slip direction defined by `aniso_angle[phase]`.

#### Scenario: Inherited δ = 4 olivine marker, aniso_angle = 0

- **WHEN** a phase has `ani_fstrain = 3`, `aniso_db = 1` (Hansen 38-point), `aniso_factor = 4.0`, `aniso_angle = 0°`
- **THEN** the marker's initial `F` SHALL produce `FS_AR_init ≈ 2.69` via the closed-form SVD
- **AND** `aniso_delta_fn(FS_AR_init)` for `aniso_db = 1` SHALL equal `4.0 ± 1e-9`
- **AND** the principal extension axis of `F_init` (eigenvector of `F_init^T·F_init` with the larger eigenvalue) SHALL be parallel to +x

#### Scenario: Inherited δ = 4 olivine marker, aniso_angle = 80°

- **WHEN** the same phase has `aniso_angle = 80°`
- **THEN** the marker's initial `F` SHALL produce `FS_AR_init ≈ 2.69` (unchanged from the 0° case)
- **AND** the principal extension axis of `F_init` SHALL be rotated 80° counter-clockwise from +x

#### Scenario: Isotropic init is a no-op

- **WHEN** a phase has `ani_fstrain >= 1` but `aniso_factor = 1.0`
- **THEN** the marker's initial `F` SHALL be the identity tensor
- **AND** `FS_AR_init` SHALL equal `1.0`

### Requirement: `aniso_delta` and `aniso_delta_fs_prev` are init-consistent with `F`

When a phase activates the init-from-finite-strain pathway, the per-marker state SHALL satisfy `aniso_delta[k] = aniso_delta_fs_prev[k] = aniso_delta_fn(FS_AR(F_init[k])) = aniso_factor[phase]` at simulation start. The three fields SHALL NOT be initialised inconsistently with each other.

#### Scenario: Three-field consistency at step 0

- **WHEN** an `ani_fstrain = 3` phase with `aniso_factor = 4` is initialised
- **THEN** `aniso_delta_fs_prev[k]` SHALL equal `aniso_factor` (= 4.0) for every marker of that phase
- **AND** `aniso_delta[k]` SHALL also equal `aniso_factor` (= 4.0)
- **AND** the FS_AR computed from `F_init[k]` SHALL produce `aniso_delta_fn(FS_AR_init) = 4.0`

### Requirement: Cold-limit operator-split telescope is exact

With the consistent initialisation above, the cold-limit telescope of the `ani_fstrain = 3` operator split (`δ_prod = aniso_delta + (δ_FS_new − aniso_delta_fs_prev)`, `δ_new = 1 + (δ_prod − 1)·exp(−Δt / τ_eff)` for `τ_eff → ∞`) SHALL satisfy `aniso_delta_N → aniso_delta_fn(FS_AR_N)` independently of the initial `aniso_factor` value, with zero +offset.

#### Scenario: Cold-limit cell never exceeds Hansen ceiling

- **WHEN** a marker on a phase with `aniso_factor = 4`, `aniso_db = 1`, `ani_fstrain = 3` is held in a cold regime (Boneh τ ≫ run time) and accumulates strain γ from past + simulation totalling γ_total
- **THEN** at all times `aniso_delta[k] ≤ aniso_delta_fn(FS_AR(γ_total)) ≤ 14.13`
- **AND** when γ_total is large enough that the Hansen fit saturates, `aniso_delta[k] → 14.13` exactly (no `aniso_factor − 1` overshoot)

#### Scenario: Hot-quiescent cell relaxes from inherited toward isotropic

- **WHEN** a marker on the same phase is held at a hot, quiescent temperature with `aniso_delta_fs_prev = aniso_factor` and no further strain accumulation (constant `F` → constant `δ_FS = aniso_factor`)
- **THEN** `aniso_delta[k]` SHALL decay analytically as `1 + (aniso_factor − 1)·exp(−Σ·Δt / τ)` toward 1
- **AND** at the half-life `Σ·Δt = τ·ln 2`, `aniso_delta` SHALL equal `1 + (aniso_factor − 1)/2`

### Requirement: Out-of-range `aniso_factor` is clamped, not crashed

If the user configures `aniso_factor[phase]` above the analytic saturation of the chosen `aniso_db` (e.g. `aniso_factor = 20` for `aniso_db = 1` where `δ_max = 14.13`), the initialiser SHALL clamp the target δ to `δ_max − ε` (ε = 1e-9), emit `LOG_WARN` naming the phase and the clamped value, and continue with the saturated `γ_eff`.

#### Scenario: Over-saturated request

- **WHEN** a phase has `aniso_factor = 20` and `aniso_db = 1` (Hansen, δ_max = 14.13)
- **THEN** the init SHALL log a warning
- **AND** `FS_AR_init` SHALL correspond to `δ_target = 14.13 − 1e-9`
- **AND** simulation SHALL continue without aborting

### Requirement: Per-`aniso_db` analytic inverse coverage

For every `aniso_db` case currently dispatched by `aniso_delta_fn` in `FlowLaws.c::ReadDataAnisotropy` (cases `{1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15}` per the audited enumeration), a corresponding analytic-inverse helper SHALL be available. For saturating-exponential cases `{1, 2, 3, 4, 5, 6, 7, 8, 9, 13}` the inverse SHALL return `γ_eff = −γ_e · ln(1 − (δ − 1)/c2)` where `c2, γ_e` are the mineral-specific calibration constants from the forward fit. For case 15 (decay form `δ(γ) = c2·exp(−γ/γ_e) + 1`) the inverse SHALL return `γ_eff = −γ_e · ln((δ − 1)/c2)` with reversed monotonicity semantics ("high `aniso_factor` ↔ low past strain"). A bisection fallback SHALL handle any future `aniso_db` whose forward has no closed-form inverse.

#### Scenario: All currently-implemented `aniso_db` have analytic inverses

- **WHEN** the init is invoked on a phase with `aniso_db ∈ {1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15}` and `aniso_factor > 1`
- **THEN** the helper SHALL return γ_eff in O(1) time
- **AND** the forward `aniso_delta_fn(FS_AR(γ_eff))` SHALL reproduce the target δ to within 1e-9

#### Scenario: Decay-form case 15 emits a semantics warning

- **WHEN** the init is invoked on a phase with `aniso_db == 15` and `aniso_factor > 1`
- **THEN** the helper SHALL emit a one-shot `LOG_WARN` clarifying that for case 15 a HIGH `aniso_factor` corresponds to LOW past strain (opposite to all other cases)
- **AND** the F-init shall proceed with the decay-form inverse

### Requirement: Reseeded and inflow markers inherit consistent state via existing G2P

The init-from-finite-strain pathway SHALL run once at simulation start, before the timestep loop. New markers created during the run by reseeding, free-surface insertion, or boundary inflow SHALL inherit their `F`, `aniso_delta`, and `aniso_delta_fs_prev` values via the existing grid-to-particle interpolation at insertion (cell-mean values), with no special-case re-initialisation.

#### Scenario: Reseeded marker continuity

- **WHEN** a reseeded marker is created mid-run in a cell whose existing markers have evolved `F` and `aniso_delta` to specific values
- **THEN** the new marker's `F`, `aniso_delta`, and `aniso_delta_fs_prev` SHALL each match the cell-mean of its existing-marker neighbours to within G2P interpolation precision
- **AND** no recall of the original `aniso_factor[phase]` SHALL be made for that new marker

