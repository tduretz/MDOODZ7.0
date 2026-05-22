## Why

The `ani_fstrain == 3` viscous-anisotropy delta-relaxation (archived change `ani-fstrain3-delta-relaxation`, capability `aniso-delta-relaxation`) carries a known initialisation inconsistency that breaks the operator split's claimed cold-limit recovery to the Hansen 38-point fit. With the standard "inherited cratonic fabric" initial condition (`aniso_factor = 4`, `aniso_angle = 80°`), the current code initialises three quantities mutually inconsistently:

- `particles->aniso_delta[k] = aniso_factor[phase] = 4` (set by `InitialiseAnisoDeltaParticles`),
- `particles->aniso_delta_fs_prev[k] = 1.0` (PartInit default),
- `particles->Fxx[k] = Fzz[k] = 1`, `Fxz[k] = Fzx[k] = 0` (PartInit identity), giving `FS_AR_init = 1` and `aniso_delta_fn(FS_AR_init) = 1` for `aniso_db = 1`.

The cold-limit telescope of the operator-split formula `δ_prod = aniso_delta + (δ_FS − aniso_delta_fs_prev)`, `δ_new = 1 + (δ_prod − 1)·exp(−Δt/τ_eff)` then yields `δ_N → δ_FS(FS_AR_N) + (aniso_delta_0 − aniso_delta_fs_prev_0) = δ_FS + (aniso_factor − 1)`, i.e. a permanent **`+offset = aniso_factor − 1`** above the Hansen ceiling — empirically observed as `δ_max = 16.04` in the v2 rifting run at 5 Ma when the Hansen 38-point ceiling is `24.5·0.536 + 1 = 14.13`. The current archived spec `aniso-delta-relaxation` claims that the cold limit "recovers the `ani_fstrain == 2` form" (Requirement "Cold-limit behaviour recovers the `ani_fstrain == 2` form"); with the current initialisation this claim is false.

The +offset is a real-physics overshoot: lab-saturated olivine cannot exceed δ ≈ 14.13 regardless of past history (Hansen 12/14/16; Tasaka 16; Kumamoto 19; Bystricky 00). Whatever inherited fabric a cratonic substrate carries is itself a representation of past strain, not an additive magnitude on top of new strain. The fix is to represent inherited fabric the way it physically arose — as a finite-strain history — by initialising the marker's `F`-tensor consistent with the desired `aniso_factor` and `aniso_angle`.

## What Changes

- **New capability** `aniso-init-from-finite-strain` — formalises the init-time mapping `(aniso_factor, aniso_angle) → F_init` for any phase with `ani_fstrain != 0` whose `aniso_db` defines an analytic `aniso_delta_fn`, and the consistency invariants between `aniso_delta`, `aniso_delta_fs_prev`, and `F` that follow from it.
- **Modified capability** `aniso-delta-relaxation` — the "Cold-limit behaviour recovers the `ani_fstrain == 2` form" requirement is **upgraded from "approximate / additive" to "exact"**: with the new init-from-F, the cold-limit telescope yields `aniso_delta → aniso_delta_fn(FS_AR)` exactly, with no `+offset`, because `aniso_delta_0 − aniso_delta_fs_prev_0 = aniso_factor − aniso_delta_fn(FS_AR_init) = 0`.
- **Per-marker `F`-tensor initialisation** for `ani_fstrain >= 1` phases with `aniso_factor > 1`: invert the per-phase `aniso_delta_fn` to find the equivalent γ_eq, build the corresponding pure-shear or simple-shear `F_canonical(γ_eq)` at the principal-direction frame, rotate by `aniso_angle` to align with the configured easy-slip direction.
- **Per-marker `aniso_delta_fs_prev` initialisation** SHALL be set to `aniso_delta_fn(FS_AR_init) = aniso_factor` (no longer 1.0) on `ani_fstrain == 3` markers, completing the cold-limit-telescoping fix.
- **`δ ≤ aniso_fac_max` and `δ ≤ aniso_delta_fn(FS_AR_CAP)`** hold by construction without further clamping in the cold limit (still clamped defensively at `[1, ani_fac_max]`).
- **BREAKING** for any existing setup that depends on the legacy `F = I` initial state on a phase with `aniso_factor > 1` and `ani_fstrain == 3` — those runs were producing the +offset behaviour (now considered a bug); behaviour after this change is the strict-Hansen one and will differ. Documented in tasks. (Scope narrowed from "`ani_fstrain >= 1`" per design Open-Question OQ1: existing `ani_fstrain ∈ {1, 2}` calibration SETs are explicitly preserved.)
- Backwards-compatible defaults: when a phase has `aniso_factor == 1` or `ani_fstrain != 3`, `F_init` SHALL remain the identity (no change to existing isotropic-init / `ani_fstrain == 2` workflows). When a phase has `ani_fstrain == 3` AND `aniso_factor > 1` AND a populated `aniso_delta_fn_inv`, the init-from-finite-strain path activates automatically.
- HDF5 output, restart, and reseeding pathways SHALL be audited to confirm the new non-identity `F_init` survives intact (no implicit assumption of `F = I` at step 0 anywhere downstream).

## Capabilities

### New Capabilities

- `aniso-init-from-finite-strain`: per-marker initialisation of the deformation-gradient tensor `F` and `aniso_delta_fs_prev` consistent with the configured `aniso_factor` (magnitude) and `aniso_angle` (direction), so that inherited fabric is encoded as past finite strain rather than as a magnitude offset.

### Modified Capabilities

- `aniso-delta-relaxation`: the cold-limit-recovery requirement is upgraded from approximate (with documented +offset caveat) to exact (`aniso_delta → aniso_delta_fn(FS_AR)` independently of `aniso_factor` initial value); the initialisation invariants subsection is extended to require `aniso_delta_fs_prev_0 = aniso_delta_fn(FS_AR_init) = aniso_factor`.

## Impact

- **Code** (no commits in this change; design + tasks only):
  - `MDLIB/RheologyParticles.c::InitialiseAnisoDeltaParticles` — extend to set per-marker `F` and `aniso_delta_fs_prev` consistent with `aniso_factor` and `aniso_angle`.
  - `MDLIB/ParticleRoutines.c::PartInit` — `aniso_delta_fs_prev = 1` default unchanged (it is overwritten by `InitialiseAnisoDeltaParticles` for relevant phases); document the order-of-operations invariant.
  - `MDLIB/AnisotropyRoutines.c` — provide a small inverse-lookup helper `InverseAnisoDeltaFn(double delta, int phase, mat_prop *)` that returns γ_eq from a target δ for the active `aniso_db` form (Hansen / calcite / quartz / etc.). For `aniso_db = 1` (Hansen olivine, 38-point fit) the analytic inverse exists: `γ_eq = −3.96·ln(1 − (δ−1)/(24.5·0.536))`. For other `aniso_db` cases the inverse is a one-shot bisection (range `[0, γ_cap]`).
  - `MDLIB/RheologyParticles.c::FiniteStrainAspectRatio` — comment update only (the cold-limit telescope now truly yields `aniso_delta_fn(FS_AR)`).
- **Tests**: CI test demonstrating the cold-limit-telescope invariant at `Nt = 1` for `aniso_factor ∈ {1, 4, 8}` × `aniso_angle ∈ {0°, 45°, 80°, 90°}`.
- **HDF5 / I/O**: `F` already written in the existing centroid output (`Fxx, Fxz, Fzx, Fzz`); the change is the *value* at step 0, not the schema.
- **Restart**: per-marker `F` and `aniso_delta` already persist across breakpoints; this change requires no new restart-state.
- **Reseeding / inflow / Free-Surface**: new markers created mid-run inherit local cell-mean `F` (already true). No additional change needed.
- **Documentation**: openspec change directory + this design doc; the user-visible note is "if `aniso_factor > 1` and `ani_fstrain >= 1`, your initial F-tensor is no longer identity."
- **Scientific implication**: the v2 rifting run's `δ_max = 16.04` overshoot would have been ≤ 14.13 with this change. The strain-rate-gated CPO build-up in shear bands is unchanged in pattern; only the absolute ceiling moves down by `~aniso_factor − 1`. Inherited fabric in cold static regions is preserved at exactly `δ = aniso_factor` (no longer drifts).
