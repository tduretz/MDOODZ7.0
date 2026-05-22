## ADDED Requirements

### Requirement: Temperature-dependent relaxation of the anisotropy factor δ

When a phase has `ani_fstrain == 3`, the viscous-anisotropy factor δ SHALL relax toward the isotropic limit (δ → 1) on the kinetic timescale `τ_relax = L_relax / V`, where the grain-boundary-migration velocity is `V = M₀ · exp(−Q/(R·T)) · μ · b² · Δρ` (Boneh et al. 2021, *G3*, Eqs. 1–4: `Q = 133 ± 18 kJ/mol`, `μ = 50 GPa`, `b = 0.6 nm`, and `M₀ = 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹` ≡ `m³·N⁻¹·s⁻¹` — Boneh 2021 prints "m³/s·J", a typo corrected per the Boneh et al. 2017 companion paper's Appendix A; the corrected value reproduces Boneh's "weeks-to-years at ~1300 °C" anchor). `V` SHALL use only the strain-energy driving force `μb²Δρ`; the surface-energy term `F_b = 3γ/d` is dropped (`ΣF ≈ F_s`, valid because `F_s ≫ F_b` for mantle grain sizes ≥ 1 mm). The relaxation SHALL be temperature-gated by the Arrhenius term — hot markers relax fast, cold markers relax negligibly — and δ SHALL approach but never cross below 1.

#### Scenario: Hot marker relaxes toward isotropic

- **WHEN** a marker with `ani_fstrain == 3` and δ > 1 is held at high temperature with no active deformation for a time much greater than τ_relax
- **THEN** δ SHALL decay monotonically toward 1
- **AND** δ SHALL NOT go below 1

#### Scenario: Cold marker is effectively frozen

- **WHEN** a marker with `ani_fstrain == 3` and δ > 1 is held at a temperature low enough that τ_relax is much greater than the model run time
- **THEN** δ SHALL remain unchanged to within floating-point precision over the run

#### Scenario: Relaxation rate follows the Arrhenius law

- **WHEN** the same marker state is evaluated at two temperatures T₁ < T₂
- **THEN** τ_relax(T₂) SHALL be shorter than τ_relax(T₁) by the factor `exp(Q/R · (1/T₁ − 1/T₂))`

### Requirement: δ is a per-marker integrated state variable

The relaxing anisotropy factor SHALL be carried as a per-marker state field (δ, or the equivalent CPO-strength index `M`), allocated for the lifetime of the run, advected with the marker, copied on reseeding and inflow, and persisted across restart. It SHALL NOT be a memoryless function of the instantaneous `FS_AR`.

#### Scenario: History-dependence

- **WHEN** two markers reach the same current `FS_AR` but one spent time hot (relaxing) and the other stayed cold
- **THEN** their stored δ values SHALL differ

#### Scenario: Restart round-trip

- **WHEN** a run is checkpointed and restarted
- **THEN** the per-marker δ state SHALL be bit-identical across the restart

### Requirement: Analytic exponential integration scheme

The per-step update of δ SHALL use the analytic exponential form `δ_new = δ_eq + (δ_old − δ_eq) · exp(−Δt / τ_relax)` with `δ_eq = 1`. This scheme SHALL be unconditionally stable for any timestep Δt, SHALL NOT overshoot past `δ_eq`, and SHALL NOT depend on the (anisotropy-unaware) Courant timestep limiter.

#### Scenario: No overshoot at large timestep

- **WHEN** the update is applied with Δt much greater than τ_relax
- **THEN** δ_new SHALL satisfy `1 ≤ δ_new ≤ δ_old`

#### Scenario: Exponential half-life

- **WHEN** δ relaxes under constant T, constant Δρ, and constant L_relax
- **THEN** δ(t) SHALL follow `1 + (δ₀ − 1)·exp(−t/τ_relax)`
- **AND** SHALL reach `1 + (δ₀ − 1)/2` at `t = τ_relax · ln 2`

### Requirement: L_relax is a per-phase length-scale parameter, GSE not required

The relaxation length scale `L_relax` SHALL be a per-phase quantity representing the characteristic DiSRX boundary-migration distance. It SHALL be sourced from the existing per-phase grain-size input (`gs_ref`) so that grain-size evolution (the `gs` flow law) is NOT required to be active. When grain-size evolution is active, `L_relax` MAY instead track the evolving grain-size field.

#### Scenario: Runs with grain-size evolution disabled

- **WHEN** a phase has `ani_fstrain == 3` and the `gs` flow law disabled
- **THEN** the relaxation SHALL still run, using `L_relax = gs_ref` for that phase

#### Scenario: L_relax monotonicity

- **WHEN** two otherwise-identical markers have `L_relax` values L₁ < L₂
- **THEN** the marker with L₂ SHALL relax more slowly, since τ_relax is proportional to L_relax

### Requirement: Δρ uses a dislocation-density proxy from existing fields

The dislocation-density driving force `Δρ` SHALL be computed from per-marker fields that MDOODZ already tracks (accumulated dislocation-creep strain and/or deviatoric stress). No new dislocation-density state variable SHALL be introduced. The proxy SHALL be finite and non-negative for all valid marker states.

#### Scenario: No new dislocation-density state variable

- **WHEN** the relaxation is implemented
- **THEN** the `markers` struct SHALL NOT gain a dislocation-density field
- **AND** Δρ SHALL be derived on the fly from existing marker fields

### Requirement: Cold-limit behaviour recovers the `ani_fstrain == 2` form

When the relaxation term is inactive (τ_relax much longer than the run time), the δ produced under `ani_fstrain == 3` SHALL track the δ produced under `ani_fstrain == 2` for the same inputs, to floating-point precision.

#### Scenario: Cold marker under ani_fstrain == 3 tracks the ani_fstrain == 2 form

- **WHEN** a marker is run under `ani_fstrain == 3` at a temperature cold enough that τ_relax is much greater than the run time
- **THEN** its δ SHALL track `aniso_delta_fn(FS_AR)` identically to an `ani_fstrain == 2` marker with the same inputs

### Requirement: Backward compatibility for `ani_fstrain` values 0, 1, and 2

The δ-relaxation SHALL be reachable only via `ani_fstrain == 3`. The dispatch arms for `ani_fstrain == 0`, `1`, and `2` SHALL produce byte-identical output to their pre-change behaviour, and all pre-existing anisotropy CI tests SHALL pass unchanged.

#### Scenario: Existing dispatch arms unchanged

- **WHEN** the existing anisotropy CI suite is run after this change
- **THEN** every test for `ani_fstrain` in {0, 1, 2} SHALL pass to its prior tolerance with no fixture changes

### Requirement: The temperature-threshold F-freeze is removed

The previous `ani_fstrain == 3` behaviour — the upstream `T < ani_T_threshold` deformation-gradient freeze and the `ani_T_threshold` per-phase parameter — SHALL be removed. The case-7 `ani_T_threshold` default SHALL be replaced with an `L_relax`-equivalent default.

**Reason**: The 1100 °C threshold was a misattributed value (its "three convergent paths" justification does not hold under primary-source review), and the freeze was a misinterpretation of the intended kinetic-relaxation directive.

**Migration**: 14 `SETS/` files reference the old behaviour (11 `.txt` scenarios that set `ani_fstrain = 3` and/or `ani_T_threshold`, plus 2 `.c` callbacks that mention `ani_T_threshold` in comments — the `Beghein14OceanicCooling*` and `Q195*_Q3_phase1*` families). No run *breaks* — orphaned `ani_T_threshold` lines are silently ignored by the parser — but `ani_fstrain == 3` changes meaning, so these obsolete freeze-experiment scenarios SHALL be removed or migrated. The `AniFstrainT_Threshold_Freeze_Sentinel` and `AniFstrainT_Threshold_Freeze_Active` GTests are removed — they exercised the freeze and were already crashing on `aniso_db = 0`.

#### Scenario: ani_T_threshold is no longer parsed

- **WHEN** the input parser reads a material block
- **THEN** `ani_T_threshold` SHALL NOT be a recognised parameter
- **AND** the `T < ani_T_threshold` F-freeze code block SHALL NOT exist in the marker deformation-gradient update

### Requirement: Analytical-unit-test CI coverage

The CI suite SHALL include GoogleTest cases that verify the relaxation against closed-form limits on a 1-cell or homogeneous box: cold freeze (τ → ∞ gives δ constant); hot relaxation matching `1 + (δ₀ − 1)·exp(−t/τ)` on the Boneh timescale; the exponential half-life with a Δt-convergence sub-test; the cold-limit equality with `ani_fstrain == 2`; L_relax monotonicity; the `ani_fac_max` clamp honoured; and the δ = 1 isotropic fixed point.

#### Scenario: Relaxation limits are CI-gated

- **WHEN** the anisotropy CI suite runs
- **THEN** it SHALL contain passing test cases for the cold-freeze, hot-relaxation, half-life, cold-limit-equality, monotonicity, clamp, and fixed-point limits listed above

### Requirement: Modelling caveats are documented inline in the code

Every modelling approximation of the relaxation law SHALL be documented as an inline comment at the relevant source location, not only in an external notes file. At minimum: the constant-`L_relax` assumption and the mid-temperature / episodic-loading regime where its bounded (~17 %) error lands; the `Δρ` proxy and its caveats; the modelling leap that δ relaxes on the DiSRX grain-growth timescale (Boneh links DiSRX to CPO weakening qualitatively — the relaxation ODE is the model's own construction); and a dimensional check on the Boneh `M₀` units.

#### Scenario: Caveat comments present at code sites

- **WHEN** the relaxation routine, the per-phase parameter parsing, and the marker δ update are reviewed
- **THEN** each of the four caveats above SHALL appear as an inline comment at the code site it concerns

### Requirement: Research note documents the kinetics and the validation plan

The change SHALL deliver a research-grade markdown note under `misc/aniso_fstrain/notes/` documenting the Boneh et al. 2021 DiSRX kinetics used, the derivation of `τ_relax = L_relax / V` and its modelling assumptions, the constant-`L_relax` regime analysis (the frozen / instantaneous / caught-partway bands), the `Δρ` proxy choice with its caveats, and the tiered validation plan (analytical unit tests as the CI gate; lab-data reproduction against Boneh et al. 2017 and Speciale et al. 2020 as future research-grade checks).

#### Scenario: Research note exists with required content

- **WHEN** the change is complete
- **THEN** the note SHALL exist under `misc/aniso_fstrain/notes/`
- **AND** it SHALL cover the kinetics, the `τ_relax` derivation and assumptions, the constant-`L_relax` regime analysis, the `Δρ` proxy, and the tiered validation plan
