## Why

`ani_fstrain == 3` is currently a **temperature-threshold hard freeze**: when a marker's `T < ani_T_threshold`, its deformation-gradient tensor `F` is frozen ([RheologyParticles.c:603-621](MDLIB/RheologyParticles.c#L603-L621)), and case 7 (damped olivine) defaults `ani_T_threshold = 1100 °C` ([FlowLaws.c:1638-1650](MDLIB/FlowLaws.c#L1638-L1650)). A multi-round literature + code investigation established three problems:

1. **The 1100 °C value is misattributed.** Its stated "three convergent paths" justification (Karato dry ~1300 °C / HK03 wet ~1100 °C / Beghein+14 ~1100–1200 °C) fails on all three: Karato and HK03 set *no* temperature threshold — the dislocation/diffusion transition they describe is stress- and grain-size-controlled — and Beghein+14 reports a *depth* (the Gutenberg discontinuity), explicitly not a temperature. The one genuine direct citation (Debayle & Ricard 2013) uses 1100 °C only as a modelling *convention*.

2. **The freeze was a misinterpretation of the original directive.** The intent was a *kinetic relaxation law* for the anisotropy factor — the deferred "Phase 2 (kinetic-law dispatch)" already noted in [AnisotropyRoutines.c:70](MDLIB/AnisotropyRoutines.c#L70). The autonomous research loop built a step-function freeze ("Phase 1") with a retro-fitted citation chain instead.

3. **A fixed-T freeze is the wrong *kind* of model.** It is binary (on/off), has no rate, and discards the physics: when deformation slows or stops, olivine CPO does not "freeze" — it *relaxes* toward isotropic over a finite, temperature-dependent timescale via discontinuous static recrystallization (DiSRX). Boneh et al. 2021 (G3) give the grain-boundary-migration kinetics for exactly this process.

This change replaces the freeze with the relaxation law that was actually intended.

## What Changes

- **BREAKING — `ani_fstrain == 3` is redefined.** From "`ani_fstrain == 2` δ-dispatch + upstream T-threshold `F`-freeze" to "`ani_fstrain == 2` δ-dispatch + a temperature-dependent kinetic *relaxation* of the anisotropy factor δ." **14 `SETS/` files reference the old behaviour** (the `Beghein14OceanicCooling*` and `Q195*_Q3_phase1*` / `Q195Run3Phase1_smoke` families — 11 `.txt` scenarios set `ani_fstrain = 3` and/or `ani_T_threshold`, and 2 `.c` callbacks mention `ani_T_threshold` in comments); all were purpose-built around the freeze and are obsoleted by this change — they are removed or migrated (tasks group 3). The two `AniFstrainT_Threshold_Freeze_*` GTests are also removed.

- **Add a temperature-dependent δ-relaxation.** When the local dislocation-creep deformation rate is low, δ relaxes toward the isotropic limit (δ → 1) on the Boneh et al. 2021 DiSRX timescale `τ_relax = L_relax / V`, grain-boundary-migration velocity `V = M₀·exp(−Q/RT)·μ·b²·Δρ` (Boneh 2021 Eqs. 1–4: `Q = 133 ± 18 kJ/mol`, `M₀ = 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹` ≡ `m³·N⁻¹·s⁻¹` — Boneh 2021 prints "m³/s·J" but that is a typo; the Boneh et al. 2017 companion paper's Appendix A prints the dimensionally-correct `m³/(N·s)`, and the corrected value reproduces Boneh's "weeks-to-years at ~1300 °C" anchor; `μ = 50 GPa`, `b = 0.6 nm`). `V` uses only the strain-energy driving force `F_s = μb²Δρ`, dropping Boneh's surface-energy term `F_b = 3γ/d` — justified because `F_s ≫ F_b` for mantle grain sizes ≥ 1 mm (Boneh 2021 §4.2, Fig. 8).

- **δ becomes a per-marker integrated state variable.** It is currently recomputed every step as a memoryless function of `FS_AR` ([AnisoFactorEvolv](MDLIB/AnisotropyRoutines.c#L85)). The relaxation needs history, so a new per-marker field carries the relaxing δ (or the equivalent CPO-strength index `M`), updated each step by the **analytic exponential update** `δ_new = δ_eq + (δ_old − δ_eq)·exp(−Δt/τ_relax)` — unconditionally stable, cannot overshoot, no dependence on the (anisotropy-unaware) `dt` limiter.

- **`L_relax` is a per-phase length-scale parameter** (the characteristic DiSRX boundary-migration distance ≈ annealed grain size). It is sourced from the existing per-phase `gs_ref` — **grain-size evolution (`gs` flow law) is NOT required**. GSE off → `L_relax = gs_ref` (per-phase constant); GSE on → it can track the evolving `mesh->d_n`.

- **`Δρ` (dislocation-density driving force) uses a proxy** built from existing per-marker fields (accumulated dislocation-creep strain `strain_pwl` and/or deviatoric stress). MDOODZ does not track dislocation density; the proxy form and its caveats are settled in design.md.

- **Deprecate `ani_T_threshold`.** The parameter and the [RheologyParticles.c:603-621](MDLIB/RheologyParticles.c#L603-L621) F-freeze block are removed; the case-7 default in [FlowLaws.c](MDLIB/FlowLaws.c#L1638) is replaced with an `ani_relax_length` default.

- **Backward compatibility.** `ani_fstrain == 0/1/2` unchanged byte-for-byte. The relaxation is reachable only via `ani_fstrain == 3`.

- **Caveats are commented in the code.** Every modelling approximation — the constant-`L_relax` assumption and the mid-T / episodic-loading band where its ≤17 % error lands; the `Δρ` proxy (including the "absolute stored strain vs. boundary contrast" mismatch and the no-recovery monotonicity of `strain_pwl`); the "δ relaxes on the DiSRX grain-growth timescale" modelling leap; the dropped surface-energy term (`ΣF ≈ F_s`); and the resolved `M₀` units (with the Boneh-2021-typo note) — is documented *inline at the relevant code site*, not only in a notes file.

## Capabilities

### New Capabilities

- `aniso-delta-relaxation`: the temperature-dependent kinetic relaxation of the viscous-anisotropy factor δ — the `τ_relax = L_relax/V` law (Boneh et al. 2021 DiSRX kinetics), δ as a per-marker integrated state variable, the analytic-exponential integration scheme, the `L_relax` and `Δρ` inputs, and the analytical-unit-test validation suite.

### Modified Capabilities

- `ci-finite-strain-anisotropy`: the `ani_fstrain` enumeration changes — `ani_fstrain == 3` is redefined from "T-threshold F-freeze" to "`ani_fstrain == 2` + δ-relaxation", `ani_T_threshold` is removed from the parameter set, and the `AniFstrainT_Threshold_Freeze_Sentinel/Active` test requirements are removed.

## Impact

- **Code under change** (MDLIB): `AnisotropyRoutines.c` (`AnisoFactorEvolv` signature + `ani_fstrain==3` dispatch, new relaxation routine `DeltaRelaxationTau` / `DeltaRhoProxy` taking per-phase constants as parameters, the per-marker δ update inside `FiniteStrainAspectRatio`, two stale comment blocks); `RheologyParticles.c` (remove the F-freeze block ~603-621); `mdoodz-private.h` (two new `markers` struct fields — `aniso_delta`, `aniso_delta_fs_prev`; updated `DeltaRelaxationTau` prototype); `include/mdoodz.h` (`mat_prop`: add `ani_relax_length[20]` plus the per-phase relaxation-kinetics fields `ani_relax_Q/M0/mu/b/drho_min/drho_max/eps_ref[20]`; remove `ani_T_threshold[20]` and its doc comment); `InputOutput.c` (parse `ani_relax_length`, drop `ani_T_threshold`, restart I/O for the two new marker fields); `FlowLaws.c` (`ReadDataAnisotropy` — case-7 default + `LOG_INFO`, and the per-phase relaxation-kinetics constants with their caveat comments, mirroring the `ani_fstrain == 2` `aniso_delta_fn` wiring); `MemoryAllocFree.c` (alloc/free the two marker fields + two new `mesh->aniso_delta_n/_s` grid fields); `ParticleRoutines.c` (init + reseeding/inflow copy); `FreeSurface.c` (init the new fields for deposited sediment markers); `HDF5Output.c` (output `aniso_delta` as a centroid field). The new-marker-field plumbing mirrors the existing `strain_pwl` pattern; the grid-field plumbing mirrors `FS_AR_n`/`FS_AR_s`; the per-phase relaxation-kinetics database mirrors `ReadDataPowerLaw` / the `ani_fstrain == 2` δ-forms.
- **Scenario files** (`SETS/`): the 11 `ani_fstrain = 3` `.txt` scenarios and the 2 `.c` callbacks with `ani_T_threshold` comments are removed or migrated.
- **Tests**: remove `AniFstrainT_Threshold_Freeze_Sentinel` and `AniFstrainT_Threshold_Freeze_Active` from `TESTS/AnisotropyBenchmarkTests.cpp` (they test the deprecated freeze and are currently crashing on `aniso_db = 0`). Add analytical unit tests for the relaxation limits: cold freeze (`τ→∞`), hot relaxation on the Boneh timescale, exponential-decay half-life + a `dt`-convergence sub-test, production-only limit recovers `ani_fstrain == 2`, `L_relax` monotonicity, `ani_fac_max` clamp honoured, isotropic fixed point. Update `TESTS/AnalyticalSolutions.md`.
- **Documentation**: a notes file under `misc/aniso_fstrain/notes/` documenting the Boneh kinetics, the constant-`L_relax` regime analysis, and the validation plan.
- **CI runtime**: +~9 short analytical unit tests on a 1-cell/homogeneous box ≈ <1 s; net change near zero (2 tests removed).
- **Risk**: **medium**. The structural change is two new per-marker state variables plus two new grid fields (struct + alloc/free/init/reseed/sediment-init/restart-I/O/HDF5) — a well-trodden pattern (`strain_pwl` and `FS_AR_n` are the templates), but restart I/O is positional so the read and write sites must be edited together. The analytic-exponential integration removes the numerical-stability (FS_AR-class runaway) risk. The BREAKING redefinition does touch 14 `SETS/` files (handled in tasks group 3), but no run *breaks* — orphaned `ani_T_threshold` lines are silently ignored by the parser, and `ani_fstrain = 3` silently changes meaning, so the migration is a cleanliness/correctness task rather than a build-breaker. The `Δρ` proxy physics and the `M₀` units are resolved in design.md.
- **Out of scope**:
  - Modelling dislocation density `ρ` as a true state variable with its own `dρ/dt` ODE. A proxy is used; explicit ρ-evolution is future work.
  - A DiSRX-specific grain-size-evolution model for `L_relax` (modelling the annealed grain size the olivine grows to). `L_relax` is a per-phase constant from `gs_ref`; the rigorous annealed-grain-size model is future work.
  - Fabric-type (A/B/C/E-type) tracking. The relaxation acts on the scalar δ *magnitude* only — it relaxes toward isotropic, not toward a different fabric type.
  - Modifying or extending MDOODZ's GSE (`gs` flow law / paleowattmeter). The investigation found MDOODZ's GSE implements dynamic-recrystallization kinetics — the wrong process for DiSRX annealing — but fixing GSE is a separate concern. This change does NOT route through GSE.
  - Director (orientation) evolution. `ani_fstrain == 3` relaxes the anisotropy *magnitude* δ; the director angle (`nx/nz`) evolution is unchanged.
  - Re-deriving or re-validating the `ani_fstrain == 2` Hansen δ(`FS_AR`) calibration. The relaxation composes with `== 2` as-is.
