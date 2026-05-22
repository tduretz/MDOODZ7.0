## Context

`ani_fstrain == 3` currently means "the `ani_fstrain == 2` δ-dispatch **plus** an upstream `T < ani_T_threshold` deformation-gradient freeze" ([RheologyParticles.c:603-621](MDLIB/RheologyParticles.c#L603-L621), [FlowLaws.c:1644-1649](MDLIB/FlowLaws.c#L1644-L1649)). A multi-round investigation established that the freeze is the wrong model: the 1100 °C threshold is misattributed, the freeze was a misinterpretation of an intended *kinetic-relaxation* directive, and a binary T-threshold cannot represent the actual physics — when deformation slows, olivine CPO *relaxes* toward isotropic over a finite, temperature-dependent timescale via discontinuous static recrystallization (DiSRX), per Boneh et al. 2021 (*G3*).

Current state of the relevant code paths (verified by consistency review):
- δ for `ani_fstrain ∈ {1,2,3}` is **memoryless**: `AnisoFactorEvolv` ([AnisotropyRoutines.c:85](MDLIB/AnisotropyRoutines.c#L85)) recomputes δ from the instantaneous `FS_AR` every call; `==3` currently shares the `==2` arm verbatim ([AnisotropyRoutines.c:90](MDLIB/AnisotropyRoutines.c#L90)).
- `AnisoFactorEvolv` is called in **6 places** inside `UpdateAnisoFactor` ([AnisotropyRoutines.c:598](MDLIB/AnisotropyRoutines.c#L598)) — three average modes × {centroid, vertex} — each wrapped in per-phase `phase_perc` weighting and harmonic/geometric averaging. The Stokes solver consumes **both** `mesh->aniso_factor_n` (centroids) and `mesh->aniso_factor_s` (vertices).
- `FiniteStrainAspectRatio` ([RheologyParticles.c:689](MDLIB/RheologyParticles.c#L689)) computes per-marker `FS_AR` into a *local* array, P2G's it to `mesh->FS_AR_n` **and** `mesh->FS_AR_s`, then frees the local array. There is no marker `FS_AR` field.
- Timestep order (`Main_DOODZ.c`): `FiniteStrainAspectRatio` (676) → `UpdateAnisoFactor` (728) → Stokes solve (~805–1000) → `DeformationGradient` / F-update (1293) → `AccumulatedStrainII` updates `strain_pwl` (1231). So at the line-676 anisotropy pass, `F`/`FS_AR` and `strain_pwl` reflect the **end of the previous timestep**; `particles->T` is current.
- The `markers` struct ([mdoodz-private.h:30-42](MDLIB/mdoodz-private.h#L30)) carries `Fxx/Fxz/Fzx/Fzz`, `strain_pwl`, `d`, `T`, `sxxd/szzd/sxz` — no anisotropy state, no dislocation density.
- MDOODZ's "GSE" is a steady-state paleowattmeter, not a transient grain-growth ODE — the *wrong* kinetics for DiSRX annealing, so this change does **not** route through it.

Constraints: must not change `ani_fstrain ∈ {0,1,2}` behaviour byte-for-byte; must be numerically robust (the FS_AR path has a documented unbounded-runaway history); the new per-marker state must survive advection, reseeding, sediment deposition, and restart.

## Goals / Non-Goals

**Goals:**
- Replace the T-threshold freeze with a temperature-dependent kinetic relaxation of δ toward the isotropic limit, on the Boneh et al. 2021 DiSRX timescale.
- Make δ a history-carrying per-marker state variable, integrated with an unconditionally-stable analytic scheme.
- Require **no** dependence on grain-size evolution (`gs` flow law) and **no** new dislocation-density state variable.
- Keep `ani_fstrain ∈ {0,1,2}` byte-identical.
- Document every modelling caveat inline at the code site.

**Non-Goals:**
- Modelling dislocation density `ρ` as a true state variable with its own ODE (a proxy is used).
- A DiSRX-specific grain-size-evolution model for `L_relax` (it is a per-phase length parameter).
- Fabric-type (A/B/C/E) tracking; director-angle evolution; re-validating the `ani_fstrain == 2` Hansen calibration.
- Fixing or extending MDOODZ's paleowattmeter GSE.
- Supporting cells that mix `ani_fstrain == 3` markers with `ani_fstrain ∈ {0,1,2}` markers — single-`ani_fstrain`-per-phase is assumed (the existing `UpdateAnisoFactor` per-phase machinery is reused, so mixed-phase cells average per-phase contributions, but no special handling is added).

## Decisions

### D1 — δ relaxation = two-field operator-split with an analytic exponential update

Per `ani_fstrain == 3` marker, carry **two** new scalar fields: `aniso_delta` (the relaxed δ state) and `aniso_delta_fs_prev` (the previous step's finite-strain δ). Each timestep:

```
delta_FS  = aniso_delta_fn(FS_AR)                          // the ani_fstrain==2 value (production)
delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev) // inherit the ==2 increment
tau_relax = L_relax / V                                    // V = M0*exp(-Q/RT)*mu*b^2*d_rho  (Fs only)
delta_new = 1.0 + (delta_prod - 1.0) * exp(-dt / tau_relax) // analytic relaxation toward 1
delta_new = fmin(delta_new, ani_fac_max)                   // clamp
delta_new = fmax(delta_new, 1.0)                           // isotropic floor
aniso_delta = delta_new ;  aniso_delta_fs_prev = delta_FS   // store for next step
```

**Why this scheme:**
- **Production = exactly the calibrated `ani_fstrain == 2` form, with no new parameter.** δ inherits the *increment* of `aniso_delta_fn(FS_AR)`, so while deforming, δ grows exactly as `ani_fstrain == 2` would.
- **Cold limit recovers `ani_fstrain == 2` exactly.** As `τ_relax → ∞`, `exp(−dt/τ) → 1`, so `delta_new = delta_prod`; the increments telescope to `aniso_delta = aniso_delta_fn(FS_AR)`.
- **Quiescent-hot limit is pure relaxation.** When `delta_FS == aniso_delta_fs_prev`, `delta_prod = aniso_delta` and `delta_new = 1 + (aniso_delta − 1)·exp(−dt/τ)`.
- **Analytic exponential** — unconditionally stable for any `dt`, cannot overshoot below `δ_eq = 1`, no dependence on the (anisotropy-unaware) Courant limiter.

**One-step lag (note, not a defect):** at the line-676 anisotropy pass, `FS_AR` and `strain_pwl` are end-of-previous-step values (see Context). The relaxation is a slow process, and this is consistent with how `ani_fstrain == 2` already consumes `FS_AR`, so the lag is acceptable — but it is stated as an inline comment.

**Alternative considered — single-field two-sided relaxation** `dδ/dt = (δ_FS − δ)·k_prod − (δ − 1)·k_relax`. Rejected: it needs a *production rate* `k_prod` (a new per-phase modelling parameter), trading one scalar field for a new free parameter and a less-obvious cold-limit equality. The two-field scheme keeps production *identically* the already-calibrated `ani_fstrain == 2` physics.

### D2 — Store δ directly (not the CPO index `M`)

`AnisoFactorEvolv` and `aniso_delta_fn` already work in δ; the relaxation target is `δ_eq = 1`. Storing δ avoids an extra `M ↔ δ` mapping and avoids coupling the field to the per-mineral `slope`.

### D3 — `Δρ` proxy = accumulated dislocation-creep strain (`strain_pwl`)

`Δρ` (the dislocation-density driving force in `V`) is proxied by the per-marker `strain_pwl` field, normalised by a reference strain (O2). Rationale: MDOODZ tracks no dislocation density; `strain_pwl` is monotonic, already per-marker, and its "sticky" character is *correct* for the post-deformation use case. **Two caveats are documented inline (D7):** (a) `strain_pwl` cannot decrease during recovery (monotonicity); (b) Boneh's `Δρ` is the dislocation-density *contrast across a migrating boundary*, whereas `strain_pwl` proxies a marker's *absolute* stored strain energy — they coincide only when the marker is significantly more strained than its neighbourhood. A `τ_II`-piezometer proxy was considered and rejected as the primary choice (it collapses exactly when deformation stops); a piezometer-cap refinement is left open.

### D4 — `L_relax` is a per-phase parameter `ani_relax_length`, defaulting to `gs_ref`

Add `ani_relax_length[20]` to `mat_prop`, parsed by `InputOutput.c` with sentinel default `−1` → "use `gs_ref[phase]`" (mirrors the `aniso_d_threshold` sentinel pattern at [InputOutput.c:1461-1462](MDLIB/InputOutput.c#L1461)). Grain-size evolution is **not** required. The case-7 (`aniso_db = 7`) default in `ReadDataAnisotropy` changes from `ani_T_threshold = 1100 °C` to `ani_relax_length = −1`.

### D5 — Boneh constants are per-phase, in the `aniso_db` material database

The relaxation kinetics constants — `Q = 133 ± 18 kJ/mol`, `μ = 50 GPa`, `b = 0.6 nm`, `M₀ = 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹` (≡ `m³·N⁻¹·s⁻¹`), and the `Δρ`-proxy parameters — are **per-phase `mat_prop` fields** (`ani_relax_Q[20]`, `ani_relax_M0[20]`, `ani_relax_mu[20]`, `ani_relax_b[20]`, `ani_relax_drho_min[20]`, `ani_relax_drho_max[20]`, `ani_relax_eps_ref[20]`), populated by `ReadDataAnisotropy` in `FlowLaws.c` — the *same function*, switched on the *same `aniso_db` index*, that sets the `ani_fstrain == 2` `aniso_delta_fn` pointer. They are genuinely mineral-specific (Boneh's `Q`, `M₀` are *olivine* grain-boundary-mobility values; the `Δρ` range is the *olivine* dislocation-density range), so the `aniso_db` database is their natural home — and a future calcite/quartz DiSRX calibration overrides them per `aniso_db` case, exactly as the δ-forms extend. `ReadDataAnisotropy` sets a Boneh-olivine default for all phases before the `switch`; the relaxation routine in `AnisotropyRoutines.c` takes the per-phase constants as parameters (`DeltaRhoProxy` / `DeltaRelaxationTau` no longer hold module-global `static const`). This keeps all anisotropy calibration in one place (`FlowLaws.c`) and follows the established `ReadDataPowerLaw` / `ReadDataGSE` convention. **`M₀` units note:** Boneh 2021 prints "m³/s·J", dimensionally inconsistent with its own Eq. 1 (`V = M_b·ΣF` would give a rate, not a velocity); the Boneh et al. 2017 companion paper Appendix A prints the correct `m³/(N·s)` ≡ `m⁴·J⁻¹·s⁻¹`, and with that unit the formula reproduces Boneh's "weeks-to-years at ~1300 °C" anchor (10 µm → 1 mm in ~0.23–23 yr for `Δρ = 10¹³–10¹¹ m⁻²`). This is committed, not open.

### D6 — Placement: per-marker δ update **inside** `FiniteStrainAspectRatio`; relaxed δ flows to the grid via dedicated grid fields

The `RheologyParticles.c:603-621` F-freeze block is **deleted**.

The new per-marker δ-relaxation update is added **inside** `FiniteStrainAspectRatio` ([RheologyParticles.c:689](MDLIB/RheologyParticles.c#L689)) — within or immediately after the existing per-marker `FS_AR` `k`-loop, **before** the local `FS_AR` array is freed (the per-marker `FS_AR` exists only there; there is no marker `FS_AR` field). For each `ani_fstrain == 3` marker it reads `materials->aniso_delta_fn[phase]`, `FS_AR`, `T`, `strain_pwl`, `L_relax`; writes `aniso_delta`, `aniso_delta_fs_prev`.

The relaxed δ reaches the Stokes solver through **new grid fields `mesh->aniso_delta_n` and `mesh->aniso_delta_s`**, populated by P2G from the marker `aniso_delta` — both centroid and vertex, mirroring exactly how `FiniteStrainAspectRatio` already P2G's `FS_AR_n`/`FS_AR_s` ([RheologyParticles.c:749-750](MDLIB/RheologyParticles.c#L749)). The grid fields are allocated/freed in `MemoryAllocFree.c` mirroring `FS_AR_n`.

`AnisoFactorEvolv` then gets **one new argument** — the P2G'd relaxed δ at that grid point — and its `ani_fstrain == 3` arm returns that value instead of recomputing `aniso_delta_fn(FS_AR)`. The prototype ([mdoodz-private.h:422](MDLIB/mdoodz-private.h#L422)) and **all 6 call sites** inside `UpdateAnisoFactor` ([AnisotropyRoutines.c:628,633,638,698,703,708](MDLIB/AnisotropyRoutines.c#L628)) are updated. This routes the relaxed δ *through* the existing per-phase `phase_perc` weighting and harmonic/geometric averaging — so mixed-phase cells, the three averaging modes, and the centroid+vertex grids all stay coherent.

**Why not "P2G `aniso_delta` straight onto `aniso_factor_n`":** that bypasses `UpdateAnisoFactor`'s per-phase weighting and would collide with the other phases' contributions to the same cell, and it ignores the vertex grid `aniso_factor_s` entirely. The dedicated-grid-field + `AnisoFactorEvolv`-argument route is the minimal change that composes with the existing grid assembly.

For `ani_fstrain ∈ {0,1,2}` the new `AnisoFactorEvolv` argument is unused and the existing path is byte-identical.

### D7 — Caveats are inline comments at the code site (hard requirement)

Per the spec, the following caveats are commented *in the source*, not only in the notes file: (a) the constant-`L_relax` assumption + the mid-T / episodic-loading band where its ≤17 % error lands — at the `L_relax` read site; (b) the `strain_pwl`-as-`Δρ` proxy — *both* the monotonicity (no-recovery) caveat *and* the absolute-vs-boundary-contrast mismatch — at the `Δρ` computation; (c) the "δ relaxes on the DiSRX grain-growth timescale" modelling leap (Boneh links DiSRX→CPO-weakening only qualitatively; the relaxation ODE is the model's own) — at the relaxation routine header; (d) the `M₀` units (resolved value + the Boneh-2021-typo / Boneh-2017-correction) and the dropped surface-energy term (`ΣF ≈ F_s`, valid because `F_s ≫ F_b` for mantle grain sizes) — at the `V` computation.

## Risks / Trade-offs

- **[New per-marker / grid-field plumbing is incomplete]** → mirror the `strain_pwl` template across struct decl, alloc/free, `PartInit`, **`FreeSurface.c` sediment-marker init**, reseeding+inflow copy, restart I/O, **and HDF5 output**; mirror the `FS_AR_n` template for the `mesh->aniso_delta_n/_s` grid fields. Restart I/O is *positional* `fread`/`fwrite` — read and write sites edited in the same commit; add a restart round-trip smoke check. `aniso_delta` is dimensionless → correctly absent from the `InputOutput.c` scaling loops (matches `strain_pwl`).
- **[`Δρ` proxy is physically approximate]** → bounded (the ρ-proxy investigation: caveats land only in the mid-T band and under episodic loading, ≤17 %-class on δ_final, never qualitative); both caveats documented inline (D7b); the piezometer-cap refinement path is left open.
- **[BREAKING redefinition of `ani_fstrain == 3`]** → 14 `SETS/` files reference the old behaviour, but **no run breaks**: orphaned `ani_T_threshold` lines are silently ignored by `ReadMatProps`, and the 2 `.c` callbacks reference `ani_T_threshold` only in comments. The redefinition silently changes the *physics* of those scenarios, several of which were purpose-built for the freeze — they are removed or migrated (tasks group 3), and the two `AniFstrainT_Threshold_Freeze_*` GTests (already crashing on `aniso_db = 0` because the iter-#784 dispatch change made their inline comments stale) are removed.
- **[Numerical runaway, as the FS_AR path once had]** → the analytic exponential update is unconditionally stable and cannot overshoot; δ is hard-clamped to `[1, ani_fac_max]`.
- **[Mixed-`ani_fstrain` cells]** → out of scope (Non-Goals); the reused `UpdateAnisoFactor` per-phase machinery averages per-phase contributions, so such cells do not crash, but the result is not specially validated.

## Migration Plan

1. Add the two marker fields + `ani_relax_length` param + the two `mesh->aniso_delta_n/_s` grid fields; wire alloc/free/init/sediment-init/reseed/restart-I/O/HDF5.
2. Add the relaxation routine + the per-marker update inside `FiniteStrainAspectRatio`; change the `AnisoFactorEvolv` signature + 6 call sites; delete the F-freeze block; remove `ani_T_threshold` from `mat_prop` (+ its doc comment), parsing, the `FlowLaws.c` case-7 default, and the two `AnisotropyRoutines.c` comment blocks.
3. Remove the two `AniFstrainT_Threshold_Freeze_*` GTests; add the analytical-unit-test cases.
4. Remove or migrate the 14 `SETS/` files; scrub `ani_T_threshold` comment refs from the 2 `.c` callbacks.
5. Build (`make build-dev`), run the full anisotropy + rheology CI suites — `ani_fstrain ∈ {0,1,2}` tests must pass unchanged.
6. Restart round-trip smoke check; manual ≥2-step run of an existing `ani_fstrain = 1` scenario to confirm no regression; `grep` confirms zero `ani_T_threshold` across `MDLIB/`, `TESTS/`, `SETS/`.
7. Write the research note (`misc/aniso_fstrain/notes/`); update `skill-anisotropy` if it documents the `ani_fstrain` enumeration or output fields.

**Rollback**: `git revert` — the change is additive except for the freeze removal, and the affected `SETS/` scenarios were obsolete.

## Open Questions

- **O1 — `M₀` / `V` units — RESOLVED.** The consistency review pinned this: `M₀ = 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹` (≡ `m³·N⁻¹·s⁻¹`); Boneh 2021's printed "m³/s·J" is a typo, corrected per Boneh et al. 2017 Appendix A; the formula then yields `V` in m/s, `τ_relax` in s, and reproduces Boneh's published "weeks-to-years at ~1300 °C" anchor. Implementation must encode this value + units + the typo note (D5, task 4.1).
- **O2 — `Δρ` normalisation.** The reference strain that normalises `strain_pwl → Δρ` is a calibration choice; pick a defensible default (tie to the Boneh dislocation-density range `10¹¹–10¹³ m⁻²` and the Hansen saturation-strain scale) and expose it as a commented constant or a per-phase parameter.
- **O3 — GSE-on path for `L_relax`.** Whether, when the `gs` flow law is active, `L_relax` should track `mesh->d_n` instead of `gs_ref`. Deferred; the `−1` sentinel makes adding this later non-breaking. No task — tracked here only.
