## 1. Per-marker and grid state plumbing

- [x] 1.1 Add `*aniso_delta` and `*aniso_delta_fs_prev` (`double*`) to the `markers` struct in `MDLIB/mdoodz-private.h` (~lines 30-42), alongside `strain_pwl`
- [x] 1.2 Allocate both fields in `PartAlloc` (~`MemoryAllocFree.c:224`) and free them in `PartFree` (~`:298`), unconditional, mirroring `strain_pwl`
- [x] 1.3 Initialise both fields in `PartInit` (`MDLIB/ParticleRoutines.c` ~line 897) — `aniso_delta = 1.0`, `aniso_delta_fs_prev = 1.0` (isotropic start)
- [x] 1.4 Copy both fields on reseeding and inflow injection — `AssignMarkerProperties` (~`ParticleRoutines.c:705`) and `AssignMarkerPropertiesInflow` (~`:776`), mirroring `strain_pwl`
- [x] 1.5 Initialise both fields for newly-deposited sediment markers in the free-surface sedimentation routine (`MDLIB/FreeSurface.c` ~line 148, alongside the `strain_pwl = 0.0` block) — `aniso_delta = 1.0`, `aniso_delta_fs_prev = 1.0` (fresh sediment is isotropic)
- [x] 1.6 Add both fields to the restart `fread` (`LoadBreakpointParticles` ~`InputOutput.c:225`) and `fwrite` (`MakeBreakpointParticles` ~`:736`) — edit read and write sites together at matching positions; both are dimensionless so they get NO entry in the `InputOutput.c` scaling loops (matches `strain_pwl`)
- [x] 1.7 Add `aniso_delta` to the HDF5 centroid output (`MDLIB/HDF5Output.c` — ~6 edit sites mirroring `strain_pwl`: 2 buffer decls, `DoodzCalloc`+`P2Mastah`, `DoodzMalloc`+`DoubleToFloat`, `AddFieldToGroup`, 2 frees)
- [x] 1.8 Add grid fields `mesh->aniso_delta_n` (centroid) and `mesh->aniso_delta_s` (vertex); allocate/free in `MemoryAllocFree.c` mirroring `FS_AR_n` / `FS_AR_s`

## 2. Per-phase parameter (`ani_relax_length`, remove `ani_T_threshold`)

- [x] 2.1 Add `ani_relax_length[20]` to `mat_prop` in `MDLIB/include/mdoodz.h`; remove `ani_T_threshold[20]` AND its doc-comment block (~lines 195-203)
- [x] 2.2 In `MDLIB/InputOutput.c`, parse `ani_relax_length` via `ReadMatProps` with sentinel default `-1.0` (→ "use `gs_ref`"), divided by `scaling.L`; remove the `ani_T_threshold` parse (~lines 1471-1472)
- [x] 2.3 In `MDLIB/FlowLaws.c` `ReadDataAnisotropy`, replace the case-7 `ani_T_threshold = 1100 °C` default (~lines 1644-1649) with `ani_relax_length = -1` (inherit `gs_ref`); update the `LOG_INFO` line accordingly
- [x] 2.4 Resolve the effective `L_relax` per phase at the use site: `ani_relax_length >= 0 ? ani_relax_length : gs_ref[phase]`

## 3. Remove the temperature-threshold freeze and obsolete scenarios

- [x] 3.1 Delete the `T < ani_T_threshold` F-freeze block in `MDLIB/RheologyParticles.c` (~lines 603-621) — comment + `if` + `continue`; confirm the only block-local variable (`phase_k`) is used nowhere else
- [x] 3.2 Remove all remaining `ani_T_threshold` references in `MDLIB/` and update **both** stale comment blocks in `MDLIB/AnisotropyRoutines.c` — the dispatch-doc block (~lines 61-70) AND the second block inside `AnisoFactorEvolv` (~lines 90-97, which describes the now-removed freeze and misattributes it to `ParticleRoutines.c`)
- [x] 3.3 Remove or migrate the 11 `ani_fstrain = 3` `SETS/*.txt` scenarios (`Beghein14OceanicCooling*`, `Q195v3_Q3_phase1_*`, `Q195v4_Q3_phase1_smoke`, `Q195Run3Phase1_smoke`); strip dead `ani_T_threshold` lines; scrub the `ani_T_threshold` comment references from `SETS/Beghein14OceanicCooling.c` and `SETS/Q195v3_smoke.c`

## 4. Relaxation routine and per-marker δ update

- [x] 4.1 Add the Boneh et al. 2021 relaxation constants as **per-phase `mat_prop` fields** (`ani_relax_Q/M0/mu/b/drho_min/drho_max/eps_ref[20]` in `MDLIB/include/mdoodz.h`), populated by `ReadDataAnisotropy` in `MDLIB/FlowLaws.c` — Boneh-olivine default for all phases before the `switch`, overridable per `aniso_db` case, mirroring the `ani_fstrain == 2` `aniso_delta_fn` wiring and the `ReadDataPowerLaw` convention. Values: `Q = 133 ± 18 kJ/mol`, `μ = 50 GPa`, `b = 0.6 nm`, `M0 = 2.0e-11` units `m^4 J^-1 s^-1` (≡ `m^3 N^-1 s^-1`) — inline-comment that Boneh 2021's printed "m³/s·J" is a typo, corrected per Boneh et al. 2017 Appendix A; the corrected value reproduces Boneh's "weeks-to-years at ~1300 °C" anchor. `DeltaRhoProxy` / `DeltaRelaxationTau` in `AnisotropyRoutines.c` take the per-phase constants as parameters (no module-global `static const`)
- [x] 4.2 Implement the `Δρ` proxy from per-marker `strain_pwl`, normalised by a reference strain (design O2 — pick a defensible default tied to the Boneh `Δρ ≈ 10¹¹–10¹³ m⁻²` range); finite and non-negative
- [x] 4.3 Implement the relaxation routine: `V = M0*exp(-Q/RT)*mu*b^2*Δρ` (strain-energy term `F_s` only — `ΣF ≈ F_s`, valid because `F_s ≫ F_b` for mantle grain sizes), `tau_relax = L_relax / V`
- [x] 4.4 Implement the per-marker δ update **inside `FiniteStrainAspectRatio`** (`MDLIB/RheologyParticles.c` ~line 689), in/after the existing per-marker `FS_AR` `k`-loop and before the local `FS_AR` array is freed, gated to `ani_fstrain == 3` phases: `delta_FS = aniso_delta_fn(FS_AR)`; `delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev)`; `delta_new = 1 + (delta_prod - 1)*exp(-dt/tau_relax)`; clamp to `[1, ani_fac_max]`; store `aniso_delta` and `aniso_delta_fs_prev`
- [x] 4.5 P2G the per-marker `aniso_delta` onto **both** `mesh->aniso_delta_n` (centroid) and `mesh->aniso_delta_s` (vertex), mirroring the two `P2Mastah` calls that `FiniteStrainAspectRatio` already makes for `FS_AR_n` / `FS_AR_s` (~`RheologyParticles.c:749-750`)
- [x] 4.6 Add one new argument to `AnisoFactorEvolv` (the P2G'd relaxed δ at the grid point); update the prototype (`mdoodz-private.h` ~line 422) and **all 6 call sites** in `UpdateAnisoFactor` (`AnisotropyRoutines.c` ~lines 628, 633, 638, 698, 703, 708); the `ani_fstrain == 3` arm returns the passed-in relaxed δ instead of recomputing `aniso_delta_fn(FS_AR)` — so it flows through the existing per-phase `phase_perc` weighting and harmonic/geometric averaging. `ani_fstrain ∈ {0,1,2}` ignore the new argument and stay byte-identical
- [x] 4.7 Require `aniso_db > 0` for `ani_fstrain == 3`; extend the `ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1` auto-default (both sites: `InputOutput.c` ~1522-1523 and `Main_DOODZ.c` ~119-120) to also cover `== 3`

## 5. Inline caveat comments (design D7 — hard requirement)

- [x] 5.1 Comment the constant-`L_relax` assumption + the mid-T / episodic-loading band where its bounded (~17 %) error lands — at the `L_relax` resolution site (task 2.4)
- [x] 5.2 Comment the `strain_pwl`-as-`Δρ` proxy — **both** the monotonicity (no-recovery) caveat AND the "absolute stored strain vs. boundary contrast" mismatch (Boneh's `Δρ` is a contrast across a migrating boundary) — at the `Δρ` computation (task 4.2)
- [x] 5.3 Comment the "δ relaxes on the DiSRX grain-growth timescale" modelling leap (Boneh links DiSRX→CPO-weakening only qualitatively; the relaxation ODE is the model's own construction) — at the relaxation routine header
- [x] 5.4 Comment, at the `V` computation: the resolved `M0` units (`m⁴·J⁻¹·s⁻¹`, with the Boneh-2021-typo / Boneh-2017-correction note), AND the dropped surface-energy term (`ΣF ≈ F_s`, valid because `F_s ≫ F_b` for mantle grain sizes ≥ 1 mm)

## 6. Tests

- [x] 6.1 Remove `AniFstrainT_Threshold_Freeze_Sentinel` (~`AnisotropyBenchmarkTests.cpp:1748`) and `AniFstrainT_Threshold_Freeze_Active` (~`:1797`), plus the `anaDeltaTThresholdSentinel` anonymous-namespace helper (~lines 1738-1746) and its comment header (~1714-1737). They share the `AnisotropyBenchmark` fixture / `mutateAniFstrainInput` / `AniFstrainEvolution.txt` with the other tests — there is NO test-specific fixture, `.txt`, or per-test CMake entry to remove
- [x] 6.2 Add `AniFstrainRelaxColdFreeze` — cold marker, `τ_relax ≫ run time`, δ constant to FP precision
- [x] 6.3 Add `AniFstrainRelaxHot` — hot marker, no production, δ(t) matches `1 + (δ₀−1)·exp(−t/τ_relax)` on the Boneh timescale; include a sub-assertion that `τ_relax(T₂)/τ_relax(T₁) = exp(Q/R·(1/T₁−1/T₂))` (the Arrhenius ratio)
- [x] 6.4 Add `AniFstrainRelaxHalfLife` — constant `T`, `Δρ`, `L_relax`: δ reaches `1 + (δ₀−1)/2` at `t = τ_relax·ln2`; include a `Δt`-convergence sub-assertion AND a large-`Δt` no-overshoot sub-assertion (`Δt ≫ τ` ⇒ `1 ≤ δ_new ≤ δ_old`)
- [x] 6.5 Add `AniFstrainRelaxColdLimitEqualsFstrain2` — `ani_fstrain == 3` at cold T tracks `aniso_delta_fn(FS_AR)` identically to `ani_fstrain == 2`
- [x] 6.6 Add `AniFstrainRelaxLengthMonotonic` — larger `L_relax` ⇒ slower relaxation
- [x] 6.7 Add `AniFstrainRelaxClampAndFixedPoint` — `ani_fac_max` clamp honoured; `δ = 1` is a fixed point
- [x] 6.8 Add `AniFstrainRelaxHistoryDependence` — two markers reaching the same `FS_AR` via different thermal histories (one hot, one cold) have different stored `aniso_delta`
- [x] 6.9 Reuse the shared `AniFstrainEvolution.txt` base fixture via `MutateInput` where possible; confirm the analytical tests run with the `gs` flow law OFF (covers the spec's "runs with GSE disabled" scenario) — make that explicit in the test setup
- [x] 6.10 Update `TESTS/AnalyticalSolutions.md`: document the `ani_fstrain == 3` enumeration value and the new test cases; remove the T-threshold-freeze rows

## 7. Documentation

- [x] 7.1 Write the research note under `misc/aniso_fstrain/notes/` — the Boneh et al. 2021 DiSRX kinetics, the `τ_relax = L_relax/V` derivation and its modelling assumptions (the `τ = L/V` construction, the `ΣF ≈ F_s` approximation, the `Δρ` proxy), the constant-`L_relax` regime analysis (frozen / instantaneous / caught-partway bands), and the tiered validation plan (analytical unit tests as the CI gate; Boneh 2017 / Speciale 2020 lab-data reproduction as future research-grade checks)
- [x] 7.2 If `skill-anisotropy` (or any visualisation field list) documents the `ani_fstrain` enumeration or the HDF5 output fields, update it for the redefined `ani_fstrain == 3` and the new `aniso_delta` centroid field — otherwise note explicitly that it does not

## 8. Verification

- [x] 8.1 Build with `make build-dev`; resolve all warnings on the touched files
- [x] 8.2 Run the full `AnisotropyBenchmarkTests` and `RheologyCreepTests` suites — every `ani_fstrain ∈ {0,1,2}` test passes unchanged to prior tolerance
- [x] 8.3 Run the new `AniFstrainRelax*` tests — all pass
- [x] 8.4 Restart round-trip smoke check — `aniso_delta` / `aniso_delta_fs_prev` are bit-identical across a checkpoint+restart
- [x] 8.5 Manual ≥2-step run of an existing `ani_fstrain = 1` anisotropy scenario — confirm no regression vs a pre-change run
- [x] 8.6 `grep` confirms zero remaining `ani_T_threshold` references in `MDLIB/`, `TESTS/`, AND `SETS/`
