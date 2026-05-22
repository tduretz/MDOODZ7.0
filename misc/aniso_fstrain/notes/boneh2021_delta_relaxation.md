# `ani_fstrain == 3` — temperature-dependent δ-relaxation (Boneh et al. 2021 DiSRX kinetics)

This note documents the kinetic relaxation law implemented for
`ani_fstrain == 3` by the `ani-fstrain3-delta-relaxation` change: the
Boneh et al. 2021 discontinuous-static-recrystallization (DiSRX)
kinetics it draws on, the derivation of `τ_relax = L_relax / V` and its
modelling assumptions, the constant-`L_relax` regime analysis, the `Δρ`
proxy, and the tiered validation plan.

It is the research-grade companion to the inline code comments in
`MDLIB/AnisotropyRoutines.c` (`DeltaRelaxationTau`, `DeltaRhoProxy`,
the Boneh constants) and `MDLIB/RheologyParticles.c`
(`FiniteStrainAspectRatio`, the per-marker δ update).

---

## 1. What `ani_fstrain == 3` was, and what it is now

**Before.** `ani_fstrain == 3` meant "the `ani_fstrain == 2` δ-dispatch
**plus** an upstream `T < ani_T_threshold` deformation-gradient freeze":
when a marker's temperature dropped below a per-phase threshold
(`ani_T_threshold`, defaulting to 1100 °C for the damped-olivine case
7), its deformation-gradient tensor `F` stopped updating, so `FS_AR`
and δ froze at whatever value they had when the marker last crossed the
threshold.

A multi-round literature + code investigation established three
problems with the freeze:

1. **The 1100 °C value was misattributed.** Its stated "three
   convergent paths" justification (Karato dry-olivine ~1300 °C / HK03
   wet-olivine ~1100 °C / Beghein et al. 2014 ~1100–1200 °C) fails on
   all three under primary-source review: Karato and Hirth & Kohlstedt
   (2003) set *no* temperature threshold — the dislocation/diffusion
   transition they describe is stress- and grain-size-controlled — and
   Beghein et al. 2014 reports a *depth* (the Gutenberg discontinuity),
   explicitly not a temperature. The one genuine direct citation
   (Debayle & Ricard 2013) uses 1100 °C only as a modelling
   *convention*.

2. **The freeze was a misinterpretation of the original directive.**
   The intent was a *kinetic relaxation law* for the anisotropy factor
   δ — a deferred "Phase 2 (kinetic-law dispatch)" was already noted in
   the code. The earlier autonomous work built a step-function freeze
   ("Phase 1") with a retro-fitted citation chain instead.

3. **A fixed-`T` freeze is the wrong *kind* of model.** It is binary
   (on/off), has no rate, and discards the physics: when deformation
   slows or stops, olivine CPO does not "freeze" — it *relaxes* toward
   isotropic over a finite, temperature-dependent timescale via
   discontinuous static recrystallization (DiSRX).

**After.** `ani_fstrain == 3` is the `ani_fstrain == 2` δ-dispatch
**plus** the temperature-dependent kinetic *relaxation* of δ that was
actually intended. `ani_T_threshold` and the `F`-freeze block are
removed.

---

## 2. The Boneh et al. 2021 DiSRX kinetics

Boneh, Wallis, Hansen, Krawczynski & Skemer (2021), "The effect of
discontinuous static recrystallization on olivine
crystallographic preferred orientation", *Geochemistry, Geophysics,
Geosystems* (G3), establish the **kinetics** of discontinuous static
recrystallization in olivine — the grain-boundary-migration (GBM)
process by which strain-free grains nucleate and grow at the expense of
deformed, dislocation-rich grains once active deformation slows.

Their grain-boundary-migration velocity (Boneh 2021, Eqs. 1–4) is

```
V = M_b · ΣF
M_b = M0 · exp(−Q / (R·T))
ΣF  = F_s + F_b ,   F_s = μ·b²·Δρ ,   F_b = 3γ_surf / d
```

with the constants

| symbol | value | meaning |
|---|---|---|
| `Q`  | 133 ± 18 kJ/mol | activation energy for GBM |
| `μ`  | 50 GPa          | shear modulus |
| `b`  | 0.6 nm          | Burgers vector |
| `M0` | 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹ | mobility prefactor |
| `Δρ` | 10¹¹–10¹³ m⁻²   | dislocation-density contrast across a migrating boundary |

`F_s` is the **strain-energy driving force** (stored dislocation
energy); `F_b` is the **surface-energy driving force** (boundary
curvature).

### 2.1 The `M0` units — resolved (Boneh-2021 typo)

Boneh et al. 2021 prints the `M0` units as **"m³/s·J"**, which is
dimensionally inconsistent with its own Eq. 1: `V = M_b · ΣF` with
`ΣF` in `J·m⁻³` (a force per unit volume) would then give

```
[m³ s⁻¹ J⁻¹] · [J m⁻³] = [s⁻¹]
```

— a *rate*, not a *velocity*. The Boneh, Skemer & Wallis (2017)
companion paper ("Dynamic stability of microstructure during static
recrystallization …") prints the dimensionally-correct
**m³/(N·s) ≡ m⁴·J⁻¹·s⁻¹** in its Appendix A. With that unit,

```
[m⁴ J⁻¹ s⁻¹] · [J m⁻³] = [m s⁻¹]
```

— a velocity, as Eq. 1 requires. And the corrected value reproduces
Boneh's published "weeks-to-years at ~1300 °C" anchor: a 10 µm grain
growing to ~1 mm takes ~0.23–23 yr for `Δρ = 10¹³–10¹¹ m⁻²`. We
therefore use `M0 = 2×10⁻¹¹ m⁴·J⁻¹·s⁻¹`. This is committed, not open;
the typo note is also an inline comment at the constant definition in
`MDLIB/AnisotropyRoutines.c`.

### 2.2 Dropping the surface-energy term — `ΣF ≈ F_s`

The implementation uses **only** the strain-energy driving force:
`ΣF ≈ F_s = μ·b²·Δρ`. The surface-energy term `F_b = 3γ_surf/d` is
dropped. This is justified because `F_s ≫ F_b` for mantle grain sizes
≥ 1 mm: Boneh 2021 §4.2 / Fig. 8 shows that for `d ≳ 1 mm` and the
dislocation densities of interest, the strain-energy term dominates the
total driving force by one to several orders of magnitude. The
approximation is documented inline at the `V` computation.

---

## 3. The relaxation law `τ_relax = L_relax / V`

### 3.1 The `τ = L / V` construction (a modelling leap)

Boneh et al. 2021 give the **kinetics** of DiSRX (the GBM velocity `V`,
and the time for an annealed grain to grow a given distance) and show
*qualitatively* that DiSRX weakens CPO. They do **not** give a
closed-form δ(t) relaxation law.

The construction used here — "δ relaxes toward the isotropic limit
(δ → 1) on the DiSRX grain-growth timescale `τ_relax = L_relax / V`" —
is **this model's own**, not Boneh's. The reasoning: when grain
boundaries sweep through a marker faster (larger `V`), the inherited,
deformation-built CPO is overwritten by strain-free, randomly-oriented
new grains faster. The characteristic distance a boundary must migrate
to substantially reset a marker's fabric is `L_relax` (≈ the annealed
grain size); the time to do so is `L_relax / V`. δ is then taken to
decay exponentially toward 1 on that timescale.

This is a defensible first-order coupling, but it is a modelling
choice: the *quantitative* δ(t) curve is not lab-validated. This caveat
is documented inline at the relaxation routine header in
`MDLIB/AnisotropyRoutines.c`.

### 3.2 The two-field operator split + analytic exponential update

δ for `ani_fstrain == 3` is a **per-marker integrated state variable**
(it carries history) — not a memoryless function of the instantaneous
`FS_AR`. Each marker carries two new scalar fields:

- `aniso_delta` — the relaxing δ state;
- `aniso_delta_fs_prev` — the previous step's finite-strain δ.

Each timestep, inside `FiniteStrainAspectRatio`:

```
delta_FS   = aniso_delta_fn(FS_AR)                          // the ani_fstrain==2 value (production)
delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev)  // inherit the ==2 increment
tau_relax  = L_relax / V                                     // V = M0·exp(-Q/RT)·μ·b²·Δρ  (Fs only)
delta_new  = 1 + (delta_prod - 1)·exp(-dt / tau_relax)       // analytic relaxation toward 1
delta_new  = fmin(delta_new, ani_fac_max)                    // clamp to the per-phase ceiling
delta_new  = fmax(delta_new, 1.0)                            // isotropic floor
aniso_delta = delta_new ;  aniso_delta_fs_prev = delta_FS     // store for next step
```

Why this scheme:

- **Production is *exactly* the calibrated `ani_fstrain == 2` form**,
  with no new free parameter. δ inherits the *increment* of
  `aniso_delta_fn(FS_AR)`, so while deforming, δ grows exactly as
  `ani_fstrain == 2` would.
- **Cold limit recovers `ani_fstrain == 2` exactly.** As
  `τ_relax → ∞`, `exp(−dt/τ) → 1`, so `delta_new = delta_prod`; the
  increments telescope to `aniso_delta = aniso_delta_fn(FS_AR)`.
- **Quiescent-hot limit is pure relaxation.** When
  `delta_FS == aniso_delta_fs_prev`, `delta_prod = aniso_delta`, so
  `delta_new = 1 + (aniso_delta − 1)·exp(−dt/τ)`.
- **Analytic exponential** — unconditionally stable for any `dt`,
  cannot overshoot below `δ_eq = 1`, and has no dependence on the
  (anisotropy-unaware) Courant timestep limiter. The FS_AR path has a
  documented unbounded-runaway history; an explicit forward-Euler
  update would reintroduce that risk, the analytic form does not.

An alternative single-field two-sided relaxation
`dδ/dt = (δ_FS − δ)·k_prod − (δ − 1)·k_relax` was considered and
rejected: it needs a *production rate* `k_prod` (a new per-phase
modelling parameter), trading one scalar field for a new free parameter
and a less-obvious cold-limit equality.

**One-step lag** (a note, not a defect): at the line-676 anisotropy
pass in `Main_DOODZ.c`, `F`/`FS_AR` and `strain_pwl` are
end-of-previous-step values, while `particles->T` is current. The
relaxation is a slow process and this lag is consistent with how
`ani_fstrain == 2` already consumes `FS_AR`, so it is acceptable; it is
stated as an inline comment.

### 3.3 The δ → grid → solver path

The relaxed per-marker `aniso_delta` reaches the Stokes solver through
two **new grid fields** `mesh->aniso_delta_n` (centroid) and
`mesh->aniso_delta_s` (vertex), populated by P2G inside
`FiniteStrainAspectRatio` exactly as `FS_AR_n`/`FS_AR_s` are.
`AnisoFactorEvolv` then takes one new argument — the P2G'd relaxed δ at
that grid point — and its `ani_fstrain == 3` arm returns that value
(clamped to `ani_fac_max`) instead of recomputing
`aniso_delta_fn(FS_AR)`. This routes the relaxed δ *through* the
existing per-phase `phase_perc` weighting and harmonic/geometric
averaging in `UpdateAnisoFactor`, so mixed-phase cells, the three
averaging modes, and the centroid+vertex grids all stay coherent. For
`ani_fstrain ∈ {0,1,2}` the new argument is unused and those arms are
byte-identical to their pre-change behaviour.

---

## 4. `L_relax` — the per-phase length scale

`L_relax` is the characteristic DiSRX boundary-migration distance — the
distance a grain boundary must sweep to substantially reset a marker's
fabric, ≈ the annealed grain size. It is the per-phase parameter
`ani_relax_length` (added to `mat_prop`):

- `ani_relax_length >= 0` → an explicit per-phase length, in metres
  (parsed and divided by `scaling.L`);
- `ani_relax_length < 0` (sentinel −1, the parser default and the
  case-7 `ReadDataAnisotropy` default) → "use `gs_ref[phase]`",
  resolved at the use site inside `FiniteStrainAspectRatio`.

**Grain-size evolution (the `gs` flow law) is NOT required.** With GSE
off, `L_relax = gs_ref` is a per-phase constant. With GSE on, `L_relax`
*could* track the evolving `mesh->d_n` — that GSE-on path is deferred
(the −1 sentinel keeps it non-breaking to add later). Note that
MDOODZ's GSE is a steady-state paleowattmeter — it implements
*dynamic*-recrystallization kinetics, the *wrong* process for DiSRX
annealing — so this change deliberately does **not** route through it.

### 4.1 Constant-`L_relax` regime analysis

With GSE off, `L_relax` is held fixed for the run. The annealed grain
size that DiSRX actually grows to is itself temperature- and
time-dependent, so a constant `L_relax` introduces an error. The
regime analysis distinguishes three bands by where a marker sits
relative to its relaxation timescale over the run:

- **Frozen band** (cold, `τ_relax ≫ run time`): δ does not move.
  `L_relax` is irrelevant — any constant value gives the same answer
  (no relaxation). **Error ≈ 0.**
- **Instantaneous band** (hot, `τ_relax ≪ run time` and `≪` the
  loading-pause duration): δ relaxes essentially to its fixed point
  (δ → 1, or → the production value while deforming) well within the
  run. The *path* depends on `L_relax` but the *endpoint* does not.
  **Error ≈ 0 on δ_final.**
- **Caught-partway band** (mid-temperature, or episodic loading where
  the quiescent intervals are comparable to `τ_relax`): δ is mid-decay
  when the run ends or deformation resumes, so the *value* depends
  directly on `τ_relax` and hence on the constant-`L_relax`
  assumption. **This is where the error lands**; the regime analysis
  bounds it at the ≤17 %-class on δ_final, and it is never qualitative
  (the sign of the trend — relax toward isotropic — is always right).

The constant-`L_relax` caveat, with the mid-T / episodic-loading band
called out, is documented inline at the `L_relax` resolution site in
`MDLIB/RheologyParticles.c`.

---

## 5. The `Δρ` proxy

`Δρ` is the dislocation-density driving force in `V`. **MDOODZ tracks
no dislocation density.** Rather than introduce a new per-marker state
variable with its own `dρ/dt` ODE (explicit ρ-evolution is future
work), `Δρ` is proxied on the fly from the per-marker accumulated
dislocation-creep strain `strain_pwl`, a field MDOODZ already carries.

The mapping (`DeltaRhoProxy` in `MDLIB/AnisotropyRoutines.c`) is a
saturating exponential tied to the Boneh dislocation-density range:

```
Δρ = Δρ_min + (Δρ_max − Δρ_min)·(1 − exp(−strain_pwl / ε_ref))
   clamped to [Δρ_min, Δρ_max]
```

with `Δρ_min = 10¹¹ m⁻²`, `Δρ_max = 10¹³ m⁻²` (the Boneh range), and
`ε_ref = 4.0` — a reference saturation strain on the scale of the
Hansen olivine δ-form saturation strain (`γ_e ≈ 3.96`). The result is
always finite and in `[Δρ_min, Δρ_max]`; a marker with no stored strain
gets `Δρ_min` (so `τ_relax` stays finite — never division by zero), a
heavily-strained marker saturates at `Δρ_max`. The `Δρ_min`/`Δρ_max`/
`ε_ref` choice is a calibration decision (design O2): it is exposed as
commented module-level constants in the source so it can be tuned.

### 5.1 Two caveats of the `strain_pwl` proxy

Both are documented inline at the `Δρ` computation in
`MDLIB/AnisotropyRoutines.c`:

**(a) Monotonicity / no-recovery.** `strain_pwl` is *accumulated*
dislocation-creep strain — it can only increase. A real dislocation
density *recovers* (drops) once deformation slows. The proxy cannot
represent that, so its `Δρ` is a "sticky" upper-ish estimate. For the
post-deformation relaxation use case this stickiness is acceptable, and
arguably even correct — the stored strain energy that drives DiSRX is
itself slow to recover — but a marker that deformed long ago and then
sat quiescent will still report a high `Δρ` and therefore relax
*faster* than a true ρ-evolution model would predict.

**(b) Absolute stored strain vs. boundary contrast.** Boneh's `Δρ` is
the dislocation-density *contrast across a migrating grain boundary* —
the thermodynamic driving force for that specific boundary.
`strain_pwl` is a marker's *absolute* accumulated stored strain energy.
The two coincide only when the marker is significantly more strained
than its neighbourhood (a strong gradient drives boundary migration
into it); for a uniformly-strained region the true driving contrast is
smaller than this proxy implies.

A `τ_II`-piezometer proxy (`Δρ` from the deviatoric stress invariant
via a piezometer) was considered and rejected as the *primary* choice:
it collapses to zero exactly when deformation stops, which is wrong for
the post-deformation use case (relaxation should continue after
loading ceases). A piezometer-cap refinement — using a piezometric
estimate as an upper bound on the `strain_pwl` proxy — is left open.

---

## 6. Tiered validation plan

### Tier 0 — analytical unit tests (the CI gate)

Implemented in `TESTS/AnisotropyBenchmarkTests.cpp`, suite
`AnisotropyBenchmark`, prefix `AniFstrainRelax*`. They exercise
`DeltaRelaxationTau` directly together with the analytic exponential
update — a deterministic per-marker calculation on a 1-marker /
homogeneous "box", with an identity scale and the `gs` flow law OFF.
Each gates a closed-form limit:

| Test | Closed-form limit gated |
|---|---|
| `AniFstrainRelaxColdFreeze` | `τ_relax → ∞` ⇒ δ constant to FP precision |
| `AniFstrainRelaxHot` | `δ(t) = 1 + (δ₀−1)·exp(−t/τ_relax)`; Arrhenius ratio `τ(T₁)/τ(T₂) = exp(Q/R·(1/T₂−1/T₁))` |
| `AniFstrainRelaxHalfLife` | half-life at `t = τ_relax·ln2`; Δt-convergence; large-Δt no-overshoot |
| `AniFstrainRelaxColdLimitEqualsFstrain2` | cold-T `ani_fstrain == 3` ≡ `ani_fstrain == 2` |
| `AniFstrainRelaxLengthMonotonic` | `τ_relax ∝ L_relax` |
| `AniFstrainRelaxClampAndFixedPoint` | `ani_fac_max` clamp; δ = 1 fixed point |
| `AniFstrainRelaxHistoryDependence` | same `FS_AR`, different thermal history ⇒ different stored δ |

These are *implementation gates*, not science validation — they verify
the code computes the law it claims to, not that the law matches
nature.

### Tier 1 — lab-data reproduction (future, research-grade)

Reproduce the published DiSRX annealing observations:

- **Boneh, Skemer & Wallis (2017)** — static-recrystallization
  microstructure stability experiments; check that the implemented `V`
  reproduces their measured grain-growth timescales at the experimental
  `T`, `d`, and stress.
- **Speciale, Behr, Hirth & Tokle (2020)** (and companion
  microstructural studies) — natural and experimental olivine
  static-recrystallization CPO weakening; check that the δ(t) decay
  rate is consistent with the observed fabric-strength loss during
  annealing.

These require building dedicated comparison scenarios under
`misc/aniso_fstrain/mdoodz_compare/` analogous to the existing
`mdoodz_vs_olivine.py` lab-overlay workflow, and are out of scope for
the CI gate.

### Tier 2 — geodynamic sanity (future)

Confirm that in a full lithospheric-cooling or post-orogenic-relaxation
scenario, `ani_fstrain == 3` produces the expected qualitative
behaviour: anisotropy *built* during active deformation, then *slowly
lost* in regions that cool or go quiescent — with the loss rate
correctly faster in hotter material.

---

## 7. Out of scope (future work)

- Modelling dislocation density `ρ` as a true state variable with its
  own `dρ/dt` ODE — a proxy is used.
- A DiSRX-specific grain-size-evolution model for `L_relax` (the
  annealed grain size DiSRX grows to) — `L_relax` is a per-phase
  constant from `gs_ref`.
- Fabric-type (A/B/C/E-type) tracking — the relaxation acts on the
  scalar δ *magnitude* only; it relaxes toward isotropic, not toward a
  different fabric type.
- Director (orientation) evolution — `ani_fstrain == 3` relaxes the
  anisotropy *magnitude* δ; the director angle (`nx`/`nz`) evolution is
  unchanged.
- Modifying or extending MDOODZ's GSE / paleowattmeter.
- The GSE-on `L_relax = mesh->d_n` path (design O3) — deferred; the −1
  sentinel makes it non-breaking to add later.

---

## 8. References

- Boneh, Y., Wallis, D., Hansen, L. N., Krawczynski, M. J., & Skemer,
  P. (2021). The effect of discontinuous static recrystallization on
  the crystallographic preferred orientation of olivine.
  *Geochemistry, Geophysics, Geosystems* (G3).
- Boneh, Y., Skemer, P., & Wallis, D. (2017). Dynamic stability of
  microstructure during static recrystallization. *(companion paper;
  Appendix A gives the dimensionally-correct `M0` units.)*
- Speciale, P. A., Behr, W. M., Hirth, G., & Tokle, L. (2020).
  Microstructural evolution during olivine static recrystallization.
- Hirth, G., & Kohlstedt, D. (2003). Rheology of the upper mantle and
  the mantle wedge: A view from the experimentalists. *(HK03 — cited
  here only to record that it sets no temperature threshold.)*
- Beghein, C., Yuan, K., Schmerr, N., & Xing, Z. (2014). Changes in
  seismic anisotropy shed light on the nature of the Gutenberg
  discontinuity. *Science*. *(reports a depth, not a temperature.)*
- Debayle, E., & Ricard, Y. (2013). Seismic observations of large-scale
  deformation at the bottom of fast-moving plates. *EPSL*. *(uses
  1100 °C as a modelling convention.)*
