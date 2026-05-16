# Visual tests — `ani_fstrain` machinery

Six explanatory plots covering the MDOODZ `ani_fstrain` (per-marker
finite-strain anisotropy) machinery, including the Phase 1 init-from-
finite-strain helper (`aniso-init-from-finite-strain` change) and the
Boneh-DiSRX `ani_fstrain=3` δ-relaxation kinetics.

All output PNGs live at `VISUAL_TESTS/img/aniso_fstrain/0*.png`.
Rendering scripts and the Python drivers live in
`VISUAL_TESTS/aniso_fstrain/`.

## Architecture (May 2026 refactor)

The suite previously consisted of **29 sister SETs** in `SETS/`
(`AniFstrainOlivineCompare_F1`, `…_F2`, `…_F3`, six `FInit_a*`, nine
`FactorSweep_f*`, four `InheritSweep_f*`, six `TSweep_T*`, plus the
single-binary 42-cell `RegimeCell` driver). Each pair of `.c` files
differed only by the default `.txt` path; each `.txt` differed only by
a handful of parameter values.

Now **one** shared C source (`SETS/AnisoFstrainBox.c`) plus **one**
shared `.txt` template (`templates/base.txt`) plus **six** Python
runners drive the whole suite. The C source accepts `argv[1]` as the
path to the `.txt`, so each runner generates a per-cell `.txt` from the
template and invokes the binary directly.

```
SETS/
  AnisoFstrainBox.c            shared kinematic callback (single binary)
  AnisoFstrainBox.txt          minimal default for CLI-only use

VISUAL_TESTS/aniso_fstrain/
  INDEX.md                     this file
  templates/
    base.txt                   single .txt with placeholders for all sweep dims
  runners/
    _lib.py                    template-substitution + binary invocation
    run_01_fstrain_compare.py  3 cells (ani_fstrain ∈ {1,2,3})
    run_02_finit_angle.py      6 cells (aniso_angle ∈ {0,30,60,80,90,120}°)
    run_03_finit_factor.py     9 cells (aniso_factor ∈ {1..14})
    run_04_inheritance.py      4 cells (aniso_factor ∈ {1,2,4,8})
    run_05_T_sweep.py          6 cells (T ∈ {500..1500} K)
    run_06_regime_map.py       42 cells (T × ε̇ regime map)
  plot_01_fstrain_compare.py
  plot_02_finit_angle_panel.py
  plot_03_finit_factor_roundtrip.py
  plot_04_inheritance_reset.py
  plot_05_T_sweep_decay.py
  plot_06_regime_map.py
  runs/                        scratch outputs (gitignored)
    01_fstrain_compare/<cell>/{input.txt, output/Output*.gzip.h5, stdout.log}
    02_finit_angle/<cell>/{...}
    03_finit_factor/<cell>/{...}
    04_inheritance/<cell>/{...}
    05_T_sweep/<cell>/{...}
    06_regime_map/<cell>/{...}, results.json
```

`SETS/CMakeLists.txt` now contains exactly **one** `add_set(AnisoFstrainBox)`
entry where there used to be 29 `add_set(AniFstrainOlivine*)` lines.
`AniFstrainOlivineBoneh_InitFS` (the Phase-1 regression SET) is
deliberately NOT part of the consolidation — separate scope.

## Reproducing

```bash
# Build the unified binary (one-time, after a cmake-build reconfigure).
cd cmake-build && cmake -USET .. && make -j 8 AnisoFstrainBox

# Run all six tests. Total wall time on an M-series mac ≈ 80 s.
cd ../VISUAL_TESTS/aniso_fstrain
python3 runners/run_01_fstrain_compare.py
python3 runners/run_02_finit_angle.py
python3 runners/run_03_finit_factor.py
python3 runners/run_04_inheritance.py
python3 runners/run_05_T_sweep.py
python3 runners/run_06_regime_map.py

# Render the six PNGs.
python3 plot_01_fstrain_compare.py
python3 plot_02_finit_angle_panel.py
python3 plot_03_finit_factor_roundtrip.py
python3 plot_04_inheritance_reset.py
python3 plot_05_T_sweep_decay.py
python3 plot_06_regime_map.py
```

To add or modify a sweep dimension, edit the corresponding
`runners/run_0X_*.py`. To add a new placeholder (currently 18 cover
SCALES, grid, time, BCs, anisotropy, and phase rheology), add a
`__NEWNAME__` token to `templates/base.txt` and a key to
`_lib.DEFAULTS`.

## Plots

### 01 — `01_fstrain_compare.png`
**ani_fstrain 1 vs 2 vs 3 comparison.** Olivine simple shear γ̇=2 to γ≈16,
T=300 K. The three modes show clearly:

- `ani_fstrain=1` (default MDOODZ): δ = FS_AR capped at `ani_fac_max=100`.
- `ani_fstrain=2` (Hansen-strict): δ = `aniso_delta_fn(FS_AR)`,
  saturates at the Hansen ceiling 14.13.
- `ani_fstrain=3` (Hansen + Boneh DiSRX): overlays `ani_fstrain=2`
  exactly under cold T (strain-rate-limited relaxation is suppressed —
  operator split collapses to F2).

Right panel shows the shared FS_AR(γ) curve (kinematics independent of
δ-mode).

### 02 — `02_finit_angle_panel.png`
**F-init geometry across `aniso_angle`.** Six panels for `aniso_angle ∈
{0°, 30°, 60°, 80°, 90°, 120°}` with `aniso_factor=4`, `ani_fstrain=3`
inheritance turned on. Black bars mark the F-major principal axes,
gray the F-minor; red dashed line is the director (perpendicular to the
fabric); blue solid line is the foliation (at `aniso_angle − 90°`).

**Acceptance check:** black bars lie on the blue foliation line for all
six angles → F-init is being constructed with the design-correct
orientation. Cell-uniform `aniso_delta = 4.0` background confirms the
round-trip (`aniso_delta_inv_hansen(4) → γ_eff → F_init → aniso_delta_fn
→ 4`) is exact.

### 03 — `03_finit_factor_roundtrip.png`
**δ round-trip across `aniso_factor`.** Nine SETs with `aniso_factor ∈
{1, 2, 3, 4, 6, 8, 10, 12, 14}`.

- Top: extracted `Centers/aniso_delta` at step 0 vs input
  `aniso_factor` — points lie exactly on the y=x identity line.
- Bottom: Hansen analytic curve δ(γ_eff) = 13.132·(1 − exp(−γ/3.96)) + 1
  with the nine SET runs overlaid at their predicted (γ_eff, δ)
  positions. Twin x-axis shows the corresponding simple-shear FS_AR.

The f=14 case sits just below the Hansen ceiling (γ_eff=18.2,
FS_AR≈400) — the design's "clamp to 14.13 − 1e-9" safety net is
visible.

### 04 — `04_inheritance_reset.png`
**Inheritance reset-strain sweep.** Four SETs with `aniso_factor ∈
{1, 2, 4, 8}` initial fabric, γ̇=2 simple shear to γ=20, T=300 K,
`ani_relax_eps_max=−1` (gate OFF — isolate F-init effect).

- δ(γ) curves colour-coded by initial aniso_factor.
- Black dashed: Hansen analytic δ(γ) — the f=1 case overlays it exactly.
- Grey band at γ ∈ [1, 5]: the Warren+08 / Boneh+15 "reset strain"
  window where inherited fabric typically gets rewritten.
- Inset: log-scale residual |δ_initN − δ_init1| vs γ — exponential
  collapse → the inheritance bonus decays with reset strain.

Note the **non-monotonic early-γ dip** for f=4 and f=8: with
`aniso_angle=80°` the F-init principal stretch is nearly orthogonal to
the imposed shear, so the first few strain increments rotate the major
axis through a temporary minimum before climbing.

### 05 — `05_T_sweep_decay.png`
**T-sweep pure δ-decay (Arrhenius staircase).** Six SETs at T ∈ {500,
700, 900, 1100, 1300, 1500} K, `mechanical=1` with tiny
`bkg_strain_rate=1e-20` (F effectively frozen), `ani_relax_eps_max=−1`,
δ_init=4. Total time ~1.27 Myr (Nt=400, dt=1e11 s).

The relaxation kinetics

τ_relax = L_relax / V,  V = M₀·exp(−Q/RT)·μ·b²·Δρ

produce a clean Arrhenius cascade:
- 500 K: frozen (τ ≈ 1.4 Tyr ≫ Myr time-window)
- 700 K: barely starting
- 900 K: mid-decay (τ ≈ 1 Myr)
- 1100 K: nearly complete decay
- 1300 K / 1500 K: instantaneous

**Caveat:** with `mechanical=1 + tiny BSR`,
`strain_pwl ≈ 0` → DeltaRhoProxy returns `drho_min = 1e10` (not
`drho_max = 1e13`), shifting the staircase UP in T by ~3 decades vs the
hot-deforming-rock case. The half-decay point sits at T ≈ 1000 K (not
850 K). Still a clean kinetic visualisation.

### 06 — `06_regime_map.png`
**2D T × ε̇ regime map.** Forty-two cells covering T ∈ {500, 700, 900,
1100, 1300, 1500} K × ε̇ ∈ {1e-16, …, 1e-10} s⁻¹. Each cell is a tiny
11×11 simple-shear run integrated to γ=5 with `pwlv=40` olivine.

Final δ at γ=5 is plotted as a heatmap. Four regimes are visible:
- **FROZEN / SATURATED** (bottom + right): δ → Hansen ceiling 12.5 (the
  γ=5 value of the Hansen forward law).
- **ERASED** (top-left): hot + quiescent → relaxation outraces FS_AR
  growth, δ falls toward the inherited δ_init=4 (and below).
- **COMPETITION**: the rifting-relevant transition zone (where the
  reference point at T=1100 K, ε̇=3e-15 sits — marked with a red star).

Overlays:
- Dashed black: kinetic boundary τ_relax(T)·ε̇ = 1 (using drho_max).
- Red dotted: gate threshold ε̇_max = 1e-13 s⁻¹ (informational here —
  the gate is OFF in this map).
- Inset top-left: 1D cut at T=1100 K showing the δ-vs-ε̇ sigmoid.

## Test sweep grids (defined in `runners/run_0X_*.py`)

| Test | Sweep dim | Values | Cells |
|------|-----------|--------|-------|
| 1 | `ani_fstrain` | 1, 2, 3 | 3 |
| 2 | `aniso_angle` | 0, 30, 60, 80, 90, 120 (degrees) | 6 |
| 3 | `aniso_factor` | 1, 2, 3, 4, 6, 8, 10, 12, 14 | 9 |
| 4 | `aniso_factor` (init) | 1, 2, 4, 8 | 4 |
| 5 | `bkg_temperature` | 500, 700, 900, 1100, 1300, 1500 (K) | 6 |
| 6 | `(T, ε̇)` 2D | 6 × 7 grid | 42 |

All cells share the same C binary; runtime differences come purely
from the per-cell `.txt`.

## Wall times (M-series mac, ARM64)

| Test | Cells | Wall (parallel) |
|------|-------|-----------------|
| 1 | 3 | ~3 s |
| 2 | 6 | ~1.3 s |
| 3 | 9 | ~1.4 s |
| 4 | 4 | ~2 s |
| 5 | 6 | ~53 s (Nt=400) |
| 6 | 42 | ~14 s |
| **total** | **70** | **~75 s** |

## Deviations / notes from the original design

1. **Single shared C binary**: avoids 29 redundant compilations and
   keeps the test surface a `.txt` problem, not a `.c` problem. The
   unified C always passes `SetTemperature = SetT` returning
   `bkg_temperature`; harmless for Tests 1-4 since `cooling=0` and
   `thermal=0` everywhere in the suite.

2. **Test 1 `ani_relax_eps_max`**: the old `Compare_F*.txt` files did
   NOT set this parameter, relying on the MDOODZ default of `1e-13`
   s⁻¹ (gate ON). The new template makes that explicit:
   `run_01_fstrain_compare.py` writes `ANI_RELAX_EPS_MAX="1e-13"` so
   the per-cell `.txt` is self-documenting.

3. **Test 1 `dt_min`**: the old `Compare_F*.txt` files used
   `dt_min=1e-10`; Tests 2-6 used `dt_min=1e-30`. The new template
   uses `1e-30` uniformly (the actual dt=0.5 is far above either
   floor — outputs are bit-identical).

4. **Test 5 dt unit**: as before, `.txt: dt = 1e11` (SI seconds) →
   `dt_scaled = 0.01` via `scaling.t = L/V = 1e13`. Test 6 computes
   `dt_SI = γ_total / (2·ε̇·Nt)` and writes it directly in SI seconds.

5. **Test 5 `mechanical = 1` with tiny BSR**: with `mechanical=0` and
   `BSR=0` the initial viscosity computation produces NaN. The new
   runner uses `mechanical=1` + `BSR=1e-20` — preserving the
   relaxation-only effect (total strain ≈ 1e-7 over 1.27 Myr).

6. **PNG verification**: every regenerated PNG is **byte-identical**
   to the pre-refactor version (`max_channel_diff = 0` over all six
   images, by `Pillow` pixel-diff).
