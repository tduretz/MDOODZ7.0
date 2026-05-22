## ADDED Requirements

### Requirement: Calcite-calibrated δ-form is available as `aniso_db = 2`

`ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c) SHALL provide a calcite-calibrated entry, dispatched when the per-phase parameter `aniso_db == 2` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( a · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with

- `γ_eff(FS_AR) = √FS_AR − 1/√FS_AR` (analytical inversion of the simple-shear identity, mineral-independent),
- `M(γ) = M_inf · (1 − exp(−γ / γ_e))` (saturating exponential),
- `M_inf = 0.41` and `γ_e = 3.71` (free LSQ fit to the combined 18-point Pieri+01 + Bruijn+11 + Barnhoorn+04 calcite torsion dataset, J→M-converted via Skemer+05; Barnhoorn high-J samples filtered for J ≤ 15 where the linear conversion is valid). **Post-archive update**: original commit was `(0.50, 3.50, δ_∞=4.0)` matching empirical PinchSwell; Barnhoorn+04 expansion brought LSQ-best to `(0.41, 3.71, δ_∞=3.46)`.
- `a = 6.0` (calcite δ-vs-M slope; calibrated to yield asymptote `δ_∞ = 4` matching the empirical `ani_fac_max = 4` calibration used by `SETS/PinchSwell*` scenarios; bracketed by Tommasi+09-style VPSC bounds for calcite saturated CPO).

The δ-vs-M slope `a = 6.0` is mineral-specific (calcite-only), distinct from the olivine `a = 24.5` taken from Hansen+12 Fig 3b. The slope reflects calcite's lower viscous-anisotropy-per-fabric-strength ratio (Barnhoorn+04: CPO contributes ~1/3 of calcite weakening, the rest from grain-size reduction — already captured by MDOODZ's wattmeter).

The closed-form δ-asymptote (before the `ani_fac_max` cap) is `δ_∞ = a · M_inf + 1 = 4.0`, significantly lower than olivine's 14.13.

The form is **opt-in** via explicit `aniso_db = 2` per phase. The auto-default (`ani_fstrain == 2 && aniso_db == 0 → aniso_db = 1`) remains **olivine-specific** — flipping it on existing scenarios would be a silent semantic change. Calcite scenarios MUST set `aniso_db = 2` explicitly. Existing `ani_fstrain ∈ {0, 1, 3}` and `aniso_db = 1` (olivine) dispatch arms remain byte-identical.

#### Scenario: aniso_db = 2 returns the calibrated formula at representative strains

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer (set by `ReadDataAnisotropy(case 2)`) with `FS_AR` set to the analytical simple-shear FS_AR at γ = 5.0 (`FS_AR ≈ 27.299`), `ani_fstrain = 2`, and `ani_fac_max = 100` (uncapped)
- **THEN** the returned δ SHALL be within `1e-6` of `6.0 · 0.50 · (1 − exp(−5.0 / 3.50)) + 1 ≈ 3.26`
- **AND** when called with `FS_AR` at γ = 10.0 (`FS_AR ≈ 102.005`), the returned δ SHALL be within `1e-6` of `6.0 · 0.50 · (1 − exp(−10.0 / 3.50)) + 1 ≈ 3.83`

#### Scenario: aniso_db = 2 saturates at δ_∞ ≈ 4.0 for large strain

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR = 1e6` (γ_eff ≈ 1000), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-3` of `6.0 · M_inf + 1 = 4.0`

#### Scenario: aniso_db = 2 reduces to δ ≈ 1 at zero strain

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR = 1.0` (γ_eff = 0, isotropic), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-9` of `1.0` — `M(0) = 0`, so δ = 6·0 + 1 = 1

#### Scenario: ani_fac_max cap is honoured under aniso_db = 2

- **WHEN** `AnisoFactorEvolv` is called via the calcite function pointer with `FS_AR` chosen such that the uncapped calcite formula yields δ > `ani_fac_max` (e.g. `FS_AR = 1e6` with `ani_fac_max = 2`)
- **THEN** the returned δ SHALL equal `ani_fac_max` exactly (within `1e-12` absolute), matching the cap semantics of all other dispatch arms

#### Scenario: aniso_db = 1 (olivine) dispatch arm is unaffected by the calcite addition

- **WHEN** `AnisoFactorEvolv` is called via the olivine function pointer (set by `ReadDataAnisotropy(case 1)`) for any value of `FS_AR ∈ [1, 1e6]`, `ani_fstrain = 2`, `ani_fac_max ∈ [1, 1e6]`
- **THEN** the returned δ SHALL equal the value the previous (olivine-only) implementation would have returned, to bitwise precision — adding case 2 SHALL NOT alter case 1's IEEE-754 output

#### Scenario: Auto-default for ani_fstrain = 2 without aniso_db remains olivine

- **WHEN** an input file (or a `MutateInput` callback) sets `ani_fstrain[k] = 2` without setting `aniso_db[k]` (left at default 0)
- **THEN** the auto-default in [MDLIB/InputOutput.c](../../MDLIB/InputOutput.c) at parse time AND in [MDLIB/Main_DOODZ.c](../../MDLIB/Main_DOODZ.c) post-`MutateInput` SHALL set `aniso_db[k] = 1` (Hansen olivine), NOT `aniso_db[k] = 2`
- **AND** the calcite calibration SHALL only be selected when the user OR the `MutateInput` callback explicitly sets `aniso_db[k] = 2`

### Requirement: CI test `AniFstrainCalcite` exercises the calcite form

The CI test suite SHALL include a GoogleTest case `AnisotropyBenchmark.AniFstrainCalcite` that runs a homogeneous prescribed simple-shear simulation with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 2` and verifies that the per-cell `aniso_factor_n` (`Centers/ani_fac` in HDF5) matches the analytical reference (the calibrated calcite formula evaluated symbolically) to FP precision over the simulated strain range.

The test SHALL share the existing base fixture [TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt](../../TESTS/AnisotropyBenchmark/AniFstrainEvolution.txt) via `MutateInput`, consistent with the existing `AniFstrainSimpleShear`, `AniFstrainPureShear`, `AniFstrainSaturation`, and `AniFstrainHansenOlivine` tests in [TESTS/AnisotropyBenchmarkTests.cpp](../../TESTS/AnisotropyBenchmarkTests.cpp). No new fixture file SHALL be added.

The reference δ at output step `k` is evaluated from the analytical FS_AR at `γ = γ_dot · (k-1) · dt`, then mapped through `γ_eff → M → δ` using the calibrated constants `M_inf = 0.41`, `γ_e = 3.71`, `slope = 6.0` (post-archive update from initial `0.50, 3.50`) from the static helper `anisoDelta_Calcite` in [MDLIB/FlowLaws.c](../../MDLIB/FlowLaws.c).

The test runs to γ = 10 (Nt = 11, dt = 0.5, γ_dot = 2) — sufficient to reach 95% of the asymptote (δ ≈ 3.29 at γ = 10 vs δ_∞ = 3.46). Running further provides no additional verification.

#### Scenario: AniFstrainCalcite matches the calibrated formula per step

- **WHEN** the test runs simple-shear with `anisotropy = 1`, `ani_fstrain = 2`, `aniso_db = 2`, `ani_fac_max = 100` (uncapped), `bkg_strain_rate = 1.0` so `γ_dot = 2`, `Nt = 11`, `dt = 0.5` (γ ∈ {0, 1, 2, ..., 10} at output steps 1..11), single phase
- **THEN** for every output step `k`, the relative error between the **mean** of `Centers/ani_fac` over interior cells and the analytical reference `δ_ana(γ) = 6.0 · 0.50 · (1 − exp(−γ_eff(FS_AR(γ)) / 3.50)) + 1` evaluated at `γ = γ_dot · (k-1) · dt` SHALL be less than `1e-6`
- **AND** the relative L2 error between the spatial field and a constant analytical vector SHALL be less than `1e-4`

#### Scenario: Spatial homogeneity is preserved under calcite calibration

- **WHEN** the same simulation completes
- **THEN** at each output step, `(max(aniso_factor) / min(aniso_factor)) - 1` SHALL be less than `1e-4`

#### Scenario: AniFstrainCalcite reuses the shared base fixture and explicitly opts in to aniso_db = 2

- **WHEN** the test source under [TESTS/AnisotropyBenchmark/](../../TESTS/AnisotropyBenchmark/) is inspected
- **THEN** no new `.txt` fixture file SHALL be added for this test
- **AND** the test SHALL call `RunMDOODZ` on `AniFstrainEvolution.txt` and rely on a `MutateInput` callback that sets `ani_fstrain[0] = 2` AND `aniso_db[0] = 2` (the latter is required because the auto-default routes to olivine)

### Requirement: Research note `SETS/AnisoFstrainResearch/calcite/calcite_calibration.md`

The change SHALL deliver a research-grade markdown record at `SETS/AnisoFstrainResearch/calcite/calcite_calibration.md` documenting the calibration of the `aniso_db = 2` calcite form. The structure SHALL mirror §7 of [SETS/AnisoFstrainResearch/hansen2012_comparison.md](../../SETS/AnisoFstrainResearch/hansen2012_comparison.md) (the Hansen olivine research note) for consistency.

The note SHALL contain:

- (a) the functional form (Eq + constants `M_inf = 0.41`, `γ_e = 3.71`, `a = 6.0`; post-archive update from initial `0.50, 3.50`),
- (b) the calibration procedure (least-squares fit on Pieri+01 + Bruijn+11 J-index data, J→M conversion via Skemer+05),
- (c) the **asymptote uncertainty**: `δ_∞ = a · M_inf + 1 = 4.0` is calibrated to match the empirical `ani_fac_max = 4` PinchSwell calibration; bracketed by Tommasi+09-style VPSC bounds (calcite δ_saturated ≈ 2-5). NOT a clean Hansen+12 Fig 3b-equivalent direct lab fit. This is the weakest link in the calibration.
- (d) the comparison plot at `SETS/AnisoFstrainResearch/calcite/mdoodz_vs_calcite_data.png`, showing the new calibrated curve overlaid with Pieri+01 + Bruijn+11 J-index lab data (J→M-converted) and the existing min(FS_AR, 4) form for comparison,
- (e) the per-datapoint verification table (γ, J_obs, M_obs, M_fit, δ_fit) listing both Pieri+01 (γ ∈ {0, 1, 2, 5, 11}) and Bruijn+11 (γ ∈ {1, 2.6, 5}) datapoints,
- (f) a discussion of the **Pieri+01 vs Bruijn+11 disagreement at γ=5** (P087 J=9.76 vs Bruijn J≈4) and how the fit smooths through it,
- (g) a discussion of the **Pieri+01 γ=11 non-monotonic outlier** (J=8.92 < J=9.76 at γ=5; recrystallization randomization) and the structural decision to smooth through it,
- (h) the simple-shear `γ_eff(FS_AR)` caveat, identical to the olivine note,
- (i) explicit non-applicability to **HPT calcite (Schuster+18)** at confining pressures > 1.6 GPa where the calcite → CaCO₃-II phase transition kicks in,
- (j) full citations to Pieri+01, Bruijn+11, Schuster+18, Skemer+05, Barnhoorn+04, Tommasi+09 (or whichever VPSC-bound references underpin the slope choice).

The accompanying calibration script SHALL be checked in at `SETS/AnisoFstrainResearch/calcite/calibrate_calcite.py` (Python, no scipy required — uses 2D grid search to fit `M_inf` and `γ_e` with `slope` fixed at 6.0).

#### Scenario: Calcite research note exists with required content

- **WHEN** `SETS/AnisoFstrainResearch/calcite/calcite_calibration.md` is read
- **THEN** it SHALL contain sections covering each of (a)-(j) above
- **AND** it SHALL reference the calibration script `calibrate_calcite.py` and the comparison plot `mdoodz_vs_calcite_data.png` via relative paths
- **AND** it SHALL NOT be referenced from any `TESTS/` `.cpp` test file or CI script (research-grade, not regression gate)

#### Scenario: Calcite calibration script and plot are checked in

- **WHEN** `SETS/AnisoFstrainResearch/calcite/` is inspected
- **THEN** `calibrate_calcite.py` SHALL exist as an executable Python script that, when run, prints the fit constants `M_inf, γ_e`, the asymptote `δ_∞`, the verification table at all Pieri+01 + Bruijn+11 datapoints, and writes `mdoodz_vs_calcite_data.png` to the same directory
- **AND** the script SHALL NOT depend on scipy (numpy + matplotlib only)

### Requirement: Quartz-calibrated δ-form is available as `aniso_db = 3`

> **Post-archive update**: this requirement was added when a free LSQ fit
> showed the Pennacchioni+10 dataset is well-fit by the existing
> saturating-exponential framework (RMS = 0.29 vs best min-form RMS = 0.52,
> 45% improvement). The earlier "non-decision" requirement was superseded.

`ReadDataAnisotropy()` in `MDLIB/FlowLaws.c` SHALL provide a quartz-calibrated entry, dispatched when the per-phase parameter `aniso_db == 3` (and `ani_fstrain == 2`):

```
δ(FS_AR) = min( 3.0 · M(γ_eff(FS_AR)) + 1 ,  ani_fac_max )
```

with `M_inf = 0.70`, `γ_e = 2.00`, `slope = 3.0`, asymptote `δ_∞ = 3.1`.

Calibrated by free LSQ fit to the 17-point Pennacchioni+10 natural-shear-zone dataset (T ~ 500°C, prism-a slip). The slope = 3.0 is fixed by VPSC-style Heilbronner-Tullis bounds. **For T < 500°C scenarios, fall back to `ani_fstrain = 1` with `ani_fac_max = 3`** since basal-a / rhomb-a slip changes the fabric kinematics and the prism-a calibration may not apply.

#### Scenario: aniso_db = 3 returns the calibrated formula at γ = 4 and γ = 8

- **WHEN** `AnisoFactorEvolv` is called via the quartz function pointer with `FS_AR ≈ 17.944` (γ=4), `ani_fstrain = 2`, `ani_fac_max = 100`
- **THEN** the returned δ SHALL be within `1e-6` of `2.82`
- **AND** at `FS_AR ≈ 66.0` (γ=8), the returned δ SHALL be within `1e-6` of `3.06`

#### Scenario: aniso_db = 3 saturates at δ_∞ ≈ 3.1

- **WHEN** `AnisoFactorEvolv` is called via the quartz function pointer with `FS_AR = 1e6`
- **THEN** the returned δ SHALL be within `1e-3` of `3.1`

### Requirement: CI test `AniFstrainQuartz` exercises the quartz form

The CI test suite SHALL include `AnisotropyBenchmark.AniFstrainQuartz` running simple-shear at γ ∈ {0, 1, ..., 8} with `ani_fstrain = 2`, `aniso_db = 3`, asserting per-step δ matches the analytical reference to `1e-6` relative.

#### Scenario: AniFstrainQuartz matches the calibrated formula per step

- **WHEN** the test runs simple-shear at `bkg_strain_rate = 1.0`, `Nt = 9`, `dt = 0.5`, single phase with `ani_fstrain = 2`, `aniso_db = 3`, `ani_fac_max = 100`
- **THEN** for every output step `k`, `EXPECT_NEAR(s.mean / δ_ana, 1.0, 1e-6)` SHALL hold

### Requirement: Research note `SETS/AnisoFstrainResearch/quartz/quartz_calibration.md`

> **Post-archive update**: this requirement was rewritten when quartz `case 3`
> was added. The earlier version recommended the min-form fallback; the
> updated version documents the calibrated case 3 form.

The change SHALL deliver a research-grade markdown record at `SETS/AnisoFstrainResearch/quartz/quartz_calibration.md` documenting the calibration of the `aniso_db = 3` quartz form, mirroring the structure of the calcite and Hansen olivine notes.

The note SHALL contain:

- (a) the functional form (Eq + constants `M_inf = 0.70`, `γ_e = 2.00`, `slope = 3.0`),
- (b) the calibration procedure (free LSQ fit on Pennacchioni+10 17 datapoints, J→M Skemer+05),
- (c) the asymptote uncertainty (VPSC-style bounds, NOT direct lab fit),
- (d) the comparison plot at `mdoodz_vs_quartz_data.png`,
- (e) the per-datapoint verification table (regime-coded WDV/MDV/SDV),
- (f) the **T-range caveat** (~500°C, prism-a slip),
- (g) a **history note** that an earlier version recommended the min-form fallback, superseded by the LSQ fit showing 45% RMS improvement,
- (h) full citations.

The accompanying scripts at `calibrate_quartz.py` and `mdoodz_vs_quartz_comparison.py` SHALL be checked in (pure stdlib + numpy + matplotlib, no scipy).

#### Scenario: Quartz research note exists with required content

- **WHEN** `SETS/AnisoFstrainResearch/quartz/quartz_calibration.md` is read
- **THEN** it SHALL contain sections covering each of (a)-(h) above
- **AND** it SHALL document the calibrated case 3 form (M_inf=0.70, γ_e=2.00, slope=3.0, δ_∞=3.1)
- **AND** it SHALL NOT be referenced from any `TESTS/` `.cpp` test file or CI script
