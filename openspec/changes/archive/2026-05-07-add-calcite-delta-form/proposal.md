## Why

The just-archived [add-hansen-olivine-delta-form](../archive/2026-05-07-add-hansen-olivine-delta-form/) change delivered a **mineral-agnostic framework** for finite-strain viscous-anisotropy dispatch in MDOODZ:

- per-phase function pointer `materials->aniso_delta_fn[k]` populated by `ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c)
- `AnisotropyRoutines.c` is fully mineral-agnostic — adding a new mineral = adding a new static helper + new case in `FlowLaws.c`, zero changes to dispatch code

So far only **case 1 = Hansen olivine** is shipped. MDOODZ scenarios using non-olivine phases (calcite — `SETS/PinchSwellGSE*.txt`, `SETS/PinchSwellAniso.txt`; quartz veins — Schmalholz & Duretz pinch-and-swell) currently fall back to `ani_fstrain = 1` (the `min(FS_AR, ani_fac_max)` form).

**Quartz** initially shipped with documentation-only ("min form is fit-for-purpose with `ani_fac_max ≈ 3`"), but a **post-archive review** showed the min-form transient overshoot is systematically wrong in the γ ∈ [1, 4] regime where MDV-class CPO is developing. A free LSQ fit to the Pennacchioni+10 17-point dataset using the existing saturating-exponential framework gives RMS = 0.29 vs the best min-form RMS = 0.52 — a **45% improvement**. Quartz was therefore added as **`case 3` of `ReadDataAnisotropy`** with constants `M_inf = 0.70, γ_e = 2.00, slope = 3.0` (asymptote `δ_∞ = 3.1`). See [SETS/AnisoFstrainResearch/quartz/quartz_calibration.md](SETS/AnisoFstrainResearch/quartz/quartz_calibration.md) for the history note.

**Calcite needs a custom δ(FS_AR) form**. With `ani_fac_max ≈ 4` (the existing PinchSwell calibration), the min form saturates at γ ≈ 1.6 while published torsion experiments show calcite CPO saturates at γ ≈ 5 (Pieri+01 J-index data: J ≈ 1.66 at γ=1, J ≈ 2.02 at γ=2, J ≈ 9.76 at γ=5; non-monotonic decline to J=8.92 at γ=11 due to recrystallization randomization). The transient overshoot of up to factor of 2 over γ ∈ [1.6, 5] matters for pinch-and-swell scenarios where the layer mostly experiences γ ∈ [0.5, 5]. This change adds a calibrated calcite form (`case 2`) using the saturating-exponential framework with calcite-specific constants.

## What Changes

- **Add `case 2 = Calcite (Pieri+01 / Bruijn+11)`** to `ReadDataAnisotropy()` in [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c). New static helper `anisoDelta_Calcite(double FS_AR)` encapsulating calcite-specific calibration constants (M-equivalent asymptote, γ_e, δ-vs-fabric slope).
- **Auto-default for ani_fstrain = 2 stays olivine** — flipping this on existing `SETS/RiftingComprehensive*.txt` and similar would be a silent semantic change. Calcite scenarios opt in by setting `aniso_db = 2` explicitly per phase.
- **Two new SETS subfolders** mirroring `SETS/AnisoFstrainResearch/`:
  - `SETS/AnisoFstrainResearch/calcite/` — calibration script `calibrate_calcite.py`, plot `mdoodz_vs_calcite_data.png` (M(γ) panel + δ(γ) panel with min-form / new case-2 form / Pieri+01 data), comparison markdown `calcite_calibration.md`.
  - `SETS/AnisoFstrainResearch/quartz/` — calibration script `calibrate_quartz.py`, plot `mdoodz_vs_quartz_data.png` (showing the existing min form with `ani_fac_max ∈ {2, 3}` against Pennacchioni+10 17-point dataset), comparison markdown `quartz_calibration.md` documenting the recommendation that **quartz scenarios use `ani_fstrain = 1` + a phase-tuned `ani_fac_max`** with explicit citations and rationale.
- **One new CI test** — `AnisotropyBenchmark.AniFstrainCalcite` — that runs a homogeneous simple-shear simulation with `ani_fstrain = 2`, `aniso_db = 2` and verifies `Centers/ani_fac` matches the analytical reference (the calibrated formula evaluated symbolically) to FP precision over γ ∈ [0, 10].
- **No quartz CI test** — quartz uses the existing `ani_fstrain = 1` form which is already covered by `AniFstrainSimpleShear` and `AniFstrainSaturation`. The quartz markdown documents the recommendation but doesn't gate it.
- **Backward compatibility (HARD)**: `case 1 = Hansen olivine` unchanged; `ani_fstrain ∈ {0, 1, 3}` unchanged; the auto-default (`ani_fstrain = 2 + aniso_db = 0` → olivine case 1) unchanged. Existing CI tests pass byte-identically. The new calcite calibration is **opt-in only** via explicit `aniso_db = 2`.

## Capabilities

### New Capabilities

(none — extending the existing capability)

### Modified Capabilities

- `ci-finite-strain-anisotropy`: extend the `aniso_db` enumeration documentation to include `case 2 = calcite`. Add a new requirement documenting the calcite calibration source, the analytical reference, and the new CI test case. Add a research-note requirement for `SETS/AnisoFstrainResearch/quartz/` documenting the min-form recommendation. Add a research-note requirement for `SETS/AnisoFstrainResearch/calcite/` analogous to the existing Hansen olivine research note.

## Impact

- **Code under change**: MDLIB additions confined to [MDLIB/FlowLaws.c](MDLIB/FlowLaws.c) — 1 new static helper (`anisoDelta_Calcite`) + 1 new switch case (~25 lines total). No changes to [AnisotropyRoutines.c](MDLIB/AnisotropyRoutines.c) (the function-pointer architecture means new minerals plug in without touching dispatch code).
- **Tests**: 1 new GTest case in [TESTS/AnisotropyBenchmarkTests.cpp](TESTS/AnisotropyBenchmarkTests.cpp). Reuses the existing `AniFstrainEvolution.txt` base fixture via `MutateInput` — no new fixture file. ~50 lines.
- **Documentation**:
  - [SETS/AnisoFstrainResearch/calcite/calcite_calibration.md](SETS/AnisoFstrainResearch/calcite/calcite_calibration.md) — calibration procedure, source breakdown (Pieri+01 + Bruijn+11), J-index → M-index conversion approach (Skemer+05), constants, comparison plot, RMS error, structural caveats.
  - [SETS/AnisoFstrainResearch/quartz/quartz_calibration.md](SETS/AnisoFstrainResearch/quartz/quartz_calibration.md) — Pennacchioni+10 17-point dataset analysis, demonstration that the existing min form fits-for-purpose with `ani_fac_max ∈ {2, 3}` and explicit T-range caveat (~500°C; basal-a / rhomb-a / prism-a slip-system regime changes outside this window).
- **Plots**: two new PNGs analogous to `hansen_olivine_calibration.png`. Calcite plot overlays MDOODZ HDF5 output (from new CI test) on Pieri+01 data. Quartz plot overlays MDOODZ output of the existing `AniFstrainSimpleShear` test (with chosen `ani_fac_max = 3`) on Pennacchioni+10 data.
- **CI runtime**: +1 test × short timesteps × small homogeneous box → ~+0.8 s. Within budget.
- **Risk**: **medium**. Calcite asymptote validation is weaker than olivine — we do not have a Hansen+16 Eq. 3-4-style independent stress-aware δ_∞ measurement. The δ-vs-fabric slope comes from published VPSC modeling (Tommasi+09 calcite Hill tensors) or J-to-M conversion (Skemer+05), not a clean linear lab fit like Hansen+12 Fig 3b. We document this as a calibration uncertainty in the calcite research note.
- **Out of scope**:
  - **Quartz case 3**: not justified — the min form with phase-tuned `ani_fac_max` is fit-for-purpose given quartz's modest δ_max. Documented as a research-note recommendation, not as MDLIB code.
  - **Quartz T-aware dispatch**: the saturation-strain question is moot since we're not adding a quartz case. If a user later finds the min-form transient is wrong for their scenario, a follow-up change can add a quartz case 3 with T-aware dispatch (`double aniso_delta_fn(double FS_AR, double T)`).
  - **High-pressure calcite** (Schuster+18): HPT data goes to γ = 80 at 1–4 GPa, but pressure-induced phase transition (calcite → CaCO₃-II at 1.6 GPa) makes it a different material. Documented for context, not used in fit.
  - **Mica, halite, anhydrite**: published torsion + M-index data is sparse or non-existent in directly-usable form. Defer.
  - **Stress-orientation-aware viscous anisotropy** (Hansen+16 Part 1 Eq. 2 / Part 2 Fig. 10): MDOODZ's δ remains a scalar magnitude. Future Hill-tensor work as a separate change.
