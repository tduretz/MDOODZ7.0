## Context

This design builds a new lithospheric-scale rifting scenario (`RiftingCombinedYield`) that combines the MLPS_Ellipses multi-layer ellipse geometry with the Popov et al. 2025 combined mode-I/mode-II yield function (`plast=2`) and partial melting. The reference implementation draws from:

- **MLPS_Ellipses** — geometry, phase layout, ellipse helpers, boundary conditions
- **RiftingMelting** — melting configuration, melt weakening, density model
- **PressurizedMagmaChamber** — `plast=2` usage, `T_st` per phase, solver tuning for combined yield

All required MDLIB features already exist (commit `8aecb92`). This is a scenario-only change — no library modifications.

## Goals / Non-Goals

**Goals:**
- Create a working rifting scenario that exercises combined mode-I/mode-II yield with partial melting in a heterogeneous lithosphere
- Provide a low-resolution "verify" configuration for fast smoke-testing (~minutes, not hours)
- Use realistic multi-layer structure with elliptical heterogeneities adapted from MLPS_Ellipses
- Assign per-phase melting models (crustal: `melt=11`, mantle: `melt=41`) and per-phase tensile strength (`T_st`) and visco-plastic viscosity (`eta_vp`)

**Non-Goals:**
- No anisotropy — intentionally omitted to isolate yield + melting interaction
- No dike injection — that's a separate capability (see TCMagmaticSystem)
- No new flow laws or library-level code changes
- Not a benchmark — this is a geodynamic scenario, not a convergence test

## Decisions

### 1. Phase structure: Adapt MLPS_Ellipses with melting assignments

**Choice**: Keep the MLPS_Ellipses 9-phase structure (IDs 0–8) and add melting + combined yield parameters per phase.

**Rationale**: The MLPS_Ellipses layout models a realistic continental lithosphere (crust with weak/strong ellipses, lithospheric mantle with weak lenses, asthenosphere). Reusing it avoids reinventing geometry while adding new physics.

**Phase-to-melting mapping:**

| ID | Role (from MLPS_Ellipses) | `melt` | `plast` | `T_st` | `eta_vp` |
|----|---------------------------|--------|---------|--------|----------|
| 0 | Continental crust (anorthite) | 11 | 2 | -10e6 | 2.5e20 |
| 1 | Lithospheric mantle (dry olivine) | 41 | 2 | -10e6 | 2.5e20 |
| 2 | Lower mantle (wet olivine) | 41 | 2 | -10e6 | 2.5e20 |
| 3 | Asthenosphere (same as 2) | 41 | 2 | -10e6 | 2.5e20 |
| 4 | Weak crust (wet quartzite ellipses) | 11 | 2 | -10e6 | 2.5e20 |
| 5 | Weak mantle (wet olivine ellipses) | 41 | 2 | -10e6 | 2.5e20 |
| 6 | Strong crust (Maryland diabase ellipses) | 11 | 2 | -10e6 | 2.5e20 |
| 7 | Sediments 1 | 11 | 2 | -10e6 | 2.5e20 |
| 8 | Sediments 2 | 11 | 2 | -10e6 | 2.5e20 |

`T_st = -10e6` (code default) is uniform across all phases, matching the MDLIB default value.

`eta_vp = 2.5e20` (visco-plastic viscosity) is set on all phases, matching RiftingMelting for lithospheric-scale stability.

**Alternatives considered**: Simplify to 4 phases (like RiftingMelting) — rejected because the multi-layer ellipse geometry is a key requirement.

### 2. Solver configuration: Picard-to-Newton with moderate iterations

**Choice**: Use `Newton=0`, `Picard2Newton=1`, `Picard2Newton_tol=1e-1`, `nit_max=20`, `line_search=1`.

**Rationale**: PressurizedMagmaChamber uses `nit_max=50` with `Picard2Newton_tol=2e-5` for a small-scale model. For a lithospheric-scale model with more degrees of freedom, start with moderate iterations (20) and a relaxed Picard-to-Newton threshold (1e-1, matching RiftingMelting). The standard `line_search=1` is sufficient — `tensile_line_search` is not needed for this scenario.

**Alternatives considered**: Pure Picard (`nit_max=1` like MLPS_Ellipses) — rejected because the combined yield function benefits from Newton iterations for convergence.

### 3. Melting configuration: Enable globally with melt weakening

**Choice**: `melting=1`, `melt_weak=30.0`, `force_melt_weak=1` (matching RiftingMelting).

**Rationale**: RiftingMelting's proven configuration. A melt weakening factor of 30 produces realistic strain localization at partially molten zones.

### 4. Two-tier verification: Quick smoke-test then medium run

**Choice**: Two separate verify configs, both sharing the same `.c` file:

1. **Quick** (`RiftingCombinedYield_quick.txt`):
   - Resolution: `Nx=101`, `Nz=81` (~3 km horizontal, ~1.5 km vertical)
   - Time steps: `Nt=10`
   - Purpose: Does it initialize and survive the first few steps? (~minutes)

2. **Medium** (`RiftingCombinedYield_medium.txt`):
   - Resolution: `Nx=101`, `Nz=81` (same coarse grid)
   - Time steps: `Nt=100`
   - Purpose: Does it remain stable over longer evolution? Catches slow divergence, thermal runaway, or melt-related instabilities that only appear after tens of steps.

**Workflow**: Run quick first. If it passes, run medium. If medium passes, proceed to production resolution.

**Rationale**: 10 steps catches initialization failures and immediate divergence but can miss instabilities that develop over time (e.g., melt feedback loops, thermal runaway). 100 steps at coarse resolution is still fast but exercises longer-term stability.

**Alternatives considered**: (a) Single verify config — rejected because 10 steps misses slow-onset issues; (b) Medium at higher resolution — rejected because the goal is speed, not accuracy.

### 5. Boundary conditions: Pure shear extension from MLPS_Ellipses

**Choice**: Reuse `SetPureShearBCVx` / `SetPureShearBCVz` with `bkg_strain_rate = -1.056e-15` (1 cm/yr extension rate), `pure_shear_ALE=1`.

**Rationale**: Identical to MLPS_Ellipses. Pure shear extension is the standard rifting mode. The built-in pure shear BC functions handle velocity boundary conditions correctly.

### 6. Thermal boundary conditions: Same as MLPS_Ellipses

**Choice**: Top = 0°C (surface/free surface), Bottom = 1330°C, lateral = zero heat flux. Shear heating on, adiabatic heating off.

**Rationale**: Standard continental geotherm boundaries. The `SetBCT` callback from MLPS_Ellipses is directly reusable.

### 7. C file structure: Minimal adaptation of MLPS_Ellipses.c

**Choice**: Copy MLPS_Ellipses.c → RiftingCombinedYield.c with these changes:
- Update `RunMDOODZ` call to reference `"RiftingCombinedYield.txt"`
- Keep all ellipse geometry functions unchanged
- Keep `SetPhase`, `SetTemperature`, `SetDensity`, `SetBCT` callbacks as-is

**Rationale**: The .c file defines geometry and callbacks. The new physics (plast=2, melting) are configured entirely through the .txt parameter file — no C-level changes needed beyond the filename.

## Risks / Trade-offs

**[Convergence with combined yield at lithospheric scale]** → Start with moderate `nit_max=20` and Picard-to-Newton. If divergence occurs, increase `nit_max` or tighten `Picard2Newton_tol`. The verify config enables fast iteration on solver parameters.

**[Melting + combined yield interaction untested at this scale]** → No existing scenario combines both at lithospheric scale. The quick-verify workflow de-risks this by catching explosions early. If interaction causes issues, melting can be disabled temporarily (`melting=0`) to isolate the problem.

**[Phase 6 high conductivity (k=30) may interact unexpectedly with melting]** → Phase 6 (impregnated asthenosphere with diabase ellipses) has anomalously high thermal conductivity. Combined with mantle melting, this could create thermal artifacts. Monitor in verification runs.

**[Tensile strength at default value]** → Uniform `T_st = -10e6` (code default) is used for all phases. Can be differentiated per phase later if needed.
