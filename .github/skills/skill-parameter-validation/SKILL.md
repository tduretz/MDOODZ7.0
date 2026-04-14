---
name: skill-parameter-validation
description: Validate MDOODZ .txt parameter files — check scaling ranges, domain/resolution bounds, per-phase material property limits, flow law index consistency, switch compatibility, Courant/CFL constraints, solver tolerances, and anisotropy settings against physical limits.
---

# Parameter Validation

Use this checklist to validate a `.txt` parameter file before running a simulation.
Values in `.txt` files are in **SI units** (not scaled).

---

## 1. Scaling Parameters

| Parameter | Field | Valid range | Typical |
|-----------|-------|-------------|---------|
| Reference viscosity | `eta` | 1e15–1e28 Pa·s | 1e20–1e22 |
| Reference length | `L` | 1e2–1e7 m | 1e4–1e6 |
| Reference velocity | `V` | 1e-15–1e-5 m/s | 1e-10–1e-8 |
| Reference temperature | `T` | 1–1e4 K | 1–1000 |

**Check**: All four must be positive. The ratio `η·V/L` gives the stress scale — verify it yields Pa-range stresses (1e4–1e9).

---

## 2. Domain and Resolution

| Parameter | Field | Valid range | Notes |
|-----------|-------|-------------|-------|
| Grid cells x | `Nx` | ≥ 20 | Typical 100–500 |
| Grid cells z | `Nz` | ≥ 20 | Typical 100–500 |
| x bounds | `xmin`, `xmax` | finite, xmin < xmax | metres |
| z bounds | `zmin`, `zmax` | finite, zmin < zmax | metres |

**Check**: Cell aspect ratio `dx/dz` should be 0.5–2.0 for solver stability. If `anisotropy = 1`, use at least Nx, Nz ≥ 50.

---

## 3. Time Stepping and CFL

| Parameter | Field | Valid range | Notes |
|-----------|-------|-------------|-------|
| Number of steps | `Nt` | ≥ 1 | |
| Initial time step | `dt` | > 0 | SI seconds (before scaling) |
| Courant number | `Courant` | 0.1–0.9 | Default 0.5, must be < 1 |
| RK order | `RK` | 1, 2, or 4 | 4 recommended |
| Max time step | `dt_max` | ≥ dt | |
| Constant dt | `constant_dt` | 0 or 1 | 0 = adaptive (recommended) |

**CFL constraint**: `Courant × dx / V_max < dt`. If `constant_dt = 1`, verify manually.

---

## 4. Per-Phase Material Properties

For each phase `k` (0 to `Nb_phases - 1`):

### Density and Thermodynamics

| Property | Field | Valid range | Units |
|----------|-------|-------------|-------|
| Density | `rho` | 800–4000 | kg/m³ |
| Shear modulus | `G` | 1e9–1e11 | Pa |
| Heat capacity | `Cp` | 500–2000 | J/kg/K |
| Thermal conductivity | `k` | 0.1–20 | W/m/K |
| Radiogenic heat | `Qr` | 0–1e-5 | W/m³ |
| Thermal expansivity | `alp` | 1e-7–1e-3 | K⁻¹ |
| Compressibility | `bet` | 1e-14–1e-8 | Pa⁻¹ |

### Strength

| Property | Field | Valid range | Units |
|----------|-------|-------------|-------|
| Cohesion | `C` | 1e5–1e9 | Pa |
| Friction angle | `phi` | 0–45 | degrees in .txt |
| Dilation angle | `psi` | 0–30 | degrees in .txt |
| Stress limiter | `Slim` | ≥ C | Pa |
| Cohesion softening | `coh_soft` | 0 or 1 | 0 = off (default is 1!) |
| Softened cohesion | `Ce` | ≤ C | Pa (only if `coh_soft=1`) |
| Softened friction | `phie` | ≤ phi | degrees (only if `phi_soft=1`) |
| Softened dilation | `psie` | ≤ psi | degrees (only if `psi_soft=1`) |

**Check**: `psi ≤ phi` always. `Slim ≥ C` always. Softened end-values ≤ initial values.

> **Trap**: `coh_soft` defaults to 1 if omitted from `.txt`. When `coh_soft=1` and `C == Ce` (or `Ce` is not set, so it defaults to `C`), MDOODZ calls `exit(122)` with "Please set a difference in cohesion". Always set `coh_soft = 0` explicitly when softening is not intended.

> **Naming**: Softened end-values use abbreviated names in `.txt` files: `Ce` (not `C_end`), `phie` (not `phi_end`), `psie` (not `psi_end`).

### Flow Laws

| Flag | Field | Valid indices | Mechanism |
|------|-------|--------------|-----------|
| Constant | `cstv` | 0 or ≥ 1 | Fixed η (uses `eta0`) |
| Power law | `pwlv` | 0, 1, 10–46 | Dislocation creep |
| Linear | `linv` | 0, 40–46 | Diffusion creep |
| GBS | `gbsv` | 0, 40–46 | Grain boundary sliding |
| Exponential | `expv` | 0, 40–46 | Peierls creep |

> **Critical**: `linv=1`, `expv=1`, `gbsv=1` do NOT exist in `FlowLaws.c` and will call `exit(12)`. Only `pwlv=1` is valid as user-defined. For built-in flow laws, use index 40+ (e.g. 40 = Dry Olivine).

**Check**: At least one mechanism must be nonzero for every phase. Common indices:
- 10 = Wet Quartzite, 11/12 = Maryland Diabase, 40 = Dry Olivine, 41 = Wet Olivine.

### Anisotropy (per phase, if `anisotropy = 1`)

| Property | Field | Valid range | Notes |
|----------|-------|-------------|-------|
| Aniso factor | `aniso_factor` | 0.01–100 | 1.0 = isotropic |
| Initial angle | `aniso_angle` | 0–360 | degrees in .txt |
| Max factor | `ani_fac_max` | ≥ aniso_factor | |
| Evolve with strain | `ani_fstrain` | 0 or 1 | |

---

## 5. Switch Compatibility

These combinations must be consistent:

| If this is set... | Then require... | Reason |
|-------------------|----------------|--------|
| `anisotropy = 1` | `finite_strain = 1` | Fabric evolution needs F tensor |
| `shear_style = 1` (simple shear) | `periodic_x = 1` | Periodic BCs needed |
| `shear_style = 0` (pure shear) | `periodic_x = 0` | No periodicity |
| `Newton = 1` | `line_search = 1` | Newton needs line search for stability |
| `anisotropy = 1` | `lin_solver = 2` | Direct solver required |
| `psi > 0` + `psi_soft = 1` | `compressible = 1` | Dilation needs compressibility |
| `elastic = 1` | `free_surface = 1` (recommended) | Stress builds need surface |
| `melting = 1` | Valid `phase_diagram` | Melt model needs lookup table |
| `surface_processes > 0` | `free_surface = 1` | Surface processes need free surface |
| No-slip velocity BCs (type 0 everywhere) | `pure_shear_ALE = 0` | Conflicts with kinematic ALE BCs → singular matrix |
| `plast = 1` + `elastic = 1` | `bkg_strain_rate` large enough that `2·G·dt·ε̇ > C·cos(φ)` | Otherwise plasticity never activates |
| `plast = 2` (combined yield) | `T_st` set per phase (negative, e.g. -10e6) | Tensile cap requires tensile strength |
| `plast = 2` (combined yield) | `eta_vp > 0` per phase | Numerical stability; use ≥ 2.5e20 at lithospheric scale with melting |
| `coh_soft = 1` | `plast = 1` or `plast = 2` + sufficient trial stress | Softening requires plastic strain accumulation |

---

## 6. Solver Parameters

| Parameter | Field | Valid range | Typical |
|-----------|-------|-------------|---------|
| Penalty | `penalty` | 1e2–1e6 | 1e3–1e5 |
| Min viscosity | `min_eta` | > 0 | 1e18 Pa·s (scaled) |
| Max viscosity | `max_eta` | > min_eta | 1e24 Pa·s (scaled) |
| Max η/min η ratio | — | ≤ 1e10 | 1e6–1e9 |
| Nonlinear iterations | `nit` | 1–100 | 5–20 |
| Linear solver | `lin_solver` | 1 or 2 | 2 (CHOLMOD direct) |
| Abs div tolerance | `lin_abs_div` | 1e-12–1e-6 | 1e-9 |
| Abs mom tolerance | `lin_abs_mom` | 1e-12–1e-6 | 1e-9 |
| Rel div tolerance | `lin_rel_div` | 1e-8–1e-3 | 1e-5 |
| Rel mom tolerance | `lin_rel_mom` | 1e-8–1e-3 | 1e-5 |

---

## 7. Boundary Conditions

| BC code | Meaning |
|---------|---------|
| 0 | Dirichlet (physical) |
| 11 | Dirichlet (non-physical) |
| 2 | Neumann (non-physical) |
| 13 | Neumann (physical) |
| -2 | Periodic |
| 30 | Air (inactive) |

**Check**: If `periodic_x = 1`, left and right BCs must use code -2 for Vx. Free surface (`free_surface = 1`) requires top boundary flagged appropriately.

---

## Quick Validation Checklist

Run through these numbered checks before launching a simulation:

1. **Scales positive**: η > 0, L > 0, V > 0, T > 0
2. **Domain valid**: xmin < xmax, zmin < zmax, Nx ≥ 20, Nz ≥ 20
3. **Aspect ratio**: 0.5 ≤ dx/dz ≤ 2.0
4. **Time step**: dt > 0, Courant < 1.0, RK ∈ {1, 2, 4}
5. **Phase count**: Nb_phases matches number of phase blocks in .txt
6. **Each phase has flow law**: at least one of cstv/pwlv/linv/gbsv/expv nonzero
7. **Flow law indices valid**: pwlv 0–46, linv/gbsv/expv 0 or 40–46 (NOT 1–39!)
8. **Density reasonable**: 800 < ρ < 4000 kg/m³
9. **Strength hierarchy**: Slim ≥ C, psi ≤ phi, softened ≤ initial
9b. **Cohesion softening**: If `coh_soft=1` (or omitted!), verify `Ce < C`. Otherwise set `coh_soft = 0`
10. **Switch compatibility**: anisotropy→finite_strain, Newton→line_search, shear_style↔periodic_x
10b. **Combined yield (`plast=2`)**: Verify `T_st` is set (negative) and `eta_vp > 0` per phase. With `melting=1`, use `eta_vp ≥ 2.5e20` at lithospheric scale
11. **Viscosity bounds**: min_eta < max_eta, ratio ≤ 1e10
12. **Thermal consistency**: If thermal = 1, check Cp > 0, k > 0 for all phases
13. **Anisotropy consistency**: If anisotropy = 1, check ani_fac_max ≥ aniso_factor per phase
14. **Output frequency**: writer_step < Nt (otherwise no output)
15. **writer_subfolder**: Set to a unique name for each simulation (especially in tests) to avoid output file collisions
16. **GBS + elastic**: If using `gbsv > 0`, set `elastic = 1` to avoid NaN from internal G override
