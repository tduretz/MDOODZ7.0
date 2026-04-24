---
name: skill-rheology
description: MDOODZ rheological framework — flow laws, viscosity computation, dislocation creep, diffusion creep, Peierls mechanism, grain boundary sliding, elasticity (Maxwell model), plasticity (Drucker-Prager, strain softening), material properties, and the mat_prop structure.
---

# Rheology and Material Properties

## Viscosity Computation Pipeline

For each grid point, the effective viscosity is computed by:

1. **Identify active mechanisms** per phase from `.txt` switches: `cstv`, `pwlv`, `linv`, `gbsv`, `expv`
2. **Compute each mechanism's viscosity** from T, P, ε̇, grain size, composition
3. **Combine via harmonic averaging**: 1/η_eff = 1/η_pwl + 1/η_lin + 1/η_gbs + 1/η_exp + 1/η_cst (weaker mechanism dominates)
4. **Apply elastic correction** (if `elastic = 1`): η_ve = η·G·Δt / (η + G·Δt) (Maxwell model)
5. **Apply plastic correction** (if yield stress exceeded): Drucker-Prager yield
6. **Clamp**: η_final ∈ [min_eta, max_eta]

The main viscosity function is `ViscosityConciseAniso` in `AnisotropyRoutines.c` (or `ViscosityConcise` for isotropic cases in `RheologyDensity.c`).

## Creep Mechanisms

### Power Law — Dislocation Creep (`pwlv`)

**Equation**: η_pwl = A^(-1/n) · f^(-r/n) · exp((Q + P·V)/(n·R·T)) · ε̇_II^(1/n - 1)

| Parameter | Symbol | Code field |
|-----------|--------|------------|
| Pre-exponential factor | A | `Apwl[phase]` |
| Stress exponent | n | `npwl[phase]` |
| Activation energy | Q | `Qpwl[phase]` |
| Activation volume | V | `Vpwl[phase]` |
| Grain size exponent | m | `mpwl[phase]` |
| Fugacity exponent | r | `rpwl[phase]` |
| Fugacity | f | `fpwl[phase]` |

**Flow Law Database** (selected entries from `FlowLaws.c`):

| `pwlv` | Name | Source | Typical use |
|--------|------|--------|-------------|
| 1 | User-defined (η₀, n, Q) | Custom | Any |
| 2 | User-defined (τ₀ constrained) | Custom | Any |
| 3 | User-defined (A, n, Q) | Custom | Any |
| 10 | Wet Quartzite | Ranalli 1995 | Upper crust |
| 11 | Maryland Diabase | Ranalli 1995 | Lower crust |
| 12 | Dry Granulite | Ranalli 1995 | Lower crust |
| 19 | Westerly Granite | Hansen & Carter 1983 | Upper crust |
| 25 | Plagioclase An75 | Ranalli 1995 | Crust |
| 40 | Dry Olivine | Hirth & Kohlstedt 2003 | Lithospheric mantle |
| 41 | Wet Olivine | Hirth & Kohlstedt 2003 | Asthenosphere |

### Linear — Diffusion Creep (`linv`)

**Equation**: η_lin = A^(-1/n) · d^(m/n) · f^(-r/n) · exp((Q + P·V)/(n·R·T)) · ε̇_II^(1/n - 1)

Grain-size sensitive (m > 0). Requires grain size initialisation via `SetGrainSize` callback.

| `linv` | Name | Source |
|--------|------|--------|
| 40 | Dry Olivine diffusion | Hirth & Kohlstedt 2003 |

When both `pwlv = 40` and `linv = 40` are set (e.g., for mantle), both dislocation and diffusion creep are active. They combine harmonically — the weaker mechanism dominates at each P-T-ε̇ condition.

### Exponential — Peierls Creep (`expv`)

Low-temperature, high-stress plasticity. Activated by `expv` switch. Only effective below T_max ≈ 1200°C (capped at `TmaxPeierls = 1200 + 273.15 K`).

| `expv` | Name | Source |
|--------|------|--------|
| 40 | Dry Olivine Peierls | Hirth & Kohlstedt 2003 |

### Grain Boundary Sliding (`gbsv`)

Intermediate mechanism between dislocation and diffusion creep. Grain-size and stress dependent.

### Constant Viscosity (`cstv`)

When `cstv > 0`, uses a constant viscosity η = η₀ (set via `eta0[phase]` in the `.txt` file). Overrides all creep mechanisms.

## Elasticity

**Maxwell visco-elastic model** (enabled by `elastic = 1`):

η_ve = η · G · Δt / (η + G · Δt)

Stress update (objective rate):
τ_new = 2·η_ve · (ε̇ + τ⁰_rotated / (2·G·Δt))

where τ⁰_rotated is the old stress rotated according to `stress_rotation` mode.

| `.txt` parameter | Description |
|-------------------|-------------|
| `elastic = 1` | Enable elasticity globally |
| `G` (per phase) | Shear modulus [Pa], typically ~1e10 |
| `stress_rotation` | Stress objectivity: `0`=none, `1`=Jaumann (default), `2`=analytical rotation |

The `stress_rotation` parameter controls how old deviatoric stresses are rotated before the visco-elastic update (`RheologyParticles.c`). Mode 0 disables rotation (valid only for irrotational flow). Mode 1 (Jaumann) uses the vorticity tensor. Mode 2 uses an analytical rotation matrix.

## Plasticity

### Plasticity Mode Switch (`plast`)

| `plast` | Model | Description |
|---------|-------|-------------|
| 0 | None | Plasticity disabled |
| 1 | Drucker-Prager | Standard shear yield (default) |
| 2 | Combined mode-I/mode-II | Popov et al. (2025) — shear + tensile yield |

### Drucker-Prager Yield (`plast = 1`)

**Yield criterion**:

τ_II ≤ C·cos(φ) + P·sin(φ)

If yield is exceeded, viscosity is reduced to enforce the yield stress:
η_pl = τ_yield / (2·ε̇_II)

### Combined Mode-I/Mode-II Yield (`plast = 2`)

Introduced in commit 8aecb92 (PR #159, Popov et al. 2025, doi:10.5194/gmd-18-7035-2025). Combines a standard Drucker-Prager shear yield surface with a circular tensile (mode-I) cap joined through a smooth delimiter point; implemented as a 3×3 implicit Newton local solve on `(τ_II, p, λ̇)` in `MDLIB/RheologyDensity.c:685-899`.

**Yield-surface geometry** (paper Eqs. 13–17) from four material parameters (`φ`, `ψ`, `c_MC`, `p_T`):

```
k    = sin(φ)              kq   = sin(ψ)               c = c_MC · cos(φ)
a    = sqrt(1 + k²)        b    = sqrt(1 + kq²)
p_y  = (p_T + c/a) / (1 − k/a)       R_y  = p_y − p_T
p_d  = p_y − R_y · k/a               τ_d  = k · p_d + c
```

- DP segment active when `τ_II·(p_y − p_d) ≥ τ_d·(p_y − p)`, otherwise cap
- DP envelope: `τ_II = k·p + c`
- Cap circle: centre `(p_y, 0)`, radius `R_y`; `τ_II² + (p − p_y)² = R_y²`

**Additional per-phase parameter**:

| Parameter | Field | Description | Typical value |
|-----------|-------|-------------|---------------|
| Tensile strength | `T_st` | Mode-I tensile cutoff [Pa] | -10e6 (must be negative) |

**Global parameter**:

| Parameter | Field | Description |
|-----------|-------|-------------|
| `tensile_line_search` | Global switch | 0 = raw Newton step; 1 = Armijo-damped line search (see caution below) |

**Requirements**:
- `T_st` set per phase (code default: -10e6 Pa)
- Typically `eta_vp > 0` (Perzyna regularisation) for numerical stability under aggressive loading, but `eta_vp = 0` (rate-independent / perfect plasticity, as in the paper's 0D tests) is supported
- Works with both Picard and Picard2Newton solvers

**Caution — `tensile_line_search = 1` can abort the process**: If the 3×3 local Newton's Armijo line search fails to find a sufficient step within 100 iterations (α shrinks below 0.01), `RheologyDensity.c:886` calls `exit(0)` — the simulation terminates silently mid-run. For setups with large per-step elastic overshoot (`2·G·ε̇·dt ≫ yield`) consider either (a) `tensile_line_search = 0` + undamped Newton, (b) smaller `dt`, or (c) a small `eta_vp ≥ 1e16` for Perzyna damping.

**Trial-pressure scheme (paper §3.3 / Eq. 31)**: The two-field mixed formulation uses `p*` = global discretised pressure (stored in `mesh->p_in`, written to HDF5 `/Centers/P`) and `p` = integration-point yield-clamped pressure. They are related by:

```
p = p* + K · θ̇_vp · dt               (Eq. 31)
```

Empirically in MDOODZ:
- **DP segment** (shear, mixed after delimiter crossing): the global Stokes solver converges → `p* ≈ p` (within 0.1 %). HDF5 `/Centers/P` already sits on the yield surface.
- **Tensile cap tip** (pure volumetric extension with `ψ > 0`): `p*` is offset from `p` by `K·θ̇_vp·dt` (≈ 8 % of τ_d for typical 0D parameters). Reconstruct via Eq. 31 using `/Centers/divu_pl`, or use a looser tolerance.

For the extensional cap case, `θ̇_vp > 0` (plastic dilation), so `p > p*` (`p_local` is less tensile than the global trial).

**0D test setup matching Popov 2025 Table 1 "Fig. 5 a, b"** (canonical verification of the local stress update):
- `G = 10¹⁰ Pa`, `K = 2·10¹¹ Pa`, `φ = 30°`, `ψ = 10°`, `c_MC = 10⁶ Pa`, `p_T = −5·10⁵ Pa`
- `η^vp = 0`, `Δt = 2 yr`
- Strain rates: volumetric `ε̇_xx = ε̇_zz = 2.333·10⁻¹⁵`; shear `ε̇_xy = 7·10⁻¹⁴`; mixed superposition
- BC convention: `Vx_W = -user2, Vx_E = +user2` → `ε̇_xx = 2·user2/Lx`, so `user2 = ε̇_xx·Lx/2` (**watch the factor of 2**)

**Example scenarios**: `Popov2025_Pureshear_VEVP`, `Popov2025_Tensile_VEVP`, `Popov2025_0DIntegration.{VolumetricExtension,DeviatoricShear,MixedStrain}` (CI), `RiftingCombinedYield`, `PressurizedMagmaChamber`, `TCMagmaticSystem`

### Strain Softening

| Parameter | Description |
|-----------|-------------|
| `C` | Initial cohesion [Pa] |
| `C_end` | Final cohesion after softening [Pa] |
| `phi` | Initial friction angle [degrees] |
| `phi_end` | Final friction angle [degrees] |
| `psi` | Dilation angle [degrees] |
| `psi_end` | Final dilation angle [degrees] |
| `pls_start` | Accumulated plastic strain at which softening begins |
| `pls_end` | Accumulated plastic strain at which softening is complete |
| `Slim` | Stress limiter [Pa] — maximum deviatoric stress |
| `sig_tens` | Tensile strength [Pa] |
| `eta_vp` | Viscoplastic regularisation viscosity [Pa·s] (see below) |
| `n_vp` | Viscoplastic stress exponent (must be 1.0 in MD7) |

### Viscoplastic Regularisation (`eta_vp`)

**Purpose**: `eta_vp` adds a rate-dependent (viscous) overstress to the plastic yield surface, preventing the stress from dropping to exactly the yield stress. This controls shear band width and improves numerical convergence.

**Modified yield condition**:

τ_II = C·cos(φ) + P·sin(φ) + η_vp · λ̇

where λ̇ is the plastic multiplier rate (related to the magnitude of plastic strain rate).

**Return mapping** (`RheologyDensity.c`, line ~630): A local Newton iteration finds λ̇ such that:

F = τ_trial - η_ve·λ̇ - τ_yield - η_vp·λ̇ = 0

The viscoplastic overstress is then `OverS = η_vp · λ̇`.

**Physical effect**: Higher `eta_vp` → wider shear bands, smoother localisation. Lower → thinner bands, more mesh-sensitive. Setting `eta_vp = 0` recovers ideal (rate-independent) plasticity.

**Typical values by model scale**:

| Application scale | Recommended `eta_vp` | Example scenario |
|-------------------|----------------------|------------------|
| Lithospheric rifting / collision | 2.5e20 Pa·s | `RiftingMelting`, `RiftingCombinedYield` |
| Crustal shear zones | 1e19–1e20 Pa·s | `CrustalShearBandingAnisotropy` |
| Localised magma chambers | 1e19 Pa·s | `PressurizedMagmaChamber`, `TCMagmaticSystem` |
| General range | 1e18–1e22 Pa·s | — |

> **Stability warning**: Using `eta_vp = 1e19` in lithospheric-scale rifting models (especially with `melting = 1`) can cause instability after ~100 time steps — negative viscosities and values ~1e60 in the output. Use `eta_vp ≥ 2.5e20` for lithospheric-scale models with melting.

**Note**: `n_vp` (viscoplastic exponent) must be 1.0 in MD7 — power-law viscoplasticity (`n_vp > 1`) is not implemented and will `exit(1)`.

## The mat_prop Structure

Holds per-phase material properties (up to 20 phases), populated from the `PHASE PROPERTIES` section of the `.txt` file. Arrays indexed by phase ID.

### Key Fields

| Field | Physical quantity | Units |
|-------|-------------------|-------|
| `rho[phase]` | Reference density | kg/m³ |
| `G[phase]` | Shear modulus | Pa |
| `Cp[phase]` | Heat capacity | J/kg/K |
| `k[phase]` | Thermal conductivity | W/m/K |
| `Qr[phase]` | Radiogenic heat production | W/m³ |
| `alp[phase]` | Thermal expansivity α | K⁻¹ |
| `bet[phase]` | Compressibility β | Pa⁻¹ |
| `C[phase]` | Cohesion | Pa |
| `phi[phase]` | Friction angle | degrees (converted to radians internally) |
| `psi[phase]` | Dilation angle | degrees |
| `eta0[phase]` | Reference viscosity | Pa·s |
| `n[phase]` | Stress exponent | - |
| `aniso_factor[phase]` | Anisotropy factor | - (1.0 = isotropic) |
| `aniso_angle[phase]` | Initial anisotropy angle | degrees |

## Density Model

**Equation of state**:

ρ = ρ₀ · (1 - α(T - T₀)) · (1 + β(P - P₀))

| Parameter | Symbol | Typical value |
|-----------|--------|---------------|
| `rho` | ρ₀ | 2700 (crust), 3300 (mantle) kg/m³ |
| `alp` | α | 3.2e-5 K⁻¹ |
| `bet` | β | 1e-11 Pa⁻¹ |

The `density_model` switch allows alternative formulations including phase-diagram-based density.

## Marker-to-Node Viscosity Averaging (`eta_average`)

When multiple material phases share a grid cell, their marker viscosities are averaged to the node. The `eta_average` .txt parameter selects the method:

| Value | Method | Best for |
|-------|--------|----------|
| 0 | Arithmetic (default) | Best convergence orders (Vx 1.02, P 0.75 in SolVi 41→81) |
| 1 | Harmonic | Lowest absolute Vx error (1.0e-3 at 51×51) but worst convergence orders (0.45) |
| 2 | Geometric | Best absolute P error at fixed resolution (18% lower than arithmetic) |

This averaging operates on the phase mixture **within** each cell. It is distinct from cell-face viscosity interpolation in the FD stencil (which is not currently implemented as harmonic).

## Key Source Files

- Flow law database: `MDLIB/FlowLaws.c`
- Viscosity evaluation: `MDLIB/RheologyDensity.c`, `MDLIB/AnisotropyRoutines.c`
- Rheology interface: `MDLIB/RheologyDensity.h`
- Particle rheology: `MDLIB/RheologyParticles.c`
