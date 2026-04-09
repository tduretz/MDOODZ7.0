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

Introduced in commit 8aecb92 (PR #159, Popov et al. 2025). Combines a standard Drucker-Prager shear yield surface with a tensile (mode-I) cap.

**Additional per-phase parameter**:

| Parameter | Field | Description | Typical value |
|-----------|-------|-------------|---------------|
| Tensile strength | `T_st` | Mode-I tensile cutoff [Pa] | -10e6 (must be negative) |

**Requirements**:
- `eta_vp > 0` for all phases (numerical stability)
- `T_st` set per phase (code default: -10e6 Pa)
- Works with both Picard and Picard2Newton solvers
- Optional: `tensile_line_search = 1` for additional solver stability in tensile-dominated problems

**Example scenarios**: `Popov2025_Pureshear_VEVP`, `Popov2025_Tensile_VEVP`, `RiftingCombinedYield`

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
