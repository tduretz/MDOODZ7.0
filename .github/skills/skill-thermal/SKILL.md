# Skill — MDOODZ Thermal Solver

## Overview

MDOODZ solves the energy equation coupled with the Stokes equations.  Temperature is **carried on markers** (particles) and interpolated to/from the staggered grid.  The thermal solve is implicit (backward-Euler) and uses a CHOLMOD direct solver.

---

## 1. Energy Equation

The discretised energy equation solved on cell centres is:

$$
\rho C_p \frac{T^{n+1} - T^n}{\Delta t} = \nabla \cdot (k \nabla T) + H_r + H_s + H_a
$$

where:
- $\rho$ — density (cell centres, `mesh->rho_n`)
- $C_p$ — specific heat capacity (`mesh->Cp`)
- $k$ — thermal conductivity (`mesh->kx`, `mesh->kz` on Vx/Vz nodes)
- $H_r$ — radiogenic heat production (`mesh->Qr`)
- $H_s$ — shear (dissipative) heating (`mesh->Wdiss`)
- $H_a$ — adiabatic heating

### Source-file location

| File | Function | Purpose |
|------|----------|---------|
| `MDLIB/ThermalRoutines.c` | `EnergyDirectSolve()` | Matrix assembly + CHOLMOD direct solve |
| `MDLIB/ThermalRoutines.c` | `ThermalSteps()` | Initial-cooling sub-stepping |
| `MDLIB/ThermalRoutines.c` | `SetThermalPert()` | Apply initial thermal perturbation |
| `MDLIB/RheologyParticles.c` | `UpdateParticleEnergy()` | Grid → particle T update (with optional subgrid diffusion) |

---

## 2. Thermal Boundary Conditions

Thermal BCs are set via the **`SetBCT`** callback registered in `MdoodzSetup`.  The BC array uses an **extended grid** of size `(Ncx+2) × (Ncz+2)` stored in `mesh->BCT_exp`.

### BCT type codes

| `type` | Meaning | Notes |
|--------|---------|-------|
| `0` | **Neumann** (zero heat flux) | Default for sidewalls |
| `1` | **Dirichlet** (fixed temperature) | Set `val` in the callback; common for top/bottom |
| `30` | **Not calculated** | Cell excluded from solve (e.g., air) |
| `-2` | **Periodic** | Links east↔west boundaries |

### SetBCT callback pattern

```c
void SetBCT(MdoodzInput *input, POSITION position,
            double *value, int *type) {
  if (position == S || position == SE || position == SW) {
    *type  = 1;          // Dirichlet
    *value = 1273.15;    // Fixed T in Kelvin
  }
  if (position == N || position == NE || position == NW) {
    *type  = 1;
    *value = 273.15;
  }
  if (position == W || position == E) {
    *type  = 0;          // Neumann (zero flux)
    *value = 0.0;
  }
}
```

---

## 3. Per-Phase Thermal Material Properties

Set in the `.txt` parameter file under each phase block:

| Parameter | `.txt` key | SI default | SI unit | Scaling divisor |
|-----------|-----------|------------|---------|-----------------|
| Thermal conductivity | `k` | `1.0e-6` | W m⁻¹ K⁻¹ | `scaling.k` |
| Specific heat capacity | `Cp` | `1.0e3` | J kg⁻¹ K⁻¹ | `scaling.Cp` |
| Thermal expansivity | `alp` | `0.0` | K⁻¹ | `1/scaling.T` |
| Compressibility | `bet` | `1.0e-40` | Pa⁻¹ | `1/scaling.S` |
| Radiogenic heat production | `Qr` | `1.0e-30` | W m⁻³ | `scaling.W / scaling.L³` |

The internal variable `materials.k_eff[phase]` is initialised from `materials.k[phase]` and is what gets interpolated to the grid via `P2Mastah`.  The **crazyConductivity** mechanism can override `k_eff` for specific phases.

---

## 4. Global Thermal Switches

Set in the `.txt` file (not per-phase):

| Parameter | `.txt` key | Default | Meaning |
|-----------|-----------|---------|---------|
| Enable thermal solver | `thermal` | `0` | `1` = solve energy equation each time step |
| Shear heating | `shear_heating` | `1` | `1` = add $H_s = \sigma'_{ij} \dot\varepsilon'_{ij}$ |
| Adiabatic heating | `adiab_heating` | `0` | `0` = off, `1` = lithostatic, `2` = full dP/dt |
| Initial cooling | `initial_cooling` | `0` | `1` = run `ThermalSteps()` before time loop |
| Cooling duration | `cooling_duration` | 1 Ga | Duration of initial cooling sub-stepping |
| Thermal perturbation | `therm_perturb` | `0` | `1` = apply localised T anomaly at start |
| Perturbation centre | `therm_perturb_x0`, `therm_perturb_z0` | `0.0` | SI position (scaled by L) |
| Perturbation radii | `therm_perturb_rad_x`, `therm_perturb_rad_z` | `0.0` | SI (scaled by L) |
| Perturbation amplitude | `therm_perturb_dT` | `0.0` | SI Kelvin (scaled by T) |

---

## 5. Heating Modes in Detail

### 5.1 Shear Heating (`shear_heating = 1`)

Dissipation power is computed during the rheological update and stored in `mesh->Wdiss`.  In the RHS:

```c
b[eqn] += mesh->Wdiss[c2];
```

### 5.2 Adiabatic Heating

**Mode 1** (`adiab_heating = 1`) — lithostatic pressure assumption:

$$
H_a = T \, \alpha \, \rho \, g \, V_z
$$

```c
Ha[c2] = mesh->T[c2] * mesh->alp[c2]
       * 0.5*(mesh->v_in[c3+1] + mesh->v_in[c3+1+nxvz])
       * model.gz * mesh->rho_n[c2];
```

**Mode 2** (`adiab_heating = 2`) — full material derivative of pressure:

$$
H_a = T \, \alpha \, \frac{P - P^{n-1}}{\Delta t}
$$

```c
Ha[c2] = mesh->T[c2] * mesh->alp[c2]
       * (mesh->p_in[c2] - mesh->p0_n[c2]) / model.dt;
```

### 5.3 Radiogenic Heat Production

Always included when `thermal = 1`:

```c
b[eqn] += Hr * mesh->Qr[c2];   // Hr = 1.0
```

---

## 6. Density — Boussinesq Coupling

The function `EvaluateDensity()` in `MDLIB/RheologyDensity.c` computes $\rho(T, P)$ with several models:

| `density_model` | Formula |
|-----------------|---------|
| `0` | $\rho = \rho_{\text{ref}}$ (constant) |
| `1` | $\rho = \rho_{\text{ref}} (1 - \alpha (T - T_0))(1 + \beta (P - P_0))$ — **Standard Boussinesq EOS** |
| `2` | Phase-diagram lookup $\rho(T, P)$ |
| `3` | $\rho = \rho_{\text{ref}} \exp(\beta P - \alpha T)$ |

For most thermal convection setups, use `density_model = 1` with appropriate `alp`, `bet`, and `T0`.

---

## 7. Temperature Advection via Markers

Temperature is **stored on markers** (`particles.T`) and advected with the flow automatically.  The grid–particle transfer cycle each time step is:

1. **Markers → Grid** (`P2Mastah`): `particles.T` → `mesh->T0_n` (old T on grid)
2. **Thermal solve** (`EnergyDirectSolve`): computes new `mesh->T` on the grid
3. **Grid → Markers** (`UpdateParticleEnergy`): incremental update $\Delta T = T^{n+1}_{\text{grid}} - T^n_{\text{grid}}$, then `particles.T += ΔT`

### Subgrid Diffusion

When `subgrid_diffusion >= 1`, the temperature increment is split:

- **Subgrid part** ($\Delta T_s$): analytical diffusion within each cell based on marker-vs-cell-average temperature difference and a diffusion timescale $\tau = \rho C_p / (k (1/\Delta x^2 + 1/\Delta z^2))$
- **Remaining part** ($\Delta T_r = \Delta T - \Delta T_s$): interpolated from grid to markers

This reduces numerical diffusion from marker–grid interpolation.

---

## 8. Matrix Assembly Notes

- The system is assembled in CSC format and solved with **CHOLMOD** (SuiteSparse).
- The stencil is 5-point (W, E, S, N, centre) on cell centres.
- Conductivity is staggered: `kx` on Vx-nodes, `kz` on Vz-nodes.
- Near free-surface cells (`type == 30` neighbours), a ghost-node Dirichlet correction is applied using averaged conductivity $k_s = \frac{1}{4}(k_W + k_E + k_N + k_S)$.
- Periodic thermal BCs (`type == -2`) wrap the east boundary equation to couple with the west-most DOF.

---

## 9. Practical Setup Guide

### Minimal thermally-coupled simulation

1. Set `thermal = 1` in the `.txt` file.
2. Provide per-phase `k`, `Cp`, and `alp` (thermal expansivity).
3. Implement `SetBCT` callback (Dirichlet top/bottom, Neumann or periodic sides).
4. Implement `SetTemperature` callback to initialise temperature on markers.
5. Set `density_model = 1` for Boussinesq buoyancy coupling.

### Checklist for thermal convection benchmarks

| Item | Recommendation |
|------|---------------|
| Scaling | For Blankenbach-style: η_ref=1, L=H, V=κ/H, T=ΔT |
| Resolution | ≥41×41 for Ra=10⁴; increase for higher Ra |
| Time stepping | Enough steps for thermal equilibrium (monitor Nu, Vrms) |
| `shear_heating` | Usually `0` for convection benchmarks |
| `adiab_heating` | Usually `0` for Boussinesq benchmarks |
| `free_surface` | `0` for enclosed-box convection |
| Gravity | Set `gz` such that $Ra = \rho_0 \alpha \Delta T g H^3 / (\kappa \eta)$ matches target |

### Common issues

| Symptom | Likely cause |
|---------|-------------|
| Temperature stays constant | `thermal = 0` (default!) — set to `1` |
| No convection develops | `density_model = 0` (constant ρ) — switch to `1` |
| T exceeds boundary values | BC callback sets wrong `type` or `value` |
| Very slow convergence | Conductivity too low → CFL-limited time step |
| Checkerboard pattern in T | Resolution too low or marker density insufficient |
