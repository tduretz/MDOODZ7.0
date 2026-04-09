---
name: skill-melting
description: MDOODZ partial melting module — melt model indices, solidus/liquidus parameterizations, melt fraction computation, melt weakening, per-phase melt configuration, kinetic relaxation, and melting scenario examples.
---

# Partial Melting

## Overview

MDOODZ supports pressure-dependent partial melting via parameterized solidus/liquidus curves. Melting is computed per marker using the marker's temperature and pressure. The melt fraction field (`X`) is stored at cell centres and written to HDF5 as `Centers/X`.

## Enabling Melting

### Global Switches (in `.txt` file)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `melting` | 0 | Enable melting globally (0=off, 1=on) |
| `force_melt_weak` | 0 | Enable melt-induced viscosity weakening (0=off, 1=on) |
| `melt_weak` | 0.0 | Activation volume override for weakening [J/mol] — typical: 30.0 |

### Per-Phase Parameter (in `.txt` PHASE PROPERTIES blocks)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `melt` | 0 | Melt model index (0=off, see table below) |
| `tau_kin` | 0.0 | Kinetic relaxation time constant [s] — 0.0 = instantaneous equilibrium |

## Melt Model Indices

| `melt` | Name | Solidus | Liquidus | Latent heat | Reference |
|--------|------|---------|----------|-------------|-----------|
| 1 | Simple T-dependent | 600°C + 273 K (fixed) | 1200°C + 273 K (fixed) | 320 kJ/kg | Stüwe (1995) |
| 10 | Crustal/sediment (P-dependent) | P-dependent formula | 1423 + 0.105·P (MPa) | 380 kJ/kg | — |
| 11 | **Felsic/crustal** (P-dependent) | P-dependent formula | 1262 + 0.09·P (MPa) | 300 kJ/kg | Schenker et al. (2012) |
| 40 | Mantle peridotite (P-dependent) | P-dependent formula | 2073 + 0.114·P (MPa) | 400 kJ/kg | — |
| 41 | **Mantle wet olivine** (P-dependent) | P-dependent formula | 2073 + 0.114·P (MPa) | 300 kJ/kg | Schenker et al. (2012) |

### Choosing the Right Model

| Phase type | Recommended `melt` | Rationale |
|------------|---------------------|-----------|
| Continental crust (e.g., quartzite, anorthite, granite) | 11 | Felsic solidus, Schenker 2012 |
| Sediments | 10 or 11 | Lower solidus than mantle |
| Lithospheric mantle (dry olivine) | 41 | P-dependent peridotite |
| Asthenosphere (wet olivine) | 41 | Same model, wet conditions |
| Diabase / mafic lower crust | 11 | Treated as crustal |

## Melt Fraction Computation

Source: `MDLIB/MeltingRoutines.c`, function `PartialMelting()`

**Equilibrium melt fraction**:

φ_eq = (exp(α·T) − exp(α·T_s)) / (exp(α·T_l) − exp(α·T_s))

where T_s = solidus temperature, T_l = liquidus temperature, both P-dependent.

**With kinetic relaxation** (`tau_kin > 0`):

φ = (dt / (τ_k + dt)) · φ_eq + (τ_k / (τ_k + dt)) · φ_0

where φ_0 is the previous melt fraction and τ_k = `tau_kin`.

**Clamped** to [0, 1].

## Melt Weakening

When `force_melt_weak = 1` and `melt_weak > 0`, the code overrides the dislocation creep activation volume for melt-bearing phases:

```c
if (model->force_melt_weak == 1) mat->apwl[k] = model->melt_weak;
```

This effectively reduces viscosity in regions with high melt fraction. The `melt_weak` value replaces the power-law activation volume (`Vpwl`), making melt-bearing regions weaker.

**Typical value**: `melt_weak = 30.0` (used in `RiftingMelting`, `RiftingCombinedYield`, and other scenarios).

## Interaction with Other Features

### With `eta_vp` (viscoplastic regularisation)

When combining melting with plasticity (`plast = 1` or `plast = 2`), `eta_vp` must be large enough to stabilise the weakened rheology.

| Model scale | Without melting | With melting |
|-------------|-----------------|--------------|
| Lithospheric rifting | eta_vp ≥ 1e19 | **eta_vp ≥ 2.5e20** |
| Localised chamber | eta_vp ≥ 1e18 | eta_vp ≥ 1e19 |

> **Stability warning**: `eta_vp = 1e19` with `melting = 1` at lithospheric scale can produce instabilities (negative viscosity, values ~1e60) after ~100 steps. Use `eta_vp = 2.5e20` or higher.

### With combined yield (`plast = 2`)

No special interaction — `plast = 2` works with melting. Ensure `T_st` and `eta_vp` are set per phase.

## Example Configuration

```
### Global (in header of .txt)
melting           = 1
force_melt_weak   = 1
melt_weak         = 30.0

### Per phase
ID   = 0           // Continental crust
melt = 11          // Schenker 2012 felsic

ID   = 1           // Lithospheric mantle
melt = 41          // Schenker 2012 mantle

ID   = 2           // Asthenosphere
melt = 41          // Schenker 2012 mantle
```

## Scenarios Using Melting

| Scenario | Melt models used | Application |
|----------|-----------------|-------------|
| `RiftingMelting` | 11 (crust), 41 (mantle) | Extension with decompression melting |
| `RiftingCombinedYield` | 11 (crust), 41 (mantle) | Rifting + combined mode-I/II yield |
| `SubductionMelting` | Fluid-released melting | Slab-edge weakening |
| `TCMagmaticSystem` | 10, 11 | Thermo-chemical magmatic evolution |
| `PressurizedMagmaChamber` | — (thermal only) | Pressurised intrusion |
| `MeltingOverpressure` | — | Melting-induced overpressure |
| `ChristmasTree` | Multi-branch | Dike geometry |

## Key Source Files

- Melting computation: `MDLIB/MeltingRoutines.c`
- Melt weakening: `MDLIB/FlowLaws.c` (activation volume override)
- Parameter parsing: `MDLIB/InputOutput.c`
- HDF5 output: `Centers/X` (melt fraction field)
