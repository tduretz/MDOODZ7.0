---
name: skill-anisotropy
description: MDOODZ mechanical anisotropy module — director-based fabric tracking, anisotropy factors (δ_n, δ_s), viscous transverse isotropy, ViscosityConciseAniso, finite strain evolution, anisotropic rifting methodology, and anisotropy scenario examples.
---

# Mechanical Anisotropy

## Overview

MDOODZ models mechanical anisotropy via **director-based transverse isotropy**. Each material particle carries a director vector (n₁, n₂) representing the orientation of the weak plane (foliation). The viscous response depends on the angle between the director and the applied stress — deformation is easier along the foliation than across it.

This implementation follows the methodology described in the anisotropic rifting paper, where inherited crustal fabric controls rift symmetry and strain localisation.

## Enabling Anisotropy

In the `.txt` file:
```
anisotropy = 1
```

**Important**: Setting `anisotropy = 1` automatically forces Newton iteration (`model.Newton = 1`) because the anisotropic constitutive law requires consistent linearisation for convergence.

## Director Field

Each particle carries an anisotropy angle θ stored in `markers.aniso_angle`. The director is:

(n₁, n₂) = (cos θ, sin θ)

### Initialisation

Set initial fabric orientation per phase in `.txt`:
```
aniso_angle = 45    // degrees, per phase
```

Or use the `SetAnisoAngle` callback for spatially variable orientation:
```c
double SetAnisoAngle(MdoodzInput *input, Coordinates coords, int phase, double predefined) {
    return predefined;  // or compute custom angle
}
```

### Evolution

The director rotates with material flow (advected by the vorticity tensor). When `ani_fstrain = 1`, the anisotropy factor also evolves with accumulated finite strain.

## Anisotropy Factor

The anisotropy factor controls the degree of viscous anisotropy:

| Parameter | Description | Typical range |
|-----------|-------------|---------------|
| `aniso_factor` | Initial ratio of directional viscosity contrast | 1.0–10.0 |
| `ani_fac_max` | Maximum anisotropy factor (saturation limit) | Same as or > `aniso_factor` |
| `ani_fstrain` | Evolve factor with finite strain (0/1) | 0 or 1 |

- **aniso_factor = 1.0**: Isotropic (no anisotropy effect)
- **aniso_factor > 1.0**: Increasing anisotropy — deformation is easier along the foliation

### Finite Strain Evolution

When `ani_fstrain = 1`:
```c
double AnisoFactorEvolv(double FS_AR, double aniso_fac_max) {
    return MINV(FS_AR, aniso_fac_max);
}
```

The anisotropy factor equals the finite strain aspect ratio (`FS_AR`), clamped at `ani_fac_max`.

## Anisotropic Invariants

The code modifies the stress invariant to account for anisotropy:

**Isotropic**: τ_II² = ½(τ_xx² + τ_zz² + τ_yy²) + τ_xz²

**Anisotropic** (`Y2` function): τ_II² = ½(τ_xx² + τ_zz² + τ_yy²) + (δ·τ_xz)²

where δ is the anisotropy factor. The strain rate invariant `I2` remains isotropic.

## Viscosity Computation

The anisotropic viscosity is computed in `ViscosityConciseAniso`. The key steps:

1. Rotate stress and strain rate tensors into the director frame (angle θ - π/2)
2. Apply the viscosity law in the rotated frame (where the anisotropy factor modifies shear components)
3. The effective viscosity depends on the angle between director and principal stress

## Grid-Level Fields

Anisotropy data interpolated from particles to the grid:

| Field | Location | Description |
|-------|----------|-------------|
| `aniso_factor_n` | Cell centres | Anisotropy factor |
| `aniso_factor_s` | Vertices | Anisotropy factor |
| `d1_n`, `d2_n` | Cell centres | Director components |
| `d1_s`, `d2_s` | Vertices | Director components |
| `angle_n` | Cell centres | Director angle |
| `angle_s` | Vertices | Director angle |
| `FS_AR_n` | Cell centres | Finite strain aspect ratio |
| `FS_AR_s` | Vertices | Finite strain aspect ratio |

## Anisotropy Scenarios

### Rifting with Anisotropy
| Scenario | Description |
|----------|-------------|
| `RiftingAnisotropy` | Continental rift with anisotropic crust |
| `RiftingCheninAniso` | Multi-layer rift (Chenin rheology) + anisotropy |

### Simple Shear Tests
| Scenario | Description |
|----------|-------------|
| `SimpleShearAnisoHomo` | Homogeneous simple shear with anisotropy |
| `ShearTemplateAniso` | Shear template with anisotropic material |
| `ShearInclusionAniso` | Shear around anisotropic inclusion |

### Viscous Anisotropy Tests
| Scenario | Description |
|----------|-------------|
| `AnisoViscTest_MWE` | Minimal working example |
| `AnisoViscTest_HomoRandomAngle` | Random initial fabric orientations |
| `AnisoViscTest_evolv` | Evolving anisotropy factor |
| `AnisoViscTest_Ellipses` | Elliptical inclusions |
| `AnisoViscTest_Debug` | Debug configuration |

### VE-VP Anisotropy
| Scenario | Description |
|----------|-------------|
| `AnisoVEVP` | Visco-elasto-visco-plastic anisotropy |
| `AnisoHomoVEVP` | Homogeneous VEVP anisotropy |

### Benchmarks and Others
| Scenario | Description |
|----------|-------------|
| `AnisotropyDabrowski` | Published benchmark (Dabrowski et al.) |
| `AnisotropyLayer` | Layered anisotropy benchmark |
| `TestAnisotropy` | Basic verification test |
| `CollisionAnisotropy` | Collision with anisotropic rheology |
| `CollisionPolarCartesianAniso` | Polar-Cartesian collision + anisotropy |
| `CrustalShearBandingAnisotropy` | Crustal shear band formation |
| `WedgeAnisotropy` | Accretionary wedge |

## Scientific Background

The anisotropy implementation follows the **transverse isotropy** approach where:
- Inherited crustal fabric (foliation) creates directional weakness
- The anisotropy factor controls how much easier deformation is along vs. across foliation
- Initial fabric orientation controls rift symmetry: fabric parallel to rift axis promotes symmetric rifting; oblique fabric promotes asymmetry
- Finite strain evolution allows fabric to strengthen with progressive deformation

See the anisotropic rifting paper for detailed analysis of how anisotropy controls rift architecture, including the role of δ in governing strain localisation patterns.

## Key Source Files

- Anisotropy routines: `MDLIB/AnisotropyRoutines.c`
- Grid anisotropy fields: `MDLIB/mdoodz-private.h` (grid struct)
- Marker anisotropy: `MDLIB/mdoodz-private.h` (markers struct)
- Anisotropy tests: `AnistropyUnitTests/`
- CI benchmark tests: `TESTS/AnisotropyBenchmarkTests.cpp`

## Quantitative CI Benchmarks

The anisotropy module has 3 CI benchmark tests in `TESTS/AnisotropyBenchmarkTests.cpp`:

| Test | Analytical Solution | Measured | Threshold |
|------|---------------------|----------|-----------|
| `DirectorEvolution` | θ(t) = arctan(tan(θ₀) − γ̇·t) | L2(θ) = 2.34e-3 rad | < 5e-3 |
| `DirectorDtConvergence` | Forward Euler order ~1.0 | order = 0.93–0.99 (5 dt values) | ≥ 0.8 |
| `StressAnisotropy` | Rotation+scaling formula | relErr(τ_II) = 1.1e-8 | < 0.01% |

The `DtConvergence` test extracts real MDOODZ θ(t) trajectories from HDF5 at every step (writer_step=1) and writes `director_trajectory_dt*.dat` files for gnuplot. Per-step error grows monotonically due to forward Euler overshoot — from +0.003°/step to +0.004°/step at dt=0.0125, and final errors halve with each dt halving (first order).

The director evolution ODE dθ/dt = −γ̇·cos²(θ) is documented in `TESTS/AnalyticalSolutions.md` §3.
The stress comparison uses the same rotation formula as `ViscosityConciseAniso` (AnisotropyRoutines.c line 80).
