---
name: skill-scenario-gallery
description: Catalogue of 79 predefined MDOODZ scenarios in SETS/ — grouped by geodynamic application (rifting, subduction, collision, shear zones, anisotropy, necking, viscous flow, thermal, magmatic, phase transitions, topography, benchmarks). Each entry lists the .c file, purpose, and recommended use.
---

# Scenario Gallery

Every simulation is defined by a `.c` / `.txt` file pair in `SETS/`.
The `.c` file contains setup callbacks; the `.txt` file contains model parameters.

---

## Rifting

| File | Description | Key features |
|------|------------|--------------|
| `RiftingBasic.c` | Lithosphere-scale extension with crustal and mantle layers | Weak seed zone for strain localisation |
| `RiftingChenin.c` | Rheologically controlled rifting | Imposed topography perturbation, multi-layer lithosphere |
| `RiftingRoman.c` | Rifting with mantle diapirism | Perturbed Moho geometry, noise filtering |
| `RiftingMelting.c` | Extension with partial melting | Lower-crust and decompression mantle melting |
| `RiftingAnisotropy.c` | Rifting with crystal anisotropy | Olivine/OPX fabric evolution during extension |
| `RiftingCheninAniso.c` | Combined topography perturbation + anisotropy | Full aniso rifting workflow |
| `RiftingCombinedYield.c` | Rifting with combined mode-I/mode-II yield + melting | `plast=2`, 9 phases, VEP + melt weakening (Popov 2025) |

**When to use**: Start with `RiftingBasic` for a first lithospheric extension run. Use `RiftingChenin` for multi-layer setups. Add `RiftingAnisotropy` or `RiftingCheninAniso` to study fabric effects on rifting. Use `RiftingCombinedYield` for combined tensile+shear failure with melting.

---

## Subduction and Oceanic Cooling

| File | Description | Key features |
|------|------------|--------------|
| `AnneloreSubduction.c` | Plunging lithospheric slab | Crustal/mantle structure, thermal evolution |
| `DoubleSubduction.c` | Two-sided subduction | Dual oceanic plates, polar/Cartesian geometry |
| `SubductionMelting.c` | Subduction with fluid-released melting | Slab-edge rheological weakening |
| `SubductionBenchmarkCase2.c` | Benchmark subduction | Comparison with analytical solutions |
| `ThanushikaSubduction.c` | Simplified subduction | Solver performance testing |
| `OceanicCooling.c` | Plate cooling away from ridge | Thermal half-space model |
| `ThermalStructureSubduction.c` | Thermal tracking of slab | Mantle/crust temperature transitions |

**When to use**: `SubductionBenchmarkCase2` for validation. `AnneloreSubduction` for realistic slab geometry. `SubductionMelting` when fluid/melt coupling matters.

---

## Collision and Mountain Building

| File | Description | Key features |
|------|------------|--------------|
| `CollisionZone.c` | Symmetric continental collision | Two cratons + oceanic basin |
| `CollisionIra.c` | Asymmetric collision | West-central-east lithosphere configuration |
| `CollisionPolarCartesian.c` | Collision in polar coordinates | Curved Earth-like geometry |
| `CollisionPolarCartesianAniso.c` | Polar collision + anisotropy | Fabric tracking in curved geometry |
| `CollisionAnisotropy.c` | Shortening with crystal anisotropy | Mountain building + fabric evolution |

**When to use**: `CollisionZone` for basic orogenesis. Add `CollisionAnisotropy` for fabric-controlled deformation during convergence.

---

## Shear Zones and Simple Shear

| File | Description | Key features |
|------|------------|--------------|
| `ShearTemplate.c` | Simple shear with circular inclusion | Viscosity contrast testing |
| `ShearTemplateAniso.c` | Simple shear + anisotropic viscosity | Director evolution in shear |
| `ShearBoxLaeti.c` | Shear box with pre-stress | Weak inclusion, loading conditions |
| `SimpleShearAnisoHomo.c` | Homogeneous simple shear | Anisotropic deformation law testing |
| `ShearInclusionAniso.c` | Weak ellipse in shear flow | Anisotropy development around inclusion |
| `ShearHeatingDuretz14.c` | Shear heating benchmark | Viscous dissipation (Duretz et al. 2014) |
| `CrustalShearBandingAnisotropy.c` | Crustal shear zone with mica layers | Elliptical weak inclusions promoting strain bands |
| `StrainLocalization_SH_GSE.c` | Shear heating + grain size evolution | Coupled localisation mechanisms |

**When to use**: `ShearTemplate` as a starting point for shear experiments. `ShearHeatingDuretz14` for thermal-mechanical coupling. `CrustalShearBandingAnisotropy` for crustal fabric studies.

---

## Anisotropy Testing

| File | Description | Key features |
|------|------------|--------------|
| `AnisotropyDabrowski.c` | Analytical benchmark: flow around inclusion | Dani Schmid solution, pure/simple shear |
| `AnisotropyLayer.c` | Layered anisotropic material | Fabric tracking verification |
| `AnisoHomoVEVP.c` | Homogeneous VEVP + anisotropy | Full rheology coupling test |
| `AnisoVEVP.c` | VEVP deformation + anisotropy | Parallel variant |
| `AnisoViscTest_Ellipses.c` | 299 elliptical weak inclusions | Stress/strain field testing |
| `AnisoViscTest_HomoRandomAngle.c` | Random director orientation field | Statistical anisotropy |
| `AnisoViscTest_evolv.c` | Anisotropy evolution in pure shear | Time-dependent fabric growth |
| `AnisoViscTest_MWE.c` | Minimal working example (46 ellipses) | Quick debugging/verification |
| `AnisoViscTest_Debug.c` | Debug configuration | Solver diagnostics |
| `TestAnisotropy.c` | General test suite (269 elements) | Comprehensive aniso testing |
| `WedgeAnisotropy.c` | Subduction wedge + anisotropy | Mica-dominated fabric |

**When to use**: `AnisotropyDabrowski` for code validation against analytics. `AnisoViscTest_MWE` for fast debugging. `AnisoViscTest_evolv` for fabric evolution studies. `WedgeAnisotropy` for subduction-related fabric.

---

## Necking and Strain Localisation

| File | Description | Key features |
|------|------------|--------------|
| `NeckingRoman.c` | Lithospheric necking | Perturbation-driven localisation |
| `NeckingReview.c` | Systematic necking study | High-resolution strain rate tracking |
| `PinchSwellGSE.c` | Pinch-and-swell instability | Grain size evolution |
| `FaultSmearing.c` | Fault initiation and smearing | Plastic strain accumulation |
| `Shrinking.c` | Shrinking viscous inclusion | Negative buoyancy test |
| `Shrinking_New.c` | Updated shrinking inclusion | Revised parameters |

**When to use**: `NeckingRoman` for rift necking studies. `PinchSwellGSE` when grain size evolution controls folding.

---

## Viscous and Multiphase Flow

| File | Description | Key features |
|------|------------|--------------|
| `ViscousInclusion.c` | Weak circular inclusion in matrix | Analytical benchmark |
| `ViscousInclusionTestSuite.c` | Extended inclusion tests | Multiple configurations |
| `ViscousSinker.c` | Buoyant sphere / Stokes sinker | Convection benchmark |
| `MLPS_Ellipses.c` | Multilayer passive sink ellipses | Multilayer particle structure |
| `Multilayers.c` | Passive multilayer contrasts | No anisotropy |
| `Anais2Shrink.c` | Layered anisotropic + shrinking | Combined test |
| `Anais2Shrink_preload.c` | Pre-loaded aniso shrinking | Pre-stress variant |
| `AnaisReference.c` | Reference anisotropic layers | Baseline comparison |

**When to use**: `ViscousInclusion` and `ViscousSinker` for fundamental Stokes solver validation.

---

## Thermal

| File | Description | Key features |
|------|------------|--------------|
| `ThermalDiffusion.c` | Pure heat diffusion | Homogeneous medium |
| `ThermoElastic.c` | Thermo-elastic coupling | Temperature-dependent stiffness |
| `CoolingChamber.c` | Lithospheric cooling chamber | Magma chamber analogue |
| `TCMagmaticSystem.c` | Thermo-chemical magmatic evolution | Multi-process, `plast=2`, melting, dike injection |
| `MeltingOverpressure.c` | Melting-induced overpressure | Magma chamber pressurisation |
| `BlankenBench.c` | Blankenbach Case 1a convection | Ra=10⁴, isoviscous, free-slip, thermal benchmark |

**When to use**: `ThermalDiffusion` to verify thermal solver. `BlankenBench` for thermal convection benchmark (Boussinesq coupling). `TCMagmaticSystem` for coupled magmatic studies.

---

## Magmatic and Chemical

| File | Description | Key features |
|------|------------|--------------|
| `PressurizedMagmaChamber.c` | Pressurised circular chamber | Deformation around intrusion, `plast=2`, `T_st` |
| `ChristmasTree.c` | Complex dike geometry | Multi-branch intrusion |
| `ChemicalDiffusion.c` | Chemical element transport | Diffusion-driven |
| `ChemicalProduction.c` | In-situ phase production | Metamorphic reaction tracking |

**When to use**: `PressurizedMagmaChamber` for studying deformation driven by magmatic overpressure.

---

## Phase Transitions and Metamorphism

| File | Description | Key features |
|------|------------|--------------|
| `QuartzCoesite.c` | Quartz → coesite transition | Phase diagram coupling |
| `QuartzCoesite_Simple.c` | Simplified variant | Reduced complexity |
| `QuartzCoesite_paper.c` | Paper variant | Resolutions: 100, 150, 200 |
| `CrystallisationStuewe1995.c` | Crystallisation kinetics | Stüwe (1995) framework |

**When to use**: `QuartzCoesite` to test phase transition coupling. `CrystallisationStuewe1995` for metamorphic reaction kinetics.

---

## Topography and Flexure

| File | Description | Key features |
|------|------------|--------------|
| `FreeSurfaceBenchmarkFEM.c` | Free surface evolution benchmark | Initial perturbation decay |
| `TopoBenchCase1.c` | Topographic benchmark | Cosine perturbation |
| `ElasticFlexureLineLoad.c` | Elastic plate flexure | Line load on plate |
| `StressBC_FreeSurf.c` | Free surface + stress BCs | Combined boundary types |

**When to use**: `TopoBenchCase1` for free surface validation. `ElasticFlexureLineLoad` for flexure studies.

---

## Benchmarks and VEP Tests

| File | Description | Key features |
|------|------------|--------------|
| `VEP_Duretz18.c` | Visco-elasto-plastic benchmark | Duretz et al. (2018) |
| `Popov2025_Pureshear_VEVP.c` | Pure shear VEVP test | Combined mode-I/mode-II yield benchmark (`plast=2`) |
| `Popov2025_Tensile_VEVP.c` | Tensile/necking VEVP test | Tensile failure with `plast=2`, `T_st` |
| `PlasticityLayers.c` | Plasticity in layered systems | Layer-controlled localisation |
| `StressBC_ShearTemplate.c` | Shear with stress BCs | Non-velocity-driven loading |

**When to use**: `VEP_Duretz18` for standard VEP benchmarking. `Popov2025_Pureshear_VEVP` and `Popov2025_Tensile_VEVP` for testing/validating the combined mode-I/mode-II yield function (`plast=2`).

---

## Regional Studies

| File | Description | Key features |
|------|------------|--------------|
| `RiverTom.c` | Large-scale lithospheric model | Regional study |
| `RiverTomLowerCrust.c` | Lower crustal rheology variant | Focus on deep crust |

---

## Quick-Start Recommendations

| Goal | Start with |
|------|-----------|
| First simulation ever | `ViscousInclusion.c` (simple, has analytical solution) |
| Basic rifting | `RiftingBasic.c` |
| Anisotropy introduction | `AnisoViscTest_MWE.c` (fast) or `AnisotropyDabrowski.c` (validated) |
| Full-featured rifting | `RiftingCheninAniso.c` |
| Subduction | `SubductionBenchmarkCase2.c` |
| Solver verification | `ViscousSinker.c` or `TopoBenchCase1.c` |
| Thermal solver only | `ThermalDiffusion.c` |
| Thermal convection benchmark | `BlankenBench.c` |
| VEP rheology testing | `VEP_Duretz18.c` |
| Combined yield (plast=2) | `Popov2025_Pureshear_VEVP.c` or `Popov2025_Tensile_VEVP.c` |
| Rifting + melting + combined yield | `RiftingCombinedYield.c` |
| Melting + melt weakening | `RiftingMelting.c` |
