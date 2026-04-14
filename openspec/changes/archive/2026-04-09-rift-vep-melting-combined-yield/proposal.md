## Why

No existing rifting scenario combines the new Popov et al. 2025 combined mode-I/mode-II yield function (`plast=2`) with partial melting in a multi-layer lithosphere with embedded heterogeneities. The combined yield criterion (commit `8aecb92`) enables tensile failure transitioning to frictional shear — critical for realistic fault nucleation during rifting — but is currently only demonstrated in small-scale benchmarks (Popov2025_Pureshear/Tensile_VEVP) and magmatic chamber models (PressurizedMagmaChamber, TCMagmaticSystem). A lithospheric-scale rifting scenario is needed to exercise these new capabilities together.

## What Changes

- **New scenario files**: `SETS/RiftingCombinedYield.c` and `SETS/RiftingCombinedYield.txt` — a rifting model based on the MLPS_Ellipses setup
- **Rheology**: Visco-Elastic-Plastic with the Popov combined mode-I/mode-II yield function (`plast=2`, per-phase tensile strength `T_st`)
- **Melting**: Partial melting enabled for crustal and mantle phases (melt models 10/11)
- **Structure**: Multi-layer lithosphere (crust, lithospheric mantle, asthenosphere) with embedded elliptical heterogeneities of varying rheology — similar to MLPS_Ellipses but adapted for the new yield function
- **Thermal**: Full thermal coupling with shear heating and radiogenic heat production
- **No anisotropy**: This model intentionally omits anisotropy to isolate the combined yield + melting interaction
- **Verification mode**: A low-resolution configuration (coarse grid, fewer time steps) for quick smoke-test validation that the model runs without diverging, before committing to high-resolution production runs

### Comparison with existing scenarios

| Feature | MLPS_Ellipses | RiftingMelting | RiftingBasic | **New scenario** |
|---|---|---|---|---|
| Elasticity | ✓ | ✓ | ✓ | ✓ |
| Thermal | ✓ | ✓ | ✓ | ✓ |
| Melting | ✗ | ✓ | ✗ | **✓** |
| Multi-layer + ellipses | ✓ (7 phases) | ✗ | ✓ (9 phases) | **✓** |
| Combined yield (plast=2) | ✗ | ✗ | ✗ | **✓** |
| Tensile strength (T_st) | ✗ | ✗ | ✗ | **✓** |
| Anisotropy | ✗ | ✗ | ✗ | ✗ |
| Quick-verify config | ✗ | ✗ | ✗ | **✓** |

### What existing scenarios lack

- **MLPS_Ellipses**: Has the multi-layer ellipse structure we want but uses only standard Drucker-Prager (`plast=1` implied), no melting, no tensile failure
- **RiftingMelting**: Has melting but simple 4-phase layered setup, no elliptical heterogeneities, no combined yield
- **Popov2025_* benchmarks**: Demonstrate the new yield function but are small-scale mechanical tests, not lithospheric rifting
- **PressurizedMagmaChamber / TCMagmaticSystem**: Use new yield + melting/dike injection but model crustal-scale magmatic processes, not continental rifting

## Capabilities

### New Capabilities
- `scenario-rift-combined-yield`: The new rifting scenario (.c/.txt pair) with multi-layer ellipse geometry, combined mode-I/mode-II plasticity, partial melting, and thermal coupling
- `scenario-quick-verify`: A low-resolution verification configuration and workflow for fast smoke-testing that a scenario runs without diverging before high-res runs

### Modified Capabilities
_(none — this change adds new scenario files without modifying existing capabilities)_

## Impact

- **New files**: `SETS/RiftingCombinedYield.c`, `SETS/RiftingCombinedYield.txt`, `SETS/RiftingCombinedYield_quick.txt` (10 steps, coarse), `SETS/RiftingCombinedYield_medium.txt` (100 steps, coarse)
- **CMake**: `SETS/AddSet.cmake` updated to register the new scenario
- **No library changes**: All required code features (`plast=2`, `T_st`, melting models) already exist in MDLIB from commit `8aecb92`
- **Dependencies**: Same as existing scenarios (SuiteSparse, HDF5, BLAS/LAPACK)
