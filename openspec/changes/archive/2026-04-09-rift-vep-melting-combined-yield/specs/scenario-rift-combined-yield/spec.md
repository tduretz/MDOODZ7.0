## ADDED Requirements

### Requirement: Scenario files exist and build
The scenario SHALL consist of `SETS/RiftingCombinedYield.c` and `SETS/RiftingCombinedYield.txt`. The `.c` file SHALL be registered in `SETS/AddSet.cmake` so that CMake generates a build target. The executable SHALL compile without errors and link against the `mdoodz` library.

#### Scenario: CMake target builds successfully
- **WHEN** the user runs `make build` (or the CMake build)
- **THEN** the `RiftingCombinedYield` executable is produced in `cmake-exec/RiftingCombinedYield/` alongside a copy of `RiftingCombinedYield.txt`

### Requirement: Multi-layer lithosphere with elliptical heterogeneities
The scenario SHALL define a 9-phase lithospheric structure adapted from MLPS_Ellipses:
- Phase 0: Continental crust (anorthite flow law, `pwlv=19`)
- Phase 1: Lithospheric mantle (dry olivine, `pwlv=40`, `linv=40`, `expv=40`)
- Phase 2: Lower mantle (wet olivine, `pwlv=40`, `linv=40`, `expv=40`)
- Phase 3: Asthenosphere (same rheology as phase 2)
- Phase 4: Weak crustal ellipses (wet quartzite, `pwlv=13`)
- Phase 5: Weak mantle ellipses (wet olivine, `pwlv=13`)
- Phase 6: Strong crustal ellipses (Maryland diabase, `pwlv=11`)
- Phase 7: Sediments 1
- Phase 8: Sediments 2

The `SetPhase` callback SHALL assign phases based on depth layers and embedded ellipses matching the MLPS_Ellipses geometry (crust 30 km, lithospheric mantle 40 km, with elliptical bodies of varying size).

#### Scenario: Phase assignment matches MLPS_Ellipses geometry
- **WHEN** the model initialises particles
- **THEN** the phase field reproduces the MLPS_Ellipses multi-layer structure with embedded ellipses at the same positions and sizes

### Requirement: Combined mode-I/mode-II yield function enabled
All phases SHALL use `plast=2` (Popov et al. 2025 combined yield criterion) with a uniform tensile strength `T_st = -10e6` Pa and visco-plastic viscosity `eta_vp = 2.5e20` Pa.s.

#### Scenario: Combined yield parameters present in .txt
- **WHEN** the `.txt` file is parsed
- **THEN** every phase block contains `plast = 2`, `T_st = -10e6`, and `eta_vp = 2.5e20`

### Requirement: Partial melting enabled per phase
The scenario SHALL enable partial melting globally (`melting=1`) with per-phase melt models:
- Crustal phases (0, 4, 6, 7, 8): `melt=11` (Schenker et al. 2012 crustal solidus/liquidus)
- Mantle phases (1, 2, 3, 5): `melt=41` (mantle peridotite solidus/liquidus)

Melt weakening SHALL be active: `melt_weak=30.0`, `force_melt_weak=1`.

#### Scenario: Melting configuration in .txt
- **WHEN** the `.txt` file is parsed
- **THEN** `melting=1` is set globally, each crustal phase has `melt=11`, each mantle phase has `melt=41`, and `melt_weak=30.0` is set

### Requirement: Visco-elastic-plastic rheology with thermal coupling
The scenario SHALL enable elasticity (`elastic=1`), thermal solving (`thermal=1`), and shear heating (`shear_heating=1`). Anisotropy SHALL NOT be enabled.

#### Scenario: Rheology switches in .txt
- **WHEN** the `.txt` file is parsed
- **THEN** `elastic=1`, `thermal=1`, `shear_heating=1` are set and no anisotropy parameters are present

### Requirement: Pure shear extension boundary conditions
The scenario SHALL apply pure shear extension (`pure_shear_ALE=1`) at a background strain rate of `bkg_strain_rate = -1.056e-15` (1 cm/yr full extension rate). Velocity BCs SHALL use the built-in `SetPureShearBCVx` and `SetPureShearBCVz` functions.

#### Scenario: Extension BCs applied
- **WHEN** the model runs
- **THEN** lateral and bottom boundaries apply velocities consistent with pure shear at the specified strain rate

### Requirement: Thermal boundary conditions
The `.c` file SHALL define a `SetBCT` callback with:
- South boundary: constant temperature 1330°C + 273.15 K
- North / free surface: constant temperature 0°C + 273.15 K
- West / East: zero heat flux

#### Scenario: Thermal BCs match continental geotherm
- **WHEN** the thermal solver evaluates boundary conditions
- **THEN** bottom is at mantle temperature (1330°C), top is at surface temperature (0°C), and lateral boundaries are insulating

### Requirement: Solver configured for combined yield convergence
The `.txt` file SHALL use Picard-to-Newton iteration: `Newton=0`, `Picard2Newton=1`, `Picard2Newton_tol=1e-1`, `nit_max=20`, `line_search=1`.

#### Scenario: Solver parameters in .txt
- **WHEN** the `.txt` file is parsed
- **THEN** `Newton=0`, `Picard2Newton=1`, `nit_max=20`, `line_search=1` are set

### Requirement: Domain and resolution
The production `.txt` SHALL define the domain as:
- `xmin=-150e3`, `xmax=150e3` (300 km wide)
- `zmin=-110e3`, `zmax=10e3` (120 km deep)
- `Nx=200`, `Nz=160`

#### Scenario: Domain matches MLPS_Ellipses
- **WHEN** the `.txt` file is parsed
- **THEN** domain bounds and grid resolution match the specified values
