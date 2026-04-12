## MODIFIED Requirements

### Requirement: Multi-layer lithosphere with elliptical heterogeneities
The scenario SHALL define a 9-phase lithospheric structure adapted from MLPS_Ellipses:
- Phase 0: Continental crust (anorthite flow law, `pwlv=19`) — isotropic
- Phase 1: Lithospheric mantle (dry olivine, `pwlv=40`, `linv=40`, `expv=40`) — anisotropic, GSE enabled
- Phase 2: Lower mantle (wet olivine, `pwlv=40`, `linv=40`, `expv=40`) — anisotropic
- Phase 3: Asthenosphere (same rheology as phase 2) — anisotropic
- Phase 4: Weak crustal ellipses (wet quartzite, `pwlv=13`) — isotropic
- Phase 5: Weak mantle ellipses (wet olivine, `pwlv=13`) — anisotropic
- Phase 6: Strong crustal ellipses (Maryland diabase, `pwlv=11`) — isotropic
- Phase 7: Sediments 1 — isotropic
- Phase 8: Sediments 2 — isotropic

The `SetPhase` callback SHALL assign phases based on depth layers and embedded ellipses matching the MLPS_Ellipses geometry (crust 30 km, lithospheric mantle 40 km, with elliptical bodies of varying size).

The `RiftingComprehensive` scenario SHALL extend this phase structure with anisotropy parameters on mantle phases (1, 2, 3, 5) and grain-size evolution on the dry olivine mantle phase (1). All other material properties SHALL remain identical to `RiftingCombinedYield`.

#### Scenario: Phase assignment matches MLPS_Ellipses geometry
- **WHEN** the model initialises particles
- **THEN** the phase field reproduces the MLPS_Ellipses multi-layer structure with embedded ellipses at the same positions and sizes

#### Scenario: Anisotropy parameters on mantle phases
- **WHEN** the `RiftingComprehensive` `.txt` file is parsed
- **THEN** phases 1, 2, 3, 5 include anisotropy material parameters and phases 0, 4, 6, 7, 8 do not

#### Scenario: GSE on dry olivine
- **WHEN** the `RiftingComprehensive` `.txt` file is parsed
- **THEN** phase 1 has `gsel = 1` with valid grain-size parameters

### Requirement: Visco-elastic-plastic rheology with thermal coupling
The scenario SHALL enable elasticity (`elastic=1`), thermal solving (`thermal=1`), and shear heating (`shear_heating=1`). The `RiftingComprehensive` variant SHALL additionally enable anisotropy (`aniso=1`, `aniso_fstrain=1`) globally while `RiftingCombinedYield` SHALL NOT have anisotropy enabled.

#### Scenario: RiftingCombinedYield rheology switches
- **WHEN** the `RiftingCombinedYield` `.txt` file is parsed
- **THEN** `elastic=1`, `thermal=1`, `shear_heating=1` are set and no anisotropy parameters are present

#### Scenario: RiftingComprehensive adds anisotropy
- **WHEN** the `RiftingComprehensive` `.txt` file is parsed
- **THEN** `elastic=1`, `thermal=1`, `shear_heating=1`, `aniso=1`, `aniso_fstrain=1` are all set
