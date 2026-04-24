### Requirement: 0D stress integration tests for combined yield surface
The CI test suite SHALL include GTest cases that verify the local stress update of the combined tensile-cap + Drucker-Prager yield surface (`plast = 2`, Popov et al. 2025) for the three canonical 0D loading paths described in paper Fig. 5 / Table 1 "0D Fig. 5 a, b".

All three cases run a homogeneous 21×15 grid with custom pure-shear / outward-normal BCs, and use material parameters matching the paper exactly: `G = 10¹⁰ Pa`, `K = 2·10¹¹ Pa`, `φ = 30°`, `ψ = 10°`, `c_MC = 10⁶ Pa`, `p_T = −5·10⁵ Pa`, `η^vp = 0` (perfect plasticity), `Δt = 2 yr`.

Assertions are made directly on yield-surface membership at the final step (rather than full-trajectory L2) to remain robust to the one-step trial-pressure offset that MDOODZ exposes through the HDF5 global `/Centers/P` field (see paper Eq. 31).

#### Scenario: Volumetric extension saturates on the tensile cap
- **WHEN** the test runs `Popov2025/Popov0D_VolumetricExtension.txt` with `ε̇_xx = ε̇_zz = 2.333·10⁻¹⁵ s⁻¹` (outward normal BCs `user2 = user3 = 1.167·10⁻¹⁵ m/s`), `plast = 2`, `η^vp = 0`
- **THEN** the centre-cell `/Centers/divu_pl` SHALL be strictly positive (dilatant flow active on the cap), the reconstructed `sII` SHALL be at most 15 % of `τ_d`, and the trajectory point `(sII, P)` SHALL lie within 30 % of the tensile cap circle radius (`|R̂_y − R_y| / R_y < 0.30`)

#### Scenario: Deviatoric shear saturates on the Drucker-Prager envelope
- **WHEN** the test runs `Popov2025/Popov0D_DeviatoricShear.txt` with `bkg_strain_rate = 7·10⁻¹⁴ s⁻¹` (standard pure-shear BCs), `plast = 2`, `η^vp = 0`
- **THEN** the centre-cell `/Centers/eII_pl` SHALL be strictly positive (plastic flow active) and the reconstructed `sII` SHALL satisfy `|sII − (k·P + c)| / (k·P + c) < 0.05` (DP envelope membership within 5 %), with `|P| < 10·c_MC` (bounded dilatant-coupling growth)

#### Scenario: Mixed loading final state lies on the combined yield surface
- **WHEN** the test runs `Popov2025/Popov0D_MixedStrain.txt` with asymmetric outward normal BCs (`user2 = −3.325·10⁻¹⁴`, `user3 = +2.573·10⁻¹⁴`) that produce paper's `trace(ε̇) = 7·10⁻¹⁵` and `ε̇_II_dev = 7·10⁻¹⁴` simultaneously, `plast = 2`, `η^vp = 0`
- **THEN** at least one of `/Centers/eII_pl` or `/Centers/divu_pl` SHALL be non-zero at the final step, AND the smaller of the cap-circle residual (`|R̂_y − R_y|`) and the DP-envelope residual (`|sII − k·P − c|`) SHALL be less than 30 % of `τ_d`

### Requirement: Corrected SETS/Popov2025 demo scenarios
The two shipped demo scenarios `SETS/Popov2025_Pureshear_VEVP.{c,txt}` and `SETS/Popov2025_Tensile_VEVP.{c,txt}` SHALL be rewritten so that running them actually exercises the combined yield surface of Popov et al. (2025), with parameters matching Table 1 column "Regularization" of the paper.

#### Scenario: Both phases use plast = 2
- **WHEN** either `Popov2025_Pureshear_VEVP.txt` or `Popov2025_Tensile_VEVP.txt` is read by MDOODZ
- **THEN** every phase block SHALL set `plast = 2`, the global block SHALL set `tensile_line_search = 1`, and the obsolete `tensile = 0/1` key SHALL not appear

#### Scenario: Setup parameters align with paper Table 1
- **WHEN** either file is read by MDOODZ
- **THEN** the material properties (φ, ψ, `c_MC_init`, `c_MC_min`, Hc, `T_st`, `η^vp`, G) and domain size SHALL match the "Regularization" column of paper Table 1, and a comment block at the top of each `.txt` file SHALL reference the paper section and figure number

### Requirement: Popov 2025 test documentation
The test documentation SHALL describe the new 0D suite in both `TESTS/README.md` (operational description, how to run) and `TESTS/AnalyticalSolutions.md` (reference solutions).

#### Scenario: README.md includes a new suite entry
- **WHEN** a developer reads `TESTS/README.md`
- **THEN** a section "Suite 20 — Popov2025 (Combined Tensile Cap / DP Yield)" SHALL appear in the table of contents and body, listing each of the three test cases, their wall-time budget, and the `.txt` fixture files used

#### Scenario: AnalyticalSolutions.md includes the reference derivation
- **WHEN** a developer reads `TESTS/AnalyticalSolutions.md`
- **THEN** a section "Popov et al. (2025) Combined Yield Surface" SHALL document: (a) the yield-surface geometry (paper Eqs. 13–17), (b) the embedded reference plot rendered from the CI test output, (c) the test-case-to-assertion mapping, and (d) an explicit list of the 2D-localisation tests that were deferred from this change

### Requirement: Paper Fig. 5 reproduction script
The repository SHALL ship a gnuplot post-processing script that reproduces paper Fig. 5 (0D stress integration) from the HDF5 output of the CI tests.

The script lives under `TESTS/Popov2025/plots/` alongside an HDF5 → ASCII extractor `extract_popov2025.cpp` (modeled on `TESTS/BlankenBench/extract_fields.cpp`) registered as a standalone CMake target. The script is a developer aid; it is NOT invoked by CI and NOT wired into the `VISUAL_TESTS/` visual-regression pipeline.

#### Scenario: fig5 script runs end-to-end against test HDF5 output
- **WHEN** a developer builds the `extract_popov2025` and `Popov2025Tests` targets, runs the test binary, invokes `./extract_popov2025` on each of the three 0D output directories, and then runs `gnuplot fig5_0d_stress.gp` from `TESTS/Popov2025/plots/`
- **THEN** a PNG SHALL land at `TESTS/Popov2025/plots/out/fig5_0d_stress.png` showing the four panels `(a) volumetric extension`, `(b) deviatoric shear`, `(c) mixed strain`, `(d) meridional P–τ_II trajectories`, with pressure reconstructed as `p_local = p* + K·θ̇_vp·dt` (paper Eq. 31) so that saturation values sit on the yield surface

### Requirement: CMake target registration
The new test binary SHALL be wired into the project build system and run under `make run-tests`.

#### Scenario: Popov2025Tests target exists and is part of run-tests
- **WHEN** the project is configured with `cmake -B cmake-build -DTEST=ON` and built
- **THEN** a `Popov2025Tests` binary SHALL be produced under `cmake-build/TESTS/`, linked against `mdoodz` and `GTest::gtest_main`, registered with `add_test(NAME Popov2025Tests COMMAND Popov2025Tests)` and `LABELS "popov2025"`, and invoked by `make run-tests`
