## Why

MDOODZ phase-1 testing established 26 test cases across 6 suites, but coverage remains shallow — most assertions only check "did it run?" (iteration count, η > 0, stress < 1e40) rather than "did it compute correctly?". Critical physics modules have **zero test coverage**: free surface topography evolution, marker advection, shear heating, visco-elastic stress memory, density/thermal expansion, velocity field verification, compressibility, finite strain tracking, composition transport, and stress (Neumann) boundary conditions. Several existing tests also lack quantitative assertions on the very fields they are designed to exercise (e.g., strain-rate partitioning fields eII_pwl/eII_lin/eII_exp/eII_gbs are written to HDF5 but never read).

This leaves the codebase vulnerable to silent regressions during refactoring. The goal is to reach a coverage level where any future code change that breaks a physical mechanism will be caught by at least one failing assertion — giving researchers confidence that refactored code still satisfies the governing equations.

## What Changes

**New test suites** (each a new GTest executable + parameter files + README documentation):

1. **ViscoElasticTests** — Maxwell stress accumulation: verify σ grows over multiple elastic time steps, and that elastic strain rate eII_el is non-zero when `elastic=1`. Stress memory test (load → unload → check residual stress).
2. **PowerLawCreepTests** — Isolated dislocation creep (pwlv=40): verify eII_pwl > 0 and viscosity drops with increasing strain rate (n=3 shear-thinning). Also verify strain-rate partitioning: for a composite run, eII_pwl + eII_lin ≈ eII_total.
3. **DensityTests** — Thermal expansion: verify ρ(T_hot) < ρ(T_cold) for α > 0. Hydrostatic pressure: verify P ≈ ρgz for a static column with no imposed strain rate.
4. **ShearHeatingTests** — Verify temperature increases when `shear_heating=1` under high strain rate. Quantitative check: ΔT ≈ 2η·ε̇²·Δt / (ρCp).
5. **FreeSurfaceTests** — Verify topography evolves: a dense block sinks and surface deflects. Read Topo/z_grid from HDF5. Verify free-surface marker chain remains consistent.
6. **AdvectionTests** — Verify markers advect correctly: a rigid-body rotation should preserve marker positions (within RK error). Also test periodic wrapping: marker crossing x-boundary reappears on the other side.
7. **VelocityFieldTests** — Read Vx/Vz from HDF5. For pure shear with known ε̇, verify Vx ≈ ε̇·x at boundaries. For simple shear, verify Vx ≈ γ̇·z. Ensure no-slip test actually has V=0 at boundaries.
8. **CompressibilityTests** — `compressible=1` with dilatant plasticity: verify div(v) ≠ 0. Compare against `compressible=0` where div(v) ≈ 0.
9. **FiniteStrainTests** — Run 5 steps with `finite_strain=1`, read Fxx/Fxz/Fzx/Fzz. Verify F departs from identity under deformation. For pure shear, Fxx > 1 and Fzz < 1.
10. **NeumannBCTests** — Stress boundary conditions (type 2/13): apply known traction on one boundary, verify resulting stress field is consistent.
11. **StrainRatePartitioningTests** — For composite rheology runs, read eII_pwl, eII_lin, eII_exp, eII_gbs, eII_el, eII_pl from HDF5. Verify they sum to roughly eII_total and that each mechanism contributes when enabled.
12. **ConvergenceRateTests** — Compare Newton (quadratic) vs Picard (linear) convergence: run same problem with both, verify Newton uses fewer iterations. Read LogIsNewtonStep to validate.

**Strengthened assertions in existing tests** (via modified specs):

13. Existing `ci-rheology-tests`: add eII_lin > 0 assertion for diffusion creep, verify grain-size effect on eII quantitatively.
14. Existing `ci-boundary-condition-tests`: add Vx/Vz field checks (pure shear Vx≈ε̇·x, no-slip V≈0 at boundary).
15. Existing `ci-thermal-tests`: add steady-state linearity check (max gradient ≈ ΔT/Lz), add heat flux conservation.

**Documentation**: Every new test will be documented in TESTS/README.md with the same scientific rigour as existing entries (governing equations, cause→effect, parameter tables, code assertions).

**Validation**: All tests must be built and verified in WSL (native Linux filesystem for speed). The apply phase must not terminate until `ctest --output-on-failure` passes 100%.

## Capabilities

### New Capabilities
- `ci-viscoelastic-tests`: Maxwell visco-elastic stress accumulation, elastic strain partitioning, stress memory
- `ci-powerlaw-creep-tests`: Isolated power-law dislocation creep, shear-thinning verification, strain-rate partitioning
- `ci-density-tests`: Thermal expansion ρ(T), hydrostatic pressure P=ρgz validation
- `ci-shear-heating-tests`: Viscous dissipation H_s=2ηε̇², temperature increase from mechanical work
- `ci-free-surface-tests`: Topography evolution, surface deflection under gravity, marker chain consistency
- `ci-advection-tests`: RK2/RK4 marker transport, periodic wrapping, position conservation
- `ci-velocity-field-tests`: Vx/Vz HDF5 field validation, BC enforcement verification
- `ci-compressibility-tests`: Dilatant flow, div(v)≠0, volumetric plastic strain
- `ci-finite-strain-tests`: Deformation gradient tensor F evolution, anisotropy fabric tracking
- `ci-neumann-bc-tests`: Stress (traction) boundary conditions, type 2/13 enforcement
- `ci-strain-partitioning-tests`: Mechanism-specific strain rate (eII_pwl, eII_lin, eII_exp, eII_gbs, eII_el, eII_pl) validation
- `ci-convergence-rate-tests`: Newton vs Picard iteration count comparison, LogIsNewtonStep verification

### Modified Capabilities
- `ci-rheology-tests`: Strengthen assertions — add eII_lin > 0, grain-size quantitative check
- `ci-boundary-condition-tests`: Add Vx/Vz field verification to existing BC tests
- `ci-thermal-tests`: Add steady-state linearity check, heat flux conservation

## Impact

- **Code**: New files in `TESTS/` — 12 new `.cpp` test files, ~40 new `.txt` parameter files, updates to `TESTS/CMakeLists.txt`
- **Existing tests**: 3 existing test files modified (additional assertions, no removal)
- **Documentation**: `TESTS/README.md` expanded with 12 new suite sections
- **CI**: 12 new CTest executables registered, total test count ~60+ cases
- **Build time**: ~15-20 additional seconds for compilation, ~30-60 seconds for test execution
- **Dependencies**: No new external dependencies (same GTest + HDF5 + SuiteSparse stack)
- **Skills**: `skill-testing-guide` updated with new test patterns
