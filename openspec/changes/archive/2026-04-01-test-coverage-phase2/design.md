## Context

MDOODZ 7.0 has 26 CI test cases across 6 suites from phase 1. These tests verify basic solver completion and field boundedness but lack quantitative validation of physical mechanisms. Fifteen coverage gaps were identified: free surface, marker advection, density/buoyancy, shear heating, visco-elastic stress memory, isolated power-law creep, velocity field validation, compressibility, finite strain, Neumann BCs, strain-rate partitioning, convergence rate comparison, and strengthening of 3 existing suites.

All HDF5 fields needed for assertions are confirmed available: Vx/Vz (VxNodes, VzNodes), rho_n, eII_pwl/eII_lin/eII_exp/eII_gbs/eII_el/eII_pl, divu/divu_el/divu_pl, Fxx/Fxz/Fzx/Fzz (requires `finite_strain=1`), Topo/z_grid (requires `free_surface=1`), LogIsNewtonStep, X (composition), and the Particles group (marker positions). The existing `TestHelpers.h` functions (`getMaxFieldValue`, `getMinFieldValue`, `getMeanFieldValue`, `getStepsCount`, `getFinalResidual`) are fully generic and accept any HDF5 group name.

The build/test environment is WSL2 Ubuntu 22.04 with native filesystem (not Windows-mounted), gcc 11.4, 8 cores, all dependencies pre-installed. Build+test cycle is ~10 seconds.

## Goals / Non-Goals

**Goals:**

- Cover every major physics module in MDLIB/ with at least one quantitative assertion (not just "did it run?")
- Verify cause→effect relationships: when a mechanism is enabled, its signature field changes in the expected direction
- Catch regressions that would silently produce wrong physics (e.g., zero plastic strain when yielding should occur)
- Document each test with sufficient scientific context for geodynamic researchers
- Keep total test execution time reasonable (target under 60s on 8-core WSL, but correctness takes priority over speed)
- Update skills (testing-guide, physics-theory, rheology, model-setup, parameter-validation) with lessons learned

**Non-Goals:**

- Manufactured solutions or analytical convergence studies (these require much finer grids and are research-grade benchmarks, not CI tests)
- Multi-physics coupled tests (e.g., thermal-mechanical feedback loops) — each test isolates one mechanism
- Performance benchmarks or profiling
- Testing the HDF5 I/O layer itself
- Visual regression tests (phase 1 handles those separately)

## Decisions

### D1: New test executables vs. adding to existing

**Decision:** Create 9 new GTest executables. Merge 3 closely related suites into existing files as additional test cases.

**Rationale:** The 12 new capabilities from the proposal group naturally:
- **New executables (9):** ViscoElasticTests, DensityTests, ShearHeatingTests, FreeSurfaceTests, VelocityFieldTests, CompressibilityTests, FiniteStrainTests, NeumannBCTests, ConvergenceRateTests
- **Added to existing (3):** PowerLawCreep and StrainRatePartitioning into RheologyCreepTests (same fixture, same domain); AdvectionTests is dropped as separate suite — marker advection is tested indirectly via composition field in a new CompositionalAdvection test inside FiniteStrainTests

**Why not one monolith?** CTest runs executables in parallel. Separate executables = faster CI on multi-core machines. Each executable owns its parameter files in a named subdirectory.

### D2: Marker advection testing strategy

**Decision:** Test advection via the grid-interpolated composition field `X` (Centers/X) and finite-strain tensor `F`, rather than reading raw marker positions from Particles group.

**Rationale:** 
- The Particles group contains marker positions but the array size varies (markers are reseeded dynamically), making comparison across steps fragile.
- The composition field `X` is the grid-interpolated view of marker phases — if advection is broken, X will smear or develop artifacts.
- The finite-strain tensor F directly records accumulated deformation from marker transport.
- This gives robust, grid-based assertions without needing to track individual markers.

### D3: Free surface test grid size

**Decision:** Use 41×41 grid (not 21×21) for free surface tests. Run 3 steps only.

**Rationale:** Free surface stability requires sufficient resolution for the surface stabilisation algorithm (`free_surface_stab=1`). A 21×21 grid is at the numerical stability boundary. 41×41 is the smallest safe resolution while keeping runtime under 5 seconds. Only 3 steps are needed to see initial topographic deflection.

### D4: Parameter file conventions

**Decision:** Follow phase-1 conventions exactly:
- Each test suite has a named subdirectory under `TESTS/` (e.g., `TESTS/ViscoElastic/`)
- One `.txt` parameter file per test case
- Domain sizes: 1×1 non-dimensional for mechanics tests, 10km×10km for rheology, 100km×100km for thermal/free-surface
- Grid: 21×21 or 31×31 for fast execution, 41×41 only when physics demands it
- Steps: 1–5 for most tests, 3 for free surface, 5 for stress memory
- Plasticity: always set `coh_soft = 0` explicitly unless testing softening
- `Slim` ≥ `C` always (phase-1 lesson)

### D5: Assertion tolerance strategy

**Decision:** Use relative-order assertions (>, <, factor-of-N) rather than exact numerical tolerances.

**Rationale:** MDOODZ is a non-linear solver — exact field values depend on grid resolution, iteration count, and convergence tolerance. Tests should verify:
- Direction: "η decreased when grain size decreased" → `EXPECT_LT(eta_small_d, eta_large_d)`
- Mechanism activation: "eII_pwl > 0 when pwlv is enabled" → `EXPECT_GT(maxEiiPwl, 0.0)`
- Sign/order: "Fxx > 1 under extension" → `EXPECT_GT(maxFxx, 1.0)`
- Bound: "P within 50% of ρgH" → `EXPECT_NEAR(maxP, rho*g*H, 0.5*rho*g*H)`

Exact floating-point comparisons are avoided. Phase-1 lesson: `EXPECT_LT(res, 1e-9)` flaked because Picard converges linearly. Use `1e-8` floor.

### D6: TestHelpers.h extensions

**Decision:** Add two new helper functions:
1. `getFieldValueAt(filename, group, dataset, index)` — read a single value at a known grid index (for boundary checks)
2. `getModelParam(filename, index)` — read from Model/Params array (for extracting dt, Lx, Lz, dx, dz)

**Rationale:** Velocity BC checks need values at specific boundary nodes. Model parameters (dt, domain size) are needed for quantitative assertions like ΔT ≈ 2ηε̇²Δt/(ρCp). The existing `getModelParam` is already present, but `getFieldValueAt` is new.

### D7: Skills updates

**Decision:** After implementation, update these skills with new lessons:
- `skill-testing-guide`: New test patterns (velocity field, free surface, composition), new TestHelpers functions, extended pitfalls
- `skill-physics-theory`: Shear heating formula, hydrostatic pressure, Maxwell stress accumulation equations as verified by tests
- `skill-rheology`: Power-law shear-thinning test pattern, strain-rate partitioning validation approach
- `skill-model-setup`: Free surface parameter requirements, Neumann BC setup pattern
- `skill-parameter-validation`: New switch compatibility rows (free_surface requirements, compressible+dilatant)

### D8: Grouping to limit total executables

**Decision:** Consolidate from 12 proposal capabilities to 9 new executables:

| Executable | Capabilities covered |
|---|---|
| ViscoElasticTests | ci-viscoelastic-tests |
| DensityTests | ci-density-tests |
| ShearHeatingTests | ci-shear-heating-tests |
| FreeSurfaceTests | ci-free-surface-tests |
| VelocityFieldTests | ci-velocity-field-tests |
| CompressibilityTests | ci-compressibility-tests |
| FiniteStrainTests | ci-finite-strain-tests + ci-advection-tests (via composition) |
| NeumannBCTests | ci-neumann-bc-tests |
| ConvergenceRateTests | ci-convergence-rate-tests |

Plus modifications to 3 existing executables:
- RheologyCreepTests ← ci-powerlaw-creep-tests + ci-strain-partitioning-tests
- BoundaryConditionTests ← ci-boundary-condition-tests modifications
- ThermalTests ← ci-thermal-tests modifications

## Risks / Trade-offs

**[Risk] Free surface test may be flaky on different platforms** → Mitigation: Use 41×41 grid and only 3 steps. Assert topographic deflection direction (z_grid changed), not exact magnitude. Add 10% tolerance on expected deflection.

**[Risk] Shear heating ΔT depends on viscosity, which is non-linear** → Mitigation: Use constant viscosity (cstv=1) for the shear heating test so ΔT = 2η₀·ε̇²·Δt/(ρCp) is exact. Only check order of magnitude: `EXPECT_GT(deltaT, 0.5 * expected)` and `EXPECT_LT(deltaT, 2.0 * expected)`.

**[Risk] Hydrostatic pressure test requires no imposed deformation** → Mitigation: Set `bkg_strain_rate = 0`, `pure_shear_ALE = 0`, `gz = -10`. Use density contrast only. Phase-1 lesson: no-slip + pure_shear_ALE = singular matrix.

**[Risk] Adding 9 executables increases CI build time** → Mitigation: Each test is 1 .cpp file compiled in parallel. On 8 cores, incremental compile is ~2s. Total test runtime stays under 60s because each case runs in <5s.

**[Risk] Neumann BC test may not produce clean results on coarse grid** → Mitigation: Use a 31×31 grid with constant viscosity and a simple applied traction. Check that boundary stress matches applied value within 20%.

**[Trade-off] Marker advection not tested directly** → Accepted: Composition field X and finite-strain tensor F give grid-based coverage of the advection+interpolation pipeline without fragile marker-index tracking.

**[Trade-off] No multi-physics coupling tests** → Accepted: Each test isolates one mechanism. Coupling bugs are harder to diagnose from a single failing test. Phase 3 could add coupled tests.
