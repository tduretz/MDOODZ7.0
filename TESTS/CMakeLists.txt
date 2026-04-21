if (VALGRIND)
    find_program( MEMORYCHECK_COMMAND valgrind )
    set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
    message("Valgrind is ON for testing")
endif()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ShearTemplate DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(NewtownIterationConvergence_test NewtonIterationConvergence.cpp)
target_link_libraries(NewtownIterationConvergence_test LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(NewtownIterationConvergence NewtownIterationConvergence_test)

# Rheology Creep Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/RheologyCreep DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(RheologyCreepTests RheologyCreepTests.cpp)
target_link_libraries(RheologyCreepTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(RheologyCreepTests RheologyCreepTests)

# Plasticity Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/Plasticity DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(PlasticityTests PlasticityTests.cpp)
target_link_libraries(PlasticityTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(PlasticityTests PlasticityTests)

# Thermal Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/Thermal DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(ThermalTests ThermalTests.cpp)
target_link_libraries(ThermalTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(ThermalTests ThermalTests)

# Thermal Verification Tests (oceanic cooling, thermo-elastic coupling)
add_executable(ThermalVerificationTests ThermalVerificationTests.cpp)
target_link_libraries(ThermalVerificationTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(ThermalVerificationTests ThermalVerificationTests)

# Boundary Condition Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/BoundaryCondition DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(BoundaryConditionTests BoundaryConditionTests.cpp)
target_link_libraries(BoundaryConditionTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(BoundaryConditionTests BoundaryConditionTests)

# Solver Mode Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/SolverMode DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(SolverModeTests SolverModeTests.cpp)
target_link_libraries(SolverModeTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(SolverModeTests SolverModeTests)

# Visco-Elastic Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ViscoElastic DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(ViscoElasticTests ViscoElasticTests.cpp)
target_link_libraries(ViscoElasticTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(ViscoElasticTests ViscoElasticTests)

# Density Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/Density DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(DensityTests DensityTests.cpp)
target_link_libraries(DensityTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(DensityTests DensityTests)

# Shear Heating Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ShearHeating DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(ShearHeatingTests ShearHeatingTests.cpp)
target_link_libraries(ShearHeatingTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(ShearHeatingTests ShearHeatingTests)

# Free Surface Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/FreeSurface DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(FreeSurfaceTests FreeSurfaceTests.cpp)
target_link_libraries(FreeSurfaceTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(FreeSurfaceTests FreeSurfaceTests)

# Velocity Field Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/VelocityField DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(VelocityFieldTests VelocityFieldTests.cpp)
target_link_libraries(VelocityFieldTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(VelocityFieldTests VelocityFieldTests)

# Compressibility Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/Compressibility DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(CompressibilityTests CompressibilityTests.cpp)
target_link_libraries(CompressibilityTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(CompressibilityTests CompressibilityTests)

# Interp Equivalence Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/InterpEquiv DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(InterpEquivalenceTests InterpEquivalenceTests.cpp)
target_link_libraries(InterpEquivalenceTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(InterpEquivalenceTests InterpEquivalenceTests)

# Finite Strain Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/FiniteStrain DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(FiniteStrainTests FiniteStrainTests.cpp)
target_link_libraries(FiniteStrainTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(FiniteStrainTests FiniteStrainTests)

# Neumann BC Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/NeumannBC DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(NeumannBCTests NeumannBCTests.cpp)
target_link_libraries(NeumannBCTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(NeumannBCTests NeumannBCTests)

# Convergence Rate Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ConvergenceRate DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(ConvergenceRateTests ConvergenceRateTests.cpp)
target_link_libraries(ConvergenceRateTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(ConvergenceRateTests ConvergenceRateTests)

# SolVi Benchmark Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/SolViBenchmark DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(SolViBenchmarkTests SolViBenchmarkTests.cpp)
target_link_libraries(SolViBenchmarkTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(SolViBenchmarkTests SolViBenchmarkTests)

# Anisotropy Benchmark Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/AnisotropyBenchmark DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(AnisotropyBenchmarkTests AnisotropyBenchmarkTests.cpp)
target_link_libraries(AnisotropyBenchmarkTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(AnisotropyBenchmarkTests AnisotropyBenchmarkTests)

# TopoBench Tests
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/TopoBench DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(TopoBenchTests TopoBenchTests.cpp)
target_link_libraries(TopoBenchTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(TopoBenchTests TopoBenchTests)

# Blankenbach Thermal Convection Tests (manual only — too slow for CI, ~30 min at 41×41)
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/BlankenBench DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(BlankenBenchTests BlankenBenchTests.cpp)
target_link_libraries(BlankenBenchTests LINK_PUBLIC mdoodz GTest::gtest_main)
if (OMP AND TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(BlankenBenchTests LINK_PUBLIC OpenMP::OpenMP_CXX)
endif()
add_test(NAME BlankenBenchTests COMMAND BlankenBenchTests)
set_tests_properties(BlankenBenchTests PROPERTIES LABELS "manual" TIMEOUT 3600)
set_property(TEST BlankenBenchTests PROPERTY DISABLED TRUE)

# Blankenbach HDF5-to-gnuplot extraction utility (no GTest, just HDF5)
add_executable(extract_blankenbach BlankenBench/extract_fields.cpp)
target_link_libraries(extract_blankenbach LINK_PUBLIC ${HDF5_C_LIBRARIES})
target_include_directories(extract_blankenbach PRIVATE ${HDF5_INCLUDE_DIRS})

# SolVi Marker Comparison (experimental — not in CI, run manually)
add_executable(SolViMarkerComparison SolViMarkerComparison.cpp)
target_link_libraries(SolViMarkerComparison LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(NAME SolViMarkerComparison COMMAND SolViMarkerComparison)
set_tests_properties(SolViMarkerComparison PROPERTIES LABELS "experimental" TIMEOUT 300)

# Rotation Advection Test (reseed mode comparison)
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/RotationAdvection DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
add_executable(RotationAdvectionTests RotationAdvectionTests.cpp)
target_link_libraries(RotationAdvectionTests LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(NAME RotationAdvectionTests COMMAND RotationAdvectionTests)
set_tests_properties(RotationAdvectionTests PROPERTIES LABELS "experimental" TIMEOUT 600)

# Geometric multigrid unit tests (add-gmg-stokes-solver change)
add_executable(MultigridTests MultigridTests.cpp)
target_link_libraries(MultigridTests LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(MultigridTests PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB)
add_test(NAME MultigridTests COMMAND MultigridTests)

# GMG Stokes solver (matvec, Vanka, V-cycle, FGMRES) unit tests
add_executable(MultigridStokesTests MultigridStokesTests.cpp)
target_link_libraries(MultigridStokesTests LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(MultigridStokesTests PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB)
add_test(NAME MultigridStokesTests COMMAND MultigridStokesTests)
set_tests_properties(MultigridStokesTests PROPERTIES TIMEOUT 300)

# GMG-vs-CHOLMOD Stokes solve equivalence tests (requires lin_solver=3 wiring
# with full MDOODZ stencil matching; marked experimental until that work lands).
add_executable(GmgStokesEquivalence GmgStokesEquivalence.cpp)
target_link_libraries(GmgStokesEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(GmgStokesEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB)
add_test(NAME GmgStokesEquivalence COMMAND GmgStokesEquivalence)
set_tests_properties(GmgStokesEquivalence PROPERTIES LABELS "experimental" TIMEOUT 1800)

# Stokes matvec golden cross-check (design D11, task 15.9–15.13): verifies the
# fine-level ApplyStokesOperatorMDOODZ reproduces BuildStokesOperatorDecoupled's
# assembled A·x to 1e-12 on a small representative mesh.
add_executable(StokesMatvecEquivalence StokesMatvecEquivalence.cpp)
target_link_libraries(StokesMatvecEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(StokesMatvecEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB)
add_test(NAME StokesMatvecEquivalence COMMAND StokesMatvecEquivalence)
set_tests_properties(StokesMatvecEquivalence PROPERTIES TIMEOUT 60)

# AnalyticalBench fixtures (SolKz GMG↔CHOLMOD equivalence, §4.2; the
# analytic-SolKz benchmark test binary is SolKzGmgEquivalence below).
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/AnalyticalBench DESTINATION ${PROJECT_BINARY_DIR}/TESTS)

# -----------------------------------------------------------------------
# add-gmg-stokes-defence §4 — easy GMG↔CHOLMOD equivalence fixtures.
# Each test runs a GMG twin (`lin_solver = 3`) and a CHOLMOD twin
# (`lin_solver = 0`) end-to-end and asserts the two linear solvers
# produce the same Stokes solution to 1e-8 relative L2. Live tests
# (not DISABLED_); the "DisabledBenchmarks" folder name refers to the
# upstream DISABLED_* fixture moniker that §3 cleared.
# -----------------------------------------------------------------------

# TopoRelax — free-surface relaxation (§4.1). Fixtures live in TopoBench/.
add_executable(TopoRelaxGmgEquivalence DisabledBenchmarks/TopoRelaxGmgEquivalence.cpp)
target_link_libraries(TopoRelaxGmgEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(TopoRelaxGmgEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME TopoRelaxGmgEquivalence COMMAND TopoRelaxGmgEquivalence)
set_tests_properties(TopoRelaxGmgEquivalence PROPERTIES TIMEOUT 900)

# SolKz high-contrast Stokes (§4.2). Fixtures live in AnalyticalBench/.
add_executable(SolKzGmgEquivalence DisabledBenchmarks/SolKzGmgEquivalence.cpp)
target_link_libraries(SolKzGmgEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(SolKzGmgEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME SolKzGmgEquivalence COMMAND SolKzGmgEquivalence)
set_tests_properties(SolKzGmgEquivalence PROPERTIES TIMEOUT 900)

# -----------------------------------------------------------------------
# add-gmg-stokes-defence §7 — anisotropic + shear-band equivalence fixtures.
# Fixtures live in TESTS/SolViBenchmark/ and TESTS/ShearBand/ (copied above).
# Each test runs a GMG twin (`lin_solver = 3`) and a CHOLMOD twin
# (`lin_solver = 0` → KillerSolver under Newton/anisotropy, or
#  `lin_solver = -1` → DirectStokesDecoupledComp where machine-precision
# reference is required) end-to-end and asserts the §7.x acceptance
# criteria for that fixture.
# -----------------------------------------------------------------------

# AnisotropyShearBand — anisotropic SolVi (§7.1). Fixtures in SolViBenchmark/.
add_executable(AnisotropyShearBandGmgEquivalence
    DisabledBenchmarks/AnisotropyShearBandGmgEquivalence.cpp)
target_link_libraries(AnisotropyShearBandGmgEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(AnisotropyShearBandGmgEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME AnisotropyShearBandGmgEquivalence COMMAND AnisotropyShearBandGmgEquivalence)
set_tests_properties(AnisotropyShearBandGmgEquivalence PROPERTIES TIMEOUT 900)

# Shear-band fixtures live in TESTS/ShearBand/ — copy them into the build tree.
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ShearBand DESTINATION ${PROJECT_BINARY_DIR}/TESTS)

# SingleShearBand — Newton + Drucker-Prager (§7.2). Fixtures in ShearBand/.
add_executable(SingleShearBandGmgEquivalence
    DisabledBenchmarks/SingleShearBandGmgEquivalence.cpp)
target_link_libraries(SingleShearBandGmgEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(SingleShearBandGmgEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME SingleShearBandGmgEquivalence COMMAND SingleShearBandGmgEquivalence)
set_tests_properties(SingleShearBandGmgEquivalence PROPERTIES TIMEOUT 1800)

# ConjugateShearBands — two crossing bands (§7.3). Fixtures in ShearBand/.
add_executable(ConjugateShearBandsGmgEquivalence
    DisabledBenchmarks/ConjugateShearBandsGmgEquivalence.cpp)
target_link_libraries(ConjugateShearBandsGmgEquivalence LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(ConjugateShearBandsGmgEquivalence PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME ConjugateShearBandsGmgEquivalence COMMAND ConjugateShearBandsGmgEquivalence)
set_tests_properties(ConjugateShearBandsGmgEquivalence PROPERTIES TIMEOUT 1800)

# -----------------------------------------------------------------------
# add-gmg-stokes-defence §6 — V-cycle animation dump generator.
# Runs the 41×41 SolCx-analog fixture with `gmg_dump_vcycle = 1` and
# emits HDF5 snapshots under SolViRes41_gmg_vcycle_anim/vcycle_dump/
# consumed by benchmarks/vcycle_vis/make_gif.py. Tagged experimental
# so default ctest skips the regeneration step — only re-run when
# updating defence figures.
# -----------------------------------------------------------------------
add_executable(VcycleAnimGenerate VcycleAnimGenerate.cpp)
target_link_libraries(VcycleAnimGenerate LINK_PUBLIC mdoodz GTest::gtest_main)
target_include_directories(VcycleAnimGenerate PRIVATE ${PROJECT_SOURCE_DIR}/MDLIB ${PROJECT_SOURCE_DIR}/TESTS)
add_test(NAME VcycleAnimGenerate COMMAND VcycleAnimGenerate)
set_tests_properties(VcycleAnimGenerate PROPERTIES LABELS "experimental" TIMEOUT 60)
