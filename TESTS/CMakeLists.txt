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

# AnalyticalBench fixtures (SolKz and related analytic benchmarks).
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/AnalyticalBench DESTINATION ${PROJECT_BINARY_DIR}/TESTS)

# ShearBand fixtures (shear-zone benchmarks).
file(COPY ${PROJECT_SOURCE_DIR}/TESTS/ShearBand DESTINATION ${PROJECT_BINARY_DIR}/TESTS)
