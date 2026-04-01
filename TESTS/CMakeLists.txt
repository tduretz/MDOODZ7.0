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
