cmake_minimum_required(VERSION 3.16)
project(MDOODZ C)

set(CMAKE_C_STANDARD 11)


if (NOT MODEL)
    set(MODEL ShearTemplate)
    message("MODEL is not set, setting a default '${MODEL}' model")
    message("To set a model use -DMODEL=ModelName option")
endif ()

include_directories(SOURCE)
set(SOURCE_FILES
        SOURCE/AdvectionRoutines.c
        SOURCE/BC.c
        SOURCE/ChemicalRoutines.c
        SOURCE/FD_Jacobian.c
        SOURCE/FlowLaws.c
        SOURCE/free_surf.c
        SOURCE/GridRoutines.c
        SOURCE/HDF5Output.c
        SOURCE/header_MDOODZ.h
        SOURCE/InputOutput.c
        SOURCE/Main_DOODZ.c
        SOURCE/MemoryAllocFree.c
        SOURCE/MiscFunctions.c
        SOURCE/ParticleReseeding.c
        SOURCE/ParticleRoutines.c
        SOURCE/RheologyDensity.c
        SOURCE/RheologyParticles.c
        SOURCE/Routines_BEN.c
        SOURCE/set_${MODEL}.c
        SOURCE/Solvers.c
        SOURCE/SparseTools.c
        SOURCE/StokesAssemblyCoupled.c
        SOURCE/StokesAssemblyDecoupled.c
        SOURCE/StokesRoutines.c
        SOURCE/ThermalRoutines.c
        SOURCE/ThermalSolver.c)

set(SuiteSparse_DIR ${CMAKE_SOURCE_DIR}/lib/suitesparse/install/lib/cmake/suitesparse-5.4.0)

include(FindHDF5)
find_package(HDF5 COMPONENTS C)
find_package(SuiteSparse CONFIG REQUIRED)
configure_file("SOURCE/${MODEL}.txt" "${CMAKE_BINARY_DIR}/setup.txt" COPYONLY)


set(DEPENDENCIES
        SuiteSparse::suitesparseconfig
        SuiteSparse::amd
        SuiteSparse::btf
        SuiteSparse::camd
        SuiteSparse::ccolamd
        SuiteSparse::colamd
        SuiteSparse::cholmod
        SuiteSparse::cxsparse
        SuiteSparse::umfpack
        hdf5::hdf5)

add_executable(MDOODZ SOURCE/main.c ${SOURCE_FILES})
target_link_libraries(MDOODZ PUBLIC ${DEPENDENCIES})

message(HDF5_C_DEFINITIONS=${HDF5_C_DEFINITIONS})
message(HDF5_LIBRARIES=${HDF5_LIBRARIES})
message(HDF5_C_INCLUDE_DIRS=${HDF5_C_INCLUDE_DIRS})


# TESTING

include_directories(TESTS)

add_executable(ShearTemplate_test TESTS/ShearTemplate_test.c ${SOURCE_FILES})
target_link_libraries(ShearTemplate_test PRIVATE ${DEPENDENCIES})
target_include_directories(ShearTemplate_test PUBLIC ${HDF5_INCLUDE})

enable_testing()

add_test(shearTemplate ShearTemplate_test)
