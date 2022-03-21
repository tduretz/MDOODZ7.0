cmake_minimum_required(VERSION 3.16)

project(MDOODZ C)

set(CMAKE_C_STANDARD 11)

if (NOT MODEL)
    set(MODEL LithoScale)
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
        SOURCE/Solvers.c
        SOURCE/SparseTools.c
        SOURCE/StokesAssemblyCoupled.c
        SOURCE/StokesAssemblyDecoupled.c
        SOURCE/StokesRoutines.c
        SOURCE/ThermalRoutines.c
        SOURCE/ThermalSolver.c)

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/suitesparse/install/lib/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/lib/suitesparse/install/lib/cmake/suitesparse-5.4.0)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/lib/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/lib/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
else()
    message(FATAL_ERROR "Please set up a path to the suitesparse-config.cmake file")
endif()

set(HDF5_DIR ${PROJECT_SOURCE_DIR}/lib/hdf5/install/share/cmake)
FIND_PACKAGE(SuiteSparse CONFIG)
FIND_PACKAGE(HDF5 CONFIG)

if (EXISTS ${PROJECT_SOURCE_DIR}/SOURCE/${MODEL}.txt)
    configure_file("${PROJECT_SOURCE_DIR}/SOURCE/${MODEL}.txt" "${CMAKE_BINARY_DIR}/setup.txt" COPYONLY)
    message("Using setup file ${MODEL}.txt")
else()
    message("Setup file ${MODEL}.txt isn't found, using default defaultSetup.txt")
    configure_file("${PROJECT_SOURCE_DIR}/SOURCE/defaultSetup.txt" "${CMAKE_BINARY_DIR}/setup.txt" COPYONLY)
endif()


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
        hdf5-static)

add_executable(MDOODZ SOURCE/main.c ${SOURCE_FILES} SOURCE/set_${MODEL}.c)
target_link_libraries(MDOODZ PRIVATE ${DEPENDENCIES})

option(OMP "Set this for using OpenMP" OFF)
option(OPT "Set this for using optimisation" OFF)

if (OPT)
    message("Using optimization")
    set(CMAKE_C_FLAGS "-O3 -ftree-vectorize -funroll-loops")
else()
    message("Not using optimization")
    set(CMAKE_C_FLAGS "-g -Wall -O0 -fno-inline -fno-omit-frame-pointer")
endif()

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    target_link_libraries(MDOODZ PUBLIC OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "-fopenmp -D _OMP_")
    message("Using OpenMP")
endif()

# TESTING

include(Tests.cmake)
