cmake_minimum_required(VERSION 3.21)
project(MDOODZ C)

set(CMAKE_C_STANDARD 11)

include_directories(SOURCE)

if (NOT MODEL)
    set(MODEL ShearTemplate)
    message("MODEL is not set, setting a default '${MODEL}' model")
    message("To set a model use -DMODEL=ModelName option")
endif ()

add_executable(MDOODZ
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

FIND_PACKAGE(SuiteSparse)
include(FindHDF5)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C)

configure_file("SOURCE/${MODEL}.txt" "${CMAKE_BINARY_DIR}/setup.txt" COPYONLY)

target_link_libraries(MDOODZ PRIVATE
        SuiteSparse::suitesparseconfig
        SuiteSparse::amd
        SuiteSparse::btf
        SuiteSparse::camd
        SuiteSparse::ccolamd
        SuiteSparse::colamd
        SuiteSparse::cholmod
        SuiteSparse::cxsparse
        SuiteSparse::umfpack
        ${HDF5_LIBRARIES})
