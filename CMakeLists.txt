cmake_minimum_required(VERSION 3.21)
project(MDOODZ7 C)

set(CMAKE_C_STANDARD 11)

include_directories(SOURCE)

add_executable(MDOODZ7
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
        SOURCE/set_ShearTemplate.c
        SOURCE/Solvers.c
        SOURCE/SparseTools.c
        SOURCE/StokesAssemblyCoupled.c
        SOURCE/StokesAssemblyDecoupled.c
        SOURCE/StokesRoutines.c
        SOURCE/ThermalRoutines.c
        SOURCE/ThermalSolver.c)

include(FindHDF5)


FIND_PACKAGE(SuiteSparse)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C)

target_link_libraries(MDOODZ7 PRIVATE SuiteSparse::cholmod)

target_link_libraries(MDOODZ7 PRIVATE HDF5)
