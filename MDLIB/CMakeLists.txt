add_library(MDLIB AdvectionRoutines.c
        BC.c
        ChemicalRoutines.c
        FD_Jacobian.c
        FlowLaws.c
        free_surf.c
        GridRoutines.c
        HDF5Output.c
        InputOutput.c
        Main_DOODZ.c
        MemoryAllocFree.c
        MiscFunctions.c
        ParticleReseeding.c
        ParticleRoutines.c
        RheologyDensity.c
        RheologyParticles.c
        Routines_BEN.c
        Solvers.c
        SparseTools.c
        StokesAssemblyCoupled.c
        StokesAssemblyDecoupled.c
        StokesRoutines.c
        ThermalRoutines.c
        ThermalSolver.c)
target_link_libraries(MDLIB ${DEPENDENCIES})
target_include_directories(MDLIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})