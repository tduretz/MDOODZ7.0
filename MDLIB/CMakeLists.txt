add_library(mdoodz SHARED
        AdvectionRoutines.c
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
target_link_libraries(mdoodz ${DEPENDENCIES})
target_include_directories(mdoodz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(mdoodz PROPERTIES
        PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/mdoodz.h)

install(TARGETS mdoodz
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)