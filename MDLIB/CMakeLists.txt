add_library(mdoodz SHARED
        AdvectionRoutines.c
        AnisotropyRoutines.c
        ChemicalRoutines.c
        FD_Jacobian.c
        FlowLaws.c
        FreeSurface.c
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
        Solvers.c
        SparseTools.c
        StokesAssemblyDecoupled.c
        StokesRoutines.c
        ThermalRoutines.c
        ThermalSolver.c
        Setup.c
        RheologyDensity.h)
target_link_libraries(mdoodz ${DEPENDENCIES})
target_include_directories(mdoodz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(mdoodz PROPERTIES
        PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/mdoodz.h)

install(TARGETS mdoodz
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
