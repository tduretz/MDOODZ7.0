find_package(suitesparse REQUIRED)
set(DEPENDENCIES
        SuiteSparse::suitesparseconfig
        SuiteSparse::amd
        SuiteSparse::btf
        SuiteSparse::camd
        SuiteSparse::ccolamd
        SuiteSparse::colamd
        SuiteSparse::cholmod
        SuiteSparse::cxsparse
        SuiteSparse::umfpack)
find_package(hdf5 CONFIG REQUIRED)
set(DEPENDENCIES ${DEPENDENCIES} hdf5::hdf5-static)

link_directories(C:/Users/User/CLionProjects/MDOODZ7.0/cmake-build-debug/MDLIB)

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    set(DEPENDENCIES ${DEPENDENCIES} OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -D _OMP_")
    message("Using OpenMP")
endif()
