if (EXISTS ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib/cmake/suitesparse-5.4.0)
elseif(EXISTS ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
else()
    message(FATAL_ERROR "Please set up a path to the suitesparse-config.cmake file")
endif()

set(HDF5_DIR ${PROJECT_SOURCE_DIR}/deps/hdf5/install/share/cmake)
FIND_PACKAGE(SuiteSparse CONFIG)
FIND_PACKAGE(HDF5 CONFIG)

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

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    set(DEPENDENCIES ${DEPENDENCIES} OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "-fopenmp -D _OMP_")
    message("Using OpenMP")
endif()