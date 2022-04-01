if (EXISTS ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib/cmake/suitesparse-5.4.0)
    FIND_PACKAGE(SuiteSparse CONFIG)
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
elseif(EXISTS ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
    set(SuiteSparse_DIR ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib64/cmake/suitesparse-5.4.0)
    FIND_PACKAGE(SuiteSparse CONFIG)
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
else()
    include(FindSuiteSparse.cmake)
    message(SuiteSparse_LIBRARIES=${SuiteSparse_LIBRARIES})
    message(SuiteSparse_INCLUDE_DIRS=${SuiteSparse_INCLUDE_DIRS})
    set(DEPENDENCIES ${SuiteSparse_LIBRARIES})
    include_directories(${SuiteSparse_INCLUDE_DIRS})
endif()

include(FindHDF5)
if (HDF5_FOUND)
    find_package(HDF5 COMPONENTS C)
    set(DEPENDENCIES ${DEPENDENCIES} ${HDF5_LIBRARIES})
    include_directories(${HDF5_INCLUDE_DIR})
    message("HDF5 include dir ${HDF5_INCLUDE_DIR}")
else()
    set(HDF5_DIR ${PROJECT_SOURCE_DIR}/deps/hdf5/install/share/cmake)
    FIND_PACKAGE(SuiteSparse CONFIG)
    set(DEPENDENCIES ${DEPENDENCIES} hdf5-static)
endif()

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    set(DEPENDENCIES ${DEPENDENCIES} OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "-fopenmp -D _OMP_")
    message("Using OpenMP")
endif()