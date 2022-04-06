if (JULIA AND CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(DEPENDENCIES
            "${CMAKE_INSTALL_PREFIX}/bin/libsuitesparseconfig.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libamd.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libbtf.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libcamd.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libccolamd.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libcolamd.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libcholmod.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libcxsparse.dll"
            "${CMAKE_INSTALL_PREFIX}/bin/libcumfpack.dll")
elseif (EXISTS ${PROJECT_SOURCE_DIR}/deps/suitesparse/install/lib/cmake/suitesparse-5.4.0)
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

if (JULIA AND CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(DEPENDENCIES ${DEPENDENCIES} "${CMAKE_INSTALL_PREFIX}/bin/libhdf5-0.dll")
else()
    include(FindHDF5)
    find_package(HDF5 COMPONENTS C)
    set(DEPENDENCIES ${DEPENDENCIES} ${HDF5_LIBRARIES})
    include_directories(${HDF5_INCLUDE_DIR})
    message("HDF5 include dir ${HDF5_INCLUDE_DIR}")
endif()

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    set(DEPENDENCIES ${DEPENDENCIES} OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "-fopenmp -D _OMP_")
    message("Using OpenMP")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm")