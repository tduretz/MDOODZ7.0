include(FindSuiteSparse.cmake)
message(SuiteSparse_LIBRARIES=${SuiteSparse_LIBRARIES})
message(SuiteSparse_INCLUDE_DIRS=${SuiteSparse_INCLUDE_DIRS})
set(DEPENDENCIES ${SuiteSparse_LIBRARIES})
include_directories(${SuiteSparse_INCLUDE_DIRS})

include(FindHDF5)
find_package(HDF5 COMPONENTS C)

set(DEPENDENCIES ${DEPENDENCIES} ${HDF5_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIR})
message("HDF5 include dir ${HDF5_INCLUDE_DIR}")

if (OMP)
    set(OpenMP_lgomp_LIBRARY "lgomp")
    find_package(OpenMP)
    set(DEPENDENCIES ${DEPENDENCIES} OpenMP::OpenMP_C)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -D _OMP_")
    message("Using OpenMP")
endif()

# Add the following to fetch pcg-c-basic
include(FetchContent)

FetchContent_Declare(
        pcg_c_basic
        GIT_REPOSITORY https://github.com/imneme/pcg-c-basic.git
        GIT_TAG        master  # Optionally specify a commit, tag, or branch
)

# !!! THE COMMENTED SECTION BELOW IS DEPRECATED !!!
# FetchContent_GetProperties(pcg_c_basic)
# if(NOT pcg_c_basic_POPULATED)
#     FetchContent_Populate(pcg_c_basic)
#     add_library(pcg_c_basic STATIC ${pcg_c_basic_SOURCE_DIR}/pcg_basic.c)
#     target_include_directories(pcg_c_basic PUBLIC ${pcg_c_basic_SOURCE_DIR})
# endif()

# !!! BELOW IS THE MODERN IMPLEMENTATION !!!!
# After a prior FetchContent_Declare(pcg_c_basic ...)

# Populate (modern API; no warning) — if the project had a CMakeLists.txt,
# this would also add_subdirectory(). For a single .c file, it just populates.
FetchContent_MakeAvailable(pcg_c_basic)

# Make sure the fetched content is populated and we have the variables
FetchContent_GetProperties(pcg_c_basic)
# pcg_c_basic_SOURCE_DIR is defined after MakeAvailable()

# Define our own target from the fetched source (guard against double definition)
if (NOT TARGET pcg_c_basic)
  add_library(pcg_c_basic STATIC
    ${pcg_c_basic_SOURCE_DIR}/pcg_basic.c
  )
  target_include_directories(pcg_c_basic PUBLIC
    ${pcg_c_basic_SOURCE_DIR}
  )
endif()

# Add pcg_c_basic to your dependencies
set(DEPENDENCIES ${DEPENDENCIES} pcg_c_basic)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm")