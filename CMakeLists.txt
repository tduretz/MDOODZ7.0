cmake_minimum_required(VERSION 3.21)
project(MDOODZ7 C)

set(CMAKE_C_STANDARD 11)

include_directories(SOURCE)

add_executable(MDOODZ7 SOURCE/Main_DOODZ.c)

FIND_PACKAGE(SuiteSparse NAMES "SuiteSparse"
        PATHS ${PROJECT_SOURCE_DIR}/lib/SuiteSparse/_build/CMakeFiles/Export/lib/cmake/suitesparse-5.4.0)

target_link_libraries(MDOODZ7 PRIVATE SuiteSparse::cholmod)