if (VALGRIND)
    find_program( MEMORYCHECK_COMMAND valgrind )
    set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
    message("Valgrind is ON for testing")
endif()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

configure_file("${PROJECT_SOURCE_DIR}/TESTS/ShearTemplate.txt" "${PROJECT_BINARY_DIR}/TESTS/ShearTemplate.txt" COPYONLY)
add_executable(NewtownIterationConvergence_test NewtonIterationConvergence.cpp)
target_link_libraries(NewtownIterationConvergence_test LINK_PUBLIC mdoodz GTest::gtest_main)
add_test(NewtownIterationConvergence NewtownIterationConvergence_test)
