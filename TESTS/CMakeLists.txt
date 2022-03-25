if (VALGRIND)
    find_program( MEMORYCHECK_COMMAND valgrind )
    set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
endif()

add_executable(ShearTemplate_test ShearTemplate_test.c ShearTemplate.txt)
target_link_libraries(ShearTemplate_test LINK_PUBLIC MDLIB)
add_test(shearTemplate ShearTemplate_test)
