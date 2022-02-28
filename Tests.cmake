include_directories(TESTS)

configure_file("TESTS/ShearTemplate.txt" "${CMAKE_BINARY_DIR}/ShearTemplate.txt" COPYONLY)
add_executable(ShearTemplate_test TESTS/ShearTemplate_test.c ${SOURCE_FILES} SOURCE/set_ShearTemplate.c)
target_link_libraries(ShearTemplate_test PRIVATE ${DEPENDENCIES})

enable_testing()

add_test(shearTemplate ShearTemplate_test)
