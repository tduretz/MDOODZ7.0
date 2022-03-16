include_directories(TESTS)

configure_file("TESTS/ShearTemplate.txt" "${CMAKE_BINARY_DIR}/ShearTemplate.txt" COPYONLY)
add_executable(ShearTemplate_test TESTS/ShearTemplate_test.c ${SOURCE_FILES} SOURCE/set_ShearTemplate.c)
target_link_libraries(ShearTemplate_test PRIVATE ${DEPENDENCIES})

configure_file("TESTS/Setup03_strongCircleAnisoConst.txt" "${CMAKE_BINARY_DIR}/Setup03_strongCircleAnisoConst.txt" COPYONLY)
add_executable(shearPwl TESTS/Shear_pwl_test.c ${SOURCE_FILES} SOURCE/set_Shear_pwl.c)
target_link_libraries(shearPwl PRIVATE ${DEPENDENCIES})

enable_testing()

add_test(shearTemplate ShearTemplate_test)
add_test(shearPwl Shear_pwl_test)
