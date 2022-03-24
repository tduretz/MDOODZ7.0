add_executable(ShearTemplate_test ShearTemplate_test.c)
configure_file("TESTS/ShearTemplate.txt" "${CMAKE_BINARY_DIR}/ShearTemplate.txt" COPYONLY)
target_link_libraries(ShearTemplate_test LINK_PUBLIC MDLIB)

add_executable(shearPwl_test TESTS/Shear_pwl_test.c ${SOURCE_FILES} SOURCE/set_Shear_pwl.c)
configure_file("TESTS/Setup03_strongCircleAnisoConst.txt" "${CMAKE_BINARY_DIR}/Setup03_strongCircleAnisoConst.txt" COPYONLY)
target_link_libraries(shearPwl_test LINK_PUBLIC MDLIB)

enable_testing()

add_test(shearTemplate ShearTemplate_test)
add_test(shearPwl Shear_pwl_test)
