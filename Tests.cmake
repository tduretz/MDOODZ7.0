include_directories(TESTS)

configure_file("TESTS/ShearTemplate.txt" "${CMAKE_BINARY_DIR}/ShearTemplate.txt" COPYONLY)
add_executable(ShearTemplate_test TESTS/ShearTemplate_test.c ${SOURCE_FILES} SOURCE/set_ShearTemplate.c)
target_link_libraries(ShearTemplate_test PRIVATE ${DEPENDENCIES})

configure_file("TESTS/Setup03_strongCircleAnisoConst.txt" "${CMAKE_BINARY_DIR}/Setup03_strongCircleAnisoConst.txt" COPYONLY)
add_executable(shearPwl_test TESTS/Shear_pwl_test.c ${SOURCE_FILES} SOURCE/set_Shear_pwl.c)
target_link_libraries(shearPwl_test PRIVATE ${DEPENDENCIES})

add_executable(matrices_test
        TESTS/matrices_test.c
        UTILS/matrices.c
        UTILS/matrices.h)

add_executable(hdf5read_test
        TESTS/hdf5read_test.c
        UTILS/hdf5read.c
        UTILS/hdf5read.h)
target_link_libraries(hdf5read_test hdf5-static)

enable_testing()

add_test(shearTemplate ShearTemplate_test)
add_test(shearPwl Shear_pwl_test)
add_test(matrices matrices_test)
add_test(hdf5read hdf5read_test)

