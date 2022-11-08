set (CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
include(FindHDF5)
find_package(HDF5 COMPONENTS CXX)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

add_executable(visualtests
        VisualTests.cpp
        HDF5pp.h
        visual-tests.h
        RiftingChenin.cpp
        TopoBenchCase1.cpp
        ShearTemplate.cpp
        ShearHeatingDuretz14.cpp
        VEP.cpp
        GSE.cpp)

target_link_libraries(visualtests LINK_PUBLIC mdoodz Eigen3::Eigen ${HDF5_LIBRARIES} GTest::gtest_main)

add_test(GnuplotTests visualtests)

message(CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR})

configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingChenin.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/RiftingChenin.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingChenin.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/RiftingChenin.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingCheninReference.h5" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/RiftingCheninReference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearTemplate.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearTemplate.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplateReference.h5" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearTemplateReference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate1Reference.h5" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearTemplate1Reference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/TopoBenchCase1.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/TopoBenchCase1.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/TopoBenchCase1.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/TopoBenchCase1.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14Reference.h5" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/ShearHeatingDuretz14Reference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/inferno.pal" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/inferno.pal" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/PinchSwellGSE.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/PinchSwellGSE.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/gse.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/gse.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/turbo.pal" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/turbo.pal" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/VEP.gnu" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/VEP.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/VEP_Duretz18.txt" "${CMAKE_BINARY_DIR}/VISUAL_TESTS/VEP_Duretz18.txt" COPYONLY)
file(COPY ${PROJECT_SOURCE_DIR}/VISUAL_TESTS/GSE_REF DESTINATION ${CMAKE_BINARY_DIR}/VISUAL_TESTS/)
file(COPY ${PROJECT_SOURCE_DIR}/VISUAL_TESTS/VEP_REF DESTINATION ${CMAKE_BINARY_DIR}/VISUAL_TESTS/)
