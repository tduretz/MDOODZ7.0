set (CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
include(FindHDF5)
find_package(HDF5 COMPONENTS CXX)

add_executable(visualtests
        main.cpp
        HDF5pp.h
        visual-tests.h
        RiftingPauline.cpp
        TopoBenchCase1.cpp
        ShearTemplate.cpp
        ShearHeatingDuretz14.cpp
        GSE.cpp)

target_link_libraries(visualtests LINK_PUBLIC mdoodz Eigen3::Eigen ${HDF5_LIBRARIES})

set_target_properties(visualtests
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/visualtests-out)

configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingPauline.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/RiftingPauline.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingPauline.gnu" "${PROJECT_SOURCE_DIR}/visualtests-out/RiftingPauline.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/RiftingPaulineReference.h5" "${PROJECT_SOURCE_DIR}/visualtests-out/RiftingPaulineReference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearTemplate.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate1.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearTemplate1.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate.gnu" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearTemplate.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplateReference.h5" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearTemplateReference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearTemplate1Reference.h5" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearTemplate1Reference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/TopoBenchCase1.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/TopoBenchCase1.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/TopoBenchCase1.gnu" "${PROJECT_SOURCE_DIR}/visualtests-out/TopoBenchCase1.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.gnu" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearHeatingDuretz14.gnu" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearHeatingDuretz14.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/ShearHeatingDuretz14Reference.h5" "${PROJECT_SOURCE_DIR}/visualtests-out/ShearHeatingDuretz14Reference.h5" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/inferno.pal" "${PROJECT_SOURCE_DIR}/visualtests-out/inferno.pal" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/PinchSwellGSE.txt" "${PROJECT_SOURCE_DIR}/visualtests-out/PinchSwellGSE.txt" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/VISUAL_TESTS/gse.gnu" "${PROJECT_SOURCE_DIR}/visualtests-out/gse.gnu" COPYONLY)
file(COPY ${PROJECT_SOURCE_DIR}/VISUAL_TESTS/GSE_REF DESTINATION ${PROJECT_SOURCE_DIR}/visualtests-out/)