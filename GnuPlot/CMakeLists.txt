set (CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
include(FindHDF5)
find_package(HDF5 COMPONENTS CXX)

add_executable(StrainStressGif StrainStressGif.cpp)

target_link_libraries(StrainStressGif LINK_PUBLIC Eigen3::Eigen ${HDF5_LIBRARIES})

set_target_properties(StrainStressGif
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/GnuPlot)
