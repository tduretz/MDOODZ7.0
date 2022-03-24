add_executable(TopoBenchCase1 TopoBenchCase1.c)

target_link_libraries(TopoBenchCase1 LINK_PUBLIC MDLIB ${DEPENDENCIES})

if (EXISTS ${PROJECT_SOURCE_DIR}/SETS/TopoBenchCase1.txt)
    configure_file("${PROJECT_SOURCE_DIR}/SETS/TopoBenchCase1.txt" "${CMAKE_BINARY_DIR}/SETS/setup.txt" COPYONLY)
    message("Using setup file ${MODEL}.txt")
else()
    message("Setup file ${MODEL}.txt isn't found, using default defaultSetup.txt")
    configure_file("${PROJECT_SOURCE_DIR}/SOURCE/defaultSetup.txt" "${CMAKE_BINARY_DIR}/setup.txt" COPYONLY)
endif()