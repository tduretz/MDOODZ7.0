function(add_set SETNAME)
    add_executable(${SETNAME} ${SETNAME}.c)
    target_link_libraries(${SETNAME} LINK_PUBLIC MDLIB ${DEPENDENCIES})
    set_target_properties(${SETNAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/SETS/${SETNAME})
    configure_file("${PROJECT_SOURCE_DIR}/SETS/${SETNAME}.txt" "${CMAKE_BINARY_DIR}/SETS/${SETNAME}/${SETNAME}.txt" COPYONLY)
endfunction()