function(add_set SETNAME)
    add_executable(${SETNAME} ${SETNAME}.c)
    target_link_libraries(${SETNAME} LINK_PUBLIC mdoodz)
    set_target_properties(${SETNAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME})
    if (TXT AND SET)
        if (NOT EXISTS ${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME}/${SETNAME}.txt)
            configure_file("${PROJECT_SOURCE_DIR}/SETS/${TXT}" "${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME}/${TXT}" COPYONLY)
        endif()
    endif()

    if (NOT EXISTS ${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME}/${SETNAME}.txt)
        configure_file("${PROJECT_SOURCE_DIR}/SETS/${SETNAME}.txt" "${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME}/${SETNAME}.txt" COPYONLY)
    endif()
endfunction()