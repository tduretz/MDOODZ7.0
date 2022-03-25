function(add_set SETNAME)
    add_executable(${SETNAME} ${SETNAME}.c ${SETNAME}.txt)
    target_link_libraries(${SETNAME} LINK_PUBLIC MDLIB)
    set_target_properties(${SETNAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-exec/${SETNAME})
endfunction()