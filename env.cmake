set(C_COMPILER gcc)
set(SuiteSparse_DIR /usr/lib/x86_64-linux-gnu/)

if (JULIA)
    set(SuiteSparse_VERBOSE ON)
    set(SuiteSparse_DIR ${CMAKE_INSTALL_PREFIX})
else()
    set(SuiteSparse_DIR /usr/lib/x86_64-linux-gnu/)
endif()