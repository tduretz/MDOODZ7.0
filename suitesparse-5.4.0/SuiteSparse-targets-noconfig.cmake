#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

if (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
    set(LIBRARY_PATH /usr/lib/x86_64-linux-gnu)
else()
    set(LIBRARY_PATH ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Import target "SuiteSparse::metis" for configuration ""
set_property(TARGET SuiteSparse::metis APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::metis PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "m"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libmetis.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::metis )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::metis "${LIBRARY_PATH}/libmetis.a" )

# Import target "SuiteSparse::suitesparseconfig" for configuration ""
set_property(TARGET SuiteSparse::suitesparseconfig APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::suitesparseconfig PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libsuitesparseconfig.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::suitesparseconfig )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::suitesparseconfig "${LIBRARY_PATH}/libsuitesparseconfig.a" )

# Import target "SuiteSparse::amd" for configuration ""
set_property(TARGET SuiteSparse::amd APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::amd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libamd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::amd )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::amd "${LIBRARY_PATH}/libamd.a" )

# Import target "SuiteSparse::btf" for configuration ""
set_property(TARGET SuiteSparse::btf APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::btf PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libbtf.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::btf )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::btf "${LIBRARY_PATH}/libbtf.a" )

# Import target "SuiteSparse::camd" for configuration ""
set_property(TARGET SuiteSparse::camd APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::camd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libcamd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::camd )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::camd "${LIBRARY_PATH}/libcamd.a" )

# Import target "SuiteSparse::ccolamd" for configuration ""
set_property(TARGET SuiteSparse::ccolamd APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::ccolamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libccolamd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::ccolamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::ccolamd "${LIBRARY_PATH}/libccolamd.a" )

# Import target "SuiteSparse::colamd" for configuration ""
set_property(TARGET SuiteSparse::colamd APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::colamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libcolamd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::colamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::colamd "${LIBRARY_PATH}/libcolamd.a" )

# Import target "SuiteSparse::cholmod" for configuration ""
set_property(TARGET SuiteSparse::cholmod APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::cholmod PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libcholmod.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::cholmod )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::cholmod "${LIBRARY_PATH}/libcholmod.a" )

# Import target "SuiteSparse::cxsparse" for configuration ""
set_property(TARGET SuiteSparse::cxsparse APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::cxsparse PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libcxsparse.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::cxsparse )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::cxsparse "${LIBRARY_PATH}/libcxsparse.a" )

# Import target "SuiteSparse::klu" for configuration ""
set_property(TARGET SuiteSparse::klu APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::klu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libklu.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::klu )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::klu "${LIBRARY_PATH}/libklu.a" )

# Import target "SuiteSparse::ldl" for configuration ""
set_property(TARGET SuiteSparse::ldl APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::ldl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libldl.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::ldl )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::ldl "${LIBRARY_PATH}/libldl.a" )

# Import target "SuiteSparse::umfpack" for configuration ""
set_property(TARGET SuiteSparse::umfpack APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::umfpack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libumfpack.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::umfpack )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::umfpack "${LIBRARY_PATH}/libumfpack.a" )

# Import target "SuiteSparse::spqr" for configuration ""
set_property(TARGET SuiteSparse::spqr APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SuiteSparse::spqr PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${LIBRARY_PATH}/libspqr.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SuiteSparse::spqr )
list(APPEND _IMPORT_CHECK_FILES_FOR_SuiteSparse::spqr "${LIBRARY_PATH}/libspqr.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
