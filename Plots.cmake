include_directories(PLOTS)

configure_file("PLOTS/ShearTemplatePlot.txt" "${CMAKE_BINARY_DIR}/ShearTemplatePlot.txt" COPYONLY)
configure_file("PLOTS/ShearTemplate.py" "${CMAKE_BINARY_DIR}/ShearTemplate.py" COPYONLY)
configure_file("PLOTS/ShearTemplateReference.gzip.h5" "${CMAKE_BINARY_DIR}/ShearTemplateReference.gzip.h5" COPYONLY)
add_executable(ShearTemplate_plot PLOTS/ShearTemplate.c ${SOURCE_FILES} SOURCE/set_ShearTemplate.c)
target_link_libraries(ShearTemplate_plot PRIVATE ${DEPENDENCIES})

enable_testing()

add_test(shearTemplate ShearTemplate_plot)