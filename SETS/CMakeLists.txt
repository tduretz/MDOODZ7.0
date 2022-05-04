include(AddSet.cmake)

if (SET)
    add_set(${SET})
else()
    add_set(ShearTemplate)
    add_set(RiftingPauline)
    add_set(TopoBenchCase1)
endif()

