include(AddSet.cmake)

if (SET)
    add_set(${SET})
else()
    add_set(ShearHeatingDuretz14)
    add_set(ViscousInclusion)
    add_set(PinchSwellGSE)
    add_set(ShearTemplate)
    add_set(ShearTemplateAniso)
    add_set(RiftingPauline)
    add_set(TopoBenchCase1)
    add_set(QuartzCoesite) # not yet finished
    add_set(MLPS_Ellipses) # not yet finished
    add_set(CollisionZone) # not yet finished
    add_set(AnisotropyLayer)
endif()

