module JuliaVisualisation

const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

using Printf, HDF5

include("LoadDataFiles.jl")
export ExtractField, ReadFile, ExtractData, Print2Disk

include("ContoursQuivers.jl")
export AddCountourQuivers!, PrincipalStress

end # module JuliaVisualisation
