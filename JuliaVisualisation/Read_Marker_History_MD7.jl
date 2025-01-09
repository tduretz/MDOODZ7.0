using DataFramesMeta, DelimitedFiles

function ExtractHistory(
    df     :: DataFrame,
    field  :: String,
    header :: AbstractArray,
    Idx    :: Int
)

    # Find column of queried field
    idx_col = findfirst(x -> x == field, header)
    if typeof(df[Idx, idx_col]) == Int
        return df[Idx, idx_col]
    elseif typeof(df[Idx, idx_col]) == SubString{String}
        return df[Idx, idx_col] |> x -> strip(x, [']', '[']) |> x -> split(x, ',') |> x -> parse.(Float64, x)
    end
end
function ReadMarkerHistory()
    
    # Set path
    fname = "/Users/lcandiot/Developer/MDOODZ7.0/cmake-exec/RiftingChenin/MarkerData.csv"

    # Load data
    data, header = readdlm(fname, ',', header = true)
    df = DataFrame(data, vec(header))

    Tm  = ExtractHistory(df, "T [K]", vec(header), 5)
    Pm  = ExtractHistory(df, "P [Pa]", vec(header), 5)
    xm  = ExtractHistory(df, "x [m]", vec(header), 5)
    zm  = ExtractHistory(df, "z [m]", vec(header), 5)
    phm = ExtractHistory(df, "ph [Int]", vec(header), 5)

    @show Tm
    @show Pm
    @show xm
    @show zm
    @show phm

    # Return
    return nothing
end

ReadMarkerHistory();