import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using HDF5, CairoMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)
My = 1e6*365*24*3600

function main()

    # Set the path to your files"
    path = "/users/lcandiot/nas/D2c/PhD_LCandioti_2016_2021/project3_UHProckExhumation/dataAndSourceCode/MODEL1_1x1k_Serp_25Phi_1Prefpwl_NewCode_6kmSerp_WesterlyGranite_ReducedConvRate10mmyr_ParticlesCorr/"
    markerTimeFileName = string(path, "MarkerData_xmin_50000.0_xmax_250000.0_zmin_-20000.0_zmax_10000.0_untilTime_177.922503077374.h5")
    
    # Plot type
    # field = :PT_path
    field = :vz_history
    
    # Unit conversion
    cm_yr = 100*3600*24*365.25
    Myr   = 1.0/(1e6*3600*24*365.25)
    res   = 4

    # Figure
    f = Figure(resolution = (500, 500), fontsize=25)

    # P-T- path
    if field==:PT_path
        # Read History
        P = ExtractData( markerTimeFileName, "History/Pm_time")
        T = ExtractData( markerTimeFileName, "History/Tm_time")

        # Visualize
        ax1 = Axis(f[1, 1], title = L"$P$ - $T$ - paths", xlabel = L"$T$ [C]", ylabel = L"$P$ [GPa]")
        for iPlot = 1:5000:size(P,1) 
            lines!(ax1, T[iPlot,:] .- 273, P[iPlot,:]/1e9)
        end
        Print2Disk(f, path, field, res)
    end
    
    # Vertical velocity
    if field==:vz_history
        # Read History
        vzm = ExtractData( markerTimeFileName, "History/vzm_time")
        t   = ExtractData( markerTimeFileName, "History/tm_time" )

        # Visualize
        ax1 = Axis(f[1, 1], title = L"$Vz$ History", xlabel = L"$t$ [Myr]", ylabel = L"$Vz$ [cm/yr]")
        for iPlot = 1:5000:size(vzm,1) 
            lines!(ax1, t[1,:] .* Myr, vzm[iPlot,:] .* cm_yr)
        end
        Print2Disk(f, path, field, res)
    end
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, path, field, res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*".png", f, px_per_unit = res) 
end

main()
