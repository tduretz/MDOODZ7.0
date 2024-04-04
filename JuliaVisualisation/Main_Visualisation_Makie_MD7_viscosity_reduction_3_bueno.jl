import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie")
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics

# Select your Makie of Choice: NOTE one has to quit the session when switching from CairoMakie to GLMakie or vice versa.
choice = 2 # 1 => GLMakie; 2 => CairoMakie
MakieOptions = ["import GLMakie as MoC",        # necessary for    interactive 'DataInspector', not supported on any 'headless server'
                "import CairoMakie as MoC"]     # does not support interactive 'DataInspector', runs also on any     'headless server'
eval(Meta.parse(MakieOptions[choice]))
MoC.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
MoC.Makie.inline!(false)

function plotData(NameVect, TimeVect, SxzVect, path)

    s  = 1
    d  = 24*3600*s
    y  = 365*d
    My = 1e6*y

    # Switches
    resol       = 1200

    # Scaling
    Lc = 1.
    tc = s; tMytext= "s"
    Vc = 1e-9


    #####################################

    f = MoC.Figure(resolution = (2.0*resol, resol), fontsize=35)

    ax1 = MoC.Axis(f[1, 1], title = L"\text{Reduction of effective shear viscosity }η_{\text{eff}}", xlabel = L"\text{shear deformation }γ", ylabel = L"\frac{η_{\text{eff}}}{η_{\text{0}}}")
    for i = 1:length(NameVect)

        Time    = TimeVect[i]
        sxz     = SxzVect[i]

        Time    = Time[2:end]   # Time = 0.0 no good - unphysically high stresses
        sxz     = sxz[2:end]    # Time = 0.0 no good
    
        #@show size(Time[sxz .> 0.0])
        #@show Time[1:10]
        #@show Time[end-10:end]
        #@show sxz[1:10]
        #@show sxz[end-10:end]
        #@show sxz .> 0.0

        MoC.lines!(ax1, Time[sxz .> 0.0] .* 2, sxz[sxz .> 0.0]./sxz[1], label=NameVect[i])
    end

    # define plotting limits
    xmax = 50   # 6 short, 16 normal, 50 long, 60 or more for custom extralong
    ymax = 1.35 # 1.35 normal

    MoC.limits!(0,xmax,0,ymax) 
    ax1.yticks = 0:0.25:ymax

    #MoC.Legend(f[1, 1], ax1)
    #MoC.axislegend(""; position=:rt, labelsize=20, nbanks = 2)    
    MoC.axislegend(""; halign = :right, valign = :top, labelsize=20, nbanks = 6)    
    Print2Disk( f, path, "ViscosityReduction")

end


function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, path, field; res=4)
    mkpath(path)
    MoC.save(path*"$field"*".png", f, px_per_unit = res) 
end

function collectData!(path, name, file_end_max, NameVect, TimeVect, SxzVect, autoManu)

    # find file_end
    # Check through all Output files (*.gzip.h5)
    Files  = sort(filter(x->endswith(x,".h5"), readdir(path)))
    file_end  = parse(Int64, replace(Files[end] , r"Output"=>"", r".gzip.h5"=>""))
    if file_end > file_end_max; file_end = file_end_max; end
    filename = string(path, @sprintf("Output%05d.gzip.h5", file_end))
    @show filename
        
    if autoManu == "auto"
        Time_time       = Float64.(ExtractData( filename, "/TimeSeries/Time_time")); 
        #Short_time      = Float64.(ExtractData( filename, "/TimeSeries/Short_time")); 
        #Work_time       = Float64.(ExtractData( filename, "/TimeSeries/Work_time")); 
        #Uthermal_time   = Float64.(ExtractData( filename, "/TimeSeries/Uthermal_time")); 
        #Uelastic_time   = Float64.(ExtractData( filename, "/TimeSeries/Uelastic_time")); 
        #T_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/T_mean_time")); 
        #P_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/P_mean_time")); 
        #τxx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxxd_mean_time")); 
        #τzz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/szzd_mean_time")); 
        sxz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxz_mean_time")); 
        #τII_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Tii_mean_time")); 
        #ε̇xx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exxd_mean_time")); 
        #ε̇zz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/ezzd_mean_time")); 
        #ε̇xz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exz_mean_time")); 
        #ε̇II_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Eii_mean_time")); 

    else # if autoManu == "manu"
        Time_time        = []
        sxz_mean_time    = []
        file_start = 0
        file_step = 100
        for i = file_start:file_step:file_end
            filename = string(path, @sprintf("Output%05d.gzip.h5", i))
            model  = ExtractData( filename, "/Model/Params") 
            t      = model[1]
            nvx    = Int(model[4])
            nvz    = Int(model[5])
            τxz    = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
            sxz_mean        = mean(τxz)
            Time_time       = push!(Time_time, t)
            sxz_mean_time   = push!(sxz_mean_time, sxz_mean)
            if mod(i,10*file_step) == 0; @show file_end, i; end
        end
    end

    NameVect    = push!(NameVect, name)
    TimeVect    = push!(TimeVect, Time_time)
    SxzVect     = push!(SxzVect, sxz_mean_time)
    return
end

function main()

    NameVect    = Vector{String}()
    TimeVect    = []
    SxzVect     = []

    @show NameVect

    # Collect all necessary data from desired output files
    
    # all aniso setups
    collectData!("/users/whalter1/work/aniso_fix/A2_1000/" , "A2 1000" , 100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/A3_1000/" , "A3 1000" , 100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A4_1000/" , "A4 1000" , 100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A5_1000/" , "A5 1000" , 100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/MWE_1000/", "MWE 1000", 100000, NameVect, TimeVect, SxzVect, "auto")
    
    # all iso setups
    collectData!("/users/whalter1/work/aniso_fix/A1_1000/", "A1 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A6_1000/", "A6 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A7_1000/", "A7 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B1_1000/", "B1 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B2_1000/", "B2 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B2_1500/", "B2 1500",  100000, NameVect, TimeVect, SxzVect, "auto")

    collectData!("/users/whalter1/work/aniso_fix/B3_250/" , "B3 250" ,  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_500/" , "B3 500" ,  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B3_750/" , "B3 750" ,  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B3_1000/", "B3 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_1250/", "B3 1250",  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B3_1500/", "B3 1500",  100000, NameVect, TimeVect, SxzVect, "manu")

    collectData!("/users/whalter1/work/aniso_fix/B4_1500/", "B4 1500",  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B5_1000/", "B5 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    
    collectData!("/users/whalter1/work/aniso_fix/B6_1000/", "B6 1000",  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B6_1500/", "B6 1500",  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B6_2000/", "B6 2000",  100000, NameVect, TimeVect, SxzVect, "auto")
    
    collectData!("/users/whalter1/work/aniso_fix/B7_1000/", "B7 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    
    collectData!("/users/whalter1/work/aniso_fix/B8_1000/", "B8 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B9_1000/", "B9 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    
    collectData!("/users/whalter1/work/aniso_fix/B10_1000/", "B10 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B11_1000/", "B11 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B12_1000/", "B12 1000",  100000, NameVect, TimeVect, SxzVect, "auto")

    collectData!("/users/whalter1/work/aniso_fix/B13_1000/", "B13 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B14_1000/", "B14 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B15_1000/", "B15 1000",  100000, NameVect, TimeVect, SxzVect, "auto")

    collectData!("/users/whalter1/work/aniso_fix/B16_1000/", "B16 1000",  100000, NameVect, TimeVect, SxzVect, "manu")
    collectData!("/users/whalter1/work/aniso_fix/B16_1500/", "B16 1500",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B16_2000/", "B16 2000",  100000, NameVect, TimeVect, SxzVect, "auto")

    collectData!("/users/whalter1/work/aniso_fix/B17_1000/", "B17 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B18_1000/", "B18 1000",  100000, NameVect, TimeVect, SxzVect, "auto")
    
    # Plot Data
    plotData(NameVect, TimeVect, SxzVect, "/users/whalter1/work/aniso_fix/viscosityReduction/")

end

main()


