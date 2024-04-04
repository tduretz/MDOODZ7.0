import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie")
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, LsqFit

# Select your Makie of Choice: NOTE one has to quit the session when switching from CairoMakie to GLMakie or vice versa.
choice = 2 # 1 => GLMakie; 2 => CairoMakie
MakieOptions = ["import GLMakie as MoC",        # necessary for    interactive 'DataInspector', not supported on any 'headless server'
                "import CairoMakie as MoC"]     # does not support interactive 'DataInspector', runs also on any     'headless server'
eval(Meta.parse(MakieOptions[choice]))
MoC.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
MoC.Makie.inline!(false)

function plotData(NameVect, ResVect, TimeVect, SxzVect, GammaVect, path)

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

    f = MoC.Figure(resolution = (resol, resol), fontsize=35)



    ax1 = MoC.Axis(f[1, 1], title = L"\text{Reduction of effective shear viscosity }η_{\text{eff}}", xlabel = L"\text{horizontal grid resolution nx}", ylabel = L"\frac{η_{\text{eff}}}{η_{\text{0}}}")
    
    

    ResVect2_old = 0.0
    SXZ_atG_old  = 0.0
    for j in eachindex(GammaVect)

        Time_atG   = []
        Gamma_atG  = []
        SXZ_atG    = []
        ResVect2   = []

        gamma = GammaVect[j]

        for i in eachindex(NameVect)
            Time    = TimeVect[i]
            sxz     = SxzVect[i]
    
            Time    = Time[2:end]   # Time = 0.0 no good - unphysically high stresses
            sxz     = sxz[2:end]    # Time = 0.0 no good


            time_atG   = Time[abs.(Time .- gamma/2) .== minimum(abs.(Time .- gamma/2))] # time at timestep closest to the corresponding gamma
            time_atG   = time_atG[1]                    # convert vector (with one element) to a single number
            gamma_atG  = time_atG*2    
            sxz_atG    = sxz[Time .== time_atG]/sxz[1]   # normalize by initial shear stress (to get reduction)
            sxz_atG    = sxz_atG[1]                    # convert vector (with one element) to a single number
            @show NameVect[i], gamma_atG, sxz_atG

            tolerance = 1e-2
            if abs(gamma - gamma_atG) < tolerance
                Time_atG   = push!(Time_atG, time_atG)
                Gamma_atG  = push!(Gamma_atG, gamma_atG)
                SXZ_atG    = push!(SXZ_atG, sxz_atG)
                ResVect2   = push!(ResVect2, ResVect[i])
            end

        end

        ResVect2 = Float64.(ResVect2)
        SXZ_atG = Float64.(SXZ_atG)

        if j > 1
            #ResVect2 = ResVect2[ SXZ_atG .<= SXZ_atG_old ]
            #SXZ_atG  = SXZ_atG[ SXZ_atG .<= SXZ_atG_old ]
            for k_old in eachindex(ResVect2_old)        # idx in ResVect2_old
                K_new = (ResVect2 .== ResVect2_old[k_old])  # idx vector in ResVect2
                if sum(K_new) > 0
                    if (SXZ_atG[K_new])[1] > SXZ_atG_old[k_old] && gamma > 1.0
                        # kill element
                        ResVect2 = ResVect2[.!(K_new)]
                        SXZ_atG  = SXZ_atG[.!(K_new)]
                    else
                        SXZ_atG_old[k_old] = (SXZ_atG[K_new])[1]
                    end
                end
            end
        else
            ResVect2_old = ResVect2
            SXZ_atG_old  = SXZ_atG
        end

        #MoC.lines!(ax1, ResVect2, SXZ_atG, label=L"\gamma = %$(GammaVect[j])")
        #v = j/length(GammaVect)
        #MoC.scatterlines!(ax1, ResVect2, SXZ_atG, label=L"\gamma = %$(GammaVect[j])", color = MoC.RGBf(v, 0, 1-v))
        MoC.scatterlines!(ax1, ResVect2, SXZ_atG, label=L"%$(GammaVect[j])")

        ## exponential curve fitting; idea from: https://julianlsolvers.github.io/LsqFit.jl/latest/getting_started/
        #model(t,p) = p[1] * exp.(-p[2] * t)
#
        #tdata   = ResVect2
        #ydata   = SXZ_atG
        #p0      = [1.0, 0.15]
#
        #fit = curve_fit(model, tdata, ydata, p0)
        #param = fit.param
        ##std_error = standard_error(fit)
#
        #@show param#, std_error

    end


    MoC.limits!(0,2000,0,1.3)
    ax1.xticks = 0:250:1750
    ax1.yticks = 0:0.25:1.50

    MoC.axislegend(L"\text{gamma }\gamma"; position=:rt, labelsize=20, nbanks = 2)    

    Print2Disk( f, path, "z_ViscRed_B3converg")

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

function collectData!(path, name, file_end_max, NameVect, ResVect, TimeVect, SxzVect, autoManu)

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
        file_step = 10
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

    # get x-resolution
    model  = ExtractData( filename, "/Model/Params") 
    nvx    = Int(model[4])
    ncx    = nvx-1

    NameVect    = push!(NameVect, name)
    ResVect     = push!(ResVect, ncx)
    TimeVect    = push!(TimeVect, Time_time)
    SxzVect     = push!(SxzVect, sxz_mean_time)
    return
end

function main()

    NameVect    = Vector{String}()
    ResVect     = []
    TimeVect    = []
    SxzVect     = []

    #GammaVect   = [1,2,4]
    #GammaVect   = 0.1:0.1:5.9
    GammaVect   = 0.1:0.1:8

    # Collect all necessary data from desired output files
    
    # all setups B3 at different resolutions
    collectData!("/users/whalter1/work/aniso_fix/B3_250/" , "B3 250" ,  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_500/" , "B3 500" ,  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_750/" , "B3 750" ,  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_1000/", "B3 1000",  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_1250/", "B3 1250",  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    collectData!("/users/whalter1/work/aniso_fix/B3_1500/", "B3 1500",  25000, NameVect, ResVect, TimeVect, SxzVect, "auto")
    
    # Plot Data
    plotData(NameVect, ResVect, TimeVect, SxzVect, GammaVect, "/users/whalter1/work/aniso_fix/viscosityReduction/")

end

main()


