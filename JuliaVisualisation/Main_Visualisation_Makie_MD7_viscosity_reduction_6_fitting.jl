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
        MoC.lines!(ax1, Time, sxz, label=NameVect[i])
    end

    #MoC.limits!(0,16,0,1.3) # 6 short, 16 normal, 50 long
    #ax1.yticks = 0:0.25:1.25
    #MoC.limits!(0,16,0,1.5) # 6 short, 16 normal, 50 long
    #ax1.yticks = 0:0.25:1.50
    MoC.limits!(0,50,0,1.35) # 6 short, 16 normal, 50 long
    ax1.yticks = 0:0.25:1.35

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

function collectData!(path, name, file_end_max, NameVect, TimeVect, SxzVect, autoManu, doFitData)

    # find file_end
    # Check through all Output files (*.gzip.h5)
    Files  = sort(filter(x->endswith(x,".h5"), readdir(path)))
    file_end  = parse(Int64, replace(Files[end] , r"Output"=>"", r".gzip.h5"=>""))
    if file_end > file_end_max; file_end = file_end_max; end
    filename = string(path, @sprintf("Output%05d.gzip.h5", file_end))
    @show filename
        
    if autoManu == "auto"
        Time_time       = Float64.(ExtractData( filename, "/TimeSeries/Time_time")); 
        sxz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxz_mean_time")); 

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

    Time_time       = Time_time[2:end]      # Time = 0.0 no good - unphysically high stresses
    sxz_mean_time   = sxz_mean_time[2:end]  # Time = 0.0 no good
    Time_time       = Time_time[sxz_mean_time .> 0.0] .* 2                  # i.e. gives γ (now not time anymore...)
    sxz_mean_time   = sxz_mean_time[sxz_mean_time .> 0.0]./sxz_mean_time[1] # normalization by first value

    NameVect    = push!(NameVect, name)
    TimeVect    = push!(TimeVect, Time_time)
    SxzVect     = push!(SxzVect, sxz_mean_time)

    if doFitData
        
        # Define the model and fit it to the data
        model(x, p) = func(x, p)
        #fit = curve_fit(model, x_data, y_data, [1.0, 1.0, 1.0])
        #fit = curve_fit(model, x_data, y_data, [-0.0,0.0,0.0,0.0,1.0,-1.0,0.0])
        fit = curve_fit(model, x_data, y_data, [-1.0,-0.5,0.0,0.0,-1.0,1.0,0.5])

        x_fit = collect(range(0, stop=50, length=1000))
        y_fit = model(x_fit, fit.param)

        #println("Fitted parameters: a=$(fit.param[1]), b=$(fit.param[2]), c=$(fit.param[3])")
        #println("Fitted parameters: a=$(fit.param[1]), b=$(fit.param[2]), c=$(fit.param[3]), d=$(fit.param[4]), m=$(fit.param[5]), n=$(fit.param[6]), o=$(fit.param[7])")
        println("Fitted parameters: p=$(fit.param)")

        name_fit = string(name,"_fit")

        NameVect    = push!(NameVect, name_fit)
        TimeVect    = push!(TimeVect, x_fit)
        SxzVect     = push!(SxzVect, y_fit)
    end

    return
end

# Define the function to fit
function func(x, p)
    #return p[1] * exp.(-p[2] * x) .+ p[3]

    a,b,c,d = p[1],p[2],p[3],p[4]
    m,n,o   = p[5],p[6],p[7]

    y_poly2 = a*(x.-d).^2 .+ b.*(x.-d) .+ c

    #y_expo1 = (a*exp.(b.*(x.-c)) .+ d)

    y_expo2 = (m*exp.(n.*(x.-o)) .+ 0.035)

    #y_expo3 = 0.9720593640446269*exp.(-0.6415103930237273.*(x.-1.516131406446425)) .+ 0.035

    #return y_poly2 .+ y_expo2
    return y_expo2
    #return y_expo1 .+ y_expo2
end

function main()

    NameVect    = Vector{String}()
    TimeVect    = []
    SxzVect     = []

    @show NameVect

    # Collect all necessary data from desired output files
    
    # all iso setups
    collectData!("/users/whalter1/work/aniso_fix/B5_1000/", "B5 1000",  50000, NameVect, TimeVect, SxzVect, "auto", true)
    #collectData!("/users/whalter1/work/aniso_fix/B6_1000/", "B6 1000",  100000, NameVect, TimeVect, SxzVect, "auto", true)
    
    # Plot Data
    plotData(NameVect, TimeVect, SxzVect, "/users/whalter1/work/aniso_fix/viscosityReduction_fit/")

end

main()


