import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie")
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG

# Select your Makie of Choice: NOTE one has to quit the session when switching from CairoMakie to GLMakie or vice versa.
choice = 2 # 1 => GLMakie; 2 => CairoMakie
MakieOptions = ["import GLMakie as MoC",        # necessary for    interactive 'DataInspector', not supported on any 'headless server'
                "import CairoMakie as MoC"]     # does not support interactive 'DataInspector', runs also on any     'headless server'
eval(Meta.parse(MakieOptions[choice]))
MoC.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
MoC.Makie.inline!(false)

s  = 1
d  = 24*3600*s
y  = 365*d
My = 1e6*y

function main(path)

    # Select one field to visualise
    field = :ViscosityReduction

    # Switches
    displayfig  = false     # display figures on screen (not possible on any "headless server")
    printfig    = true      # print figures to disk
    resol       = 1600



    # Scaling
    Lc = 1.
    tc = s; tMytext= "s"
    Vc = 1e-9

    # Time loop
    if displayfig==1 || printfig==1

        istep=file_end 
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model  = ExtractData( filename, "/Model/Params")
        xc     = ExtractData( filename, "/Model/xc_coord")
        zc     = ExtractData( filename, "/Model/zc_coord")
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")
        xv_hr  = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr  = ExtractData( filename, "/VizGrid/zviz_hr")
        xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        @show "Time step" istep
        #@show "Model apect ratio" Lx/Lz
        @show "Model time" t/My


        #####################################

        f = MoC.Figure(resolution = (Lx/Lz*resol, resol), fontsize=25)

        if field==:ViscosityReduction
            ax1 = MoC.Axis(f[1, 1], title = L"Reduction of effective shear viscosity at $t$ = %$(tMy) %$(tMytext)", xlabel = L"shear deformation γ", ylabel = L"η_{s eff} / η_{s init}")
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)

            ax2 = MoC.Axis(f[2, 1], xlabel = L"$x$ [m]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)

            ax3 = MoC.Axis(f[3, 1], xlabel = L"$x$ [m]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)

            @show minimum(Vx_grid)
            @show minimum(Vx_mark)
            @show maximum(Vx_grid)
            @show maximum(Vx_mark)
        end
        
        if printfig Print2Disk( f, path, string(field), istep, xshift) end

        if inspector
            DataInspector(f)
        end
        if displayfig
            display(f)
        end
        sleep(nap)
    end
end


function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, path, field, istep, xshift; res=4)
    if xshift > 0.0
        path1 = path*"/_$field/shift/"
    else
        path1 = path*"/_$field/"
    end
    mkpath(path1)
    MoC.save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

function collectData(path, name, file_end_max)

    # find file_end
    # Check through all Output files (*.gzip.h5)
    Files  = sort(filter(x->endswith(x,".h5"), readdir(path)))
    file_end  = parse(Int64, replace(Files[end] , r"Output"=>"", r".gzip.h5"=>""))
    if file_end > file_end_max; file_end = file_end_max; end
    filename = string(path, @sprintf("Output%05d.gzip.h5", file_end))
    @show filename

    Time_time       = Float64.(ExtractData( filename, "/TimeSeries/Time_time")); 
    Short_time      = Float64.(ExtractData( filename, "/TimeSeries/Short_time")); 
    Work_time       = Float64.(ExtractData( filename, "/TimeSeries/Work_time")); 
    Uthermal_time   = Float64.(ExtractData( filename, "/TimeSeries/Uthermal_time")); 
    Uelastic_time   = Float64.(ExtractData( filename, "/TimeSeries/Uelastic_time")); 
    T_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/T_mean_time")); 
    P_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/P_mean_time")); 
    τxx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxxd_mean_time")); 
    τzz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/szzd_mean_time")); 
    sxz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxz_mean_time")); 
    τII_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Tii_mean_time")); 
    ε̇xx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exxd_mean_time")); 
    ε̇zz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/ezzd_mean_time")); 
    ε̇xz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exz_mean_time")); 
    ε̇II_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Eii_mean_time")); 

    return name, Time_time, sxz_mean_time
end

# Set the path to your files
collectData("/users/whalter1/work/aniso_fix/A1_1000/", "A1 1000", 12000)
collectData("/users/whalter1/work/aniso_fix/A1_500/" , "A1 500" , 12000)
collectData("/users/whalter1/work/aniso_fix/A1_250/" , "A1 250" , 12000)

collectData("/users/whalter1/work/aniso_fix/A2_1000/", "A2 1000", 12000)
collectData("/users/whalter1/work/aniso_fix/A2_500/" , "A2 500" , 12000)
collectData("/users/whalter1/work/aniso_fix/A2_250/" , "A2 250" , 12000)

collectData("/users/whalter1/work/aniso_fix/A3_1000/", "A3 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/A3_500/" , "A3 500" ,  12000)
collectData("/users/whalter1/work/aniso_fix/A3_250/" , "A3 250" ,  12000)

collectData("/users/whalter1/work/aniso_fix/A4_1000/", "A4 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/A4_500/" , "A4 500" ,  12000)
collectData("/users/whalter1/work/aniso_fix/A4_250/" , "A4 250" ,  12000)

collectData("/users/whalter1/work/aniso_fix/A5_1000/", "A5 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/A5_500/" , "A5 500" ,  12000)


collectData("/users/whalter1/work/aniso_fix/A6_1000/", "A6 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/A6_500/" , "A6 500" ,  12000)

collectData("/users/whalter1/work/aniso_fix/B1_1000/", "B1 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/B1_500/" , "B1 500" ,  12000)
collectData("/users/whalter1/work/aniso_fix/B1_250/" , "B1 250" ,  12000) 
collectData("/users/whalter1/work/aniso_fix/B2_1000/", "B2 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/B2_500/" , "B2 500" ,  12000)
collectData("/users/whalter1/work/aniso_fix/B2_250/" , "B2 250" ,  12000)
collectData("/users/whalter1/work/aniso_fix/B3_1000/", "B3 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/B3_500/" , "B3 500" ,  12000)
collectData("/users/whalter1/work/aniso_fix/B3_250/" , "B3 250" ,  12000)


collectData("/users/whalter1/work/aniso_fix/MWE_1000/", "MWE 1000",  12000)
collectData("/users/whalter1/work/aniso_fix/MWE_500/" , "MWE 500" ,  12000)


main("/users/whalter1/work/aniso_fix/viscosityReduction/")


