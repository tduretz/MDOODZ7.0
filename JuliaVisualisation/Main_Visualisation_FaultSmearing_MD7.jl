import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
using CairoMakie, GLMakie
const Mak = GLMakie
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

function LoadField(path, istep)
    name     = @sprintf("Output%05d.gzip.h5", istep)
    @info "Reading $(name)"
    filename = string(path, name)

    model  = ExtractData( filename, "/Model/Params")
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1

    fields = (
        model    = ExtractData( filename, "/Model/Params"),
        xc       = ExtractData( filename, "/Model/xc_coord"),
        zc       = ExtractData( filename, "/Model/zc_coord"),
        xv       = ExtractData( filename, "/Model/xg_coord"),
        zv       = ExtractData( filename, "/Model/zg_coord"),
        xvz      = ExtractData( filename, "/Model/xvz_coord"),
        zvx      = ExtractData( filename, "/Model/zvx_coord"),
        xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr"),
        zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr"),
        τzz_t    = ExtractData( filename, "TimeSeries/szzd_mean_time"),
        P_t      = ExtractData( filename, "TimeSeries/P_mean_time"),
        τxz_t    = ExtractData( filename, "TimeSeries/sxz_mean_time"),
        t_t      = ExtractData( filename, "TimeSeries/Time_time"),
        # ηc       = ExtractField(filename,  "/Centers/eta_n", centroids, true, mask_air),
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz)),        
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz)),              
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15,   
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz)),              
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz)),         
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2))),
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1))),
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz)),
        τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz)),
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz)),
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz)),
        ε̇zz   = Float64.(reshape(ExtractData( filename, "/Centers/ezzd"), ncx, ncz)),
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz)),
    )
    dummy = (
        t     = fields.model[1],
        τyy   = -(fields.τzz .+ fields.τxx),
        ε̇yy   = -(fields.ε̇xx .+ fields.ε̇zz),
    )
    fields = merge(fields, dummy)
    dummy = (
        τII   = sqrt.( 0.5*(fields.τxx.^2 .+ fields.τyy.^2 .+ fields.τzz.^2 .+ 0.5*(fields.τxz[1:end-1,1:end-1].^2 .+ fields.τxz[2:end,1:end-1].^2 .+ fields.τxz[1:end-1,2:end].^2 .+ fields.τxz[2:end,2:end].^2 ) ) ),
        ε̇II   = sqrt.( 0.5*(fields.ε̇xx.^2 .+ fields.ε̇yy.^2 .+ fields.ε̇zz.^2 .+ 0.5*(fields.ε̇xz[1:end-1,1:end-1].^2 .+ fields.ε̇xz[2:end,1:end-1].^2 .+ fields.ε̇xz[1:end-1,2:end].^2 .+ fields.ε̇xz[2:end,2:end].^2 ) ) ), 
        τxzc  = 0.25*(fields.τxz[1:end-1,1:end-1] .+ fields.τxz[2:end,1:end-1] .+ fields.τxz[1:end-1,2:end] .+ fields.τxz[2:end,2:end]), 
    )
    return merge(fields, dummy)
end

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 

function main()

    # Set the path to your files
    path = (
        "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/_p10_e18_t3/",
        "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/_p10_e18_t3_nonconv/",
        "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/_p10_e18_t2_nonconv/",
        "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/_p10_e18_t1_nonconv/",
    )
    names = (
        "T3: Fully converged",
        "T3: 1 iteration",
        "T2: 1 iteration",
        "T1: 1 iteration",
    )
    n_files = length(path)

    # File numbers
    file_start = 000
    file_step  = 10
    file_end   = 2000

    # Select field to visualise
    field = :Phases
    # field = :Cohesion
    # field = :Density
    field = :Viscosity 
    # field = :PlasticStrainrate
    # field = :Stress
    # field = :StrainRate
    # field = :Pressure
    # field = :Divergence
    # field = :Temperature
    # field = :Velocity_x
    # field = :Velocity_z
    # field = :Velocity
    # field = :GrainSize
    # field = :Topography
    # field = :TimeSeries
    # field = :AnisotropyFactor
    # field = :MeltFraction
    # field = :TimeSeries
    field = :EffectiveFrictionTime

    # Switches
    printfig    = false  # print figures to disk
    printvid    = false
    framerate   = 3
    PlotOnTop = (
        ph_contours = false,  # add phase contours
        T_contours  = false,   # add temperature contours
        fabric      = false,  # add fabric quiver (normal to director)
        topo        = false,
        σ1_axis     = false,
        vel_vec     = false,
    )
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 
    vel_arrow   = 5
    vel_scale   = 300000
    nap         = 0.1    # pause for animation 
    resol       = 1000
    Lx, Lz      = 1.0, 1.0

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    Lc = 1.0
    tc = My
    Vc = 1.0

    probe = []
    for i_file=1:n_files
        probe_data = (ϕeff = Float64.([]), t  = Float64.([]))
        push!(probe, probe_data)
    end
    cm_yr = 100.0*3600.0*24.0*365.25

    # Time loop
    f = Figure(size = (Lx/Lz*resol*1.2, resol), fontsize=25)

    for istep=file_start:file_step:file_end
    
        dat = []
        for i_file=1:n_files
            data = LoadField(path[i_file], istep)
            push!(dat, data)
        end
            
        #####################################
        empty!(f)
        f = Figure(size = (Lx/Lz*resol*1.2, resol), fontsize=25)

        if field==:EffectiveFrictionTime
            for i_file=1:n_files
                ϕ_eff = -mean(dat[i_file].τxzc[:,end] ./ (dat[i_file].τzz[:,end] .- dat[i_file].P[:,end]) )
                if istep==0 ϕ_eff = 0. end
                push!(probe[i_file].ϕeff, ϕ_eff) 
                push!(probe[i_file].t, dat[i_file].t) 
            end
            if istep==file_end
                ax1 = Axis(f[1, 1], title = L"$ϕ_\mathrm{eff}$", xlabel = L"$t$", ylabel = L"$ϕ_\mathrm{eff}$")
                for i_file=1:n_files
                    lines!(ax1, Float64.(probe[i_file].t), Float64.(probe[i_file].ϕeff), label=names[i_file])
                end
                axislegend(ax1, position = :rb)
            end

        end

        if field!=:EffectiveFrictionTime || istep==file_end
            DataInspector(f)
            display(f)
            sleep(nap)
        end

    end

    yscale = Lz/Lx

    if printfig && printvid
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1080:1080*$(yscale)" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)
    end

end

function PrincipalStress(τxx, τzz, τxz, P)
    σ1   = (x=zeros(size(τxx)), z=zeros(size(τxx)) )
    τxzc = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end-0,1:end-1] .+ τxz[1:end-1,2:end-0] .+ τxz[2:end-0,2:end-0]) 
    for i in eachindex(τxzc)
        if P[i]>1e-13
            σ = [-P[i]+τxx[i] τxzc[i]; τxzc[i] -P[i]+τzz[i]]
            v = eigvecs(σ)
            σ1.x[i] = v[1,1]
            σ1.z[i] = v[2,1]
        end
    end
    return σ1
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, path, field, istep; res=4)
     path1 = path*"/_$field/"
     mkpath(path1)
     save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()