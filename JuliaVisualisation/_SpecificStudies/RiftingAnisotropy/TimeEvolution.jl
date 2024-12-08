# import Pkg
# Pkg.activate(normpath(joinpath(@__DIR__, "../..")))
using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, UnPack
using CairoMakie#, GLMakie
Mak = CairoMakie
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
# fontsize_theme = Theme(fontsize=200)
# set_theme!(fontsize_theme)
const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

@views function main()

    # Set the path to your files
    path  ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_HR/"
    
    # File numbers
    step = [500; 1100; 2000]
    step = [700; 1400; 2900]

    # Select field to visualise
    # field = :Phases
    # field = :Cohesion
    # field = :Density
    # field = :Viscosity  
    # field = :PlasticStrainrate
    # field = :Stress
    field = :StrainRate
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
    # field = :EffectiveFrictionTime
    # field = :ChristmasTree

    # # Define Tuple for enlargment window
    zoom = ( 
        xmin = -150, 
        xmax = 150,
        zmin = -70,
        zmax = 2.0,
    )

    # Switches
    PlotOnTop = (
        ph_contours   = true,  # add phase contours
        fabric        = true,  # add fabric quiver (normal to director)
        T_contours    = true,   # add temperature contours
        topo          = false,
        quiver_origin = false,
        σ1_axis       = false,
        ε̇1_axis       = false,
        vel_vec       = false,
        ϕ_contours    = false,
        PT_window     = false,
    )
    options = (
        printfig    = true,  # print figures to disk
        printvid    = false,
        framerate   = 6,
        α_heatmap   = 0.75,   # transparency of heatmap 
        vel_arrow   = 5,
        vel_scale   = 10000,
        vel_step    = 18,
        nap         = 0.1,    # pause for animation 
        resol       = 1000,
        Lx          = 1.0,
        Lz          = 1.0,
        LAB_color   = true,
        LAB_T       = 1250,
    )

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    scales = (
        Lc = 1e3,
        tc = My,
        Vc = 1.0,
        τc = 1e6,
    )

    probe = (ϕeff = Float64.([]), t  = Float64.([]))
    cm_yr = 100.0*3600.0*24.0*365.25

    #####################################

    # Color palette for phase map
    cmap    = zeros(RGB{Float64}, 7)
    cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
    cmap[2] = RGBA{Float64}(207/255, 089/255, 087/255, 1.)  
    cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
    cmap[4] = RGBA{Float64}(213/255, 213/255, 209/255, 1.) 
    cmap[5] = RGBA{Float64}(117/255, 099/255, 097/255, 1.) 
    cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
    cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
    phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)
      
    #####################################

    @unpack Lc, tc, Vc, τc = scales

    model = ReadFile(path, step[1], scales, options, PlotOnTop)
    @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, Fab, height, group_phases, Δ = model

    if @isdefined zoom
        Lx = (zoom.xmax - zoom.xmin)./scales.Lc 
        Lz = (zoom.zmax - zoom.zmin)./scales.Lc
        window = zoom
    else
        window = ( 
            xmin = minimum(xc)/scales.Lc, 
            xmax = maximum(xc)/scales.Lc,
            zmin = minimum(zc)/scales.Lc,
            zmax = maximum(zc)/scales.Lc,
        )
    end

    #####################################
    ftsz =  12*options.resol/500
    f = Figure(size = (0.8*Lx/4/Lz*options.resol*1.2, options.resol), fontsize=ftsz)

    cmap = :vik

    if field==:StrainRate

        model = ReadFile(path, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model

        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[1, 1], title = L"A) $\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy_string) Ma", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (cmap, options.α_heatmap), colorrange=(-17, -13))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, 1.5Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 170]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 170]) ./ ε̇1.x[:, 170] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        
        ##############################################################

        model = ReadFile(path, step[2], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[2, 1], title = L"B) $\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy_string) Ma", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (cmap, options.α_heatmap), colorrange=(-17, -13))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, 1.5Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[2, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)


        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 350]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 350]) ./ ε̇1.x[:, 350] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        

        # ##############################################################

        model = ReadFile(path, step[3], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[3, 1], title = L"C) $\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (cmap, options.α_heatmap), colorrange=(-17, -13))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, 1.5Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[3, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

                # lines!(ax1, xc./Lc, log10.(ε̇II[:, 700]) )
                # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 700]) ./ ε̇1.x[:, 700] )
                # xlims!(ax1, window.xmin, window.xmax)
                # ylims!(ax1, -180, 180)
                

        ##############################################################
    
        save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/TimeEvolution.png", f, px_per_unit = 4)   
    end

    display(f)
end

function PrincipalStress(τxx, τzz, τxz, P)
    σ1   = (x=zeros(size(τxx)), z=zeros(size(τxx)) )
    τxzc = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end-0,1:end-1] .+ τxz[1:end-1,2:end-0] .+ τxz[2:end-0,2:end-0]) 
    for i in eachindex(τxzc)
        if P[i]>-1e-13
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
     @show name = path1*"$field"*@sprintf("%05d", istep)*".png"
     save(name, f, px_per_unit = res) 
end

main()