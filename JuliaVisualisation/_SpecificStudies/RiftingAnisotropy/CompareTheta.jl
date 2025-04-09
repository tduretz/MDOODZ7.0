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
    path_LR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_MR/"
    path_MR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t45/"
    path_HR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t80/"


    # File numbers
    step = [1360; 860; 1660]
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

    # Define Tuple for enlargment window
    zoom = ( 
        xmin = -150, 
        xmax = 150,
        zmin = -70,
        zmax = 2,
    )

    # Switches
    PlotOnTop = (
        ph_contours   = true,  # add phase contours
        fabric        = true,  # add fabric quiver (normal to director)
        T_contours    = false,   # add temperature contours
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
        vel_step    = 20,
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

    model = ReadFile(path_LR, step[1], scales, options, PlotOnTop)
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
    f = Figure(size = (0.85*Lx/4/Lz*options.resol*1.2, options.resol), fontsize=ftsz)

    # if field==:Phases
    #     ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
    #     hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
    #     # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = :turbo)
    #     hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_dual_hr, colormap = :turbo)
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Viscosity
    #     ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Density
    #     ax1 = Axis(f[1, 1], title = L"$\rho$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap), colorrange=(2600, 3000))  
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$\rho$ [kg.m$^{-3}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Stress
    #     ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, τII./1e6, colormap = (:turbo, α_heatmap)) 
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Pressure
    #     ax1 = Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, P./1e9, colormap = (:turbo, α_heatmap)) #, colorrange=(1,1.2)1e4*365*24*3600
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label =  L"$P$ [GPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Divergence
    #     ax1 = Axis(f[1, 1], title = L"∇⋅V at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, divu, colormap = (:turbo, α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    if field==:StrainRate

        cmap = :vik

        model = ReadFile(path_LR, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model

        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[1, 1], title = L"A) $\theta_\textrm{ini} = 10^\circ$ - $t$ = %$(tMy_string) Ma", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 170]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 170]) ./ ε̇1.x[:, 170] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        
        ##############################################################

        model = ReadFile(path_MR, step[2], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[2, 1], title = L"B) $\theta_\textrm{ini} = 45^\circ$ - $t$ = %$(tMy_string) Ma", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[2, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)


        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 350]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 350]) ./ ε̇1.x[:, 350] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        

        # ##############################################################

        model = ReadFile(path_HR, step[3], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[3, 1], title = L"C) $\theta_\textrm{ini} = 80^\circ$ - $\tau_\textrm{II}$ at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        @show size(xc), size(zc), size(log10.(τII))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[3, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

                # lines!(ax1, xc./Lc, log10.(ε̇II[:, 700]) )
                # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 700]) ./ ε̇1.x[:, 700] )
                # xlims!(ax1, window.xmin, window.xmax)
                # ylims!(ax1, -180, 180)
                

        ##############################################################
    
        # save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/CompareTheta.png", f, px_per_unit = 4)     
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