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

include("ChristmasTreeAnisotropy_AsFunction.jl")

@views function main()

    # Set the path to your files
    path  ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_HR/"
    
    # File numbers
    step = [500; 1100; 2000]
    step = [00; 1400; 2900]

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
    # zoom = ( 
    #     xmin = -150, 
    #     xmax = 150,
    #     zmin = -70,
    #     zmax = 2.0,
    # )

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
    cmap[1] = RGBA{Float64}(213/255, 185/255, 148/255, 1.)  
    cmap[2] = RGBA{Float64}(1/255, 1/255, 1/255, 1.)  
    cmap[3] = RGBA{Float64}(190/255, 216/255, 172/255, 1.) 
    cmap[4] = RGBA{Float64}(156/255, 205/255, 144/255, 1.) 
    cmap[5] = RGBA{Float64}(190/255, 216/255, 172/255, 1.) 
    cmap[6] = RGBA{Float64}(255/255, 255/255, 255/255, 1.) 
    cmap[7] = RGBA{Float64}(108/255, 183/255, 166/255, 1.) 
    phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)
      
    #####################################

    @unpack Lc, tc, Vc, τc = scales

    model = ReadFile(path, step[1], scales, options, PlotOnTop)
    @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, Fab, height, all_phases, Δ = model

    zoom = ( 
            xmin = -60., 
            xmax = 60.,
            zmin = minimum(zc)/scales.Lc,
            zmax = maximum(zc)/scales.Lc,
        )

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
    ftsz =  24*options.resol/500
    f = Figure(size = (2.1*options.resol, options.resol), fontsize=ftsz)

    cmap = :vik

    if field==:StrainRate

        model = ReadFile(path, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, Δ = model

        tMy_string = @sprintf("%1.2lf", tMy)
        ax1 = Axis(f[1, 1], title = L"$$A) Initial model configuration", ylabel = L"$y$ [%$(length_unit)]", xlabel = L"$x$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[1, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )

        Colorbar(f, hm, label = L"$\tau_\textrm{II}$ [Pa]", width = 320, height = 30,
        labelsize = 40, ticklabelsize = 30, bbox=BBox(400, 1250, -450, 900),
        alignmode = Outside(20), halign = :left, ticklabelcolor = :black, labelcolor = :black, vertical=false,
        tickcolor = :black)
        xb = [-20.; 60.]
        yb = [-150.; -120.]
        heatmap!(ax1, xb, yb, zeros(2,2), color=(:broc, 0.5))

        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

        Tprofile = T[1,:]
        yprofile =  zc./Lc

        τii_rot1, yc, τii_cart1, τii_yield, τii_yield_min, T, T0, P = ChrismasTreeAniso()

        #####################

        ax2 = Axis(f[1, 2], xlabel = L"B) $T$ [$^\circ$C]", ylabel = L"$z$ [km]",  xticklabelcolor = :tomato, xaxisposition=:top, yticklabelsvisible=false) # , title = L"$$Pressure profile"
        ax1 = Axis(f[1, 2], xlabel = L"$P$ [GPa]", xticklabelcolor = :blue,  xaxisposition=:bottom, yticklabelsvisible=false) # , title = L"$$Temperature profile"
        hidespines!(ax2)
        hideydecorations!(ax2)
        
        lines!(ax2, T .- T0, -35*ones(size(T)), color=:gray, linestyle=:dash, linewidth=3)
        lines!(ax1, P./1e9, yc./1e3, color=:blue, linewidth=4 )
        lines!(ax2, Tprofile, yprofile, color=:tomato, linewidth=4)

        # lines!(ax2, T .- T0, yc./1e3, color=:blue ) 
        lines!(ax1, P./1e9, -125*ones(size(P)), color=:gray, linestyle=:dash, linewidth=3)

        ylims!(ax2,-150, 10)
        ylims!(ax1,-150, 10)

        #####################
        Mak.colgap!(f.layout, 10)

        ax3 = Axis(f[1, 3], title = L"C) $τ_\mathrm{II}$ ($\delta = 3$)", xlabel = L"$τ_\mathrm{II}$ [GPa]", yticklabelsvisible=false)
        lines!(ax3, τii_rot1./1e9, yc./1e3, label=L"$\sqrt{Y_2^{'}}$", color=:blue, linewidth=4 )
        xlims!(ax3, 0, .7)
        # lines!(ax3, τii_rot2./1e9, yc./1e3 )
        lines!(ax3, τii_cart1./1e9, yc./1e3, label=L"$\sqrt{J_2}$", color=:tomato, linewidth=4 )

        lines!(ax3, τii_yield./1e9,     yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = 0.0$)", linestyle=:dash, color=:blue, linewidth=4 )
        lines!(ax3, τii_yield_min./1e9, yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = π/4$)", linestyle=:dash, color=:tomato, linewidth=4 )

        ylims!(ax3,-150, 10)

        # lines!(ax3, τii_cart2./1e9, yc./1e3 )
        axislegend(position=:rb)
        display(f)

        Mak.colgap!(f.layout, 50)





        # ax3 = Axis(f[1, 2], title = L"B) $$Stress profile ($\delta = 3$)", xlabel = L"$τ_\mathrm{II}$ [GPa]")
        # lines!(ax3, τii_rot1./1e9, yc./1e3, label=L"$\sqrt{Y_2^{'}}$", color=:blue )
        # xlims!(ax3, 0, .7)
        # # lines!(ax3, τii_rot2./1e9, yc./1e3 )
        # lines!(ax3, τii_cart1./1e9, yc./1e3, label=L"$\sqrt{J_2}$", color=:tomato )

        # lines!(ax3, τii_yield./1e9,     yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = 0.0$)", linestyle=:dash, color=:blue )
        # lines!(ax3, τii_yield_min./1e9, yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = π/4$)", linestyle=:dash, color=:tomato )

        # # lines!(ax3, τii_cart2./1e9, yc./1e3 )
        # axislegend(position=:rb)
        # display(f)

        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 170]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 170]) ./ ε̇1.x[:, 170] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        
        ##############################################################

        # model = ReadFile(path, step[2], scales, options, PlotOnTop)
        # @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        # tMy_string = @sprintf("%1.2lf", tMy)
        # ax1 = Axis(f[2, 1], title = L"B) $\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy_string) Ma", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (cmap, options.α_heatmap), colorrange=(-17, -13))
        # AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, 1.5Δ, Mak)                
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[2, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # Mak.colgap!(f.layout, 20)
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, window.zmin, window.zmax)


        # # lines!(ax1, xc./Lc, log10.(ε̇II[:, 350]) )
        # # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 350]) ./ ε̇1.x[:, 350] )
        # # xlims!(ax1, window.xmin, window.xmax)
        # # ylims!(ax1, -180, 180)
        

        # # ##############################################################

        # model = ReadFile(path, step[3], scales, options, PlotOnTop)
        # @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, group_phases, Δ = model
        
        # tMy_string = @sprintf("%1.2lf", tMy)
        # ax1 = Axis(f[3, 1], title = L"C) $\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false,)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (cmap, options.α_heatmap), colorrange=(-17, -13))
        # AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, 1.5Δ, Mak)                
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[3, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # Mak.colgap!(f.layout, 20)
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, window.zmin, window.zmax)

        #         # lines!(ax1, xc./Lc, log10.(ε̇II[:, 700]) )
        #         # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 700]) ./ ε̇1.x[:, 700] )
        #         # xlims!(ax1, window.xmin, window.xmax)
        #         # ylims!(ax1, -180, 180)
                

        ##############################################################
    
        save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/InitialConfiguration.png", f, px_per_unit = 4)   
    end

    display(f)
end


main()