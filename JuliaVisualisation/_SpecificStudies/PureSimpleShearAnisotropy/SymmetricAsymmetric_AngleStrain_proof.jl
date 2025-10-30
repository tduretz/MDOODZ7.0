# import Pkg
# Pkg.activate(normpath(joinpath(@__DIR__, "../..")))
using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, UnPack, GeometryTypes, LaTeXStrings, MakieTeX
using CairoMakie#, GLMakie
Mak = CairoMakie
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
# fontsize_theme = Theme(fontsize=200)
# set_theme!(fontsize_theme)
const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

function Moho_LAB_offset(all_phases, T, coords, height)

    @show v,indx = findmin(height)
    @show x_moho = coords.v.x[indx] / 1e3
    𝐗_moho = [x_moho, height[indx]]

    #    # Find x position of highest moho
    #    z_moho = zeros(size(all_phases.ph, 1))
    #    for ix in axes(all_phases.ph, 1)
    #        indz = findfirst(all_phases.ph[ix,:] .== 0)
    #        z_moho[ix] = coords.c.z[indz] / 1e3
    #    end
    #    v,indx = findmax(z_moho)
   
    #    @show v, indx
    #    x_moho = coords.c.x[indx] / 1e3
    #    𝐗_moho = [x_moho, z_moho[indx]]

    # # Find x position of highest moho
    # z_moho = zeros(size(all_phases.ph_hr, 1))
    # for ix in axes(all_phases.ph_hr, 1)
    #     indz = findfirst(all_phases.ph_hr[ix,:] .== 0)
    #     z_moho[ix] = coords.c_hr.z[indz] / 1e3
    # end
    # v,indx = findmax(z_moho)

    # @show v, indx
    # @show z_moho[indx-1], z_moho[indx], @show z_moho[indx+1]
   
    # @show  sum(z_moho .== z_moho[indx])
    # x_moho = coords.c_hr.x[indx] / 1e3
    # 𝐗_moho = [x_moho, z_moho[indx]] 

    # # Find x position of highest T
    # z_lab = zeros(size(T, 1))
    # for ix in axes(T, 1)
    #     indz = findfirst(T[ix,:] .< 50)
    #     z_lab[ix] = coords.c.z[indz] / 1e3
    # end
    # _,indx = findmax(z_lab)
    # x_lab =  coords.c.x[indx] / 1e3
    # 𝐗_moho = [x_lab, z_lab[indx]]

    # Find x position of highest T
    z_lab = zeros(size(T, 1))
    for ix in axes(T, 1)
        indz = findfirst(T[ix,:] .< 1300)
        z_lab[ix] = coords.c.z[indz] / 1e3
    end
    _,indx = findmax(z_lab)
    x_lab =  coords.c.x[indx] / 1e3
    𝐗_lab = [x_lab, z_lab[indx]]

    return 𝐗_moho, 𝐗_lab 
end

@views function main()

    # Set the path to your files
    path_LR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1/"
    path_LR1 ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t0/"


    path_MR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1_5/"
    path_HR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d6/"

    path_MR1 ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t80/"
    path_HR1 ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t5/"


    # File numbers
    step = [1000; 1200; 1200]
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
        ph_contours   = false,  # add phase contours
        fabric        = false,  # add fabric quiver (normal to director)
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
        vel_step    = 20,
        nap         = 0.1,    # pause for animation 
        resol       = 1000,
        Lx          = 1.0,
        Lz          = 1.0,
        LAB_color   = true,
        LAB_T       = 1200,
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
    @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, Fab, height, all_phases, Δ = model

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
    ftsz =  18*options.resol/500
    f = Figure(size = (2*0.85*Lx/4/Lz*options.resol*1.2, options.resol), fontsize=ftsz,  fonts = (; regular = "helvetica") )

    if field==:StrainRate

        # Switches
        PlotOnTop = (
            ph_contours   = true,  # add phase contours
            fabric        = false,  # add fabric quiver (normal to director)
            T_contours    = true,   # add temperature contours
            topo          = false,
            quiver_origin = false,
            σ1_axis       = false,
            ε̇1_axis       = false,
            vel_vec       = false,
            ϕ_contours    = false,
            PT_window     = false,
        )

        cmap = :amp

        model = ReadFile(path_LR, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model

        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[1, 1], title = L"A) Isotropic model (δ$ = 1.0$, \text{\theta}$_\text{ini} = 10^\text{o}$) -  $\mathrm{\varepsilon}$ = %$(stretch_string) %", ylabel = L"$y$ (%$(length_unit))", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        ax1 = Axis(f[1, 1], title = L="A) Isotropic model (δ = 1.0, θini = 10ᵒ) -  ε = 19.7 %", ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))

        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = Reverse(:amp), colorrange=(-1, 1))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))


        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (cmap, options.α_heatmap), colorrange=(-1, 1))
        # hm = heatmap!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, all_phases.ph_dual_hr, colormap = :turbo)

        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ/2, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)

        # text!(ax1, window.zmin, window.zmax)


        # Colorbar(f, hm, label = L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 30,
        # labelsize = 40, ticklabelsize = 30, bbox=BBox(1660, 1850, 690, 980),
        # alignmode = Outside(20), halign = :right, ticklabelcolor = :black, labelcolor = :black,
        # tickcolor = :black)

        # colsize!(f.layout, 2, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[1, 3], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.Colorbar(f[1, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 50)

        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)


        # Wite box with transparency
        x1 = [65.; 150.]
        z1 = [-100; -48.]
        heatmap!(ax1, x1, z1, ones(length(x1), length(z1)), alpha= 0.7, colormap=:broc)
        # poly!(ax1, Rectangle{Float32}(125, 25, 50, 50),    color = :gray, strokewidth = 2, strokecolor = :red)
        Colorbar(f, hm, label = "log(εΙΙ) ", width = 320, height = 30,
        labelsize = 40, ticklabelsize = 30, bbox=BBox(700, 1850, 690, 900),
        alignmode = Outside(20), halign = :left, ticklabelcolor = :black, labelcolor = :black, vertical=false,
        tickcolor = :black)

        hidexdecorations!(ax1)


        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 170]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 170]) ./ ε̇1.x[:, 170] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)

        ##############################################################

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

        model = ReadFile(path_LR1, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model

        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[1, 2], title = L"B) Anisotropic model ($\delta = 4$, $\theta_\text{ini} = 0^\text{o}$) -  $\varepsilon$ = %$(stretch_string) %", ylabel = L"$y$ (%$(length_unit))", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        ax1 = Axis(f[1, 2], title = L="B) Anisotropic model: δ = 4.0 (θini = 0ᵒ) -  ε = 18.2 %", ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))

        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))
        
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)


        # # Wite box with transarency
        # x1 = [65.; 150.]
        # z1 = [-100; -48.]
        # heatmap!(ax1, x1, z1, ones(length(x1), length(z1)), alpha= 0.9, colormap=:broc)
        # # poly!(ax1, Rectangle{Float32}(125, 25, 50, 50),    color = :gray, strokewidth = 2, strokecolor = :red)
        # Colorbar(f, hm, label = L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 320, height = 30,
        # labelsize = 40, ticklabelsize = 30, bbox=BBox(1660+50, 1850, 690, 900),
        # alignmode = Outside(20), halign = :left, ticklabelcolor = :black, labelcolor = :black, vertical=false,
        # tickcolor = :black)


        # colsize!(f.layout, 2, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[1, 3], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.Colorbar(f[1, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 50)

        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

        hidexdecorations!(ax1)
        hideydecorations!(ax1)


        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 170]) )
        # lines!(ax1, xc./Lc, atand.(ε̇1.z[:, 170]) ./ ε̇1.x[:, 170] )
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, -180, 180)
        
        ##############################################################

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

        model = ReadFile(path_MR, step[2], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model
        
        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[2, 1], title = L"C) Anisotropic model: $\delta = 1.5$ ($\theta_\text{ini} = 10^\text{o}$) - $\varepsilon$ = %$(stretch_string) %", ylabel = L"$y$ (%$(length_unit))", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        ax1 = Axis(f[2, 1], title = L="C) Anisotropic model: δ = 1.5 (θini = 10ᵒ) -  ε = 22.6 %", ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))

        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (cmap, options.α_heatmap), colorrange=(-1, 1))
        # hm = heatmap!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, all_phases.ph_dual_hr, colormap = :turbo)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = Reverse(:lapaz), colorrange=(-1, 1))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))

        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # # Mak.Colorbar(f[2, 2], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.Colorbar(f[2, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)
        hidexdecorations!(ax1)


        # ##############################################################

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

        model = ReadFile(path_HR, step[3], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model
        
        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[2, 2], title = L"D) Anisotropic model: $\delta = 6.0$ ($\theta_\text{ini} = 10^\text{o}$) - $\varepsilon$ = %$(stretch_string) %", ylabel = L"$y$ (%$(length_unit))", xlabel = L"$x$ (%$(length_unit))", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        ax1 = Axis(f[2, 2], title = L="D) Anisotropic model: δ = 6.0 (θini = 10ᵒ) -  ε = 22.0 %", ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))

        @show size(xc), size(zc), size(log10.(τII))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (cmap, options.α_heatmap), colorrange=(-1, 1))
        # hm = heatmap!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, all_phases.ph_dual_hr, colormap = :turbo)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = Reverse(:lapaz), colorrange=(-1, 1))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))

        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)

        # # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # # Mak.Colorbar(f[3, 2], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)
        hideydecorations!(ax1)
        hidexdecorations!(ax1)


        ##############################################################

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

        model = ReadFile(path_MR1, step[2], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model
        
        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[3, 1], title = L"E) Anisotropic model: $\theta_\text{ini} = 80^\text{o}$ ($\delta = 4.0$) - $\varepsilon$ = %$(stretch_string) %", xgridvisible = false, ygridvisible = false, aspect=DataAspect(), xlabel = L"$x$ (%$(length_unit)) ", ylabel = L"$y$ (%$(length_unit))")
        ax1 = Axis(f[3, 1], title = L="E) Anisotropic model: θini = 80ᵒ (δ = 4.0) -  ε = 18.7 %", xlabel = "x (km)", ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (cmap, options.α_heatmap), colorrange=(-1, 1))
        # hm = heatmap!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, all_phases.ph_dual_hr, colormap = :turbo)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = Reverse(:lapaz), colorrange=(-1, 1))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))

        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        
        # # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # # Mak.Colorbar(f[2, 3], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.Colorbar(f[2, 2], hm, label =  L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)



        # ##############################################################

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

        model = ReadFile(path_HR1, step[3], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, εII, C, Δ = model
        
        𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)

        stretch = tMy * (1e-15*3600*24*365.25*1e6) *100
        stretch_string = @sprintf("%1.2lf", stretch)
        # ax1 = Axis(f[3, 2], title = L"F) Anisotropic model: $\theta_\text{ini} = 5^\text{o}$ ($\delta = 4.0$) - $\varepsilon$ = %$(stretch_string) %", xlabel = L"$x$ [%$(length_unit)]", xgridvisible = false, ygridvisible = false, aspect=DataAspect())
        ax1 = Axis(f[3, 2], title = L="F) Anisotropic model: θini = 5ᵒ (δ = 4.0) -  ε = 20.8 %", xlabel = "x (km)",  ylabel = "y (km)", titlefont = :regular, xgridvisible = false, ygridvisible = false, aspect=DataAspect(), yticks = ([0, -25, -50], ["0", "25", "50"]))
        @show size(xc), size(zc), size(log10.(τII))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(τII), colormap = (cmap, options.α_heatmap), colorrange=(6.69, 8.69))
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (cmap, options.α_heatmap), colorrange=(-1, 1))
        # hm = heatmap!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, all_phases.ph_dual_hr, colormap = :turbo)
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = Reverse(:lapaz), colorrange=(-1, 1))
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = cmap, colorrange=(-1, 1))

        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases.group_phases, Δ, Mak)                
        
        # lines!(𝐗_moho[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)
        # lines!( 𝐗_lab[1].*ones(size(zc)), zc./Lc, color=:black, linewidth=3)

        # # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # # Mak.Colorbar(f[3, 3], hm, label =  L"$\log_{10}$ $\varepsilon_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # # Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)
        hideydecorations!(ax1)
                

        ##############################################################
    
        save("/Users/tduretz/PowerFolders/_manuscripts/PureSimpleShearAnisotropy/Figures/SymmetricAsymmetric_AngleStrain_proof.png", f, px_per_unit = 4)     end

    display(f)
end

# function PrincipalStress(τxx, τzz, τxz, P)
#     σ1   = (x=zeros(size(τxx)), z=zeros(size(τxx)) )
#     τxzc = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end-0,1:end-1] .+ τxz[1:end-1,2:end-0] .+ τxz[2:end-0,2:end-0]) 
#     for i in eachindex(τxzc)
#         if P[i]>-1e-13
#             σ = [-P[i]+τxx[i] τxzc[i]; τxzc[i] -P[i]+τzz[i]]
#             v = eigvecs(σ)
#             σ1.x[i] = v[1,1]
#             σ1.z[i] = v[2,1]
#         end
#     end
#     return σ1
# end

# function ExtractData( file_path, data_path)
#     data = h5open(file_path, "r") do file
#         read(file, data_path)
#     end
#     return data
# end

# function Print2Disk( f, path, field, istep; res=4)
#      path1 = path*"/_$field/"
#      mkpath(path1)
#      @show name = path1*"$field"*@sprintf("%05d", istep)*".png"
#      save(name, f, px_per_unit = res) 
# end

main()