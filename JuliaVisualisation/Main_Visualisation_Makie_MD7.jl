import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
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
    path ="/home/larafriedrichs/repositories/MDOODZ7.0/MDLIB/"
    #path=raw"C:\Users\49176\OneDrive\Desktop\Test_c_code\\"
    path="/home/larafriedrichs/repositories/MDOODZ7.0/runs/firstmodel/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    path = "/Users/lcandiot/Developer/MDOODZ7.0/cmake-exec/ThanushikaSubduction/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1/"
    path ="/Users/tduretz/Downloads/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingMelting/"


    # File numbers
    file_start = 1000
    file_step  = 100
    file_end   = 1000

    # Select field to visualise
    field = :Phases
    # field = :Cohesion
    # field = :Density
    # field = :Viscosity  
    # field = :PlasticStrainrate
    # field = :Stress
    # field = :σxx
    # field = :σzz
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
    # field = :X
    # field = :TimeSeries
    # field = :EffectiveFrictionTime
    # field = :ChristmasTree
    field = :GPE

    # Define Tuple for enlargment window
    # zoom = ( 
    #     xmin = -200e3, 
    #     xmax = 200e3,
    #     zmin = -5e3,
    #     zmax = 1e3,
    # )

    # Switches
    printfig    = false  # print figures to disk
    printvid    = false
    framerate   = 12
    PlotOnTop = (
        ph_contours   = false,  # add phase contours
        fabric        = false,   # add fabric quiver (normal to director)
        T_contours    = false,  # add temperature contours
        topo          = false,
        quiver_origin = false,
        σ1_axis       = false,
        ε̇1_axis       = false,
        vel_vec       = false,
        ϕ_contours    = false,
        PT_window     = false,
    )
    α_heatmap   = 1.0   # transparency of heatmap 
    vel_arrow   = 5
    vel_scale   = 0.00001
    vel_step    = 10
    nap         = 0.1    # pause for animation 
    resol       = 250
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0
    LAB_color   = false
    LAB_T       = 1250

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    Lc = 1
    tc = My
    Vc = 1.0
    τc = 1

    probe = (ϕeff = Float64.([]), t  = Float64.([]))
    cm_yr = 100.0*3600.0*24.0*365.25

    # Time loop
    f = Figure(size = (Lx/Lz*resol, resol), fontsize=40, figure_padding=10)

    for istep=file_start:file_step:file_end
    
        name     = @sprintf("Output%05d.gzip.h5", istep)
        @info "Reading $(name)"
        filename = string(path, name)
        model    = ExtractData( filename, "/Model/Params")
        xc       = ExtractData( filename, "/Model/xc_coord")
        zc       = ExtractData( filename, "/Model/zc_coord")
        xv       = ExtractData( filename, "/Model/xg_coord")
        zv       = ExtractData( filename, "/Model/zg_coord")
        xvz      = ExtractData( filename, "/Model/xvz_coord")
        zvx      = ExtractData( filename, "/Model/zvx_coord")
        xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr")
        τzz_t    = ExtractData( filename, "TimeSeries/szzd_mean_time")
        P_t      = ExtractData( filename, "TimeSeries/P_mean_time")
        τxz_t    = ExtractData( filename, "TimeSeries/sxz_mean_time")
        t_t      = ExtractData( filename, "TimeSeries/Time_time")

        xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)
        coords = ( c=(x=xc, z=zc), c_hr=(x=xc_hr, z=zc_hr), v=(x=xv, z=zv))

        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        Lx>1e3 ? length_unit="km" :  length_unit="m"

        if @isdefined zoom
            Lx = (zoom.xmax - zoom.xmin)./Lc 
            Lz = (zoom.zmax - zoom.zmin)./Lc
            window = zoom
        else
            window = ( 
                xmin = minimum(xv)/Lc, 
                xmax = maximum(xv)/Lc,
                zmin = minimum(zv)/Lc,
                zmax = maximum(zv)/Lc,
            )
        end

        @info "Model info"
        @show "Model apect ratio" Lx/Lz
        @show "Model time" t/My
        @show length_unit
        centroids, vertices = (ncx, ncz), (nvx, nvz)

        ph           = ExtractField(filename,  "/VizGrid/compo", centroids, false, 0)
        mask_air     = ph .== -1.00 
        ph_hr        = ExtractField(filename,  "/VizGrid/compo_hr", (ncx_hr, ncz_hr), false, 0)
        ph_dual_hr   = ExtractField(filename,  "/VizGrid/compo_dual_hr", (ncx_hr, ncz_hr), false, 0)
        group_phases = copy(ph_hr); ph_hr[ph_hr.==-1.00] .=   NaN; ph_dual_hr[ph_dual_hr.==-1.00] .=   NaN;
        ηc           = ExtractField(filename,  "/Centers/eta_n", centroids, true, mask_air)
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz))
        τyy   = -(τzz .+ τxx)
        σzz   = -P + τzz
        σxx   = -P + τxx
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇zz   = Float64.(reshape(ExtractData( filename, "/Centers/ezzd"), ncx, ncz))
        ε̇yy   = -(ε̇xx .+ ε̇zz)
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        τII   = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
            
        τxzc  = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end,1:end-1] .+ τxz[1:end-1,2:end] .+ τxz[2:end,2:end]) 
        C     = Float64.(reshape(ExtractData( filename, "/Centers/cohesion"), ncx, ncz))
        ϕ     = ExtractField(filename, "/Centers/phi", centroids, false, 0)
        X     = ExtractField(filename, "/Centers/X", centroids, false, 0)
        divu  = ExtractField(filename, "/Centers/divu", centroids, false, 0)
       
        T_hr  = zeros(size(ph_hr)); T_hr[1:2:end-1,1:2:end-1] .= T; T_hr[2:2:end-0,2:2:end-0] .= T
        if LAB_color
            ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.<LAB_T] .= 2
            ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.>LAB_T] .= 3
            ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.<LAB_T] .= 2
            ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.>LAB_T] .= 3
            ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.<LAB_T] .= 2
            ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.>LAB_T] .= 3
        end

        Fab = 0.
        if PlotOnTop.fabric
            δani    = ExtractField(filename, "/Centers/ani_fac", centroids, false, 0)
            Nx      = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
            Nz      = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
            Fab     = (x=-Nz./Nx, z=ones(size(Nz)))
            nrm     = sqrt.(Fab.x.^2 .+ Fab.z.^2)
            Fab.x ./= nrm
            Fab.z ./= nrm
        end
        height = 0.
        if PlotOnTop.topo
            height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
            Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
            Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
            Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
            Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
            x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
            z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
        end
        Vxc   = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc   = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])
        V = (x=Vxc, z=Vzc, arrow=vel_arrow, step=vel_step, scale=vel_scale)
        σ1, ε̇1 = 0., 0.
        if PlotOnTop.σ1_axis 
            σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        end
        if PlotOnTop.ε̇1_axis 
            ε̇1 = PrincipalStress(ε̇xx, ε̇zz, ε̇xz, zeros(size(ε̇xx))) 
        end
        PT = (P.>2.2e9 .&& P.<3.0e9 .&& T.>430 .&& T.<530).*ones(size(T))
        @show extrema(P)
        @show extrema(τII)
        @show extrema(ε̇II)
     
        
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

        # # Group phases for contouring
        # group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        # group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        # group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################
        empty!(f)
        ftsz =  30*resol/500
        f = Figure(size = (Lx/Lz*resol, resol), fontsize=ftsz)

        if field==:Phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]", aspect = Lx/Lz)
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = :turbo)
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_dual_hr, colormap = :turbo)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = "Phases", labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Density
            ax1 = Axis(f[1, 1], title = L"$\rho$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap), colorrange=(2600, 5000))  
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\rho$ [kg.m$^{-3}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:σxx
            ax1 = Axis(f[1, 1], title = L"$\sigma_\textrm{xx}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, σxx./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\sigma_\textrm{xx}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:σzz
            ax1 = Axis(f[1, 1], title = L"$\sigma_\textrm{zz}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, σzz./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\sigma_\textrm{zz}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Pressure
            ax1 = Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:turbo, α_heatmap)) #, colorrange=(1,1.2)1e4*365*24*3600
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label =  L"$P$ [GPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Divergence
            ax1 = Axis(f[1, 1], title = L"∇⋅V at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, divu, colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end #if printfig Print2Disk( f, path, string(field),ε̇BG) end
        end

        if field==:PlasticStrainrate
            ε̇pl[ε̇pl.==0.0] .= 1e-30
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity
            ax1 = Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            Vx_BG = 0*xc .- 2*zc'  
            V     = sqrt.( (Vxc .- 0.0*Vx_BG).^2 + (Vzc).^2)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, V, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            # xlims!(ax1, 0., 3.e-3)
            # ylims!(ax1, 0., 3.e-3)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$V$ [m.s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_x
            ax1 = Axis(f[1, 1], title = L"$Vx$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xv, zvx[2:end-1], Vx[2:end-1,:], colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vx$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_z
            ax1 = Axis(f[1, 1], title = L"$Vz$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xvz[2:end-1], zv, Vz[:,2:end-1]*cm_yr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vz$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:GrainSize
            ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            Mak.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:AnisotropyFactor
            ax1 = Axis(f[1, 1], title = L"$δ_\textrm{ani}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, δani, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            # xminz, xmaxz = -0.4, 0.4
            # zminz, zmaxz = -0.17, 0.17
            # Lx = xmaxz - xminz
            # Lz = zmaxz - zminz
            # xlims!(ax1, -0.4, 0.4)
            # ylims!(ax1, -0.17, 0.17)
            Mak.Colorbar(f[1, 2], hm, label = L"$δ_\textrm{ani}$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:MeltFraction
            ax1 = Axis(f[1, 1], title = L"ϕ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ϕ, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$ϕ$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:X
            ax1 = Axis(f[1, 1], title = L"$X$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, X, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$ϕ$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Cohesion
            ax1 = Axis(f[1, 1], title = L"$C$ [MPa] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, C./1e6, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$C$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Temperature
            ax1 = Axis(f[1, 1], title = L"$T$ [C] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (Reverse(:bilbao), α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$T$ [C]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            rowsize!(f.layout, 1, Aspect(1, Lz/Lx))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Topography
            height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
            Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
            Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
            Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
            Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
            x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
            z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));

            ax1 = Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$h$ [km]")
            @show mean(height)
            @show mean(z_mark)
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)
            xlims!(ax1, window.xmin, window.xmax)

            ax2 = Axis(f[2, 1], xlabel = L"$x$ [km]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)
            xlims!(ax1, window.xmin, window.xmax)

            ax3 = Axis(f[3, 1], xlabel = L"$x$ [km]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)
            xlims!(ax1, window.xmin, window.xmax)
        end

        if field==:TimeSeries
            ax1 = Axis(f[1, 1], title = L"$τ_{xz}$", xlabel = L"$t$", ylabel = L"$\tau_{xz}$")
            τxz_t[1] = 0.
            @show .-τxz_t
            lines!(ax1, t_t, .-τxz_t./(.-P_t.+τzz_t))
        end

        if field==:EffectiveFrictionTime
            ϕ_eff = -mean(τxzc[:,end] ./ (τzz[:,end] .- P[:,end]) )
            if istep==0 ϕ_eff = 0. end
            push!(probe.ϕeff, ϕ_eff) 
            push!(probe.t, t) 
            if istep==file_end
                ax1 = Axis(f[1, 1], title = L"$ϕ_\mathrm{eff}$", xlabel = L"$t$", ylabel = L"$ϕ_\mathrm{eff}$")
                lines!(ax1, Float64.(probe.t), Float64.(probe.ϕeff))
            end
        end

        if field==:ChristmasTree
            ax1 = Axis(f[1, 1], title = L"Stress profile at $t$ = %$(tMy) Ma", xlabel = L"$τII$ [MPa]", ylabel = L"$z$ [km]")
            lines!(ax1, mean(τII, dims=1)[:]/τc, coords.c.z./Lc/1e3, )
            lines!(ax1, mean(T, dims=1)[:], coords.c.z./Lc/1e3, )
            ax2 = Axis(f[1, 2], title = L"Temperature profile at $t$ = %$(tMy) Ma", xlabel = L"$T$ [C]", ylabel = L"$h$ [km]")
            lines!(ax2, mean(T, dims=1)[:], coords.c.z./Lc/1e3, )
        end

        if field==:GPE

            f = Figure(size = (Lx/Lz/1.8*resol*2.5, 2.5*resol), fontsize=ftsz)

            # Here we load data without NaNs
            NaN_air = false
            τxx  = ExtractField(filename,  "/Centers/sxxd",   centroids, NaN_air, mask_air)
            P    = ExtractField(filename,  "/Centers/P",      centroids, NaN_air, mask_air)
            T    = ExtractField(filename,  "/Centers/T",      centroids, NaN_air, mask_air)
            ρc   = ExtractField(filename,  "/Centers/rho_n",  centroids, NaN_air, mask_air)
            ρv   = ExtractField(filename,  "/Vertices/rho_s", vertices,  NaN_air, mask_air)
            tagv = ExtractField(filename,  "/Flags/tag_s",    vertices,  NaN_air, mask_air)
            tagc = ExtractField(filename,  "/Flags/tag_n",    centroids, NaN_air, mask_air)

            ρc   = 0.25(ρv[1:end-1,1:end-1] .+ ρv[2:end,1:end-1] .+ ρv[1:end-1,2:end] .+ ρv[2:end,2:end])

            height  = Float64.(ExtractData( filename, "/Topo/z_grid"));

            # 
            for I in CartesianIndices(ρc)
                i, j = I[1], I[2]
                ρ = 0.0
                n = 0
                if tagv[i,j] != 30
                    ρ += ρv[i,j]
                    n += 1
                end
                if tagv[i+1,j] != 30
                    ρ += ρv[i+1,j]
                    n += 1
                end
                if tagv[i,j+1] != 30
                    ρ += ρv[i,j+1]
                    n += 1
                end
                if tagv[i+1,j+1] != 30
                    ρ += ρv[i+1,j+1]
                    n += 1
                end
                if n > 0
                    ρc[i,j] = ρ / n
                else
                    ρc[i,j] = 0
                end
            end

            weight = zeros(size(ρc))

            # 2D loop on centroids
            for I in CartesianIndices(ρc)
                i, j = I[1], I[2]
                if j<size(ρc,2) # avoid to boundary
                    if tagc[i,j] == -1
                        weight[i,j] = 1.0
                    end
                    if tagc[i,j+1] == 31 
                        h = 0.5*(height[i] + height[i+1])  # Height of surface within the cell        
                        weight[i,j] = (h - zc[j]-Δz/2)/Δz  # Volume of rocks versus air                   
                    end
                end
            end

            # Example how to extract surface density 
            ρsurf = zeros(size(ρc,1))
            index = zeros(size(ρc,1))
            for i in axes(ρc,1)
                ind      = findlast(x->x>0, ρc[i,:]) 
                index[i] = ind
                ρsurf[i] = ρc[i,ind]
            end
        
            # Relationships from Schmalholz et al. (2019) 
            gz   = -9.81
            σxx  = - P .+ τxx
            ∫σxx = sum(σxx, dims=2)*Δz
            Pl   = -(gz*reverse(cumsum(reverse(ρc.*weight, dims=2), dims=2))*Δz)
            Po   = P .- Pl
            ∫Po  = sum(Po, dims=2)*Δz
            ∫τxx = sum(τxx, dims=2)*Δz
            Fx   = ∫τxx - ∫Po
            GPE  = sum(Pl, dims=2)*Δz
            ∫σxz = sum(τxz, dims=2)*Δz
            Pb   = -diff(∫σxz[:])/Δx

            # Visualise
            Lc    = 1000

            ax1 = Axis(f[1, 1], title = L"Cell weights at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc/Lc, zc/Lc, weight, colormap = (Reverse(:bilbao), α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$w$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            ylims!(ax1, -5, 5)
 
            ax3 = Axis(f[2, 1], title = L"∫σxx dz at $t$ = %$(tMy) Ma", xlabel = L"$τII$ [MPa]", ylabel = L"$∫σxx$ [TN/m]")
            lines!(ax3, xc/Lc, (∫σxx[:] .- ∫σxx[:][1])./1e12 )
            xlims!(ax3, -200, 200)

            ax2 = Axis(f[3, 1], title = L"Force at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$F$ [TN/m]")
            lines!(ax2, xc/Lc, (Fx[:]  .-  Fx[1] )./1e12)
            lines!(ax2, xc/Lc, (GPE[:] .-  GPE[1])./1e12)
            xlims!(ax2, -200, 200)

            ax4 = Axis(f[4, 1], title = L"Basal overpressure at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$Po$ [MPa]")
            lines!(ax4, xc/Lc, (Po[:,1] )./1e6)
            lines!(ax4, xc/Lc,  Pb ./1e6)
            xlims!(ax4, -200, 200)


            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field!=:EffectiveFrictionTime || istep!=file_end
            DataInspector(f)
            display(f)
            sleep(nap)
        end

    end

    yscale = Lz/Lx

    if printfig && printvid
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale='bitand(oh*dar, 65534)':'bitand(ih/2, 65534)', setsar=1" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)

        # FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1080:1080*$(yscale)" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)
    end

end

main()