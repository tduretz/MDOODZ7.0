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
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"

    # File numbers
    file_start = 0
    file_step  = 10
    file_end   = 0

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
    field = :X
    # field = :TimeSeries
    # field = :EffectiveFrictionTime
    # field = :ChristmasTree

    # Define Tuple for enlargment window
    # zoom = ( 
    #     xmin = -500, 
    #     xmax = 500,
    #     zmin = -300,
    #     zmax = 5,
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
    resol       = 500
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0
    LAB_color   = true
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
    f = Figure(size = (Lx/Lz*resol*1.2, resol), fontsize=40)

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

        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        xc = LinRange(xmin+Δx/2, xmax-Δx/2, ncx)
        zc = LinRange(zmin+Δz/2, zmax-Δz/2, ncz)
        xv = LinRange(xmin, xmax, nvx)
        zv = LinRange(zmin, zmax, nvz)
        coords = ( c=(x=xc, z=zc), c_hr=(x=xc_hr, z=zc_hr), v=(x=xv, z=zv))

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
        ftsz = 30*resol/500
        ftsz = 18

        f = Figure(size = (1.1*Lx/Lz*resol*1.2, resol), fontsize=ftsz)

        if field==:Phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = :turbo)
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_dual_hr, colormap = :turbo)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
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
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:X

            k_chem  = 5e-10
            σ       = 0.1
            # X_exact = exp.( -(xc.^2 .+ (zc').^2)/(σ^2) )

            X_exact = 1/(1+4*t*k_chem/σ^2) * exp.( -(xc.^2 .+ (zc').^2)/(σ^2 + 4*t*k_chem) )
            # @show size(X_exact), size(X)


            ax2 = Axis(f[1, 1], title = L"$X$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            lines!(ax2, xc,       X[:, Int64(floor(size(X,2)/2))])
            lines!(ax2, xc, X_exact[:, Int64(floor(size(X,2)/2))])

            ax1 = Axis(f[2, 1], title = L"$X$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, X, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            Mak.Colorbar(f[2, 2], hm, label = L"$X$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
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
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        display(f)

    end

    yscale = Lz/Lx

    if printfig && printvid
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale='bitand(oh*dar, 65534)':'bitand(ih/2, 65534)', setsar=1" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)

        # FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1080:1080*$(yscale)" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)
    end

end

main()