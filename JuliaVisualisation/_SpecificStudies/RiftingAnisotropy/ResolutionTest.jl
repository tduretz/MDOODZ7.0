import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../..")))
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, UnPack
using CairoMakie#, GLMakie
Mak = CairoMakie
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
# fontsize_theme = Theme(fontsize=200)
# set_theme!(fontsize_theme)
const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 

function AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, Fab, height, Lc, cm_y, group_phases, Δ)
    if PlotOnTop.T_contours 
        contour!(ax1, coords.c.x./Lc, coords.c.z./Lc, T, levels=0:200:1400, linewidth = 4, color=:black )  
    end  
    if PlotOnTop.ph_contours 
        contour!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
    end
    if PlotOnTop.fabric 
        arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, Fab.x[1:V.step:end,1:V.step:end], Fab.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5)
    end 
    if PlotOnTop.σ1_axis
        arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, σ1.x[1:V.step:end,1:V.step:end], σ1.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5, color=:red)
    end 
    if PlotOnTop.ε̇1_axis
        arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, ε̇1.x[1:V.step:end,1:V.step:end], ε̇1.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5, color=:black)
    end  
    if PlotOnTop.vel_vec
        arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, V.x[1:V.step:end,1:V.step:end]*cm_y, V.z[1:V.step:end,1:V.step:end]*cm_y, arrowsize = V.arrow, lengthscale = V.scale)
    end 
    if PlotOnTop.topo
        lines!(ax1, xv./Lc, height./Lc)
    end
    if PlotOnTop.ϕ_contours 
        contour!(ax1, coords.c.x./Lc, coords.c.z./Lc, ϕ, levels=0:1:0.1, linewidth = 6, color=:white, linestyle= :dashdot )  
    end 
end

function ReadFile(path, step, scales, options, PlotOnTop)

    name     = @sprintf("Output%05d.gzip.h5", step)
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
    tMy    = round(t/scales.tc, digits=6)
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1
    xmin, xmax = xv[1], xv[end]
    zmin, zmax = zv[1], zv[end]
    Lx, Lz    = (xmax-xmin)/scales.Lc, (zv[end]-zv[1])/scales.Lc
    Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
    (xmax-xmin)>1e3 ? length_unit="km" :  length_unit="m"

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
    divu  = ExtractField(filename, "/Centers/divu", centroids, false, 0)
    T_hr  = zeros(size(ph_hr)); T_hr[1:2:end-1,1:2:end-1] .= T; T_hr[2:2:end-0,2:2:end-0] .= T
    if options.LAB_color
        ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.<options.LAB_T] .= 2
        ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.>options.LAB_T] .= 3
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.<options.LAB_T] .= 2
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.>options.LAB_T] .= 3
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.<options.LAB_T] .= 2
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.>options.LAB_T] .= 3
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
    V = (x=Vxc, z=Vzc, arrow=options.vel_arrow, step=options.vel_step, scale=options.vel_scale)
    σ1, ε̇1 = 0., 0.
    if PlotOnTop.σ1_axis 
        σ1 = PrincipalStress(τxx, τzz, τxz, P) 
    end
    if PlotOnTop.ε̇1_axis 
        ε̇1 = PrincipalStress(ε̇xx, ε̇zz, ε̇xz, -1/3*divu) 
        @show extrema(ε̇1.x)

    end

    return (tMy=tMy,length_unit=length_unit, Lx=Lx, Lz=Lz, xc=xc, zc=zc, ε̇II=ε̇II,
    coords=coords, V=V, T=T, ϕ=ϕ, σ1=σ1, ε̇1=ε̇1, Fab=Fab, height=height, group_phases=group_phases, Δ=Δ)
end


@views function main()

    # Set the path to your files
    path_LR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_LR/"
    path_MR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_MR/"
    path_HR ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_HR/"


    # File numbers
    step = [120; 200; 500]
    # step = [500; 1000; 2000]
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
        xmin = 75, 
        xmax = 200,
        zmin = -30,
        zmax = 0.e3,
    )

    # Switches
    PlotOnTop = (
        ph_contours = true,  # add phase contours
        fabric      = false,  # add fabric quiver (normal to director)
        T_contours  = false,   # add temperature contours
        topo        = false,
        σ1_axis     = true,
        ε̇1_axis     = true,
        vel_vec     = false,
        ϕ_contours  = true,
    )
    options = (
        printfig    = false,  # print figures to disk
        printvid    = false,
        framerate   = 6,
        α_heatmap   = 0.5,   # transparency of heatmap 
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
    @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, Fab, height, group_phases, Δ = model

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
    f = Figure(size = (1.1*Lx/4/Lz*options.resol*1.2, options.resol), fontsize=ftsz)

    # if field==:Phases
    #     ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
    #     hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
    #     # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = :turbo)
    #     hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_dual_hr, colormap = :turbo)
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
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
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
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
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
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
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
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
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
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
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    if field==:StrainRate

        options = (
            printfig    = false,  # print figures to disk
            printvid    = false,
            framerate   = 6,
            α_heatmap   = 0.5,   # transparency of heatmap 
            vel_arrow   = 5,
            vel_scale   = 10000,
            vel_step    = 5,
            nap         = 0.1,    # pause for animation 
            resol       = 1000,
            Lx          = 1.0,
            Lz          = 1.0,
            LAB_color   = true,
            LAB_T       = 1250,
        )

        model = ReadFile(path_LR, step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, ε̇1, Fab, height, group_phases, Δ = model

        ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma - 300 $\times$ 200 cells", ylabel = L"$y$ [%$(length_unit)]")
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, options.α_heatmap), colorrange=(-17, -13))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, Fab, height, Lc, cm_y, group_phases, Δ/2)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)
        
        ##############################################################


        options = (
        printfig    = false,  # print figures to disk
        printvid    = false,
        framerate   = 6,
        α_heatmap   = 0.5,   # transparency of heatmap 
        vel_arrow   = 5,
        vel_scale   = 10000,
        vel_step    = 10,
        nap         = 0.1,    # pause for animation 
        resol       = 1000,
        Lx          = 1.0,
        Lz          = 1.0,
        LAB_color   = true,
        LAB_T       = 1250,
    )

        model = ReadFile(path_MR, step[2], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, ε̇1, Fab, height, group_phases, Δ = model
        
        ax1 = Axis(f[2, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma - 600 $\times$ 400 cells", ylabel = L"$y$ [%$(length_unit)]")
        # hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, options.α_heatmap), colorrange=(-17, -13))
        # AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, Fab, height, Lc, cm_y, group_phases, Δ)                
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        # Mak.Colorbar(f[2, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        # Mak.colgap!(f.layout, 20)
        # xlims!(ax1, window.xmin, window.xmax)
        # ylims!(ax1, window.zmin, window.zmax)
        # lines!(ax1, xc./Lc, log10.(ε̇II[:, 350]) )
        lines!(ax1, xc./Lc, atand.(σ1.z[:, 350]) ./ σ1.x[:, 350] )

        xlims!(ax1, window.xmin, window.xmax)
        

        # ##############################################################

        options = (
            printfig    = false,  # print figures to disk
            printvid    = false,
            framerate   = 6,
            α_heatmap   = 0.5,   # transparency of heatmap 
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

        model = ReadFile(path_HR, step[3], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, coords, V, T, ϕ, σ1, ε̇1, Fab, height, group_phases, Δ = model
        
        ax1 = Axis(f[3, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma - 1200 $\times$ 800 cells", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]")
        hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, options.α_heatmap), colorrange=(-17, -13))
        AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, Fab, height, Lc, cm_y, group_phases, Δ*2)                
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[3, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
        Mak.colgap!(f.layout, 20)
        xlims!(ax1, window.xmin, window.xmax)
        ylims!(ax1, window.zmin, window.zmax)

        ##############################################################

        if options.printfig Print2Disk( f, path, string(field), istep) end #if printfig Print2Disk( f, path, string(field),ε̇BG) end
    end

    # if field==:PlasticStrainrate
    #     ε̇pl[ε̇pl.==0.0] .= 1e-30
    #     ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Velocity
    #     ax1 = Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     Vx_BG = 0*xc .- 2*zc'  
    #     V     = sqrt.( (Vxc .- 0.0*Vx_BG).^2 + (Vzc).^2)
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, V, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     # xlims!(ax1, 0., 3.e-3)
    #     # ylims!(ax1, 0., 3.e-3)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$V$ [m.s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Velocity_x
    #     ax1 = Axis(f[1, 1], title = L"$Vx$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xv./Lc, zvx[2:end-1]./Lc, Vx[2:end-1,:]*cm_yr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$Vx$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Velocity_z
    #     ax1 = Axis(f[1, 1], title = L"$Vz$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xvz[2:end-1]./Lc, zv./Lc, Vz[:,2:end-1]*cm_yr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = L"$Vz$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:GrainSize
    #     ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     xminz, xmaxz = -0.4, 0.4
    #     zminz, zmaxz = -0.17, 0.17
    #     Lx = xmaxz - xminz
    #     Lz = zmaxz - zminz
    #     xlims!(ax1, -0.4, 0.4)
    #     ylims!(ax1, -0.17, 0.17)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     Mak.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:AnisotropyFactor
    #     ax1 = Axis(f[1, 1], title = L"$δ_\textrm{ani}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, δani, colormap = (:bilbao, α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     # xminz, xmaxz = -0.4, 0.4
    #     # zminz, zmaxz = -0.17, 0.17
    #     # Lx = xmaxz - xminz
    #     # Lz = zmaxz - zminz
    #     # xlims!(ax1, -0.4, 0.4)
    #     # ylims!(ax1, -0.17, 0.17)
    #     Mak.Colorbar(f[1, 2], hm, label = L"$δ_\textrm{ani}$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:MeltFraction
    #     ax1 = Axis(f[1, 1], title = L"ϕ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, ϕ, colormap = (:bilbao, α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     Mak.Colorbar(f[1, 2], hm, label = L"$ϕ$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Cohesion
    #     ax1 = Axis(f[1, 1], title = L"$C$ [MPa] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, C./1e6, colormap = (Reverse(:bilbao), α_heatmap))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     Mak.Colorbar(f[1, 2], hm, label = L"$C$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Temperature
    #     ax1 = Axis(f[1, 1], title = L"$T$ [C] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    #     hm = heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = Reverse(:bilbao))
    #     AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, Fab, height, Lc, cm_y, group_phases, Δ)                
    #     Mak.Colorbar(f[1, 2], hm, label = L"$T$ [C]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
    #     Mak.colgap!(f.layout, 20)
    #     colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    #     xlims!(ax1, window.xmin, window.xmax)
    #     ylims!(ax1, window.zmin, window.zmax)
    #     if printfig Print2Disk( f, path, string(field), istep) end
    # end

    # if field==:Topography
    #     height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
    #     Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
    #     Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
    #     Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
    #     Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
    #     x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
    #     z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));

    #     ax1 = Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$h$ [km]")
    #     @show mean(height)
    #     @show mean(z_mark)
    #     lines!(ax1, xv./Lc, height./Lc)
    #     scatter!(ax1, x_mark./Lc, z_mark./Lc)
    #     xlims!(ax1, window.xmin, window.xmax)

    #     ax2 = Axis(f[2, 1], xlabel = L"$x$ [km]", ylabel = L"$Vx$ [km]")
    #     lines!(ax2, xv./Lc, Vx_grid./Vc)
    #     scatter!(ax2, x_mark./Lc, Vx_mark./Vc)
    #     xlims!(ax1, window.xmin, window.xmax)

    #     ax3 = Axis(f[3, 1], xlabel = L"$x$ [km]", ylabel = L"$Vz$ [km]")
    #     lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
    #     scatter!(ax3, x_mark/Lc, Vz_mark./Vc)
    #     xlims!(ax1, window.xmin, window.xmax)
    # end

    # if field==:TimeSeries
    #     ax1 = Axis(f[1, 1], title = L"$τ_{xz}$", xlabel = L"$t$", ylabel = L"$\tau_{xz}$")
    #     τxz_t[1] = 0.
    #     @show .-τxz_t
    #     lines!(ax1, t_t, .-τxz_t./(.-P_t.+τzz_t))
    # end

    # if field==:EffectiveFrictionTime
    #     ϕ_eff = -mean(τxzc[:,end] ./ (τzz[:,end] .- P[:,end]) )
    #     if istep==0 ϕ_eff = 0. end
    #     push!(probe.ϕeff, ϕ_eff) 
    #     push!(probe.t, t) 
    #     if istep==file_end
    #         ax1 = Axis(f[1, 1], title = L"$ϕ_\mathrm{eff}$", xlabel = L"$t$", ylabel = L"$ϕ_\mathrm{eff}$")
    #         lines!(ax1, Float64.(probe.t), Float64.(probe.ϕeff))
    #     end
    # end

    # if field==:ChristmasTree
    #     ax1 = Axis(f[1, 1], title = L"Stress profile at $t$ = %$(tMy) Ma", xlabel = L"$τII$ [MPa]", ylabel = L"$z$ [km]")
    #     lines!(ax1, mean(τII, dims=1)[:]/τc, coords.c.z./Lc/1e3, )
    #     lines!(ax1, mean(T, dims=1)[:], coords.c.z./Lc/1e3, )
    #     ax2 = Axis(f[1, 2], title = L"Temperature profile at $t$ = %$(tMy) Ma", xlabel = L"$T$ [C]", ylabel = L"$h$ [km]")
    #     lines!(ax2, mean(T, dims=1)[:], coords.c.z./Lc/1e3, )
    # end

    # if field!=:EffectiveFrictionTime || istep!=file_end
    #     DataInspector(f)
        display(f)
    #     sleep(nap)
    # end

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
     save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()