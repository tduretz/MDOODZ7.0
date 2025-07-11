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

ξcr  = [0.28 47.58 69.00 100.45 105.12 110.91 105.66921f0 98.70723f0 91.94]

function Thicknesses(d1, Δz)
    d1_ph_lit = d1.all_phases.group_phases.==4 .|| d1.all_phases.group_phases.==0 .|| d1.all_phases.group_phases.==2
    d1_lit_thick  = (sum(d1_ph_lit, dims=2) .*Δz)[:]
    d1_ph_crust = d1.all_phases.group_phases.==4 .|| d1.all_phases.group_phases.==0 
    d1_crust_thick  = (sum(d1_ph_crust, dims=2) .*Δz)[:]
    return (crust=d1_crust_thick, lith=d1_lit_thick)
end
 
@views function main()

    # Set the path to your files
    path = ["/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1_5/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d2/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d3/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d5/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d6/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d7/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d8/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d10/",
            "/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_MR/"]

    step = [1100, 1100, 1100, 1100, 1140, 1100, 1100, 1100, 1240, 1200]    

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
    options = (
        printfig    = true,  # print figures to disk
        printvid    = false,
        framerate   = 6,
        α_heatmap   = 0.75,   # transparency of heatmap 
        vel_arrow   = 5,
        vel_scale   = 10000,
        vel_step    = 20,
        nap         = 0.1,    # pause for animation 
        resol       = 500,
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
    Lc = 1,
    tc = My,
    Vc = 1.0,
    τc = 1,
    )

    probe = (ϕeff = Float64.([]), t  = Float64.([]))
    cm_yr = 100.0*3600.0*24.0*365.25

    
        # name     = @sprintf("Output%05d.gzip.h5", istep)
        # @info "Reading $(name)"
        # filename = string(path, name)
        # model    = ExtractData( filename, "/Model/Params")
        # xc       = ExtractData( filename, "/Model/xc_coord")
        # zc       = ExtractData( filename, "/Model/zc_coord")
        # xv       = ExtractData( filename, "/Model/xg_coord")
        # zv       = ExtractData( filename, "/Model/zg_coord")
        # xvz      = ExtractData( filename, "/Model/xvz_coord")
        # zvx      = ExtractData( filename, "/Model/zvx_coord")
        # xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr")
        # zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr")
        # τzz_t    = ExtractData( filename, "TimeSeries/szzd_mean_time")
        # P_t      = ExtractData( filename, "TimeSeries/P_mean_time")
        # τxz_t    = ExtractData( filename, "TimeSeries/sxz_mean_time")
        # t_t      = ExtractData( filename, "TimeSeries/Time_time")

        # xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        # zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        # ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)
        # coords = ( c=(x=xc, z=zc), c_hr=(x=xc_hr, z=zc_hr), v=(x=xv, z=zv))

        # t      = model[1]
        # tMy    = round(t/tc, digits=6)
        # nvx    = Int(model[4])
        # nvz    = Int(model[5])
        # ncx, ncz = nvx-1, nvz-1
        # xmin, xmax = xv[1], xv[end]
        # zmin, zmax = zv[1], zv[end]
        # Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        # Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        # Lx>1e3 ? length_unit="km" :  length_unit="m"

        # if @isdefined zoom
        #     Lx = (zoom.xmax - zoom.xmin)./Lc 
        #     Lz = (zoom.zmax - zoom.zmin)./Lc
        #     window = zoom
        # else
        #     window = ( 
        #         xmin = minimum(xv)/Lc, 
        #         xmax = maximum(xv)/Lc,
        #         zmin = minimum(zv)/Lc,
        #         zmax = maximum(zv)/Lc,
        #     )
        # end

        # @info "Model info"
        # @show "Model apect ratio" Lx/Lz
        # @show "Model time" t/My
        # @show length_unit
        # centroids, vertices = (ncx, ncz), (nvx, nvz)

        # ph           = ExtractField(filename,  "/VizGrid/compo", centroids, false, 0)
        # mask_air     = ph .== -1.00 
        # ph_hr        = ExtractField(filename,  "/VizGrid/compo_hr", (ncx_hr, ncz_hr), false, 0)
        # ph_dual_hr   = ExtractField(filename,  "/VizGrid/compo_dual_hr", (ncx_hr, ncz_hr), false, 0)
        # all_phases = copy(ph_hr); ph_hr[ph_hr.==-1.00] .=   NaN; ph_dual_hr[ph_dual_hr.==-1.00] .=   NaN;
        # ηc           = ExtractField(filename,  "/Centers/eta_n", centroids, true, mask_air)
        # ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        # P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        # T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        # d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        # ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
        # Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        # Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
        # τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        # τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz))
        # τyy   = -(τzz .+ τxx)
        # σzz   = -P + τzz
        # σxx   = -P + τxx
        # τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        # ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        # ε̇zz   = Float64.(reshape(ExtractData( filename, "/Centers/ezzd"), ncx, ncz))
        # ε̇yy   = -(ε̇xx .+ ε̇zz)
        # ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        # τII   = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        # ε̇II   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
            
        # τxzc  = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end,1:end-1] .+ τxz[1:end-1,2:end] .+ τxz[2:end,2:end]) 
        # C     = Float64.(reshape(ExtractData( filename, "/Centers/cohesion"), ncx, ncz))
        # ϕ     = ExtractField(filename, "/Centers/phi", centroids, false, 0)
        # divu  = ExtractField(filename, "/Centers/divu", centroids, false, 0)
       
        # T_hr  = zeros(size(ph_hr)); T_hr[1:2:end-1,1:2:end-1] .= T; T_hr[2:2:end-0,2:2:end-0] .= T
        # if LAB_color
        #     ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.<LAB_T] .= 2
        #     ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.>LAB_T] .= 3
        #     ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.<LAB_T] .= 2
        #     ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.>LAB_T] .= 3
        #     ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.<LAB_T] .= 2
        #     ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.>LAB_T] .= 3
        # end

        # Fab = 0.
        # if PlotOnTop.fabric
        #     δani    = ExtractField(filename, "/Centers/ani_fac", centroids, false, 0)
        #     Nx      = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
        #     Nz      = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
        #     Fab     = (x=-Nz./Nx, z=ones(size(Nz)))
        #     nrm     = sqrt.(Fab.x.^2 .+ Fab.z.^2)
        #     Fab.x ./= nrm
        #     Fab.z ./= nrm
        # end
        # height = 0.
        # if PlotOnTop.topo
        #     height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
        #     Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
        #     Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
        #     Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
        #     Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
        #     x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
        #     z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
        # end
        # Vxc   = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        # Vzc   = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])
        # V = (x=Vxc, z=Vzc, arrow=vel_arrow, step=vel_step, scale=vel_scale)
        # σ1, ε̇1 = 0., 0.
        # if PlotOnTop.σ1_axis 
        #     σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        # end
        # if PlotOnTop.ε̇1_axis 
        #     ε̇1 = PrincipalStress(ε̇xx, ε̇zz, ε̇xz, zeros(size(ε̇xx))) 
        # end
        # PT = (P.>2.2e9 .&& P.<3.0e9 .&& T.>430 .&& T.<530).*ones(size(T))
     
        d1 = ReadFile(path[1], step[1], scales, options, PlotOnTop)
        @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, Δ = d1
        d1_5 = ReadFile(path[2], step[2], scales, options, PlotOnTop)
        d2 = ReadFile(path[3], step[3], scales, options, PlotOnTop)
        d3 = ReadFile(path[4], step[4], scales, options, PlotOnTop)
        d5 = ReadFile(path[5], step[5], scales, options, PlotOnTop)
        d6 = ReadFile(path[6], step[6], scales, options, PlotOnTop)
        d7 = ReadFile(path[7], step[7], scales, options, PlotOnTop)
        d8 = ReadFile(path[8], step[8], scales, options, PlotOnTop)
        d10 = ReadFile(path[9], step[9], scales, options, PlotOnTop)
        d4_ref = ReadFile(path[10], step[10], scales, options, PlotOnTop)

        
        @show d1.tMy
        @show d5.tMy
        @show d10.tMy

        Δz = (zc[2] - zc[1])/2

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
        # all_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        # all_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        # all_phases[ ph_hr.==3 ]                                             .= 3

        #####################################
        ftsz =  20*options.resol/500
        f = Figure(size = (1.2*options.Lx/options.Lz*options.resol*1.2, options.resol), fontsize=ftsz)

        if field==:Phases
            tMy_string = @sprintf("%1.2lf", tMy)
            # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)

            d1_thick = Thicknesses(d1, Δz)
            d1_5_thick = Thicknesses(d1_5, Δz)
            d2_thick = Thicknesses(d2, Δz)
            d3_thick = Thicknesses(d3, Δz)
            d5_thick = Thicknesses(d5, Δz) 
            d6_thick = Thicknesses(d6, Δz)
            d7_thick = Thicknesses(d7, Δz)
            d8_thick = Thicknesses(d8, Δz)
            d10_thick = Thicknesses(d10, Δz)
            d4_thick = Thicknesses(d4_ref, Δz)
           
            ##################################
            ax1 = Axis(f[1, 1:2], title = L"A) Crustal thickness at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$H^\textrm{crust}$ [km]")
            # hidexdecorations!(ax1)
            lines!(ax1, coords.c_hr.x./1e3, d1_thick.crust./1e3, label=L"$\delta = 1$")  
            lines!(ax1, coords.c_hr.x./1e3, d4_thick.crust./1e3, label=L"$\delta = 4$ (ref)", color=:black, linestyle=:dash)  
            lines!(ax1, coords.c_hr.x./1e3, d5_thick.crust./1e3, label=L"$\delta = 5$")  
            lines!(ax1, coords.c_hr.x./1e3, d10_thick.crust./1e3, label=L"$\delta = 10$")  

            axislegend(framevisible=false, position=:lb)


            ax2 = Axis(f[2, 1:2], title = L"B) Lithosphere thickness at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$H^\textrm{lith}$ [km]")
            lines!(ax2, coords.c_hr.x./1e3, d1_thick.lith./1e3  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  
            lines!(ax2, coords.c_hr.x./1e3, d4_thick.lith./1e3, label=L"$\delta = 4$ (ref)", color=:black, linestyle=:dash)  
            lines!(ax2, coords.c_hr.x./1e3, d5_thick.lith./1e3  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  
            lines!(ax2, coords.c_hr.x./1e3, d10_thick.lith./1e3  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  

            axislegend(framevisible=false, position=:lb)


            # ax1 = Axis(f[1, 1:2], title = L"Crustal thickness at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$\Delta H^\textrm{crust}$ (%)")
            # hidexdecorations!(ax1)

            # lines!(ax1, coords.c_hr.x./1e3, (d1_thick.crust.-d1_thick.crust)./d1_thick.crust.*100, label=L"$\delta = 1$")  
            # lines!(ax1, coords.c_hr.x./1e3, (d5_thick.crust.-d1_thick.crust)./d1_thick.crust.*100, label=L"$\delta = 5$")  
            # lines!(ax1, coords.c_hr.x./1e3, (d10_thick.crust.-d1_thick.crust)./d1_thick.crust.*100, label=L"$\delta = 10$")  

            # axislegend(framevisible=false, position=:lt)

            # ax2 = Axis(f[2, 1:2], title = L"Lithosphere thickness at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$\Delta H^\textrm{lit}$ (%)")

            # lines!(ax2, coords.c_hr.x./1e3, (d1_thick.lith.-d1_thick.lith)./d1_thick.lith.*100  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  
            # lines!(ax2, coords.c_hr.x./1e3, (d5_thick.lith.-d1_thick.lith)./d1_thick.lith.*100  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  
            # lines!(ax2, coords.c_hr.x./1e3, (d10_thick.lith.-d1_thick.lith)./d1_thick.lith.*100  ,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  


            @info "Crust anisotropy"
            @show mean((d1_5_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d2_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d3_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d5_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d6_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d7_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d8_thick.crust.-d1_thick.crust)./d1_thick.crust)
            @show mean((d10_thick.crust.-d1_thick.crust)./d1_thick.crust)

            @info "Lithosphere anisotropy"
            @show mean((d1_5_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d2_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d3_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d5_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d6_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d7_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d8_thick.lith.-d1_thick.lith)./d1_thick.lith)
            @show mean((d10_thick.lith.-d1_thick.lith)./d1_thick.lith)



            δ = [1 1.5 2 3 5 6 7 8 10]
            # ξlit = [0.051307674f0 7.64 10.84 16.83 22.91 26.42 29.65 31.86 33.428 35.20] 
            # ξcr  = [0.28 47.58 69.00 100.45 105.12 110.91 105.66921f0 98.70723f0 91.94 81.53]
            
            ξcr  = [
                0
                0.023999203f0
                0.048868004f0
                0.07952324f0
                0.10476836f0
                0.12021628f0
                0.13891824f0
                0.14822909f0
                0.17336552f0
            ]

            ξlit = [
                0
                0.052580155f0
                0.15157382f0
                0.22192997f0
                0.3184466f0
                0.37015787f0
                0.42853257f0
                0.47614223f0
                0.5619342f0
            ]
            # ax3 = Axis(f[3, 1:2], title = L"$$C) Deviation caused by anisotropy", ylabel=L"Difference $(%)$", xlabel=L"$\delta$ [-]" )
            # lines!(ax3, δ[:], ξlit[:].*100, label="Lithosphere", color=:black)
            # lines!(ax3, δ[:], ξcr[:].*100, label="Crust", color=:gray)
            # # ylims!(ax3, 0, 120)
            # axislegend(framevisible=false, position=:lt)

            save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/ThicknessVariations_step1.png", f, px_per_unit = 4)   
            
            # # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_crust, colormap = :turbo)
            # AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            # Mak.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            # Mak.colgap!(f.layout, 20)
            # xlims!(ax1, window.xmin, window.xmax)
            # ylims!(ax1, window.zmin, window.zmax)
            # if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Density
            ax1 = Axis(f[1, 1], title = L"$\rho$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap), colorrange=(2600, 3000))  
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\rho$ [kg.m$^{-3}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end

            @show size(τII)
            @show mean((all_phases[1:599,:] .- all_phases[end:-1:600,:]))*100

            error    = 0.
            sumtot   = 0
            for i=1:299, j=1:size(τII,2)
                left  = τII[i,j]
                right = τII[end+1-i,j]
                if !isnan(left) && !isnan(right) && left>0.
                    sum   += 1
                    error += abs(left-right) / abs(left)
                end
            end
            @show error/sumtot*100
        end

        if field==:σxx
            ax1 = Axis(f[1, 1], title = L"$\sigma_\textrm{xx}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, σxx./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\sigma_\textrm{xx}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:σzz
            ax1 = Axis(f[1, 1], title = L"$\sigma_\textrm{zz}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, σzz./τc, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\sigma_\textrm{zz}$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Pressure
            ax1 = Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:turbo, α_heatmap)) #, colorrange=(1,1.2)1e4*365*24*3600
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label =  L"$P$ [GPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Divergence
            ax1 = Axis(f[1, 1], title = L"∇⋅V at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, divu, colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
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
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
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
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            # xlims!(ax1, 0., 3.e-3)
            # ylims!(ax1, 0., 3.e-3)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$V$ [m.s$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_x
            @show Vx
            ax1 = Axis(f[1, 1], title = L"$Vx$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xv, zvx[2:end-1], Vx[2:end-1,:], colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vx$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_z
            ax1 = Axis(f[1, 1], title = L"$Vz$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xvz[2:end-1], zv, Vz[:,2:end-1]*cm_yr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vz$ [cm.yr$^{-1}$]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:GrainSize
            ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:AnisotropyFactor
            ax1 = Axis(f[1, 1], title = L"$δ_\textrm{ani}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, δani, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            # xminz, xmaxz = -0.4, 0.4
            # zminz, zmaxz = -0.17, 0.17
            # Lx = xmaxz - xminz
            # Lz = zmaxz - zminz
            # xlims!(ax1, -0.4, 0.4)
            # ylims!(ax1, -0.17, 0.17)
            Mak.Colorbar(f[1, 2], hm, label = L"$δ_\textrm{ani}$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:MeltFraction
            ax1 = Axis(f[1, 1], title = L"ϕ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ϕ, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$ϕ$", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Cohesion
            ax1 = Axis(f[1, 1], title = L"$C$ [MPa] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, C./1e6, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$C$ [MPa]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            xlims!(ax1, window.xmin, window.xmax)
            ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Temperature
            ax1 = Axis(f[1, 1], title = L"$T$ [C] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (Reverse(:bilbao), α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, all_phases, Δ, Mak)                
            Mak.Colorbar(f[1, 2], hm, label = L"$T$ [C]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
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

        if field!=:EffectiveFrictionTime || istep!=file_end
            DataInspector(f)
            display(f)
        end


    yscale = Lz/Lx

   

end

main()