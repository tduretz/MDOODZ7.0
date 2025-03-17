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

ξcr  = [0.28 47.58 69.00 100.45 105.12 110.91 105.66921f0 98.70723f0 91.94]

δ = [1 1.5 2 3 4 5 6 7 8 10]
off7 = [0.1 41.33 53.03 65.56 74.74 83.92 103.95 110.63 63.08 NaN]

function Moho_LAB_offset(all_phases, T, coords, height)

    # @show v,indx = findmin(height)
    # @show x_moho = coords.v.x[indx] / 1e3
    # 𝐗_moho = [x_moho, height[indx]]

    crust_thick = sum((all_phases.ph_hr .==0), dims=2)*(coords.c_hr.z[2]-coords.c_hr.z[1])

    @show size(crust_thick)
@show size(coords.c_hr.x)

    # @show crust_thick/1e3
    @show v,indx = findmin(crust_thick)
    @show x_moho = coords.c_hr.x[indx] / 1e3
    𝐗_moho = [x_moho, x_moho]

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
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1/"

    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_MR/"

    # File numbers
    file_start = 900
    file_step  = 100
    file_end   = 900
    
    # file_start = 1000
    # file_step  = 100
    # file_end   = 1000

    offset = []

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
    printfig    = true  # print figures to disk
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
        εII_pl= Float64.(reshape(ExtractData( filename, "/Centers/strain_pl"), ncx, ncz))
        εII_pwl= Float64.(reshape(ExtractData( filename, "/Centers/strain_pwl"), ncx, ncz))
        εII   = εII_pl + εII_pwl
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
        cmap[1] = RGBA{Float64}(221/255, 205/255, 176/255, 1.)  
        cmap[2] = RGBA{Float64}(1/255, 1/255, 1/255, 1.)  
        cmap[3] = RGBA{Float64}(190/255, 216/255, 172/255, 1.) 
        cmap[4] = RGBA{Float64}(157/255, 199/255, 189/255, 1.) 
        cmap[5] = RGBA{Float64}(190/255, 216/255, 172/255, 1.) 
        cmap[6] = RGBA{Float64}(255/255, 255/255, 255/255, 1.) 
        cmap[7] = RGBA{Float64}(157/255, 199/255, 189/255, 1.) 
        phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)

        # # Group phases for contouring
        # group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        # group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        # group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################
        empty!(f)
        ftsz =  30*resol/500
        f = Figure(size = (0.75*Lx/Lz*resol*1.2, resol), fontsize=ftsz)

        if field==:Phases

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

            scales = (
                Lc = 1e3,
                tc = My,
                Vc = 1.0,
                τc = 1e6,
            )

            model = ReadFile(path, istep, scales, options, PlotOnTop)
            @unpack tMy, length_unit, Lx, Lz, xc, zc, ε̇II, τII, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, all_phases, Δ = model
            
            𝐗_moho, 𝐗_lab = Moho_LAB_offset(all_phases, T, coords, height)
            off = (abs(𝐗_moho[1] - 𝐗_lab[1]))

            tMy_string = @sprintf("%1.2lf", tMy)
            off_string = @sprintf("%1.2lf", off)

            # ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy_string) Ma - \delta = 4 - offset = %$(off_string) km", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            ax1 = Axis(f[1, 1], title = L"Isotropic model (\delta = 1)  at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")

            # ax1 = Axis(f[1, 1], title = L"Anisotropic model (\delta = 6)  at $t$ = %$(tMy_string) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc_hr./scales.Lc, zc_hr./scales.Lc, ph_hr, colormap = phase_colors)
            
            group_phases1 = copy(ph_hr)
            group_phases1[isnan.(group_phases1)] .= -1 # add air
            group_phases1[group_phases1.==1] .= 2   # remove inclusion

            Mak.contour!(ax1, coords.c_hr.x./scales.Lc, coords.c_hr.z./scales.Lc, group_phases1, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  

            # (Reverse(:bilbao), 0.6)
            cmap  = cgrad(Reverse(:bilbao), length(εII), categorical=true)
            colors = [(cmap[i],1-i/length(εII)) for i in 1:length(εII)]

            hm = heatmap!(ax1, xc./scales.Lc, zc./scales.Lc, (εII), colormap = (colors, 0.6)) #, colorrange=(0,10)


            lines!(𝐗_moho[1].*ones(size(zc)), zc./scales.Lc, color=:white, linestyle=:dash, linewidth=3)
            lines!( 𝐗_lab[1].*ones(size(zc)), zc./scales.Lc, color=:white, linestyle=:dash, linewidth=3)
            
            display(f)

            push!(offset, abs(𝐗_moho[1] - 𝐗_lab[1]))

            # ax1 = Axis(f[1, 1], title = L"Thickness ratio at $t$ = %$(tMy_string) Ma - \delta = 8", xlabel = L"$x$ [km]", ylabel = L"$H^\textrm{left} / H^\textrm{right}$ [-]")

            # ph_lit = ph_dual_hr.==4 .|| ph_dual_hr.==0 .|| ph_dual_hr.==2
            # lit_thick  = (sum(ph_lit, dims=2) .*Δz)[:]

            # ph_crust = ph_dual_hr.==4 .|| ph_dual_hr.==0 #.|| ph_dual_hr.==2
            # crust_thick  = (sum(ph_crust, dims=2) .*Δz)[:]
            
            # Nx_2 = Int64(size(ph_dual_hr,1)/2)
            # crust_left  = crust_thick[1:Nx_2]
            # crust_right = crust_thick[end:-1:end-Nx_2+1]

            # lit_left  = lit_thick[1:Nx_2]
            # lit_right = lit_thick[end:-1:end-Nx_2+1]

            # error_lit = sum(abs.(lit_left.-lit_right)./abs.(lit_right))
            # error_cr  = sum(abs.(crust_left.-crust_right)./abs.(crust_right))
            # @show error_lit
            # @show error_cr
            # # lines!(ax1, xc_hr[1:Nx_2]./1e3, crust_left ./1e3)  
            # # lines!(ax1, xc_hr[1:Nx_2]./1e3, crust_right./1e3)
            # # lines!(ax1, xc_hr[1:Nx_2]./1e3, lit_left ./1e3)  
            # # lines!(ax1, xc_hr[1:Nx_2]./1e3, lit_right./1e3)

            # lines!(ax1, xc_hr[1:Nx_2]./1e3, crust_left ./crust_right, label=L"$H_\textrm{crust}^\textrm{left} / H_\textrm{crust}^\textrm{right}$")  
            # lines!(ax1, xc_hr[1:Nx_2]./1e3, lit_left   ./lit_right,   label=L"$H_\textrm{lith}^\textrm{left} / H_\textrm{lith}^\textrm{right}$")  
            # axislegend(framevisible=false, position=:lt)

            # δ = [1 1.5 2 3 4 5 6 7 8 10]
            # ξlit = [0.051307674f0 7.64 10.84 16.83 22.91 26.42 29.65 31.86 33.428 35.20] 
            # ξcr  = [0.28 47.58 69.00 100.45 105.12 110.91 105.66921f0 98.70723f0 91.94 81.53]
            # ax2 = Axis(f[1, 2], title = L"$$Effect of $\delta$ on symmetry", ylabel=L"Deviation from symmetry $(%)$", xlabel=L"$\delta$ [-]" )
            # lines!(ax2, δ[:], ξlit[:], label="Lithosphere")
            # lines!(ax2, δ[:], ξcr[:], label="Crust")
            # ylims!(ax2, 0, 120)
            # axislegend(framevisible=false, position=:lt)

            # save("/Users/tduretz/PowerFolders/_manuscripts/PureSimpleShearAnisotropy/Figures/Systematics.png", f, px_per_unit = 4)   
            
            # # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_crust, colormap = :turbo)
            # AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)                
            # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            # Mak.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            # Mak.colgap!(f.layout, 20)
            # xlims!(ax1, window.xmin, window.xmax)
            # ylims!(ax1, window.zmin, window.zmax)
            if printfig Print2Disk( f, path, string(field), istep) end
        end
    end
    @show offset
    @show extrema(offset[offset.<200])

end

main()