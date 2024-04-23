import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, CairoMakie, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

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

function main()

    # Set the path to your files
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB//"
    # path ="/Users/tduretz/Downloads/TEST_ShearBandsHomo_SRC/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/DoubleSubduction_OMP16/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/NR00/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/qcoe_ref/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/qcoe_LR/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/qcoe_simp2/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/qcoe_simp_tau1e10/"
    # path = "/home/whalter1/github/MDOODZ7.0/cmake-exec/ShearTemplate/"
    # path = "/home/whalter1/github/MDOODZ7.0/cmake-exec/WilliamB6/"
    # path = "/home/whalter1/github/MDOODZ7.0/cmake-exec/William_Dabrowski_1incl/"
    path = "/home/whalter1/github/MDOODZ7.0/cmake-exec/ShearTemplate/"

    # File numbers
    file_start = 0
    file_step  = 10
    file_end   = 130

    # Select field to visualise
    # field = :Phases
    # field = :Cohesion
    # field = :Density
    # field = :Viscosity 
    # field = :PlasticStrainrate
    # field = :Stress
    field = :StrainRate
    # field = :Strain
    # field = :Strain_xx # normal strain
    # field = :Strain_zz # normal strain
    # field = :Strain_xz # shear strain
    # field = :FabricStrain_xx # normal strain along foliation (i.e. εxxp)
    # field = :FabricStrain_zz # normal strain along director (i.e. εzzp)
    # field = :FabricStrain_xz # shear strain along and across foliation (i.e. εxzp)
    # field = :Fxx
    # field = :Fzz
    # field = :Fxz
    # field = :Fzx
    # field = :Fxxp
    # field = :Fzzp
    # field = :Fxzp
    # field = :Fzxp
    # field = :Pressure
    # field = :Temperature
    # field = :Velocity_x
    # field = :Velocity_z
    # field = :Velocity
    # field = :GrainSize
    # field = :Topography
    # field = :TimeSeries
    # field = :AnisotropyFactor

    # Switches
    printfig    = true  # print figures to disk
    printvid    = true
    framerate   = 3
    ph_contours = false  # add phase contours
    T_contours  = false  # add temperature contours
    fabric      = true  # add fabric quiver (normal to director)
    topo        = false
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 
    σ1_axis     = false
    nap         = 0.3    # pause for animation 
    resol       = 1000
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    Lc = 1.0
    tc = My
    Vc = 1.0

    # Time loop
    for istep=file_start:file_step:file_end
    
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model  = ExtractData( filename, "/Model/Params")
        xc     = ExtractData( filename, "/Model/xc_coord")
        zc     = ExtractData( filename, "/Model/zc_coord")
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")
        xv_hr  = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr  = ExtractData( filename, "/VizGrid/zviz_hr")
        # τxz_t  = ExtractData( filename, "TimeSeries/sxz_mean_time")
        # t_t    = ExtractData( filename, "TimeSeries/Time_time")

        xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

        @show minimum(xc)
        @show maximum(xc)

        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zmax-zmin)/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        Lx>1e3 ? length_unit="km" :  length_unit="m"
        @info "Model info"
        @show "Model apect ratio" Lx/Lz
        @show "Model time" t/My
        @show length_unit
        centroids, vertices = (ncx, ncz), (nvx, nvz)

        ph           = ExtractField(filename,  "/VizGrid/compo", centroids, false, 0)
        mask_air     = ph .== -1.00 
        ph_hr        = ExtractField(filename,  "/VizGrid/compo_hr", (ncx_hr, ncz_hr), false, 0)
        ph_dual_hr   = ExtractField(filename,  "/VizGrid/compo_dual_hr", (ncx_hr, ncz_hr), false, 0)
        group_phases = copy(ph_hr); ph_hr[ph_hr.==-1.00] .= NaN
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
        εII   = Float64.(reshape(ExtractData( filename, "/Centers/strain"), ncx, ncz));          εII[mask_air] .= NaN
        εxx   = Float64.(reshape(ExtractData( filename, "/Centers/strain_xx"), ncx, ncz));       εxx[mask_air] .= NaN
        εzz   = Float64.(reshape(ExtractData( filename, "/Centers/strain_zz"), ncx, ncz));       εzz[mask_air] .= NaN
        εxz   = Float64.(reshape(ExtractData( filename, "/Centers/strain_xz"), ncx, ncz));       εxz[mask_air] .= NaN
        εxxp  = Float64.(reshape(ExtractData( filename, "/Centers/strain_xxp"), ncx, ncz));      εxxp[mask_air].= NaN
        εzzp  = Float64.(reshape(ExtractData( filename, "/Centers/strain_zzp"), ncx, ncz));      εzzp[mask_air].= NaN
        εxzp  = Float64.(reshape(ExtractData( filename, "/Centers/strain_xzp"), ncx, ncz));      εxzp[mask_air].= NaN
        Fxx   = Float64.(reshape(ExtractData( filename, "/Centers/Fxx"), ncx, ncz));             Fxx[mask_air] .= NaN
        Fzz   = Float64.(reshape(ExtractData( filename, "/Centers/Fzz"), ncx, ncz));             Fzz[mask_air] .= NaN
        Fxz   = Float64.(reshape(ExtractData( filename, "/Centers/Fxz"), ncx, ncz));             Fxz[mask_air] .= NaN
        Fzx   = Float64.(reshape(ExtractData( filename, "/Centers/Fzx"), ncx, ncz));             Fzx[mask_air] .= NaN
        Fxxp  = Float64.(reshape(ExtractData( filename, "/Centers/Fxxp"), ncx, ncz));            Fxxp[mask_air].= NaN
        Fzzp  = Float64.(reshape(ExtractData( filename, "/Centers/Fzzp"), ncx, ncz));            Fzzp[mask_air].= NaN
        Fxzp  = Float64.(reshape(ExtractData( filename, "/Centers/Fxzp"), ncx, ncz));            Fxzp[mask_air].= NaN
        Fzxp  = Float64.(reshape(ExtractData( filename, "/Centers/Fzxp"), ncx, ncz));            Fzxp[mask_air].= NaN
        τII   = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
        C     = Float64.(reshape(ExtractData( filename, "/Centers/cohesion"), ncx, ncz))
        
        if fabric
            δani  = ExtractField(filename, "/Centers/ani_fac", centroids, false, 0)
            Nx    = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
            Nz    = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
            Fab_x = -Nz./Nx
            Fab_z = ones(size(Nz))
            nrm   = sqrt.(Fab_x.^2 .+ Fab_z.^2)
            Fab_x ./= nrm
            Fab_z ./= nrm
        end
        if topo
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
        if σ1_axis 
            σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        end
        
        #####################################

        # Color palette for phase map
        cmap    = zeros(RGB{Float64}, 7)
        cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
        cmap[2] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
        cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
        cmap[4] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cmap[5] = RGBA{Float64}(217/255, 099/255, 097/255, 1.) 
        cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
        cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)

        # # Group phases for contouring
        # group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        # group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        # group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################

        f = Figure(size = (Lx/Lz*resol, resol), fontsize=25)

        if field==:Phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            # hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = :turbo)
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_dual_hr, colormap = :turbo)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 3, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 6, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end 
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Density
            ax1 = Axis(f[1, 1], title = L"$ρ$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end 
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$ρ$ [kg.m$^3$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII./1e6, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end   
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [MPa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Pressure
            ax1 = Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$P$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, εII, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            @show minimum(εII), maximum(εII)
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{II}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain_xx
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xxd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxx), maximum(εxx)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εxx, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{xxd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain_zz
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{zzd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εzz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εzz), maximum(εzz)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εzz, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{zzd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain_xz
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xzd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxz), maximum(εxz)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εxz, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{xzd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:FabricStrain_xx
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xxd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxxp), maximum(εxxp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εxxp, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{xxd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:FabricStrain_zz
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{zzd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εzzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εzzp), maximum(εzzp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εzzp, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{zzd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:FabricStrain_xz
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xzd}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxzp), maximum(εxzp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, εxzp, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{xzd}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fxx
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{xx}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxx), maximum(Fxx)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fxx, colormap = (:turbo, α_heatmap), colorrange=(0.5,1.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{xx}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fzz
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{zz}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzz), maximum(Fzz)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fzz, colormap = (:turbo, α_heatmap), colorrange=(0.5,1.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{zz}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fxz
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{xz}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxz), maximum(Fxz)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fxz, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{xz}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fzx
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{zx}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzx), maximum(Fzx)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fzx, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{zx}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fxxp
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{xxp}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxxp), maximum(Fxxp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fxxp, colormap = (:turbo, α_heatmap), colorrange=(0.5,1.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{xxp}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fzzp
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{zzp}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzzp), maximum(Fzzp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fzzp, colormap = (:turbo, α_heatmap), colorrange=(0.5,1.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{zzp}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fxzp
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{xzp}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxzp), maximum(Fxzp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fxzp, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{xzp}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Fzxp
            ax1 = Axis(f[1, 1], title = L"$\textrm{F_{zxp}}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzxp), maximum(Fzxp)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Fzxp, colormap = (:turbo, α_heatmap), colorrange=(-0.5,0.5))
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label =  L"$\textrm{F_{zxp}}$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:PlasticStrainrate
            ε̇pl[ε̇pl.==0.0] .= 1e-30
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity
            ax1 = Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            Vx_BG = 0*xc .- 2*zc'  
            V     = sqrt.( (Vxc .- 0.0*Vx_BG).^2 + (Vzc).^2)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, V, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            # xlims!(ax1, 0., 3.e-3)
            # ylims!(ax1, 0., 3.e-3)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$V$ [m.s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:GrainSize
            ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:AnisotropyFactor
            ax1 = Axis(f[1, 1], title = L"$δ_\textrm{ani}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, δani, colormap = (:bilbao, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            # xminz, xmaxz = -0.4, 0.4
            # zminz, zmaxz = -0.17, 0.17
            # Lx = xmaxz - xminz
            # Lz = zmaxz - zminz
            # xlims!(ax1, -0.4, 0.4)
            # ylims!(ax1, -0.17, 0.17)
            # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$δ_\textrm{ani}$", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Cohesion
            ax1 = Axis(f[1, 1], title = L"$C$ [MPa] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, C./1e6, colormap = (:bilbao, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$C$ [MPa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_x
            ax1 = Axis(f[1, 1], title = L"$V_{x}$ [cm/y] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, Vxc.*cm_y, colormap = (:bilbao, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$V_{x}$ [cm/y]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Temperature
            ax1 = Axis(f[1, 1], title = L"$T$ [C] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (:bilbao, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Fab_x, Fab_z, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            CairoMakie.Colorbar(f[1, 2], hm, label = L"$T$ [C]", width = 20, labelsize = 25, ticklabelsize = 14 )
            CairoMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Topography
            ax1 = Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$h$ [km]")
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)

            ax2 = Axis(f[2, 1], xlabel = L"$x$ [km]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)

            ax3 = Axis(f[3, 1], xlabel = L"$x$ [km]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)
        end

        if field==:TimeSeries
            ax1 = Axis(f[1, 1], title = L"$τ_{xz}$", xlabel = L"$t$", ylabel = L"$\tau_{xz}$")
            lines!(ax1, t_t, τxz_t)
        end

        DataInspector(f)
        display(f)
        sleep(nap)

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