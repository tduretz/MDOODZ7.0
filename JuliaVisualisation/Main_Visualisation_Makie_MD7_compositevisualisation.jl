import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG

# Select your Makie of Choice: NOTE one has to quit the session when switching from CairoMakie to GLMakie or vice versa.
choice = 2 # 1 => GLMakie; 2 => CairoMakie
MakieOptions = ["import GLMakie as MoC",        # necessary for    interactive 'DataInspector', not supported on any 'headless server'
                "import CairoMakie as MoC"]     # does not support interactive 'DataInspector', runs also on any     'headless server'
eval(Meta.parse(MakieOptions[choice]))
MoC.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
MoC.Makie.inline!(false)

My = 1e6*365*24*3600

function main()

    # Set the path to your files
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/NonLinearPureshearAnisotropic/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/NR00/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/1_NR07/"
    path ="/users/whalter1/MDOODZ_output_folders_ready_for_download/CollisionPolarCartesian_v22_polar/"
    path ="/users/whalter1/work/william/MDOODZ_output_folders_ready_for_download/CollisionPolarCartesian_v27_polar/"
    #path ="/users/whalter1/MDOODZ7.0/cmake-exec/ShearTemplate"
    path = "/Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/"
    # File numbers
    file_start = 0
    file_step  = 10
    file_end   = 200

    # Select one field to visualise
    #field = :Phases
    #field = :Density
    #field = :Viscosity 
    #field = :Stress
    field = :Stress_xx
    #field = :Stress_zz
    #field = :Stress_xz
    #field = :StrainRate
    #field = :PlasticStrainrate
    #field = :Pressure
    #field = :Temperature
    #field = :Velocity
    #field = :Velocity_x
    #field = :Velocity_z
    #field = :Strain # cumulated strain
    #field = :GrainSize # doesn't function properly yet
    #field = :Topography # doesn't function properly yet

    # Switches
    displayfig  = false     # display figures on screen (not possible on any "headless server")
    printfig    = true      # print figures to disk
    printvid    = true      # print a video (printfig must be on)
    ph_contours = false     # add phase contours
    T_contours  = false     # add temperature contours
    fabric      = false     # add fabric quiver (normal to director)
    velocity    = false     # add velocity quiver
    α_heatmap   = 0.85      # transparency of heatmap 
    σ1_axis     = false
    nap         = 0.3       # pause for animation 
    resol       = 800
    inspector   = false     # (not possible on any "headless server")
    framerate   = 2         # Frame rate for the movie
    movieName   = "$(path)/_$(field)/$(field)"  # Name of the movie

    # Scaling
    Lc = 1000.
    tc = My
    Vc = 1e-9

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
        @show "Model apect ratio" Lx/Lz
        @show "Model time" t/My

        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));                  mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr));     
        group_phases = copy( ph_hr); ph_hr[ph_hr.==-1.00] .= NaN    
        ηc    = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz));                  ηc[mask_air]  .= NaN
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));                  ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));                      P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;            T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));                      d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));                 ε̇pl[mask_air] .= NaN
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))      
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))      
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz));                   τxx[mask_air] .= NaN
        τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz));                   τzz[mask_air] .= NaN
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        τxzc  = 0.25*(τxz[1:end-1,1:end-1].+τxz[2:end,1:end-1].+τxz[1:end-1,2:end].+τxz[2:end,2:end]);  τxzc[mask_air] .= NaN
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        εII   = Float64.(reshape(ExtractData( filename, "/Centers/strain"), ncx, ncz));                 εII[mask_air] .= NaN
        τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(2*ε̇xx.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
        if fabric
            Nx    = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz)); 
            Nz    = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz));
        end
        height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
        Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
        Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
        Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
        Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
        x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
        z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
        Vxc     = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1]); Vxc[mask_air] .= NaN
        Vzc     = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0]); Vzc[mask_air] .= NaN
        Vc      = (Vxc.^2 .+ Vzc.^2).^0.5
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
        phase_colors = MoC.cgrad(cmap, length(cmap), categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################

        f = MoC.Figure(resolution = (Lx/Lz*resol, resol), fontsize=25)

        if field==:Phases
            ax1 = MoC.Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Density
            ax1 = MoC.Axis(f[1, 1], title = L"Density $ρ$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end 
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"Density $ρ$ [kg/m^3]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = MoC.Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end 
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end   
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress_xx
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xx}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxx, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end   
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xx}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress_zz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{zz}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τzz, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end   
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{zz}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress_xz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xz}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xv./Lc, zv./Lc, τxz, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end   
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xz}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:PlasticStrainrate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma, sum = %$(round(sum(ε̇pl),sigdigits=4)), max = %$(round(maximum(ε̇pl),sigdigits=4)), min = %$(round(minimum(ε̇pl),sigdigits=4))", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            ε̇pl[mask_air] .= NaN
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Pressure
            ax1 = MoC.Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$P$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Temperature
            ax1 = MoC.Axis(f[1, 1], title = L"$T$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$T$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity
            ax1 = MoC.Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Vc, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_x
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{x}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xv./Lc, zc./Lc, Vx[:,2:end-1], colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{x}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_z
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{z}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zv./Lc, Vz[2:end-1,:], colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{z}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain
            ax1 = MoC.Axis(f[1, 1], title = L"$\varepsilon_{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εII, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
            end
            if σ1_axis
                arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\varepsilon_{II}$ []", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:GrainSize
            ax1 = MoC.Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            if velocity 
                qv,qh = 6,3
                arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
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
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Topography
            ax1 = MoC.Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$h$ [km]")
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)

            ax2 = MoC.Axis(f[2, 1], xlabel = L"$x$ [km]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)

            ax3 = MoC.Axis(f[3, 1], xlabel = L"$x$ [km]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)

            @show minimum(Vx_grid)
            @show minimum(Vx_mark)
            @show maximum(Vx_grid)
            @show maximum(Vx_mark)
        end

        if inspector
            DataInspector(f)
        end
        if displayfig
            display(f)
        end
        sleep(nap)
    end
    if printfig && printvid
            FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1920:1080" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`)
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
    MoC.save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()

#run(`ffmpeg -framerate 2 -i /Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/_Phases/Phases%05d.png -c libx264 -pix_fmt yuv420p -y out_movie.mp4`)

# Working command
#ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4