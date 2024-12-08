import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
# using GLMakie
using CairoMakie
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, FFMPEG
Makie.update_theme!( fonts = ( regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(true)
My = 1e6*365*24*3600

function main()

    # Set the path to your files"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/1_NR01/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/DEBUG/MDLIB/"

    # File numbers
    file_start = 00
    file_step  = 50
    file_end   = 3500

    # Select field to visualise
    # field = :Phases
    # field = :ThinningFactor
    field = :PhasesRheology
    # field = :ChristmasTree
    # field = :Density
    # field = :Viscosity 
    # field = :PlasticStrainrate
    # field = :Stress
    # field = :StrainRate
    # field = :Strain
    # field = :Pressure
    # field = :Temperature
    # field = :Velocity_x
    # field = :Velocity_z
    # field = :GrainSize
    # field = :Topography

    # Switches
    printfig    = true  # print figures to disk
    printvid    = true
    framerate   = 3
    ph_contours = true  # add phase contours
    T_contours  = true  # add temperature contours
    fabric      = false  # add fabric quiver (normal to director)
    α_heatmap   = 0.95   # transparency of heatmap 
    nap         = 0.3    # pause for animation 
    resol       = 1000
    ftsz        = 40
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0

    # Scaling
    Lc = 1000.
    tc = My
    Vc = 1e-9

    Hcrust  = 0
    Hmant   = 0
    Hlit    = 0
    Hcrust0 = 0
    Hmant0  = 0
    Hlit0   = 0 

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

        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));          mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr)); 
        ph_hr2 = copy(ph_hr)
        group_phases = copy( ph_hr); ph_hr[ph_hr.==-1.00] .= NaN
        ηc    = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz));          ηc[mask_air]  .= NaN
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
        ε̇el   = Float64.(reshape(ExtractData( filename, "/Centers/eII_el"), ncx, ncz));         ε̇el[mask_air] .= NaN
        ε̇pwl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pwl"), ncx, ncz));         ε̇pwl[mask_air] .= NaN
        ε̇lin   = Float64.(reshape(ExtractData( filename, "/Centers/eII_lin"), ncx, ncz));         ε̇lin[mask_air] .= NaN
        ε̇exp   = Float64.(reshape(ExtractData( filename, "/Centers/eII_exp"), ncx, ncz));         ε̇exp[mask_air] .= NaN

        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        εII   = Float64.(reshape(ExtractData( filename, "/Centers/strain"), ncx, ncz));          εII[mask_air] .= NaN
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

        Vxc   = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc   = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])

        #####################################

        # Color palette for phase map
        cmap    = zeros(RGB{Float64}, 6)
        cmap[1] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) # 0  
        cmap[2] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) # 1
        cmap[3] = RGBA{Float64}(211/255, 65/255, 59/255, 1.) # 2
        cmap[4] = RGBA{Float64}(112/255, 200/255, 94/255, 1.) # 3
        cmap[5] = RGBA{Float64}(191/255, 145/255, 127/255, 1.) # 4
        cmap[6] = RGBA{Float64}(86/255, 142/255, 121/255, 1.) # 5
        phase_colors = cgrad(cmap, 6, categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        # Transparent turbo
        cbarPal= :turbo
        cmap = get(colorschemes[cbarPal], LinRange(0,1,100))
        turboα = [(cmap[i],i/100) for i in 1:75]

        #####################################

        ### 1. GROUP PHASES
        ph_hr_rheo = copy(ph_hr)

        for i in eachindex(xc_hr), j in eachindex(zc_hr)
            if ph_hr[i,j]==5 
                ph_hr_rheo[i,j] = 0 
            end
            if ph_hr[i,j]==4 
                ph_hr_rheo[i,j] = 0 
            end
            if ph_hr[i,j]==1 
                ph_hr_rheo[i,j] = 0 
            end

            if ph_hr[i,j]==2
                ph_hr_rheo[i,j] = 1
            end
            if ph_hr[i,j]==6
                ph_hr_rheo[i,j] = 4
            end

            if ph_hr[i,j]==3
                ph_hr_rheo[i,j] = 2
            end

            if ph_hr_rheo[i,j]==4
                ph_hr_rheo[i,j] = 1
            end
            if ph_hr_rheo[i,j]==2
                ph_hr_rheo[i,j] = 1
            end

        end
        
        ### 2. COLOR PHASES BASED ON RHEOLOGICAL MODES
        ic = 0
        for i in eachindex(xc_hr)
            if mod(i-1,2)==0 
                ic += 1
            end
            jc = 0
            for j in eachindex(zc_hr)
                if mod(j-1,2)==0 
                    jc += 1
                end

                # Rheological mode
                mode = :vis
                ε̇vis = ε̇exp[ic,jc] + ε̇pwl[ic,jc] + ε̇lin[ic,jc]
                if  ε̇pl[ic,jc] > ε̇vis && ε̇pl[ic,jc] > ε̇el[ic,jc]
                    mode = :pl
                end
                if  ε̇el[ic,jc] > ε̇vis && ε̇el[ic,jc] > ε̇pl[ic,jc]
                    mode = :el
                end

                if isnan(ph_hr[i,j]) == false
                    if ph_hr_rheo[i,j] == 0 && mode==:pl
                        ph_hr_rheo[i,j] = 2
                    end
                    if ph_hr_rheo[i,j] == 0 && mode==:el
                        ph_hr_rheo[i,j] = 4
                    end
                    if ph_hr_rheo[i,j] == 1 && mode==:pl
                        ph_hr_rheo[i,j] = 3
                    end
                    if ph_hr_rheo[i,j] == 1 && mode==:el
                        ph_hr_rheo[i,j] = 5
                    end
                end
            end
        end

        ### 3. COMPUTE THICKNESSES
        Hcrust = sum(group_phases.==0, dims=2)*Δz/2
        Hmant  = sum(group_phases.==1, dims=2)*Δz/2
        Hlit   = sum(group_phases.==0 .|| group_phases.==1, dims=2)*Δz/2
        if istep==0
            Hcrust0 = Hcrust
            Hmant0  = Hmant
            Hlit0   = Hlit
        end

        ##########################################################################

        f = Figure(resolution = (Lx/Lz*resol, resol), fontsize=ftsz)

        ##########################################################################
        if field==:ChristmasTree
            τII[isnan.(τII)] .= 0.
            strength = @sprintf("%2.2e", sum(τII[1,:])*Lz/1e12) 
            ax1 = Axis(f[1, 1], title =  "Integrated strength: "*strength*" TN/m",  xlabel=L"$τ_{II}$ [MPa]", ylabel=L"$z$ [km]")
            lines!(ax1, τII[1,:]./1e6, zc./Lc, linewidth=10)
            
            # @show τII[1,:]
            # @show sum(ones(size(τII))*Δz)
            xlims!(ax1,-5, 700) 
            ylims!(ax1, -160, 1)
            
            T[isnan.(T)] .= 0.
            ρc[isnan.(ρc)] .= 0.
            heat = @sprintf("%2.2e", sum(T[1,:].*ρc[1,:].*1050)*Lz/1e12) 
            ax2 = Axis(f[1, 2], title = "Integrated heat: "*heat* " TJ/m² ",  xlabel=L"$T$ [$^\mathrm{o}$C]")
            hideydecorations!(ax2) 
            lines!(ax2, T[1,:], zc./Lc, linewidth=10)
            lines!(ax2, T[1,:], -35.0*ones(size(zc)), linewidth=5, linestyle=:dash, color=:grey)
            lines!(ax2, 650*ones(size(zc)), zc./Lc, linewidth=5, linestyle=:dash, color=:grey)
            xlims!(ax2,0, 1400)
            ylims!( ax2, -160, 1)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:ThinningFactor
            ax1 = Axis(f[1, 1], title = L"Thinning at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$\beta$ [-]")
            lines!(ax1, xc_hr./Lc, 1.0 .- Hcrust[:]./Hcrust0[:], label=L"$$Crust", linewidth=10)
            lines!(ax1, xc_hr./Lc, 1.0 .- Hmant[:]./Hmant0[:], label=L"$$Mantle lithosphere", linewidth=10)
            lines!(ax1, xc_hr./Lc, 1.0 .- Hlit[:]./Hlit0[:], label=L"$$Lithosphere", linewidth=10)
            xlims!(ax1, xmin/Lc, xmax/Lc)
            ylims!(ax1, -0.05, 1.05)
            axislegend(position = :rt, framevisible = false)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            group_phases[group_phases.==-1] .= NaN
            hm1 = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, colormap = phase_colors, colorrange=(-0.5,5.5))
            # hm2 = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = turboα)
            # if ph_contours 
            #     contour!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr2, levels=-1:1:maximum(ph_hr2), linewidth = 1, color=(:white, 0.25)  )  
            # end
            # if T_contours 
            #     contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            # end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            cbar = Colorbar(f[1, 2], hm1, width = 20, labelsize=ftsz, ticklabelsize = 20 ) #, label = L"$$Phases"
            # Colorbar(f[1, 3], hm2, label = L"$\dot{\varepsilon}^\textrm{p}$", width = 20, labelsize = 25, ticklabelsize = 14 )
            cbar.ticks = ([0, 1, 2, 3, 4, 5], ["Crust viscous", "Mantle viscous", "Crust plastic", "Mantle plastic", "Crust elastic", "Mantle elastic"])
            # colgap!(f.layout, 40)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:PhasesRheology
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm1 = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr_rheo, colormap = phase_colors, colorrange=(-0.5,5.5))
            # hm2 = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = turboα)
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr2, levels=-1:1:maximum(ph_hr2), linewidth = 3, color=(:white, 0.35)  )  
            end
            # if T_contours 
            #     contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            # end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            cbar = Colorbar(f[1, 2], hm1, width = 20, labelsize=ftsz, ticklabelsize = 20 ) #, label = L"$$Phases"
            # Colorbar(f[1, 3], hm2, label = L"$\dot{\varepsilon}^\textrm{p}$", width = 20, labelsize = 25, ticklabelsize = 14 )
            cbar.ticks = ([0, 1, 2, 3, 4, 5], ["Crust viscous", "Mantle viscous", "Crust plastic", "Mantle plastic", "Crust elastic", "Mantle elastic"])
            # colgap!(f.layout, 40)
            xlims!(ax1, -75, 125)
            ylims!(ax1, -75, 2)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end           
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = ftsz, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end          
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Strain
            ax1 = Axis(f[1, 1], title = L"${\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(εII), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
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
            GLMakie.Colorbar(f[1, 2], hm, label =  L"${\varepsilon}_\textrm{II}$ [-]", width = 20, labelsize = ftsz, ticklabelsize = ftsz )
            GLMakie.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:PlasticStrainrate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if fabric 
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
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
                arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
            end
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
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

            @show minimum(Vx_grid)
            @show minimum(Vx_mark)
            @show maximum(Vx_grid)
            @show maximum(Vx_mark)
        end

        if printfig==false
            DataInspector(f)
            display(f) 
            sleep(nap)
        end
            
    end

    yscale = Lz/Lx

    if printfig && printvid
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1080:1080*$(yscale)" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)
    end

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