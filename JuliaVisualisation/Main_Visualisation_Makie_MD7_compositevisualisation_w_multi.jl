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

s  = 1
d  = 24*3600*s
y  = 365*d
My = 1e6*y

function main(path, file_start, file_end,fac, xshift, file_step)

    # File numbers
    #file_start = 0
    #file_step  = 10
    #file_end   = file_start

    # Select one field to visualise
    #field = :Phases
    #field = :Density
    #field = :Viscosity 
    #field = :Stress
    #field = :Stress_xx
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
    field = :AnisoFactor # anisotropy factor
    #field = :GrainSize # doesn't function properly yet
    #field = :Topography # doesn't function properly yet

    # Switches
    displayfig  = false     # display figures on screen (not possible on any "headless server")
    printfig    = true      # print figures to disk
    printvid    = true      # print a video (printfig must be on)
    ph_contours = false     # add phase contours
    T_contours  = false     # add temperature contours
    fabric      = true     # add fabric quiver (normal to director)
    velocity    = false      # add velocity quiver
    deviatoricV = false      # considers deviatoric velocity for plot and quiver (deviatoric velocity = velocity minus far-field velocity)
    α_heatmap   = 0.85      # transparency of heatmap 
    σ1_axis     = false
    polar       = false     # for curved setups using Earth's curvature (transforms tensors from x-z to r-phi coord system)
    nap         = 0.3       # pause for animation 
    resol       = 1600
    inspector   = false     # (not possible on any "headless server")
    framerate   = 30         # Frame rate for the movie
    movieName   = "$(path)/_$(field)/$(field)"  # Name of the movie
    #xshift      = 0.5       # periodically shifts all values to the right, fraction 0.0-1.0

    # Scaling
    Lc = 1.
    tc = s; tMytext= "s"
    Vc = 1e-9

    # Time loop
    if displayfig==1 || printfig==1
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
        @show "Time step" istep
        #@show "Model apect ratio" Lx/Lz
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
        if field==:AnisoFactor
            AniFac= Float64.(reshape(ExtractData( filename, "/Centers/ani_fac"), ncx, ncz))            
            #@show AniFac[50,50]
            @show minimum(AniFac)
            @show maximum(AniFac)
            # weird issue that sometimes aniso factor is negative infinity (e.g. -9e12 or so) where it should actually be at maximum (e.g. 1000, saturated)!
            maxifac = maximum(AniFac)
            AniFac[AniFac .< 0] .= maxifac
            @show minimum(AniFac)
            
            AniFac[mask_air] .= NaN
        end
        if fabric
            Nx    = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
            Nz    = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
            Fab_x = -Nz./Nx
            Fab_z = ones(size(Nz))
            nrm   = sqrt.(Fab_x.^2 .+ Fab_z.^2)
            Fab_x ./= nrm
            Fab_z ./= nrm
        end
        if field==:Topography
            height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
            Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
            Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
            Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
            Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
            x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
            z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
        end
        if σ1_axis 
            σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        end

        #####################################
        # x-shift
        if xshift > 0.0
            xidx    = Int(floor(ncx*xshift))+1 # x index where to cut
            ηc      = ηc[[xidx:end;1:xidx-1],:]
            ε̇II     = ε̇II[[xidx:end;1:xidx-1],:]
            τII     = τII[[xidx:end;1:xidx-1],:]
            εII     = εII[[xidx:end;1:xidx-1],:]
            if field==:AnisoFactor
                AniFac  = AniFac[[xidx:end;1:xidx-1],:]
            end
            if fabric
                Fab_x   = Fab_x[[xidx:end;1:xidx-1],:]
                Fab_z   = Fab_z[[xidx:end;1:xidx-1],:]
            end
        end

        #####################################

        # Postprocessing deviatoric velocity
        if deviatoricV
            shear_style     = 1 # 0 = pure shear; 1 = simple shear
            bkg_strain_rate = 1 # background (far-field) strain rate

            Vx_bkg = Vx
            Vz_bkg = Vz
            if shear_style == 1
                bla = 2*bkg_strain_rate*zc'
                @show size(Vx)
                @show size(bla)
                Vx_bkg = repeat(bla, outer = nvx+1)
                @show size(Vx_bkg)
            end
            Vx = Vx-Vx_bkg
        end

        # Postprocessing velocity
        Vxc     = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1]); Vxc[mask_air] .= NaN
        Vzc     = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0]); Vzc[mask_air] .= NaN
        Vc      = (Vxc.^2 .+ Vzc.^2).^0.5

        #####################################
        # Postprocessing polar setups
        if polar
            τxx, τzz, τxz, τxz = tensorCartesianToPolar( τxx, τzz, τxz, τxz, xc, zc, ncx, ncz )
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
            ax1 = MoC.Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Density
            ax1 = MoC.Axis(f[1, 1], title = L"Density $ρ$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"Density $ρ$ [kg/m^3]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Viscosity
            ax1 = MoC.Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end 
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress_xx
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xx}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxx, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xx}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress_zz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{zz}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τzz, colormap = (:turbo, α_heatmap)) 
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{zz}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress_xz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xz}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            if polar
                hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxz, colormap = (:turbo, α_heatmap)) 
            else
                hm = MoC.heatmap!(ax1, xv./Lc, zv./Lc, τxz, colormap = (:turbo, α_heatmap)) 
            end
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end  
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xz}$ [Pa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:StrainRate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap), colorrange=(-3, 1.5))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:PlasticStrainrate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) %$(tMytext), sum = %$(round(sum(ε̇pl),sigdigits=4)), max = %$(round(maximum(ε̇pl),sigdigits=4)), min = %$(round(minimum(ε̇pl),sigdigits=4))", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            ε̇pl[mask_air] .= NaN
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Pressure
            ax1 = MoC.Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$P$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Temperature
            ax1 = MoC.Axis(f[1, 1], title = L"$T$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$T$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Velocity
            ax1 = MoC.Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Vc, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Velocity_x
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{x}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xv./Lc, zc./Lc, Vx[:,2:end-1], colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{x}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Velocity_z
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{z}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zv./Lc, Vz[2:end-1,:], colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{z}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Strain
            ax1 = MoC.Axis(f[1, 1], title = L"$\varepsilon_{II}$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εII, colormap = (:turbo, α_heatmap))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\varepsilon_{II}$ []", width = 20, labelsize = 25, ticklabelsize = 14 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:AnisoFactor
            ax1 = MoC.Axis(f[1, 1], title = L"Anisotropy factor at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(AniFac), colormap = (:turbo, α_heatmap), colorrange=(0, 3))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =L"$log_{10}\left(\text{Anisotropy factor}\right)$", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:GrainSize
            ax1 = MoC.Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
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
        end

        if field==:Topography
            ax1 = MoC.Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) %$(tMytext)", xlabel = L"$x$ [m]", ylabel = L"$h$ [km]")
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)

            ax2 = MoC.Axis(f[2, 1], xlabel = L"$x$ [m]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)

            ax3 = MoC.Axis(f[3, 1], xlabel = L"$x$ [m]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)

            @show minimum(Vx_grid)
            @show minimum(Vx_mark)
            @show maximum(Vx_grid)
            @show maximum(Vx_mark)
        end
        
        if velocity 
            #qv,qh = 10,10
            #MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = Δ*1000, lengthscale=Δ*20)
            qv,qh = 6,3
            MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
        end
        if fabric 
            #fac = 1
            qv,qh = fac,fac
            MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Fab_x[1:qh:end,1:qv:end], Fab_z[1:qh:end,1:qv:end], arrowsize = 0, lengthscale=Δ*0.7*fac, align = :center)
            #MoC.arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
        end
        if printfig Print2Disk( f, path, string(field), istep, xshift) end

        if inspector
            DataInspector(f)
        end
        if displayfig
            display(f)
        end
        sleep(nap)
    end
    end
    if printvid
            #FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1920:1080" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`) 
            zres = resol
            xres = resol
            if xshift > 0.0
                FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/shift/'*'.png -vf "scale=$(xres):$(zres)" -c:v libx264 -pix_fmt yuv420p -y "$(movieName)_shift.mov"`)
            else
                FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=$(xres):$(zres)" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`)
            end
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

function Print2Disk( f, path, field, istep, xshift; res=4)
    if xshift > 0.0
        path1 = path*"/_$field/shift/"
    else
        path1 = path*"/_$field/"
    end
    mkpath(path1)
    MoC.save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

function tensorCartesianToPolar( AXX, AZZ, AXZ, AZX, xc, zc, ncx, ncz )
    # tensor coordinate transformation from x-z to r-phi coordinate system
    # A_XZ = tensor to transform
    # A_RP = transformed tensor
    XX = ones(ncx)*zc' # matrix with all x-coords
    ZZ = xc*one(ncz)' # matrix with all z-coords
    ALPHA = asin.( XX ./ ( XX.^2 .+ ZZ.^2 ).^0.5 )
    # Idea for coord transform:
    # R = [cos.(ALPHA) sin.(ALPHA)
    #     -sin.(ALPHA) cos.(ALPHA)]
    # ARP = R' * AXZ * R
    AXZc = v2c(AXZ)
    AZXc = v2c(AZX)
    #ARR = ( AXX.*cos.(ALPHA) .- AZXc.*sin.(ALPHA) ) .* cos.(ALPHA) .- ( AXZc.*cos.(ALPHA) .- AZZ.*sin.(ALPHA) ) .* sin.(ALPHA)
    #APP = ( AXX.*sin.(ALPHA) .+ AZXc.*cos.(ALPHA) ) .* sin.(ALPHA) .+ ( AZZ.*cos.(ALPHA) .+ AXZc.*sin.(ALPHA) ) .* cos.(ALPHA)
    #ARP = ( AXX.*cos.(ALPHA) .- AZXc.*sin.(ALPHA) ) .* sin.(ALPHA) .+ ( AXZc.*cos.(ALPHA) .- AZZ.*sin.(ALPHA) ) .* cos.(ALPHA)
    #APR = ( AXX.*sin.(ALPHA) .+ AZXc.*cos.(ALPHA) ) .* cos.(ALPHA) .- ( AZZ.*cos.(ALPHA) .+ AXZc.*sin.(ALPHA) ) .* sin.(ALPHA)
    ARR = ( AXX.*cos.(ALPHA) .+ AZXc.*sin.(ALPHA) ) .* cos.(ALPHA) .+ ( AXZc.*cos.(ALPHA) .+ AZZ.*sin.(ALPHA) ) .* sin.(ALPHA)
    APP = ( AZZ.*cos.(ALPHA) .- AXZc.*sin.(ALPHA) ) .* cos.(ALPHA) .- ( AZXc.*cos.(ALPHA) .- AXX.*sin.(ALPHA) ) .* sin.(ALPHA)
    ARP = ( AXZc.*cos.(ALPHA) .+ AZZ.*sin.(ALPHA) ) .* cos.(ALPHA) .- ( AXX.*cos.(ALPHA) .+ AZXc.*sin.(ALPHA) ) .* sin.(ALPHA)
    APR = ( AZXc.*cos.(ALPHA) .- AXX.*sin.(ALPHA) ) .* cos.(ALPHA) .+ ( AZZ.*cos.(ALPHA) .- AXZc.*sin.(ALPHA) ) .* sin.(ALPHA)
    check = 1
    if check == 1
        AIIXZ = (0.5*( AXX.^2 .+ AZZ.^2 .+ AXZc.^2 .+ AZXc.^2 )).^0.5
        AIIRP = (0.5*( ARR.^2 .+ APP.^2 .+ ARP.^2 .+ APR.^2 )).^0.5
        AIIdiff = AIIRP .- AIIXZ
        @show maximum(AIIdiff)
    end
    return ARR, APP, ARP, APR
end

function v2c(Av) # interpolation vertex to center
    return 0.25 * ( Av[1:end-1,1:end-1] + Av[2:end,1:end-1] + Av[2:end,1:end-1] + Av[2:end,2:end] )
end

#main()

#run(`ffmpeg -framerate 2 -i /Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/_Phases/Phases%05d.png -c libx264 -pix_fmt yuv420p -y out_movie.mp4`)

# Working command
#ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4


# Set the path to your files
#main("/users/whalter1/work/aniso_fix/A1_1000/", 2610, 2610, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A1_500/", 11410, 11410, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A1_250/", 0, 0, 1, 0.5, 10)

#main("/users/whalter1/work/aniso_fix/A2_1000/", 2790, 3100, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A2_500/", 2680, 6190, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A2_250/", 5230, 10000, 1, 0.5, 10)
#
#main("/users/whalter1/work/aniso_fix/A3_1000/", 3320, 3640, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A3_500/", 5000, 7500, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A3_250/", 0, 5000, 1, 0.5, 10)
#
#main("/users/whalter1/work/aniso_fix/A4_1000/", 2770, 3080, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A4_500/", 5000, 7500, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A4_250/", 0, 5000, 1, 0.5, 10)
#
#main("/users/whalter1/work/aniso_fix/A5_1000/", 3460, 3880, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A5_500/", 5000, 7500, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/A5_250/", 0, 5000, 1, 0.5, 10)

#main("/users/whalter1/work/aniso_fix/B1_1000/", 0, 0, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B1_500/", 0, 0, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B1_250/", 0, 0, 1, 0.5, 10)
#
#main("/users/whalter1/work/aniso_fix/B2_1000/", 0, 0, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B2_500/", 0, 0, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B2_250/", 0, 0, 1, 0.5, 10)
#
#main("/users/whalter1/work/aniso_fix/B3_1000/", 0, 0, 4, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B3_500/", 0, 0, 2, 0.5, 10)
#main("/users/whalter1/work/aniso_fix/B3_250/", 0, 0, 1, 0.5, 10)

#
#main("/users/whalter1/work/aniso_fix/A2_500_PSALE/", 200, 200, 1, 0.0, 10)
#main("/users/whalter1/work/aniso_fix/A2_250_PSALE/", 80, 80, 1, 0.0, 10)
#
#main("/users/whalter1/work/aniso_fix/A2_500_PSDirichlet/", 280, 280, 2, 0.0, 10)
#main("/users/whalter1/work/aniso_fix/A2_250_PSDirichlet/", 461, 800, 1, 0.0, 10)

#main("/users/whalter1/work/aniso_fix/MWE_100_PSALE/", 1556, 1556, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/MWE_100_PSDirichlet/", 0, 100, 1, 0.0, 1)

#main("/users/whalter1/work/aniso_fix/A3_250_PSALE/", 915, 920, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A4_250_PSALE/", 193, 198, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A5_250_PSALE/", 168, 173, 1, 0.0, 1)

#main("/users/whalter1/work/aniso_fix/A2_250_PS180/", 85, 90, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A2_200_PSALE/", 0, 94, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A2_250_PSALE_another/", 0, 79, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A2_100_PSALE_another/", 0, 314, 1, 0.0, 1)
#main("/users/whalter1/work/aniso_fix/A2_150_PSALE/", 0, 700, 1, 0.0, 1)

#main("/users/whalter1/work/aniso_fix/MWE_500/", 400, 400, 1, 0.5, 10)


