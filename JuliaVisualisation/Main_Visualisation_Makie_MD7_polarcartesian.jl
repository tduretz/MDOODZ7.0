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

    #GeoMod23 unused
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/anievo_cartes_rand/"
    
    #GeoMod23 used
    path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/anievo_polar_rand/"#; file_start = 750
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/normal_cartes/"; file_start = 840 
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/normal_polar/"; file_start = 830
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/slim100_cartes/"; file_start = 660
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/slim100_polar/"; file_start = 640
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/slim300_cartes/"; file_start = 830
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/slim300_polar/"; file_start = 840
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/surface_cartes/"; file_start = 830
    #path = "/mnt/c/Users/whalter1/Desktop/GeoMod23/polarcartesian/surface_polar/"; file_start = 840

    path = "/users/whalter1/work/MDOODZ_polar2/hres/anievo_polar_rand/"


    # File numbers
    file_start = 750
    file_step  = 10
    file_end   = 750
    #file_end   = file_start

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
    #field = :AnisoFactor # anisotropy factor
    #field = :GrainSize # doesn't function properly yet
    #field = :Topography # doesn't function properly yet

    # Switches
    displayfig  = false     # display figures on screen (not possible on any "headless server")
    printfig    = true      # print figures to disk
    printvid    = false      # print a video (printfig must be on)
    ph_contours = false     # add phase contours
    T_contours  = false     # add temperature contours
    fabric      = false      # add fabric quiver (normal to director)
    velocity    = false     # add velocity quiver
    activeV     = false     # considers the "active velocity" (the velocity minus the background velocity)
    α_heatmap   = 0.85      # transparency of heatmap 
    σ1_axis     = false
    polar       = false     # for curved setups using Earth's curvature (transforms tensors from x-z to r-phi coord system)
    nap         = 0.3       # pause for animation 
    resol       = 1600
    inspector   = false     # (not possible on any "headless server")
    framerate   = 2         # Frame rate for the movie
    movieName   = "$(path)/_$(field)/$(field)"  # Name of the movie
    xcrop       = true      # crops on both sides from lateral boundary
    xcropfract  = 0.1       # fraction to crop away on either side

    # Scaling
    Lc = 1e3 # km
    tc = My
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
        tMy    = round(t/tc, digits=2)
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

        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));                  mask_air = ph .== -1.00; ph[mask_air]  .= NaN
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

        # Postprocessing Active Velocity
        # Active Velocity = Velocity - Background Velocity
        if activeV
            shear_style     = 1 # 0 = pure shear; 1 = simple shear
            bkg_strain_rate = 1 # background (far-field) strain rate
            if shear_style == 1
                Vx_bkg = 2*bkg_strain_rate*ones(nvx)*zc'
                @show size(Vx_bkg)
                Vx[:,2:end-1] .-= Vx_bkg
            elseif shear_style == 0
                # not implemented yet
                # feel free to add :)
            end
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
        cmap    = zeros(RGB{Float64}, 3)
        cmap[1] = RGBA{Float64}(178/255,  86/255,  82/255, 1.)
        cmap[2] = RGBA{Float64}( 48/255,  20/255,  55/255, 1.)
        cmap[3] = RGBA{Float64}( 98/255, 118/255, 186/255, 1.)
        #cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
        #cmap[2] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
        #cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
        #cmap[4] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        #cmap[5] = RGBA{Float64}(217/255, 099/255, 097/255, 1.) 
        #cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
        #cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        phase_colors = MoC.cgrad(cmap, length(cmap), categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################
        
        # Cropping indices
        xidx_start = Int(floor(ncx*xcropfract))+1
        xidx_end   = Int(floor(ncx*(1-xcropfract)))+1
        zidx_start = 1
        zidx_end   = ncz

        # Crop
        xc      = xc[xidx_start:xidx_end]
        zc      = zc[zidx_start:zidx_end]
        xc_hr   = xc_hr[xidx_start:xidx_end]
        zc_hr   = zc_hr[zidx_start:zidx_end]

        ph_hr   = ph_hr[xidx_start:xidx_end,zidx_start:zidx_end]
        ph      = ph[xidx_start:xidx_end,zidx_start:zidx_end]
        ηc      = ηc[xidx_start:xidx_end,zidx_start:zidx_end]
        τII     = τII[xidx_start:xidx_end,zidx_start:zidx_end]
        τxx     = τxx[xidx_start:xidx_end,zidx_start:zidx_end]
        Vc      = Vc[xidx_start:xidx_end,zidx_start:zidx_end]
        ε̇II     = ε̇II[xidx_start:xidx_end,zidx_start:zidx_end]
        εII     = εII[xidx_start:xidx_end,zidx_start:zidx_end]
        ε̇pl     = ε̇pl[xidx_start:xidx_end,zidx_start:zidx_end]
        if field==:AnisoFactor
            AniFac     = AniFac[xidx_start:xidx_end,zidx_start:zidx_end]
        end
        if fabric
            Fab_x     = Fab_x[xidx_start:xidx_end,zidx_start:zidx_end]
            Fab_z     = Fab_z[xidx_start:xidx_end,zidx_start:zidx_end]
        end
        #@show size(τxx)
        #@show size(xc_hr)
        #@show size(zc_hr)

        #####################################
        zres = resol
        xres = Lx/Lz*resol
        f = MoC.Figure(resolution = (xres, zres), fontsize=40)

        if field==:Phases
            ax1 = MoC.Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, ph, colormap = phase_colors)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"Density $ρ$ [kg/m^3]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end            
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [Pa]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress_xx
            ax1 = MoC.Axis(f[1, 1], title = L"Deviatoric horizontal stress at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxx/1e6, colormap = (:turbo, α_heatmap), colorrange = (-400, 400)) 
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
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xx}$ [MPa]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end         
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{zz}$ [Pa]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Stress_xz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xz}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
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
            MoC.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{xz}$ [Pa]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:StrainRate
            ax1 = MoC.Axis(f[1, 1], title = L"Strain-rate at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap), colorrange = (-19, -12.5))
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
            MoC.Colorbar(f[1, 2], hm, label =  L"$log_{10}\left(\dot{\varepsilon}_\textrm{II}\right)$", ticks = (-20:2:-12), width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:PlasticStrainrate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma, sum = %$(round(sum(ε̇pl),sigdigits=4)), max = %$(round(maximum(ε̇pl),sigdigits=4)), min = %$(round(minimum(ε̇pl),sigdigits=4))", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #ε̇pl[mask_air] .= NaN
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
            MoC.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$P$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label =  L"$T$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Velocity
            ax1 = MoC.Axis(f[1, 1], title = L"Velocity at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Vc, colormap = (:turbo, α_heatmap))#, colorrange = (1e-10, 3e-9))
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
            MoC.Colorbar(f[1, 2], hm, label = L"$V$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{x}$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            if σ1_axis
                MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
            end  
            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.Colorbar(f[1, 2], hm, label = L"$V_{z}$ [s$^{-1}$]", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:Strain
            ax1 = MoC.Axis(f[1, 1], title = L"Cumulated strain at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εII, colormap = (:turbo, α_heatmap), colorrange = (0, 150))
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
            MoC.Colorbar(f[1, 2], hm, label = L"$\varepsilon_{II}$", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
        end

        if field==:AnisoFactor
            ax1 = MoC.Axis(f[1, 1], title = L"Anisotropy factor $\delta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(AniFac), colormap = (:turbo, α_heatmap), colorrange = (0.0, 3.0))
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, AniFac, colormap = (:turbo, α_heatmap))
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
            MoC.Colorbar(f[1, 2], hm, label =L"$log_{10}\delta$", ticks = (0.0:0.5:3.0), width = 20, labelsize = 40, ticklabelsize = 40 )
            #MoC.Colorbar(f[1, 2], hm, label =L"$\delta$", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
            MoC.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 40, ticklabelsize = 40 )
            MoC.colgap!(f.layout, 20)
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
        
        if velocity 
            #qv,qh = 10,10
            #MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = Δ*1000*2, lengthscale=Δ*20*2)
            qv,qh = 6,3
            MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Vxc[1:qh:end,1:qv:end], Vzc[1:qh:end,1:qv:end], arrowsize = 1e1, lengthscale=Δ/1.5*1e10)
        end
        if fabric 
            #fac = 3 # for usual
            fac = 1 # for zoom
            qv,qh = fac,fac
            MoC.arrows!(ax1, xc[1:qh:end]./Lc, zc[1:qv:end]./Lc, Fab_x[1:qh:end,1:qv:end], Fab_z[1:qh:end,1:qv:end], arrowsize = 0, lengthscale=Δ*0.7*fac)
            #MoC.arrows!(ax1, xc./Lc, zc./Lc, Nz, Nx, arrowsize = 0, lengthscale=Δ/1.5)
        end
        if printfig Print2Disk( f, path, string(field), istep) end

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

        filename = string(path, @sprintf("Output%05d.gzip.h5", file_start))
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")
        xmin, xmax = xv[1], xv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc

        zres = resol
        xres = Int(floor(Lx/Lz*resol))+1 # solve issue later: Lx, Lz not defined in outer loop...
        #xres = resol
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=$(xres):$(zres)" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`)
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
    ARR = ( AXX.*cos.(ALPHA) .- AZXc.*sin.(ALPHA) ) .* cos.(ALPHA) .- ( AXZc.*cos.(ALPHA) .- AZZ.*sin.(ALPHA) ) .* sin.(ALPHA)
    APP = ( AXX.*sin.(ALPHA) .+ AZXc.*cos.(ALPHA) ) .* sin.(ALPHA) .+ ( AZZ.*cos.(ALPHA) .+ AXZc.*sin.(ALPHA) ) .* cos.(ALPHA)
    ARP = ( AXX.*cos.(ALPHA) .- AZXc.*sin.(ALPHA) ) .* sin.(ALPHA) .+ ( AXZc.*cos.(ALPHA) .- AZZ.*sin.(ALPHA) ) .* cos.(ALPHA)
    APR = ( AXX.*sin.(ALPHA) .+ AZXc.*cos.(ALPHA) ) .* cos.(ALPHA) .- ( AZZ.*cos.(ALPHA) .+ AXZc.*sin.(ALPHA) ) .* sin.(ALPHA)
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

main()

#run(`ffmpeg -framerate 2 -i /Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/_Phases/Phases%05d.png -c libx264 -pix_fmt yuv420p -y out_movie.mp4`)

# Working command
#ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4