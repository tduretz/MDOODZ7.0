import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie")
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

function main(path, file_start, file_end_max,fac, xshift, file_step)

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
    field = :Pressure
    #field = :Temperature
    #field = :Velocity
    #field = :Velocity_x
    #field = :Velocity_z
    #field = :Strain # cumulated strain
    #field = :StrainDev # deviatoric cumulated strain
    #field = :Strain_xx # normal strain
    #field = :Strain_zz # normal strain
    #field = :Strain_xz # shear strain
    #field = :FabricStrain_xx # normal strain along foliation (i.e. εxxp)
    #field = :FabricStrain_zz # normal strain along director (i.e. εzzp)
    #field = :FabricStrain_xz # shear strain along and across foliation (i.e. εxzp)
    #field = :Fxx
    #field = :Fzz
    #field = :Fxz
    #field = :Fzx
    #field = :Fxxp
    #field = :Fzzp
    #field = :Fxzp # (can be used for aniso factor evolv)
    #field = :Fzxp
    #field = :FxzpDev # (for the others, deviatorics is (assumed/defined as) equal to non-deviatorics)
    #field = :AnisoFactor # anisotropy factor
    #field = :GrainSize # doesn't function properly yet
    #field = :Topography # doesn't function properly yet
    #field = :ViscosityReduction
    #field = :StrainLocalisation
    #field = :XVelocityLocalisation
    #field = :XYEtaLocalisation

    # Switches
    Δγ          = 0.5       # if Δγ > 0, then only creates every Δγ-th plot and skips in-between; if Δγ = 0, then creates all plots, no skipping
    displayfig  = false     # display figures on screen (not possible on any "headless server")
    printfig    = true      # print figures to disk
    printvid    = true      # print a video (printfig must be on)
    ph_contours = false     # add phase contours
    T_contours  = false     # add temperature contours
    fabric      = true      # add fabric quiver (normal to director)
    velocity    = false     # add velocity quiver
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
    TimeSeries  = true
    phaseProportions = false

    # File numbers
    #file_step  = 10
    # 1) find file_start
    # Check index of last (highest) already generated picture file (.png)
    path_to_pics = joinpath(path,"_$(field)")
    if xshift > 0.0; path_to_pics = joinpath(path_to_pics,"shift"); end
    @show path_to_pics
    if ~isdir(path_to_pics) || ~isfile("$(path_to_pics)/$(field)$(lpad(0,5,"0")).png")
        file_start = 0
    else
        firstIsGood = false
        k = 0
        while firstIsGood == false
            #path_file_start = "$(path_to_pics)/$(field)$(lpad(k,5,"0")).png"
            path_file_above = "$(path_to_pics)/$(field)$(lpad(k+file_step,5,"0")).png"
            if isfile(path_file_above)
                k += file_step
            else
                file_start = k + file_step
                firstIsGood = true
            end
        end
    end
    # 2) find file_end
    # Check through all Output files (*.gzip.h5)
    Files  = sort(filter(x->endswith(x,".h5"), readdir(path)))
    file_end  = parse(Int64, replace(Files[end] , r"Output"=>"", r".gzip.h5"=>""))
    if file_end > file_end_max; file_end = file_end_max; end
    if file_start >= file_end && Δγ == 0; file_start = file_end; printvid = false; printfig = false; end # don't make video if (probably) already done...

    @show path

    @show Δγ
    if Δγ != 0
        filename  = string(path, @sprintf("Output%05d.gzip.h5", file_end))
        Time_time = Float64.(ExtractData( filename, "/TimeSeries/Time_time"));
        Time_time = [Time_time[1]; Time_time[Time_time .> 0.0]]
        γ_time    = Time_time .* 2
        list = Int64.([0;find_indices(γ_time, Δγ)])
        list = Int64.(round.(list./10).*10)
    end

    # Scaling
    Lc = 1.
    tc = s; tMytext= "s"
    Vc = 1e-9

    # Time loop
    if displayfig==1 || printfig==1
    if Δγ != 0
        filelist=list
        @show list
    else
        filelist=[file_end;file_start:file_step:file_end]  
        @show file_start, file_end, file_step
    end
    for istep=filelist
    #for istep = 0    
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
        Δt     = model[8]
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        γ = t * 2 # here...
        γtxt = @sprintf("%.1f", γ)

        @show "Time step" istep, γtxt
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
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), nvx, (ncz+2)))      
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), nvz))      
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz));                   τxx[mask_air] .= NaN
        τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz));                   τzz[mask_air] .= NaN
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        τxzc  = 0.25*(τxz[1:end-1,1:end-1].+τxz[2:end,1:end-1].+τxz[1:end-1,2:end].+τxz[2:end,2:end]);  τxzc[mask_air] .= NaN
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇zz   = Float64.(reshape(ExtractData( filename, "/Centers/ezzd"), ncx, ncz))
        ε̇yy   = -(ε̇xx .+ ε̇zz)
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        εII   = Float64.(reshape(ExtractData( filename, "/Centers/strain"), ncx, ncz));                 εII[mask_air] .= NaN
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
        Vxc   = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc   = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])
        if σ1_axis 
            σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        end
        if TimeSeries
            Time_time       = Float64.(ExtractData( filename, "/TimeSeries/Time_time")); 
            Short_time      = Float64.(ExtractData( filename, "/TimeSeries/Short_time")); 
            Work_time       = Float64.(ExtractData( filename, "/TimeSeries/Work_time")); 
            Uthermal_time   = Float64.(ExtractData( filename, "/TimeSeries/Uthermal_time")); 
            Uelastic_time   = Float64.(ExtractData( filename, "/TimeSeries/Uelastic_time")); 
            T_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/T_mean_time")); 
            P_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/P_mean_time")); 
            τxx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxxd_mean_time")); 
            τzz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/szzd_mean_time")); 
            sxz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxz_mean_time")); 
            τII_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Tii_mean_time")); 
            ε̇xx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exxd_mean_time")); 
            ε̇zz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/ezzd_mean_time")); 
            ε̇xz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exz_mean_time")); 
            ε̇II_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Eii_mean_time")); 
        end

        u=unique(ph)
        nbph=length(u)
        if phaseProportions
            D=Dict([(i,count(x->x==i,ph)) for i in u])
            dtot = sum(values(D))
            for i in u
                println("count for ph = $(i) is $(D[i]) \t $(D[i]/dtot*100) %")
            end
            println("count for ph total is $(dtot) \t $(dtot/dtot*100) %")
        end

        #####################################
        # x-shift
        if xshift > 0.0
            # x-coord on center location
            xidx    = Int(floor(ncx*xshift))+1 # x index where to cut

            ph      = ph[[xidx:end;1:xidx-1],:]
            ph_hr   = ph_hr[[xidx:end;1:xidx-1],:]
            ηc      = ηc[[xidx:end;1:xidx-1],:]
            ρc      = ρc[[xidx:end;1:xidx-1],:]
            P       = P[[xidx:end;1:xidx-1],:]
            T       = T[[xidx:end;1:xidx-1],:]
            d       = d[[xidx:end;1:xidx-1],:]
            ε̇pl     = ε̇pl[[xidx:end;1:xidx-1],:]
            Vz      = Vz[[xidx:end;1:xidx-1],:] # not accurate ... to solve later...!
            τxx     = τxx[[xidx:end;1:xidx-1],:]
            τzz     = τzz[[xidx:end;1:xidx-1],:]
            τxzc    = τxzc[[xidx:end;1:xidx-1],:]
            ε̇xx     = ε̇xx[[xidx:end;1:xidx-1],:]
            εII     = εII[[xidx:end;1:xidx-1],:]
            ε̇II     = ε̇II[[xidx:end;1:xidx-1],:]
            τII     = τII[[xidx:end;1:xidx-1],:]
            if field==:AnisoFactor
                AniFac  = AniFac[[xidx:end;1:xidx-1],:]
            end
            if fabric
                Nx      = Nx[[xidx:end;1:xidx-1],:]
                Nz      = Nz[[xidx:end;1:xidx-1],:]
                Fab_x   = Fab_x[[xidx:end;1:xidx-1],:]
                Fab_z   = Fab_z[[xidx:end;1:xidx-1],:]
            end
            # x-coord on vertex location
            xidx    = Int(floor(nvx*xshift))+1 # x index where to cut
            Vx      = Vx[[xidx:end;1:xidx-1],:]
            τxz     = τxz[[xidx:end;1:xidx-1],:]
            ε̇xz     = ε̇xz[[xidx:end;1:xidx-1],:]
        end

        #####################################

        # Postprocessing deviatoric velocity
        if field==:XVelocityLocalisation; deviatoricV = true; end
        if deviatoricV
            shear_style     = 1 # 0 = pure shear; 1 = simple shear
            bkg_strain_rate = 1 # background (far-field) strain rate

            Vx_bkg = Vx
            Vz_bkg = Vz
            if shear_style == 1
                bla = 2*bkg_strain_rate*zc'
                bla_first = bla[1]-(bla[2]-bla[1])
                bla_last  = bla[end]+(bla[end]-bla[end-1])
                bla = [bla_first; bla'; bla_last]'
                Vx_bkg = repeat(bla, outer = nvx)
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
        cmap[1] = RGBA{Float64}(240/255, 240/255, 240/255, 1.)  
        cmap[2] = RGBA{Float64}(  0/255, 142/255, 198/255, 1.)  
        #cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
        #cmap[2] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
        cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
        cmap[4] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cmap[5] = RGBA{Float64}(217/255, 099/255, 097/255, 1.) 
        cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
        cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cmap    = cmap[1:nbph]
        phase_colors = MoC.cgrad(cmap, length(cmap), categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################

        f = MoC.Figure(resolution = (Lx/Lz*resol, resol), fontsize=25)
        
        if field==:XYEtaLocalisation

            Vx_HorizMean = mean(Vx, dims = 1)

            ax1 = MoC.Axis(f[1, 1], title = L"\eta_{xy} Localisation at $γ$ = %$(round(γ, digits=2))", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.plot!([0.0,0.0], [-0.5,0.5], color=:red , markersize = 2, linewidth=1)
            hm = MoC.plot!(Vx_HorizMean, zc./Lc , color=:blue, markersize = 2, linewidth=1)

            #MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            #MoC.xlims!(ax1, -0.5, 0.5)
            MoC.ylims!(ax1, -0.5, 0.5)

            fabric = false
        end
                
        if field==:XVelocityLocalisation

            Vx_HorizMean = mean(Vx, dims = 1)

            ax1 = MoC.Axis(f[1, 1], title = L"X Velocity Localisation at $γ$ = %$(round(γ, digits=2))", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.plot!([0.0,0.0], [-0.5,0.5], color=:red , markersize = 2, linewidth=1)
            hm = MoC.plot!(Vx_HorizMean, zc./Lc , color=:blue, markersize = 2, linewidth=1)

            #MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            #MoC.xlims!(ax1, -0.5, 0.5)
            MoC.ylims!(ax1, -0.5, 0.5)

            fabric = false
        end
        
        if field==:StrainLocalisation
            if istep==file_start
                global npoints      = (ncz+2)*1000
                global xPoints      = zeros(npoints)
                global zPoints      = Float64.(LinRange(zmin, zmax, npoints))
                global xPoints_PS   = zeros(npoints) # pure shear x coordinates at given zPoints
                global xPoints_PS0  = zeros(npoints) # xPoints_PS old
                global Δt = t
            else
                global Δt = t - t0
            end
            global t0 = t

            # get indices of closest VX and VZ points
            VX_xi = Int.(round.((xPoints .+ Lx/2        ) ./  Lx       .* (nvx-1))) .+ 1
            VX_zi = Int.(round.((zPoints .+ Lz/2 .- Δz/2) ./ (Lz + Δz) .* (ncz+1))) .+ 2
            VZ_xi = Int.(round.((xPoints .+ Lx/2 .- Δx/2) ./ (Lx + Δx) .* (ncx+1))) .+ 2
            VZ_zi = Int.(round.((zPoints .+ Lz/2        ) ./  Lz       .* (nvz-1))) .+ 1

            VX_update = Float64.(Vx[CartesianIndex.(VX_xi, VX_zi)])
            VZ_update = Float64.(Vz[CartesianIndex.(VZ_xi, VZ_zi)])

            global xPoints      = xPoints .+ VX_update * Δt
            global zPoints      = zPoints .+ VZ_update * Δt
            global xPoints_PS   = zPoints .* γ

            xPoints     = xPoints    .- (xPoints    .> Lx/2) .* Lx .+ (xPoints    .< -Lx/2) .* Lx
            zPoints     = zPoints    .- (zPoints    .> Lz/2) .* Lz .+ (zPoints    .< -Lz/2) .* Lz
            while xPoints_PS != xPoints_PS0
                xPoints_PS0 = xPoints_PS
                xPoints_PS  = xPoints_PS0 .- (xPoints_PS0 .> Lx/2) .* Lx .+ (xPoints_PS0 .< -Lx/2) .* Lx
            end

            ax1 = MoC.Axis(f[1, 1], title = L"Strain Localisation at $γ$ = %$(round(γ, digits=2))", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.plot!(xPoints_PS, zPoints, color=:red , markersize = 2, linewidth=1)
            hm = MoC.plot!(xPoints   , zPoints, color=:blue, markersize = 2, linewidth=1)

            MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
            MoC.xlims!(ax1, -0.5, 0.5)
            MoC.ylims!(ax1, -0.5, 0.5)
        end

        if field==:Phases
            if phaseProportions
                ax1 = MoC.Axis(f[1, 1], title = L"Phases at $γ$ = %$(γtxt), inclusion area = %$(D[1]/dtot*100) %", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            else
                ax1 = MoC.Axis(f[1, 1], title = L"Phases at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            end
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, ph, colormap = phase_colors)
            colorbarLabel = "Phases"
        end

        if field==:Density
            ax1 = MoC.Axis(f[1, 1], title = L"Density $ρ$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:amp, α_heatmap))
            colorbarLabel = L"Density $ρ$ [kg/m^3]"
        end

        if field==:Viscosity
            ax1 = MoC.Axis(f[1, 1], title = L"$\eta$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            colorbarLabel = L"$\eta$ [Pa.s]"
        end

        if field==:Stress
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:amp, α_heatmap), colorrange=(0.0, 4.0e3))
            colorbarLabel = L"$\tau_\textrm{II}$ [Pa]"
        end

        if field==:Stress_xx
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xx}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxx, colormap = (:balance, α_heatmap), colorrange=(-4.0e3, 4.0e3))
            colorbarLabel = L"$\tau_\textrm{xx}$ [Pa]"
        end

        if field==:Stress_zz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{zz}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τzz, colormap = (:balance, α_heatmap), colorrange=(-4.0e3, 4.0e3))
            colorbarLabel = L"$\tau_\textrm{zz}$ [Pa]"
        end

        if field==:Stress_xz
            ax1 = MoC.Axis(f[1, 1], title = L"$\tau_\textrm{xz}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            if polar
                hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, τxz, colormap = (:amp, α_heatmap), colorrange=(0.0, 4.0e3))
            else
                hm = MoC.heatmap!(ax1, xv./Lc, zv./Lc, τxz, colormap = (:amp, α_heatmap), colorrange=(0.0, 4.0e3))
            end
            colorbarLabel = L"$\tau_\textrm{xz}$ [Pa]"
        end

        if field==:StrainRate
            γtxt = @sprintf("%.1f", γ)
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap), colorrange=(-3, 1.5))
            colorbarLabel = L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]"
        end

        if field==:PlasticStrainrate
            ax1 = MoC.Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) %$(tMytext), sum = %$(round(sum(ε̇pl),sigdigits=4)), max = %$(round(maximum(ε̇pl),sigdigits=4)), min = %$(round(minimum(ε̇pl),sigdigits=4))", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            ε̇pl[mask_air] .= NaN
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            colorbarLabel = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]"
        end

        if field==:Pressure
            ax1 = MoC.Axis(f[1, 1], title = L"$P$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, P, colormap = (:balance, α_heatmap), colorrange=(-4.0e3, 4.0e3))
            colorbarLabel = L"$P$"
        end

        if field==:Temperature
            ax1 = MoC.Axis(f[1, 1], title = L"$T$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (:turbo, α_heatmap))
            colorbarLabel = L"$T$"
        end

        if field==:Velocity
            ax1 = MoC.Axis(f[1, 1], title = L"$V$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Vc, colormap = (:amp, α_heatmap), colorrange=(0.0, 0.2))
            colorbarLabel = L"$V$ [s$^{-1}$]"
        end

        if field==:Velocity_x
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{x}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xv./Lc, zc./Lc, Vx[:,2:end-1], colormap = (:balance, α_heatmap), colorrange=(-0.2, 0.2))
            colorbarLabel = L"$V_{x}$ [s$^{-1}$]"
        end

        if field==:Velocity_z
            ax1 = MoC.Axis(f[1, 1], title = L"$V_{z}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zv./Lc, Vz[2:end-1,:], colormap = (:balance, α_heatmap), colorrange=(-0.2, 0.2))
            colorbarLabel = L"$V_{z}$ [s$^{-1}$]"
        end

        if field==:Strain
            ax1 = MoC.Axis(f[1, 1], title = L"$\varepsilon_{II}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εII, colormap = (:amp, α_heatmap), colorrange=(0.0,20.0))
            colorbarLabel = L"$\varepsilon_{II}$ []"
        end

        if field==:StrainDev
            ax1 = MoC.Axis(f[1, 1], title = L"Deviatoric $\varepsilon_{II}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εII.-γ, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"Deviatoric $\varepsilon_{II}$ []"
        end

        if field==:Strain_xx
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xxd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxx), maximum(εxx)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εxx, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{xxd}$ [-]"
        end

        if field==:Strain_zz
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{zzd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εzz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εzz), maximum(εzz)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εzz, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{zzd}$ [-]"
        end

        if field==:Strain_xz
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xzd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxz), maximum(εxz)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εxz, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{xzd}$ [-]"
        end

        if field==:FabricStrain_xx
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xxd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxxp), maximum(εxxp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εxxp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{xxd}$ [-]"
        end

        if field==:FabricStrain_zz
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{zzd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εzzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εzzp), maximum(εzzp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εzzp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{zzd}$ [-]"
        end

        if field==:FabricStrain_xz
            ax1 = MoC.Axis(f[1, 1], title = L"${\varepsilon}_\textrm{xzd}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(εxzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(εxzp), maximum(εxzp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, εxzp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"${\varepsilon}_\textrm{xzd}$ [-]"
        end

        if field==:Fxx
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{xx}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxx), maximum(Fxx)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fxx, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{xx}}$ [-]"
        end

        if field==:Fzz
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{zz}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzz), maximum(Fzz)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fzz, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{zz}}$ [-]"
        end

        if field==:Fxz
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{xz}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxz), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxz), maximum(Fxz)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fxz, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{xz}}$ [-]"
        end

        if field==:Fzx
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{zx}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzx), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzx), maximum(Fzx)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fzx, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{zx}}$ [-]"
        end

        if field==:Fxxp
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{xxp}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fxxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxxp), maximum(Fxxp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fxxp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{xxp}}$ [-]"
        end

        if field==:Fzzp
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{zzp}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzzp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzzp), maximum(Fzzp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fzzp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{zzp}}$ [-]"
        end

        if field==:Fxzp
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{xzp}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fxzp), maximum(Fxzp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fxzp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{xzp}}$ [-]"
        end

        if field==:Fzxp
            ax1 = MoC.Axis(f[1, 1], title = L"$\textrm{F_{zxp}}$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            #hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(Fzxp), colormap = (:turbo, α_heatmap), colorrange=(-4.,1.))
            @show minimum(Fzxp), maximum(Fzxp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fzxp, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"$\textrm{F_{zxp}}$ [-]"
        end

        if field==:FxzpDev
            ax1 = MoC.Axis(f[1, 1], title = L"Deviatoric $γ$ along fabric at bkgd $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            @show minimum(Fxzp), maximum(Fxzp)
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, Fxzp.-γ, colormap = (:balance, α_heatmap), colorrange=(-10.0,10.0))
            colorbarLabel = L"Deviatoric $\textrm{F_{xzp}}$ [-]"
        end

        if field==:AnisoFactor
            ax1 = MoC.Axis(f[1, 1], title = L"Anisotropy factor at $γ$ = %$(γtxt)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(AniFac), colormap = (:turbo, α_heatmap), colorrange=(0, 3))
            colorbarLabel = L"$log_{10}\left(\text{Anisotropy factor}\right)$"
        end

        if field==:GrainSize
            ax1 = MoC.Axis(f[1, 1], title = L"$d$ at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = MoC.heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            colorbarLabel = "d"
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
        end

        if field==:Topography
            ax1 = MoC.Axis(f[1, 1], title = L"Topography at $γ$ = %$(γtxt)", xlabel = L"$x$ [m]", ylabel = L"$h$ [km]")
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

        if field==:ViscosityReduction
            ax1 = MoC.Axis(f[1, 1], title = L"Reduction of effective shear viscosity at $γ$ = %$(γtxt)", xlabel = L"shear deformation γ", ylabel = L"η_{s eff} / η_{s init}")
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
        
        if ph_contours 
            contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
        end
        if T_contours 
            contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white)  
        end
        if σ1_axis
            MoC.arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5)
        end  
        MoC.colsize!(f.layout, 1, MoC.Aspect(1, Lx/Lz))
        MoC.Colorbar(f[1, 2], hm, label =  colorbarLabel, width = 20, labelsize = 25, ticklabelsize = 14 )
        MoC.colgap!(f.layout, 20)
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

        if Δγ != 0; pathsave = path*"/_dg_$Δγ/"; else pathsave = path end
        if printfig Print2Disk( f, pathsave, string(field), istep, xshift) end

        if inspector
            DataInspector(f)
        end
        if displayfig
            display(f)
        end
        sleep(nap)
    end
    end
    if printvid && Δγ == 0
            #FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1920:1080" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`) 
            zres = resol
            xres = resol
            if xshift > 0.0
                FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/shift/'*'.png -vf "scale=$(xres):$(zres)" -c:v libx264 -pix_fmt yuv420p -y "$(movieName)_shift.mov"`)
            else
                FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=$(xres):$(zres)" -c:v libx264 -pix_fmt yuv420p -y "$(movieName).mov"`)
            end
        elseif printvid
            print("\nNo video printed, because Δγ != 0\n")
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

function find_indices(Time_time, Δγ)
    indices = []
    max_multiple = floor(maximum(Time_time) / Δγ)
    for multiple in 1:max_multiple
        threshold = multiple * Δγ
        index = findfirst(x -> x > threshold, Time_time)
        push!(indices, index)
    end
    return indices
end

#main()

#run(`ffmpeg -framerate 2 -i /Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/_Phases/Phases%05d.png -c libx264 -pix_fmt yuv420p -y out_movie.mp4`)

# Working command
#ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4

# Set the path to your files
localTestOnly = false
aniSetupsOnly = true
file_step = 10

if localTestOnly == false
if aniSetupsOnly == false
    # all iso setups
    
    main("/users/whalter1/work/aniso_fix/B19_2000/", 0, 100000, 8, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B19_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B19_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B6_2000/", 0, 100000, 8, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B5_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B6_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B6_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B7_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B8_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B9_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B10_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B11_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B12_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B13_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B14_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B15_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B16_2000/", 0, 100000, 8, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B16_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B16_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B17_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B18_1000/", 0, 100000, 4, 0.5, file_step)

    main("/users/whalter1/work/aniso_fix/A1_1000/", 2610, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/A6_1000/", 3460, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/A7_1000/", 3460, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B1_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B2_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B2_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B3_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B3_1250/", 0, 100000, 5, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B3_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B3_750/" , 0, 100000, 3, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B3_500/" , 0, 100000, 2, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B4_1500/", 0, 100000, 6, 0.5, file_step)
end

## benchmark setups Dabrowski 2012
main("/users/whalter1/work/aniso_fix/D1_500_XL/" , 0   , 100000, 2, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_500_XL_March/" , 0   , 100000, 2, 0.0, file_step)

main("/users/whalter1/work/aniso_fix/D1_2000/"   , 2790, 100000, 8, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_1000/"   , 2790, 100000, 4, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_1000_L/" , 2790, 100000, 4, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_1000_XL/", 2790, 100000, 4, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_500/"    , 0   , 100000, 2, 0.0, file_step)
main("/users/whalter1/work/aniso_fix/D1_500_L/"  , 0   , 100000, 2, 0.0, file_step)

## all aniso setups

#main("/users/whalter1/work/aniso_fix/C02_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C03_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C04_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C05_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/C12_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C13_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C15_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/C2_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C3_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C4_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/C5_1000/" , 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/D2_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D3_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D4_1000/" , 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D5_1000/" , 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/D02_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D03_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D04_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D05_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/D12_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D13_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/D15_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/A02_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A03_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A04_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A05_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/A12_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A13_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A15_1000/", 2790, 100000, 4, 0.5, file_step)

#########################################################################################

#main("/users/whalter1/work/aniso_fix/A02_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A03_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A04_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A05_500/" , 0   , 100000, 2, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/A12_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A13_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A14_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A15_500/" , 0   , 100000, 2, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/A10_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A00_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A0_500/"  , 0   , 100000, 2, 0.5, file_step)

#main("/users/whalter1/work/aniso_fix/A2_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A3_1000/", 3320, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A4_1000/", 2770, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A5_1000/", 3460, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/MWE_1000/", 400, 100000, 4, 0.5, file_step)

# for Strain Localization
#main("/users/whalter1/work/aniso_fix/A12_500/" , 0   , 100000, 2, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A12_1000/", 2790, 100000, 4, 0.5, file_step)

# for testing new setup...
#main("/users/whalter1/MDOODZ7.0/cmake-exec/AnisoViscTest_evolv_multi_ellipses/", 0, 18000, 2, 0.5, file_step)
#main("/users/whalter1/MDOODZ7.0/cmake-exec/Fig10_bench/", 0, 18000, 4, 0.0, 1)
#

else
# for local testing...
main("/mnt/c/Users/whalter1/Desktop/temp/aniso_fix/A14_500/" , 0   , 100000, 2, 0.5, file_step)
end
