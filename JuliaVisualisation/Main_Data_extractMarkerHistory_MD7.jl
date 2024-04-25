import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, CairoMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)
My = 1e6*365*24*3600

function main()

    # Set the path to your files
    path = "/users/lcandiot/nas/D2c/PhD_LCandioti_2016_2021/project3_UHProckExhumation/dataAndSourceCode/MODEL1_1x1k_Serp_25Phi_1Prefpwl_NewCode_6kmSerp_WesterlyGranite_ReducedConvRate10mmyr_ParticlesCorr/"
    
    # File numbers
    file_start   = 9000
    file_step    = 200
    file_end     = 13000

    # Switches
    MDOODZv4     = true      # Select MDOODZ version
    if MDOODZv4
        MDOODZv7 = false
    else
        MDOODZv7 = true
    end
    resolution   = 500
    T_contours   = true
    printfig     = false
    writeHistory = true
    visualize    = false
    # Scaling
    Lc = 1000.
    tc = My
    Vc = 1e-9

    # Fields
    field = :Phases

    # Coloring
    cmap         = zeros(RGB{Float64}, 7)
    cmap[1]      = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
    cmap[2]      = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
    cmap[3]      = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
    cmap[4]      = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
    cmap[5]      = RGBA{Float64}(217/255, 099/255, 097/255, 1.) 
    cmap[6]      = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
    cmap[7]      = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
    phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)

    # ---------------------------------------------------------- #
    #               Marker extraction routine                    #

    # Load particle data
    firstStep   = 8480
    lastStep    = 13000
    fmark_first = string(path, @sprintf("Particles%05d.gzip.h5", firstStep))
    fmark_last  = string(path, @sprintf("Particles%05d.gzip.h5", lastStep ))
    xm_first    = Float64.(ExtractData( fmark_first, "/Particles/x" ))
    zm_first    = Float64.(ExtractData( fmark_first, "/Particles/z" ))
    xm_last     = Float64.(ExtractData( fmark_last,  "/Particles/x" ))
    zm_last     = Float64.(ExtractData( fmark_last,  "/Particles/z" ))
    model_end   = Float64.(ExtractData( fmark_last,  "/Model/Params"))
    # Identify target markers based on coordinates
    targetCoord  = (50e3, 250e3, -20e3, 10e3)          # Coordinates of exhumed material
    target_ind   = Int64[ ]
    for iMarker in eachindex(xm_last)  # Loop through markers in target region ensuring that they existed at the onset of convergence
        if ( xm_last[iMarker] >= targetCoord[1] && xm_last[iMarker] <= targetCoord[2] && zm_last[iMarker] >= targetCoord[3] && zm_last[iMarker] <= targetCoord[4] && iMarker <= size(xm_first,1))
            push!(target_ind, iMarker)
        end
    end

    # Initialise storage arrays
    if writeHistory
        vxm_time   = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        vzm_time   = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        xm_time    = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        zm_time    = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        Tm_time    = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        Pm_time    = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        phm_time   = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        sxxdm_time = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        sxzm_time  = zeros(Float64, size(target_ind,1), cld(file_end - file_start, file_step) + 1)
        tm_time    = zeros(Float64, 1,                  cld(file_end - file_start, file_step) + 1)
    end

    # Storage file name 
    markerTimeFileName = string(path, "MarkerData_xmin_$(targetCoord[1])_xmax_$(targetCoord[2])_zmin_$(targetCoord[3])_zmax_$(targetCoord[4])_untilTime_$(model_end[1]/tc).h5")
    println(markerTimeFileName)
    # ---------------------------------------------------------- #
    println(Threads.nthreads())
    tic = Base.time()
    # Time loop
    Threads.@threads for iMTime in eachindex(tm_time)

        # Define table index at current time (used for marker storage)
        istep = Int64( (iMTime-1) * file_step + file_start )

        # Read grid file
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model    = ExtractData( filename, "/Model/Params")
        xc       = ExtractData( filename, "/Model/xc_coord")
        zc       = ExtractData( filename, "/Model/zc_coord")
        xv       = ExtractData( filename, "/Model/xg_coord")
        zv       = ExtractData( filename, "/Model/zg_coord")
        xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr")
        xc_hr    = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr    = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

        # Read marker file
        fmark = string(path, @sprintf("Particles%05d.gzip.h5", istep))
        xm    = Float64.(ExtractData( fmark, "/Particles/x"     ));
        zm    = Float64.(ExtractData( fmark, "/Particles/z"     ));
        vxm   = Float64.(ExtractData( fmark, "/Particles/Vx"    ));
        vzm   = Float64.(ExtractData( fmark, "/Particles/Vz"    ));
        Pm    = Float64.(ExtractData( fmark, "/Particles/P"     ));
        Tm    = Float64.(ExtractData( fmark, "/Particles/T"     ));
        phm   = Float64.(ExtractData( fmark, "/Particles/phase" ));
        sxxdm = Float64.(ExtractData( fmark, "/Particles/sxxd"  ));
        sxzm  = Float64.(ExtractData( fmark, "/Particles/sxz"   ));

        # Model properties
        t         = model[1]
        tMy       = round(t/tc, digits=6)
        nvx       = Int(model[4])
        nvz       = Int(model[5])
        ncx, ncz  = nvx-1, nvz-1
        Lx, Lz    = model[2], model[3]
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        # @show "Model apect ratio" Lx/Lz
        # @show "Model time" t/My

        # Load unreasonable amount of data
        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));          mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr)); 
        ηc    = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz));          ηc[mask_air]  .= NaN
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(2*ε̇xx.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN

        # Version dependent variable names 
        if MDOODZv4
            height  = Float64.(ExtractData( filename, "/Topo/height")) 
            Vx_grid = Float64.(ExtractData( filename, "/Topo/vx"))
            Vz_grid = Float64.(ExtractData( filename, "/Topo/vz"))  
            x_mark  = Float64.(ExtractData( filename, "/Topo/x"))         # Marker chain x-coordinate
            z_mark  = Float64.(ExtractData( filename, "/Topo/z"))         # Marker chain z-coordinate
        else
            height  = Float64.(ExtractData( filename, "/Topo/z_grid")) 
            Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"))
            Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"))  
            Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"))
            Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"))
            x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"))
            z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"))
        end

        # Store marker data at current time
        if writeHistory
            xm_time[:,    iMTime] = xm[target_ind    ]
            zm_time[:,    iMTime] = zm[target_ind    ]
            vxm_time[:,   iMTime] = vxm[target_ind   ]
            vzm_time[:,   iMTime] = vzm[target_ind   ]
            Pm_time[:,    iMTime] = Pm[target_ind    ]
            Tm_time[:,    iMTime] = Tm[target_ind    ]
            phm_time[:,   iMTime] = phm[target_ind   ]
            sxxdm_time[:, iMTime] = sxxdm[target_ind ]
            sxzm_time[:,  iMTime] = sxzm[target_ind  ]
            tm_time[1,    iMTime] = t
        end

        # Figure
        if visualize
            f = Figure(resolution = (Lx/Lz*resolution, resolution), fontsize=25)
            ax1 = Axis(f[1, 1], title = L"%$(field) at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            sc1 = scatter!(ax1, xm[target_ind[1:1000:end]]./Lc, zm[target_ind[1:1000:end]]./Lc, color=:red)
            # sc2 = scatter!(ax1, xm_first[target_ind]./Lc, zm_first[target_ind]./Lc, color=:blue)     
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end           
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Colorbar(f[1, 2], hm, label = "$(field)", width = 20, labelsize = 25, ticklabelsize = 14 )
            colgap!(f.layout, 20)
        end
        # Print for cluster visualisation
        if printfig Print2Disk( f, path, string(field), istep) end
    end
    toc = Base.time() - tic
    println("Time loop took $(toc) s to execute")
    # Write marker data to disk
    if writeHistory
        WriteData(markerTimeFileName, "/History/xm_time",    xm_time   )
        WriteData(markerTimeFileName, "/History/zm_time",    zm_time   )
        WriteData(markerTimeFileName, "/History/vxm_time",   vxm_time  )
        WriteData(markerTimeFileName, "/History/vzm_time",   vzm_time  )
        WriteData(markerTimeFileName, "/History/Pm_time",    Pm_time   )
        WriteData(markerTimeFileName, "/History/Tm_time",    Tm_time   )
        WriteData(markerTimeFileName, "/History/phm_time",   phm_time  )
        WriteData(markerTimeFileName, "/History/sxxdm_time", sxxdm_time)
        WriteData(markerTimeFileName, "/History/sxzm_time",  sxzm_time )
        WriteData(markerTimeFileName, "/History/tm_time",    tm_time   )
    end
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function WriteData( file_path, data_path, var)
    h5open(file_path, "w") do file
        write(file, data_path, var)
    end
end

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$(field)/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()