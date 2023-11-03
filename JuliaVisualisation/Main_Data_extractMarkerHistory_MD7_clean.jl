import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, CairoMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)
My = 1e6*365*24*3600

function main()

    # Set the path to your files
    simDir = "/users/lcandiot/nas/D2c/PhD_LCandioti_2016_2021/project3_UHProckExhumation/dataAndSourceCode/"
    simList = readdir(simDir)
    println(simList)
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
    writeHistory = true

    # Scaling
    Lc = 1000.
    tc = My
    Vc = 1e-9

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
    targetCoord = (-500e3, 500e3, -20e3, 60e3)          # Coordinates of exhumed material
    target_ind  = Int64[ ]
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
    # ---------------------------------------------------------- #

    # Time loop
    println(Threads.nthreads())
    tic = Base.time()
    Threads.@threads for iMTime in eachindex(tm_time)

        # Define table index at current time (used for marker storage)
        istep = Int64( (iMTime-1) * file_step + file_start )

        # Read grid file
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model    = ExtractData( filename, "/Model/Params")

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
        model =          ExtractData( fmark, "/Model/Params"     );

        # Model properties
        t         = model[1]

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
    end
    toc = Base.time() - tic
    println("Time loop took $(toc) s to execute")
    # Write marker data to disk
    if writeHistory
        h5open(markerTimeFileName, "w") do file
            write(file, "/History/xm_time", xm_time)
            write(file, "/History/zm_time", zm_time)
            write(file, "/History/vxm_time", vxm_time)
            write(file, "/History/vzm_time", vzm_time)
            write(file, "/History/Pm_time", Pm_time)
            write(file, "/History/Tm_time", Tm_time)
            write(file, "/History/phm_time", phm_time)
            write(file, "/History/sxxdm_time", sxxdm_time)
            write(file, "/History/sxzm_time", sxzm_time)
            write(file, "/History/tm_time", tm_time)
        end
    end
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end


function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$(field)/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()