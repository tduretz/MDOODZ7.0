import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie")
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, CSV, Tables, DataFrames

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

function periodicshift(A,Lx,γ)
    Aaff = A
    A0    = zeros(size(A))
    if γ > 0.2
        while Aaff != A0
            A0 = Aaff
            Aaff  = A0 .- (A0 .> Lx/2) .* Lx .+ (A0 .< -Lx/2) .* Lx
        end
    end
    return Aaff
end

function main(path, file_start, file_end_max,fac, xshift, file_step)

    # Select one field to visualise
    field = :StrainLocalisationDeviatoric

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
    if file_start >= file_end; file_start = file_end; end

    @show path
    @show file_start, file_end, file_step

    # Scaling
    Lc = 1.

    # Time loop
    #for istep=[file_end;file_start:file_step:file_end]    
    file_start = 0 # here hardcoded...
    #file_step = 1000 # here hardcoded...
    #file_end = 0 # here hardcoded...
    for istep=file_start:file_step:file_end    
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model  = ExtractData( filename, "/Model/Params")
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")

        t      = model[1]
        Δt     = model[8]
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        γ = t * 2 # here...
        #@show "Time step" istep
        #@show "Model apect ratio" Lx/Lz
        #@show "Model time" t/My

        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), nvx, (ncz+2)))      
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), nvz))      

        #####################################
        # x-shift
        if xshift > 0.0
            # x-coord on center location
            xidx    = Int(floor(ncx*xshift))+1 # x index where to cut
            Vz      = Vz[[xidx:end;1:xidx-1],:] # not accurate ... to solve later...!
            # x-coord on vertex location
            xidx    = Int(floor(nvx*xshift))+1 # x index where to cut
            Vx      = Vx[[xidx:end;1:xidx-1],:]
        end

        #####################################
                
        if field==:StrainLocalisationDeviatoric
            if istep==file_start
                global γ_nextLevel = 0.
                global γ_levelIncr = 0.5 # (hardcoded)
                global npoints = (ncz+2)*1000 # (hardcoded)
                global xP      = zeros(npoints)
                global zP      = Float64.(LinRange(zmin, zmax, npoints))
                global Δt = t
            else
                global Δt = t - t0
            end
            global t0 = t

            # get indices of closest VX and VZ points
            VX_xi = Int.(round.((periodicshift(xP,Lx,γ) .+ Lx/2        ) ./  Lx       .* (nvx-1))) .+ 1
            VX_zi = Int.(round.((zP                     .+ Lz/2 .- Δz/2) ./ (Lz + Δz) .* (ncz+1))) .+ 2
            VZ_xi = Int.(round.((periodicshift(xP,Lx,γ) .+ Lx/2 .- Δx/2) ./ (Lx + Δx) .* (ncx+1))) .+ 2
            VZ_zi = Int.(round.((zP                     .+ Lz/2        ) ./  Lz       .* (nvz-1))) .+ 1

            VX_update = Float64.(Vx[CartesianIndex.(VX_xi, VX_zi)])
            VZ_update = Float64.(Vz[CartesianIndex.(VZ_xi, VZ_zi)])

            global xP      = xP .+ VX_update * Δt
            global zP      = zP .+ VZ_update * Δt

            # save these variables to csv:
            # γ
            # zP
            # xP
            if γ >= γ_nextLevel
                saveVars(path,γ,npoints,zP,xP)
                γ_nextLevel += γ_levelIncr
            end

            # postprocess later (no need to save them):
            # zP0       = zP[1]
            # xP_PS     = zP .* γ
            # xP_DevAbs = xP .- xP_PS
            # xP_DevRel = xP_DevAbs ./ (γ * Lx)

        end
        
    end

end

function saveVars(path,γ,npoints,zP,xP)

    matrixToSave = vcat(γ,npoints,zP,xP) # put all variables in one column
    saveName = (path*"localisation.csv")

    if γ > 0.
        oldTable = CSV.read(saveName, DataFrame, header=false, delim=',')    
        oldTable = Matrix(oldTable)
        #@show oldTable
        matrixToSave = hcat(oldTable, matrixToSave)
        CSV.write(saveName,  Tables.table(matrixToSave), writeheader=false)
    else
        @show(saveName)
        @show(size(matrixToSave))
        CSV.write(saveName,  Tables.table(matrixToSave), writeheader=false)
    end
    @printf("Saved to csv at γ = %.3f\n", γ)

end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function v2c(Av) # interpolation vertex to center
    return 0.25 * ( Av[1:end-1,1:end-1] + Av[2:end,1:end-1] + Av[2:end,1:end-1] + Av[2:end,2:end] )
end

#main()

#run(`ffmpeg -framerate 2 -i /Users/lcandiot/Developer/MDOODZ7.0/RUNS/RiftingChenin/_Phases/Phases%05d.png -c libx264 -pix_fmt yuv420p -y out_movie.mp4`)

# Working command
#ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4

# Set the path to your files
localTestOnly = true
aniSetupsOnly = true
file_step = 10

if localTestOnly == false
if aniSetupsOnly == false
    # all iso setups
    
    main("/users/whalter1/work/aniso_fix/B19_2000/", 0, 100000, 8, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B19_1500/", 0, 100000, 6, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B19_1000/", 0, 100000, 4, 0.5, file_step)
    main("/users/whalter1/work/aniso_fix/B6_2000/", 0, 100000, 8, 0.5, file_step)
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
    main("/users/whalter1/work/aniso_fix/B5_1000/", 0, 100000, 4, 0.5, file_step)
end

## all aniso setups
#main("/users/whalter1/work/aniso_fix/A02_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A03_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A04_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A05_1000/", 2790, 100000, 4, 0.5, file_step)
#
#main("/users/whalter1/work/aniso_fix/A12_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A13_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A15_1000/", 2790, 100000, 4, 0.5, file_step)
#
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
main("/users/whalter1/work/aniso_fix/A14_1000/", 2790, 100000, 4, 0.5, file_step)
#main("/users/whalter1/work/aniso_fix/A12_1000/", 2790, 100000, 4, 0.5, file_step)

# for testing new setup...
#main("/users/whalter1/MDOODZ7.0/cmake-exec/AnisoViscTest_evolv_multi_ellipses/", 0, 18000, 4, 0.5, file_step)
#main("/users/whalter1/MDOODZ7.0/cmake-exec/Fig10_bench/", 0, 18000, 4, 0.0, 1)
#

else
# for local testing...
main("/mnt/c/Users/whalter1/Desktop/temp/aniso_fix/A14_500/" , 0   , 100000, 2, 0.5, file_step)
end