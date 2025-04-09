import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
using DelimitedFiles

function main()
   
    filex   = "/Users/tduretz/Downloads/CrystalGenerator_test_version/out_x_setup_r201x201_fsp02.bin"
    filez   = "/Users/tduretz/Downloads/CrystalGenerator_test_version/out_z_setup_r201x201_fsp02.bin"
    filep   = "/Users/tduretz/Downloads/CrystalGenerator_test_version/out_p_setup_r201x201_fsp02.bin"
  
    @show nmark   = Int(filesize(filex)/sizeof(Float64))
    x       = zeros(nmark)
    z       = zeros(nmark)
    p       = zeros(Int32, nmark)
    read!(filex, x) # read data
    read!(filez, z) # read data
    read!(filep, p) # read data


    # io = open(filep  ) # read data
    
    # k=0
    # while !eof(io)
    #     k+=0
    # end
    # @show k
    

     # Figure
     Lx, Lz = 1., 1.
     resolution = 1000
     f = Figure(resolution = (Lx/Lz*resolution, resolution), fontsize=25)
     ax1 = Axis(f[1, 1], title = L"Phases", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
     scatter!(ax1, x[p.==1], z[p.==1])
     scatter!(ax1, x[p.==2], z[p.==2])
     scatter!(ax1, x[p.==3], z[p.==3])
     scatter!(ax1, x[p.==4], z[p.==4])
     scatter!(ax1, x[p.==5], z[p.==5])
     display(f)     

    # # @show nmark
    # # k = 1
    # # for imark=1:nmark
    # #     for i=1:3
    # #         x[imark] = markers[k]; k +=1
    # #         y[imark] = markers[k]; k +=1
    # #         p[imark] = markers[k]; k +=1
    # #     end
    # # end 

    # @show minimum(x)
    # @show maximum(x)
    # @show minimum(z)
    # @show maximum(z)
    # @show p
end

main()