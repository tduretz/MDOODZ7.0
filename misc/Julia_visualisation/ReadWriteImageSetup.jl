import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
gr()
plotlyjs()

using Plots, FileIO, ImageCore, Printf

function ReadBinFileVisualise()
    path  = "/Users/tduretz/REPO/MDOODZ7.0/IMPORT/"
    fname = string(path, "SetupCoesite2_nxy2305.bin")
    x     = read(fname)
    nx    = Int16(sqrt(length(x)))
    x2d   = reshape(x, nx, nx)
    heatmap(x2d')
end

function ReadPNGSaveBin()

    path_in  = "/Users/tduretz/REPO/MDOODZ7.0/RUNS/"
    fname_in = string(path_in, "setup_coesite2_140.png")
    img      = channelview(FileIO.load(fname_in))
    inds     = unique(img[1,:,:])
    phases   = zeros(UInt8, size(img[1,:,:]))
    
    display(sum(img[1,:,:].==inds[1]))
    display(sum(img[1,:,:].==inds[2]))
    display(sum(img[1,:,:].==inds[3]))

    phases[img[1,:,:].==inds[2]] .= 0
    phases[img[1,:,:].==inds[1]] .= 1
    phases[img[1,:,:].==inds[3]] .= 2

    nx = size(img, 2) 

    # heatmap(phases')
    path_out  = "/Users/tduretz/REPO/MDOODZ7.0/IMPORT/"
    fname_out = string(path_out, string("SetupCoesite2_nxy", @sprintf("%04d", nx), ".bin"))
    write(fname_out, phases)
end

# @time ReadBinFileVisualise()
@time ReadPNGSaveBin()

