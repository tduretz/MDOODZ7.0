import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

ky = 1e3*365*24*3600

function LoadX(path, istep)
    filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
    model  = ExtractData( filename, "/Model/Params")
    xc     = ExtractData( filename, "/Model/xc_coord")
    zc     = ExtractData( filename, "/Model/zc_coord")
    t      = model[1]
    tky    = round(t/ky, digits=4)
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1
    X      = Float64.(reshape(ExtractData( filename, "/Centers/X"), ncx, ncz));

    return (xc=xc, zc=zc, X=X, tky=tky)
end

function main()

    path_ref = "/Users/tduretz/REPO/MDOODZ7.0/VISUAL_TESTS/Shrinking_REF/"
    path_new = "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/Shrinking_NEW/"

    step_ref = 10
    step_new = 10
    printfig = true

    ref = LoadX(path_ref, step_ref)
    new = LoadX(path_new, step_new)
        
    xminz, xmaxz = -1., 1.
    zminz, zmaxz = -1., 1.
    Lx, Lz       = xmaxz - xminz, zmaxz - zminz            

    # # #####################################

    f = Figure(resolution = (1600, 800), fontsize=25)
    
    # Reference
    ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(ref.tky) ka", xlabel = "x [km]", ylabel = "y [km]", aspect=Lx/Lz)
    hm = heatmap!(ax1, ref.xc, ref.zc, ref.X, colormap = (:turbo, 0.85), colorrange=(0, 0.3))
    xlims!(ax1, xminz, xmaxz)
    ylims!(ax1, zminz, zmaxz)
    text!(-.8, .8, text="Reference", color=:white, fontsize=29)
    GLMakie.Colorbar(f[1, 2], hm, label = L"$X$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )

    # New
    ax1 = Axis(f[1, 3], title = L"$d$ at $t$ = %$(new.tky) ka", xlabel = "x [km]", ylabel = "y [km]", aspect=Lx/Lz)
    hm = heatmap!(ax1, new.xc, new.zc, new.X, colormap = (:turbo, 0.85), colorrange=(0, 0.3))
    xlims!(ax1, xminz, xmaxz)
    ylims!(ax1, zminz, zmaxz)
    text!(-.8, .8, text="New", color=:white, fontsize=29)
    GLMakie.Colorbar(f[1, 4], hm, label = L"$X$ [-]", width = 20, labelsize = 25, ticklabelsize = 14 )

    if printfig save("_VisualTests/Shrinking.png", f, px_per_unit = 4) end

    DataInspector(f)
    display(f)
        
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

main()