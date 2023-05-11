import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

My = 1e6*365*24*3600

function LoadGrainSize(path, istep)
    filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
    model  = ExtractData( filename, "/Model/Params")
    xc     = ExtractData( filename, "/Model/xc_coord")
    zc     = ExtractData( filename, "/Model/zc_coord")
    t      = model[1]
    tMy    = round(t/My, digits=2)
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1
    d      = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));
    return (xc=xc, zc=zc, d=d, tMy=tMy)
end

function main()

    path_ref = "/Users/tduretz/REPO/MDOODZ7.0/VISUAL_TESTS/GSE_REF/"
    path_new = "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/PinchSwellGSE_NEW/"

    step_ref = 100
    step_new = 110
    printfig = true

    ref = LoadGrainSize(path_ref, step_ref)
    new = LoadGrainSize(path_new, step_new)
        
    xminz, xmaxz = -0.4, 0.4
    zminz, zmaxz = -0.17, 0.17
    Lx, Lz       = xmaxz - xminz, zmaxz - zminz            
    
    # #####################################

    f = Figure(resolution = ( Lx/Lz*1000, 1000), fontsize=25, aspect = 2.0)
    
    ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(new.tMy) Ma", xlabel = "x [km]", ylabel = "y [km]")
    hm = heatmap!(ax1, new.xc, new.zc, log10.(new.d.*1e6), colormap = (:turbo, 0.85), colorrange=(1, 3))
    xlims!(ax1, xminz, xmaxz)
    ylims!(ax1, zminz, zmaxz)
    colsize!(f.layout, 1, Aspect(1, Lx/Lz))
    GLMakie.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
    GLMakie.colgap!(f.layout, 20)
    if printfig save(path*"GrainSize"*@sprintf("%05d", istep)*".png", f, px_per_unit = 4) end

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