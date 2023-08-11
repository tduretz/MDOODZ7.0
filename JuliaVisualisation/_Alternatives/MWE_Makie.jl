import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using HDF5, GLMakie, Printf

function main()
    # Mesh
    x = LinRange(-1.0, 3, 200)
    y = LinRange(-2.0, 2.0, 100)
    Lx, Ly = x[end]-x[1], y[end]-y[1]
    # One continuous field 
    f = 1e3.*exp.(-x.^2  .- (y').^2)
    # One discontinuous field 
    g = zero(f)
    g[(x.^2 .+ (y.^2)') .< 0.5^2] .= 1.0
    # Make transparent colormap
    cbarPal= :plasma
    cmap = get(colorschemes[cbarPal], LinRange(0,1,100))
    cmap2 = [(cmap[i],i/100) for i in 1:100]
    # Compose figure
    fig = Figure(resolution = (600, 600), fontsize = 22, fonts = (;regular="CMU Serif"))
    ax = fig[1, 1] = Axis(fig, xlabel = L"x", ylabel = L"y", aspect=Lx/Ly)
    
    hm = heatmap!(ax, x, y, f, colormap = cmap2)
    contour!(ax, x, y, g)
    Colorbar(fig, hm, label = L"\log_{10}[(u^2+v^2)^{1/2}]", width = 20,
        labelsize = 14, ticklabelsize = 14, bbox=ax.scene.px_area,
        alignmode = Outside(10), halign = :right, ticklabelcolor = :white, labelcolor = :white,
        tickcolor = :white)
    GLMakie.colgap!(fig.layout, 20)
    DataInspector(fig)
    display(fig)

end

main()