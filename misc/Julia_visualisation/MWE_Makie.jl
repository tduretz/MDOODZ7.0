import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf

function main()
    # Mesh
    x = LinRange(-2.0, 2.0, 200)
    y = LinRange(-1.0, 1.0, 100)
    Lx, Ly = x[end]-x[1], y[end]-y[1]
    # One continuous field 
    f = 1e3.*exp.(-x.^2  .- (y').^2)
    # One discontinuous field 
    g = zero(f)
    g[(x.^2 .+ (y.^2)') .< 0.5^2] .= 1.0
    # Compose figure
    fig = Figure(resolution = (600, 400), fontsize = 22, fonts = (;regular="CMU Serif"))
    ax = fig[1, 1] = Axis(fig, xlabel = L"x", ylabel = L"y", aspect=Lx/Ly)
    fs = GLMakie.heatmap!(ax, x, y, f, colormap = Reverse(:plasma))
    GLMakie.contour!(ax, x, y, g)
    GLMakie.Colorbar(fig[1, 2], fs, label = L"\log_{10}[(u^2+v^2)^{1/2}]", width = 20,
        labelsize = 14, ticklabelsize = 14)
        GLMakie.colgap!(fig.layout, 20)
    DataInspector(fig)
    display(fig)

end

main()