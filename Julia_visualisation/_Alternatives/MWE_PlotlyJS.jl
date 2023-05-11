import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using PlotlyJS

function main_PlotlyJS()
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
    p = Plot([
    heatmap(x=x, y=y, z=f', colorscale=colors.plasma, colorbar_thickness=24), 
    contour(x=x, y=y, z=g', mode="lines", line_width=2, colorscale=[[1, "rgb(255,255,255)"]], contours_coloring="none"),
    heatmap(x=x, y=y, z=f', colorscale=colors.plasma, colorbar_thickness=24)],
    Layout(width=600, height=600*Ly/Lx,
    xaxis = attr(;
    range = [-2,2],
    constrain = "domain",
    ),
    yaxis = attr(;
        range = [-1,1],
        scaleanchor = "x",
        scaleratio = 1,
    )))
    display(p)
end

main_PlotlyJS()