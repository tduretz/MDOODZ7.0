using Plots
plotlyjs()

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
    p = plot()
    # p = heatmap!(x, y, f')
    p = contour!(x, y, g', c=:black, linewidth=0.1) # contour does not appear on top of heatmap, colorbar is messed up, linewidth fails
    p = plot!(aspect_ratio=1) # aspect ratio fails: does not seem 1.0 on my screen
    display(p)
end

main()