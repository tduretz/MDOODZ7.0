N = 1000
r = [(rand(Point2f0, 5000) .- 0.5) .* 25 for i = 1:N]
scene = scatter(r[1], markersize = 1, limits = FRect(-25/2, -25/2, 25, 25))
s = scene[end] # last plot in scene
record(scene, "output.mp4", r) do m
    s[1] = m
end