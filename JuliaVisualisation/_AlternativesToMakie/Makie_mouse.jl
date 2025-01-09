import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using AbstractPlotting, GLMakie

img = rand(100, 100)
scene = Scene(resolution = (500, 500))
heatmap!(scene, img, scale_plot = false)
clicks = Node(Point2f0[(0,0)])
on(scene.events.mousebuttons) do buttons
   if ispressed(scene, Mouse.left)
       pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
       clicks[] = push!(clicks[], pos)
   end
   return
end
scatter!(scene, clicks, color = :red, marker = '+', markersize = 10)

# Do not execute beyond this point!

RecordEvents(scene, "output")