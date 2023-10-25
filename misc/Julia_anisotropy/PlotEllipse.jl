using Plots

δ = 0.5  # Replacing delta with Δ

pyplot()

# Parameters
τ_ΙΙ = 1.0

θ = LinRange(0, 2*π, 400)
τ_χχ_isotropic = τ_ΙΙ .* cos.(θ)
τ_χy_isotropic = τ_ΙΙ .* sin.(θ)

# For the ellipse:
α_values = τ_ΙΙ .* cos.(θ)
β = δ * τ_ΙΙ
τ_χχ_ellipse = α_values
τ_χy_ellipse = β .* sin.(θ)


indices = findall(τ_χy_ellipse.^2 .+ τ_χχ_ellipse.^2 .> τ_ΙΙ^2)
τ_χy_ellipse[indices] .*= τ_ΙΙ ./ sqrt.(τ_χχ_ellipse[indices].^2 .+ τ_χy_ellipse[indices].^2)
τ_χχ_ellipse[indices] .*= τ_ΙΙ ./ sqrt.(τ_χχ_ellipse[indices].^2 .+ τ_χy_ellipse[indices].^2)


p = Plots.plot(size=(800,800), legend=false, aspect_ratio=:equal, grid=false, framestyle=:none, xlims=(-τ_ΙΙ-0.1, τ_ΙΙ+0.1), ylims=(-τ_ΙΙ-0.1, τ_ΙΙ+0.1))
Plots.plot!(τ_χχ_isotropic, τ_χy_isotropic, linestyle=:dash, color=:black, label="Isotropic Case", linewidth=2.0)
Plots.plot!(τ_χχ_ellipse, τ_χy_ellipse, color=:red, label="Anisotropic", linewidth=2.0)


quiver_factor = 1.05
Plots.quiver!([-τ_ΙΙ], [0], quiver=([2*τ_ΙΙ*quiver_factor], [0]), color=:gray, linewidth=2.0)
Plots.quiver!([0], [-τ_ΙΙ], quiver=([0], [2*τ_ΙΙ*quiver_factor]), color=:gray, linewidth=2.0)


font_size = 14
annotate!(0.8*τ_ΙΙ, -0.05, Plots.text("\$\\tau_{χχ}'\$", :center, font_size))
annotate!(-0.15, 0.8*τ_ΙΙ, Plots.text("\$\\tau_{χy}'\$", :center, font_size))

# Display plot
p
