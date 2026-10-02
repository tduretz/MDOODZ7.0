import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
using CairoMakie, GLMakie
using CSV, Tables

Mak = GLMakie
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

const y = 365 * 24 * 3600
const My = 1e6 * y
const cm_y = y * 100.

@views function main()
    # Define path to the output file
    path = "/home/lilou/librairies/MDOODZ7.0/save_runs/SS_inclusion_n1_Pidam2_evolaniso_20261001_181918/"

    # Step number
    istep = 10

    # Load data from the file
    name = @sprintf("Output%05d.gzip.h5", istep)
    filename = string(path, name)

    # Extract model parameters and coordinates
    model = ExtractData(filename, "/Model/Params")
    xc = ExtractData(filename, "/Model/xc_coord")
    zc = ExtractData(filename, "/Model/zc_coord")
    xv = ExtractData(filename, "/Model/xg_coord")
    zv = ExtractData(filename, "/Model/zg_coord")
    nvx = Int(model[4])
    nvz = Int(model[5])
    ncx, ncz = nvx - 1, nvz - 1
    mask_air = ExtractField(filename, "/VizGrid/compo", (ncx, ncz), false, 0) .== -1.0

    # Extract fields
    τII = ExtractField(filename, "/Centers/sxxd", (ncx, ncz), true, mask_air)
    η_lin = ExtractField(filename, "/Centers/eta_n", (ncx, ncz), true, mask_air)
    δ_ani = ExtractField(filename, "/Centers/ani_fac", (ncx, ncz), false, 0)
    Vdam = ExtractField(filename, "/Vertices/Vdam_s", (nvx, nvz), false, 0)
    Nx = Float64.(reshape(ExtractData(filename, "/Centers/nx"), ncx, ncz))
    Nz = Float64.(reshape(ExtractData(filename, "/Centers/nz"), ncx, ncz))

    # Compute anisotropy angle (in degrees)
    angle_rad = atan.(Nz, Nx)  # Angle in radians
    angle_deg = rad2deg.(angle_rad)
    angle_deg = ifelse.(angle_deg .< 0, angle_deg .+ 180, angle_deg)  # Convert to [0, 180]
    angle_deg = ifelse.(angle_deg .> 90, 180 .- angle_deg, angle_deg)  # Convert to [0, 90]

    # Create a figure with 2 rows and 2 columns
    f = Figure(size = (1500, 1000), fontsize=30)


    # Plot Vdam
    ax2 = Axis(f[1, 1],
               title = L"$V_{dam}$",
               xlabel = L"$x$ [km]",
               ylabel = L"$y$ [km]",
               aspect = 1,
               titlesize = 40)
    hm2 = heatmap!(ax2, xv, zv, Vdam, colormap = (Reverse(:bilbao), 1.0))
    Colorbar(f[1, 2], hm2, label = L"$V_{dam}$", width = 10, height = Relative(0.5))

      # Plot linear viscosity
    ax1 = Axis(f[1, 3],
               title = L"$\eta$",
               xlabel = L"$x$ [km]",
               ylabel = L"$y$ [km]",
               aspect = 1,
               titlesize = 40)
    hm1 = heatmap!(ax1, xc, zc, η_lin, colormap = (:lapaz, 1.0), colorrange=(0.3,0.8))
    Colorbar(f[1, 4], hm1, label = L"$\eta$ [Pa.s]", width = 10, height = Relative(0.5), )

      

    # Plot anisotropy ratio
    ax3 = Axis(f[2, 1],
               title = L"$\delta_{ani}$",
               xlabel = L"$x$ [km]",
               ylabel = L"$y$ [km]",
               aspect = 1,
               titlesize = 40)
    hm3 = heatmap!(ax3, xc, zc, δ_ani, colormap = (:amp, 1.0))
    Colorbar(f[2, 2], hm3, label = L"$\delta_{ani}$", width = 10, height = Relative(0.5))

    # Plot anisotropy angle
    ax4 = Axis(f[2, 3],
               title = L"$θ_{ani}$ [°]",
               xlabel = L"$x$ [km]",
               ylabel = L"$y$ [km]",
               aspect = 1,
               titlesize = 40)
    hm4 = heatmap!(ax4, xc, zc, angle_deg, colormap = (:jet, 1.0))
    Colorbar(f[2, 4], hm4, label = L"$θ_{ani}$ [°]", width = 10, height = Relative(0.5))





    # Save the figure
    save("/home/lilou/Desktop/presentations-rapports/pres_sete_2026/figures/aniso_Evol_4fields.png", f)

    display(f)
end

main()