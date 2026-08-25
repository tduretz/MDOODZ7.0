using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, UnPack, JLD2
using CairoMakie#, GLMakie
Mak = CairoMakie

function main()

    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1/"
    d1 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d1_5/"
    d1_5 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d2/"
    d2 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d3/"
    d3 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/ref_d4_MR/"
    d4 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d5/"
    d5 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d6/"
    d6 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d7/"
    d7 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d8/"
    d8 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/d10/"
    d10 = load(path*"DrivingForce.jld2")

    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t0/"
    t0 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t5/"
    t5 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t20/"
    t20 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t30/"
    t30 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t45/"
    t45 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t80/"
    t80 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t90/"
    t90 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t60/"
    t60 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t70/"
    t70 = load(path*"DrivingForce.jld2")


    ftsz = 18
    f = Figure( fontsize=ftsz)

    cmap = :roma10

    ax1 = Axis(f[1, 1:2], title = L"$$A) Driving force through time: Anisotropy strength", xlabel = L"$t$ [Ma]", ylabel = L"$F_\textrm{D}$ [TN/m]")
    # hidexdecorations!(ax1)
    lines!(ax1, d1["time"], d1["F_driving"]./1e12,     label=L"$\delta = 1$"  , color= colorschemes[cmap].colors[1])  
    lines!(ax1, d1_5["time"], d1_5["F_driving"]./1e12, label=L"$\delta = 1.5$", color= colorschemes[cmap].colors[2]) 
    lines!(ax1, d2["time"], d2["F_driving"]./1e12,     label=L"$\delta = 2$"  , color= colorschemes[cmap].colors[3])  
    lines!(ax1, d3["time"], d3["F_driving"]./1e12,     label=L"$\delta = 3$"  , color= colorschemes[cmap].colors[4])  
    lines!(ax1, d4["time"], d4["F_driving"]./1e12,     label=L"$\delta = 4$"  , color= colorschemes[cmap].colors[5])  
    lines!(ax1, d5["time"], d5["F_driving"]./1e12,     label=L"$\delta = 5$"  , color= colorschemes[cmap].colors[6])  
    lines!(ax1, d6["time"], d6["F_driving"]./1e12,     label=L"$\delta = 6$"  , color= colorschemes[cmap].colors[7])  
    lines!(ax1, d7["time"], d7["F_driving"]./1e12,     label=L"$\delta = 7$"  , color= colorschemes[cmap].colors[8])  
    lines!(ax1, d8["time"], d8["F_driving"]./1e12,     label=L"$\delta = 8$"  , color= colorschemes[cmap].colors[9])  
    lines!(ax1, d10["time"], d10["F_driving"]./1e12,   label=L"$\delta = 10$" , color= colorschemes[cmap].colors[10])  
    xlims!(ax1, 0, 10)
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14, nbanks=2)
    Colorbar(f, limits = (1, 10), colormap = :roma10,
    label = L"δ", vertical = false, flipaxis = true, bbox=BBox(450, 550, -650, 1350), ticklabelsize=12)


    cmap = :vik10

    ax2 = Axis(f[2, 1:2], title = L"$$B) Driving force through time: Initial fabric angle", xlabel = L"$t$ [Ma]", ylabel = L"$F_\textrm{D}$ [TN/m]")
    lines!(ax2, t0["time"],   t0["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 0$" , color= colorschemes[cmap].colors[1])  
    lines!(ax2, t5["time"],   t5["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 5$" , color= colorschemes[cmap].colors[3])  
    lines!(ax2, d4["time"],   d4["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 10$", color= colorschemes[cmap].colors[4])  
    lines!(ax2, t20["time"], t20["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 20$", color= colorschemes[cmap].colors[6])  
    lines!(ax2, t30["time"], t30["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 30$", color= colorschemes[cmap].colors[7])  
    lines!(ax2, t45["time"], t45["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 45$", color= colorschemes[cmap].colors[8])  
    lines!(ax2, t60["time"], t60["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 45$", color= colorschemes[cmap].colors[8])  
    lines!(ax2, t70["time"], t70["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 45$", color= colorschemes[cmap].colors[8])  
    lines!(ax2, t80["time"], t80["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 80$", color= colorschemes[cmap].colors[5])  
    lines!(ax2, t90["time"], t90["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 90$", color= colorschemes[cmap].colors[2])  
    xlims!(ax1, 0, 10)

    Colorbar(f, limits = (0, 90), colormap = :vik10,
    label = L"θ_\text{ini}", vertical = false, flipaxis = true, bbox=BBox(450, 550, -650, 900), ticklabelsize=12)
    
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14, nbanks=2)

    display(f)  

    save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/DrivingForceData.png", f, px_per_unit = 4)     

end

main()