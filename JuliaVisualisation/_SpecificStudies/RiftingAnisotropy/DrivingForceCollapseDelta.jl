using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics, UnPack, JLD2
using CairoMakie#, GLMakie
Mak = CairoMakie

using Polynomials


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

    ftsz = 18
    f = Figure( fontsize=ftsz)

    δ = [1 1.5 2 3 4 5 6 7 8 10][:]
    f1 = maximum(d1["F_driving"]./1e12)
    f1_5 = maximum(d1_5["F_driving"]./1e12)
    f2 = maximum(d2["F_driving"]./1e12)
    f3 = maximum(d3["F_driving"]./1e12)
    f4 = maximum(d4["F_driving"]./1e12)
    f5 = maximum(d5["F_driving"]./1e12)
    f6 = maximum(d6["F_driving"]./1e12)
    f7 = maximum(d7["F_driving"]./1e12)
    f8 = maximum(d8["F_driving"]./1e12)
    f10 = maximum(d10["F_driving"]./1e12)

    maxF = [f1 f1_5 f2 f3 f4 f5 f6 f7 f8 f10][:]

    logF = log10.(maxF)

    n   = length(logF)
    ∑x  = sum(δ)
    ∑y  = sum(logF)
    ∑xy = sum(δ.*logF)
    ∑x2 = sum(δ.*δ)

    a = (n*∑xy - ∑y*∑x)/(n*∑x2 - ∑x^2)
    b = (∑y - a*∑x) / n

    @show δ
    @show maxF

    degree = 2
    p = fit(δ, maxF, degree)

    @show coeffs(p)

    ax1 = Axis(f[1, 1:2], title = L"$$A) Maximum driving force Vs. anisotropy strength", xlabel = L"$\delta$", ylabel = L"max. $F_\textrm{D}$ [TN/m]")
    # hidexdecorations!(ax1)
    # lines!(ax1, d1["time"], d1["F_driving"]./1e12, label=L"$\delta = 1$")  
    # lines!(ax1, d1_5["time"], d1_5["F_driving"]./1e12, label=L"$\delta = 1.5$") 
    # lines!(ax1, d2["time"], d2["F_driving"]./1e12, label=L"$\delta = 2$")  
    # lines!(ax1, d3["time"], d3["F_driving"]./1e12, label=L"$\delta = 3$")  
    # lines!(ax1, d4["time"], d4["F_driving"]./1e12, label=L"$\delta = 4$")  
    # lines!(ax1, d5["time"], d5["F_driving"]./1e12, label=L"$\delta = 5$")  
    # lines!(ax1, d6["time"], d6["F_driving"]./1e12, label=L"$\delta = 6$")  
    # lines!(ax1, d7["time"], d7["F_driving"]./1e12, label=L"$\delta = 7$")  
    # lines!(ax1, d8["time"], d8["F_driving"]./1e12, label=L"$\delta = 8$")  
    # lines!(ax1, d10["time"], d10["F_driving"]./1e12, label=L"$\delta = 10$")  
    # scatter!(ax1, δ, logF)
    # lines!(ax1, δ, a*δ .+ b)
    # lines!(ax1, δ, 10.0.^(a*δ .+ b))
    fit_force = coeffs(p)[2]*δ .+ coeffs(p)[3]*δ.^2 .+ coeffs(p)[1]
    lines!(ax1, δ, fit_force, linestyle=:dash, color=:black)
    scatter!(ax1, δ, maxF, marker=:xcross, markersize=18)
    text!(ax1, 5, 17, text=L"$F_\text{pred} = 24.05 - 2.85\delta + 0.12\delta^2$", rotation=0.)



    # ylims!(ax1, 5, 25)
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14, nbanks=2)


    ax2 = Axis(f[2, 1:2], title = L"$$B) Scaled driving force through time", xlabel = L"$t$ [Ma]", ylabel = L"$F_\textrm{D} / F_\textrm{pred}$")
    lines!(ax2, d1["time"], d1["F_driving"]./1e12 ./  fit_force[1], label=L"$\delta = 1$")  
    lines!(ax2, d1_5["time"], d1_5["F_driving"]./1e12 ./  fit_force[2], label=L"$\delta = 1.5$") 
    lines!(ax2, d2["time"], d2["F_driving"]./1e12 ./  fit_force[3], label=L"$\delta = 2$")  
    lines!(ax2, d3["time"], d3["F_driving"]./1e12 ./  fit_force[4], label=L"$\delta = 3$")  
    lines!(ax2, d4["time"], d4["F_driving"]./1e12 ./  fit_force[5], label=L"$\delta = 4$")  
    lines!(ax2, d5["time"], d5["F_driving"]./1e12 ./  fit_force[6], label=L"$\delta = 5$")  
    lines!(ax2, d6["time"], d6["F_driving"]./1e12 ./  fit_force[7] , label=L"$\delta = 6$")  
    lines!(ax2, d7["time"], d7["F_driving"]./1e12 ./  fit_force[8], label=L"$\delta = 7$")  
    lines!(ax2, d8["time"], d8["F_driving"]./1e12 ./  fit_force[9], label=L"$\delta = 8$")  
    lines!(ax2, d10["time"], d10["F_driving"]./1e12./ fit_force[10] , label=L"$\delta = 10$")  


    # ax2 = Axis(f[2, 1:2], title = L"$$Driving force versus initial anisotropy angle", xlabel = L"$t$ [Ma]", ylabel = L"$F_\textrm{D}$ [TN/m]")
    # lines!(ax2, t5["time"],   t5["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 5$")  
    # lines!(ax2, t80["time"], t80["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 80$")  
    # lines!(ax2, t20["time"], t20["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 20$")  
    # lines!(ax2, t30["time"], t30["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 30$")  
    # lines!(ax2, t45["time"], t45["F_driving"]./1e12, label=L"$\theta_\textrm{ini} = 45$")  
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14)

    display(f) 
    
    save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/DrivingForceCollapseDelta.png", f, px_per_unit = 4)     

end

main()