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
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t60/"
    t60 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t70/"
    t70 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t80/"
    t80 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t0/"
    t0 = load(path*"DrivingForce.jld2")
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiftingAnisotropy/t90/"
    t90 = load(path*"DrivingForce.jld2")

    ftsz = 18
    f = Figure( fontsize=ftsz, size=(450, 600))

    θ = [0 5. 10 20 30 45 60 70 80 90][:]
    ft0 = maximum(t0["F_driving"]./1e12)

    ft5 = maximum(t5["F_driving"]./1e12)
    ft10 = maximum(d4["F_driving"]./1e12)
    ft20 = maximum(t20["F_driving"]./1e12)
    ft30 = maximum(t30["F_driving"]./1e12)
    ft45 = maximum(t45["F_driving"]./1e12)
    ft70 = maximum(t60["F_driving"]./1e12)
    ft60 = maximum(t70["F_driving"]./1e12)
    ft80 = maximum(t80["F_driving"]./1e12)
    ft90 = maximum(t90["F_driving"]./1e12)

   
    maxF = [ft0 ft5 ft10 ft20 ft30 ft45 ft60 ft70 ft80 ft90][:]


    degree = 2
    p = fit(θ, maxF, degree)

    ax1 = Axis(f[1, 1:2], title = L"$$A) Maximum driving force Vs. initial fabric angle", xlabel = L"$\theta_\textrm{ini}$", ylabel = L"max. $F_\textrm{D}$ [TN/m]")
    θ_fit = LinRange(0, 90, 100) 
    fit_force = coeffs(p)[2]*θ_fit .+ coeffs(p)[3]*θ_fit.^2 .+ coeffs(p)[1]

    # fit_force = 12 .+ 6*cosd.((θ_fit.+10)*3.5)

    @show coeffs(p)

    fit_force_data = coeffs(p)[2]*θ .+ coeffs(p)[3]*θ.^2 .+ coeffs(p)[1]
    lines!(ax1, θ_fit, fit_force, linestyle=:dash, color=:black)
    scatter!(ax1, θ, maxF, marker=:xcross, markersize=18)
    text!(ax1, 13, 15, text=L"$F_\text{pred} = 20.12 - 0.60\theta_\text{ini} + 0.006\theta_\text{ini}^2$", rotation=0.)


    ##########################################################################

    # ylims!(ax1, 5, 25)
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14, nbanks=2)

    # Find index of maximum driving force
    ind0 = argmax(t0["F_driving"]./1e12 ./ fit_force_data[1]) 
    @show t0["time"][ind0]
    ind5 = argmax(t5["F_driving"]./1e12 ./ fit_force_data[2]) 
    @show t5["time"][ind5]
    ind10 = argmax(d4["F_driving"]./1e12 ./ fit_force_data[3]) 
    @show d4["time"][ind10]
    ind20 = argmax(t20["F_driving"]./1e12 ./ fit_force_data[4])
    @show t20["time"][ind20]
    ind30 = argmax(t30["F_driving"]./1e12 ./ fit_force_data[5])
    @show t30["time"][ind30]
    ind45 = argmax(t45["F_driving"]./1e12 ./ fit_force_data[6])
    @show t45["time"][ind45]
    ind60 = argmax(t60["F_driving"]./1e12 ./ fit_force_data[7])
    @show t60["time"][ind60]
    ind70 = argmax(t70["F_driving"]./1e12 ./ fit_force_data[8])
    @show t70["time"][ind70]
    ind80 = argmax(t80["F_driving"]./1e12 ./ fit_force_data[9])
    @show t80["time"][ind80]
    ind90 = argmax(t90["F_driving"]./1e12 ./ fit_force_data[10])
    @show t90["time"][ind90]


    # Find index of minimum time derivative of driving force
    ind = argmin(diff(t0["F_driving"]./1e12 ./ fit_force_data[1]))
    t_t0  = 0.5.*(t0["time"][1:end-1] .+ t0["time"][2:end-0])[ind] - t0["time"][ind0]*0
    ind = argmin(diff(t5["F_driving"]./1e12 ./ fit_force_data[1]))
    t_t5  = 0.5.*(t5["time"][1:end-1] .+ t5["time"][2:end-0])[ind] - t5["time"][ind5]*0
    ind = argmin(diff(d4["F_driving"]./1e12 ./ fit_force_data[2]))
    t_t10  = 0.5.*(d4["time"][1:end-1] .+ d4["time"][2:end-0])[ind] - d4["time"][ind10]*0
    ind = argmin(diff(t20["F_driving"]./1e12 ./ fit_force_data[6]))
    t_t80 = 0.5.*(t80["time"][1:end-1] .+ t80["time"][2:end-0])[ind] - t80["time"][ind80]*0
    ind = argmin(diff(t20["F_driving"]./1e12 ./ fit_force_data[3]))
    t_t20 = 0.5.*(t20["time"][1:end-1] .+ t20["time"][2:end-0])[ind] - t20["time"][ind20]*0
    ind = argmin(diff(t30["F_driving"]./1e12 ./ fit_force_data[4]))
    t_t30 = 0.5.*(t30["time"][1:end-1] .+ t30["time"][2:end-0])[ind] - t30["time"][ind30]*0
    ind = argmin(diff(t45["F_driving"]./1e12 ./ fit_force_data[5]))
    t_t45 = 0.5.*(t45["time"][1:end-1] .+ t45["time"][2:end-0])[ind] - t45["time"][ind45]*0
    ind = argmin(diff(t90["F_driving"]./1e12 ./ fit_force_data[1]))
    t_t90  = 0.5.*(t90["time"][1:end-1] .+ t90["time"][2:end-0])[ind] - t90["time"][ind90]*0
    
    # Find when force has decreased by 20%
    @show t_t0  =   t0["time"][ findfirst( t0["F_driving"][:]./1e12 .< 0.8*maxF[1] .&& (t0["time"].> t0["time"][ind0]) ) ]
    @show t_t5  =   t5["time"][ findfirst( t5["F_driving"][:]./1e12 .< 0.8*maxF[2] .&& (t5["time"].> t5["time"][ind5]) ) ]
    @show t_t10 =  d4["time"][ findfirst( d4["F_driving"][:]./1e12  .< 0.8*maxF[3] .&& (d4["time"].> d4["time"][ind10]) ) ]
    @show t_t20 = t20["time"][ findfirst( t20["F_driving"][:]./1e12 .< 0.8*maxF[4] .&& (t20["time"].> t20["time"][ind20]) ) ]
    @show t_t30 = t30["time"][ findfirst( t30["F_driving"][:]./1e12 .< 0.8*maxF[5] .&& (t30["time"].> t30["time"][ind30]) ) ]
    @show t_t45 = t45["time"][ findfirst( t45["F_driving"][:]./1e12 .< 0.8*maxF[6] .&& (t45["time"].> t45["time"][ind45]) ) ]
    @show t_t80 = t80["time"][ findfirst( t80["F_driving"][:]./1e12 .< 0.8*maxF[7] .&& (t80["time"].> t80["time"][ind80]) ) ]
    @show t_t90 = t90["time"][ findfirst( t90["F_driving"][:]./1e12 .< 0.8*maxF[7] .&& (t90["time"].> t90["time"][ind90]) ) ]


    t_t = [t_t0 t_t5 t_t10 t_t20 t_t30 t_t45 t_t80 t_t90][:]

    degree = 2
    p = fit(θ, t_t, degree)
    fit_time = coeffs(p)[2]*θ_fit .+ coeffs(p)[3]*θ_fit.^2 .+ coeffs(p)[1]
    fit_time_data = coeffs(p)[2]*θ .+ coeffs(p)[3]*θ.^2 .+ coeffs(p)[1]
   
    

    ax2 = Axis(f[2, 1:2], title = L"$$B) Characteristic time Vs. initial fabric angle", xlabel = L"$\theta_\textrm{ini}$",  ylabel = L"$t_\textrm{char}$ [Ma]")
    lines!(ax2, θ_fit, fit_time, linestyle=:dash, color=:black)
    scatter!(ax2, θ, t_t, marker=:xcross, markersize=18)
    text!(ax2, 13, 0.5, text=L"$t_\text{pred} = 0.15 - 0.16\theta_\text{ini} - 0.002\theta_\text{ini}^2$", rotation=0.)

    ylims!(ax2, 0, 4.2)

    @show coeffs(p)


    ##########################################################################

    time_t0  = t0["time"][ind5:end]  ./fit_time_data[1] 
    time_t5  = t5["time"][ind5:end]  ./fit_time_data[2] 
    time_d4  = d4["time"][ind5:end]  ./fit_time_data[3]
    time_t20 = t20["time"][ind20:end]./fit_time_data[4]
    time_t30 = t30["time"][ind30:end]./fit_time_data[5]
    time_t45 = t45["time"][ind45:end]./fit_time_data[6]
    time_t80 = t80["time"][ind80:end]./fit_time_data[7]
    time_t90 = t90["time"][ind90:end]./fit_time_data[8]

    ax2 = Axis(f[3, 1:2], title = L"$$C) Scaled driving force through scaled time", xlabel = L"$t / t_\textrm{pred}$",  ylabel = L"$F_\textrm{D} / F_\textrm{pred}$")
    # lines!(ax2, t5["time"],   t5["F_driving"]./1e12 ./ fit_force_data[1]  , label=L"$\theta_\textrm{ini} = 5$")  
    # lines!(ax2, t80["time"], t80["F_driving"]./1e12 ./ fit_force_data[5] , label=L"$\theta_\textrm{ini} = 80$")  
    # lines!(ax2, t20["time"], t20["F_driving"]./1e12 ./ fit_force_data[2] , label=L"$\theta_\textrm{ini} = 20$")  
    # lines!(ax2, t30["time"], t30["F_driving"]./1e12 ./ fit_force_data[3] , label=L"$\theta_\textrm{ini} = 30$")  
    # lines!(ax2, t45["time"] ./1, t45["F_driving"]./1e12 ./ fit_force_data[4] , label=L"$\theta_\textrm{ini} = 45$")  
    # axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14)
    lines!(ax2, time_t0 .-time_t0[1],   t0["F_driving"][ind5:end]./1e12  ./ fit_force_data[1]  , label=L"$\theta_\textrm{ini} = 0$")  
    lines!(ax2, time_t5 .-time_t5[1],   t5["F_driving"][ind5:end]./1e12  ./ fit_force_data[2]  , label=L"$\theta_\textrm{ini} = 5$")  
    lines!(ax2, time_d4 .-time_d4[1],   d4["F_driving"][ind5:end]./1e12  ./ fit_force_data[3]  , label=L"$\theta_\textrm{ini} = 10$")  
    lines!(ax2, time_t20.-time_t20[1], t20["F_driving"][ind20:end]./1e12 ./ fit_force_data[4] , label=L"$\theta_\textrm{ini} = 20$")  
    lines!(ax2, time_t30.-time_t30[1], t30["F_driving"][ind30:end]./1e12 ./ fit_force_data[5] , label=L"$\theta_\textrm{ini} = 30$")  
    lines!(ax2, time_t45.-time_t45[1], t45["F_driving"][ind45:end]./1e12 ./ fit_force_data[6] , label=L"$\theta_\textrm{ini} = 45$")  
    lines!(ax2, time_t80.-time_t80[1], t80["F_driving"][ind80:end]./1e12 ./ fit_force_data[7] , label=L"$\theta_\textrm{ini} = 80$")  
    lines!(ax2, time_t90.-time_t90[1], t90["F_driving"][ind90:end]./1e12 ./ fit_force_data[8] , label=L"$\theta_\textrm{ini} = 90$")  

    axislegend(framevisible=false, position=:rt, orientation = :horizontal, labelsize=14, nbanks=3)

    display(f)  

    save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/DrivingForceCollapseTheta.png", f, px_per_unit = 4)     

end

main()