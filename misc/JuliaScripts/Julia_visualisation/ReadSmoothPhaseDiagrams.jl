import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))

using GLMakie, FileIO, ImageCore, Printf

function ReadBinFileVisualise(fname, nx, ny)
    data = zeros(Float64, nx, ny)
    ρ = read!(fname, data)
    return ρ 
end

function Diffusion!(X, nt, Δt, K, dT, dP)
    for _=1:nt
        dρdT = diff(X[:,2:end-1], dims=1) / dT
        d2ρdT2 = diff(dρdT, dims=1) / dT
        dρdP = diff(X[2:end-1,:], dims=2) / dP
        d2ρdP2 = diff(dρdP, dims=2) / dP
        X[2:end-1,2:end-1] .+= Δt.*K.*(d2ρdT2 .+ d2ρdP2)
    end
end

function main()
    nT      = 675
    nP      = 2500
    Tmin    = 398+1e-3   # K
    Tmax    = 800+273
    Pmin    = 0.0090     # Pa
    Pmax    = 49.9890 #/10*1e9
    dT      = (Tmax-Tmin)/(nT-1)
    dP      = (Pmax-Pmin)/(nP-1)
    T       = LinRange(Tmin, Tmax, nT)
    P       = LinRange(Pmin, Pmax, nP)

    # Read in
    pathMD  = "/Users/tduretz/REPO/MDOODZ7.0/IMPORT/"

    # ρ_nsm010 = ReadBinFileVisualise(path*"SiO2_nsm010.dat", 675, 2500)

    pathTD  = "/Users/tduretz/REPO/Thermodynamics.jl/"
    ρ_HP98   = ReadBinFileVisualise(pathTD*"SiO2_julia_revised_v4_HP98.dat", 675, 2500)

    pathTD  = "/Users/tduretz/REPO/Thermodynamics.jl/"
    ρ_HP11   = ReadBinFileVisualise(pathTD*"SiO2_julia_revised_v4_HP11.dat", 675, 2500)

    # Diffusion smoothing
    K              = 1.0
    Δt             = min(dP,dT)^2/K/4.2
    nt             = 0
    ρ_HP98_nsm10 = copy(ρ_HP98)
    ρ_HP11_nsm10 = copy(ρ_HP11)
    Diffusion!(ρ_HP98_nsm10, nt, Δt, K, dT, dP)
    Diffusion!(ρ_HP11_nsm10, nt, Δt, K, dT, dP)
    
    # heatmap(T, P, ρ_nsm010' .- ρ_Gerya_smooth' ) 
    # heatmap(T, P, ρ_Gerya' .- ρ_Gerya_smooth' ) 
    # p = plot(xlabel="P [kbar]", ylabel="ρ [kg/m³]")
    # plot!( ρ_nsm010[ 400, :], label="HP")
    # plot!( ρ_Gerya_smooth[ 400, :], label="Gerya (2004)")

    f = Figure(resolution = (1200, 1000))
    iT400 = findfirst(T.>997)
    ax1 = Axis(f[1, 1],     
    title = "HP98 - ρ [kg/m³]",
    xlabel = "T [K]",
    ylabel = "P [kbar]",
    )
    GLMakie.heatmap!(ax1, T.-273.15, P, ρ_HP98_nsm10)
    GLMakie.xlims!(ax1, 400., 800.)

    ax2 = Axis(f[1, 2],     
    title = "HP11 - ρ [kg/m³]",
    xlabel = "T [C]",
    ylabel = "P [kbar]",
    )
    GLMakie.heatmap!(ax2, T, P, ρ_HP11_nsm10)
    GLMakie.xlims!(ax2, 400., 800.)

    ax = Axis(f[2, 1:2],     
        title = "@T = $(T[iT400]) K",
        ylabel = "ρ [kg/m³]",
        xlabel = "P [kbar]",
    )
    lines!(ax, P, ρ_HP11_nsm10[ iT400, :], label="HP11")
    lines!(ax, P, ρ_HP98_nsm10[ iT400, :], label="HP98")
    axislegend(framevisible = false, position = :lb)
    current_figure()
    DataInspector(f)
    display(f)

    write(pathMD*"SiO2_HP98_nsm010.dat", ρ_HP98_nsm10)
    write(pathMD*"SiO2_HP11_nsm010.dat", ρ_HP11_nsm10)
end

main()