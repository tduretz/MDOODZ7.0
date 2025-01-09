import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using Plots
plotlyjs()

function main_simple_yields()

    # Elatic/visco-elastic predictor
    ηe = 1e20
    ξe = 4e20
    # Plastic yields
    C  = 5e7
    P1 = 0.
    # Shear
    ϕs = 35.0*π/180
    ψs = 10*π/180
    # Tension
    ϕt = 85.0*π/180
    ψt = 85.0*π/180
    # Trial presure and stres 
    P    = collect(LinRange(-100e6, 300e6, 1000))
    (v, iP) = findmin(abs.(P.-P1))
    Ps      = P[iP:end]
    Pt      = P[1:iP]
    τiis    = C*cos(ϕs) .+ Ps.*sin(ϕs)
    τiit    = C*cos(ϕs) .+ Pt.*tan(ϕt)
    τii     = 1.1*C     .+ P.*sin(2*ϕs)
    # Evaluate trial yields 
    Fs = τii .- C*cos(ϕs) .- P.*sin(ϕs)
    Ft = τii .- C*cos(ϕs) .- P.*tan(ϕt)
    # Quick and dirty mode selection
    Y = zeros(Int64, size(P))
    Y[Ft.>0.0 .&& P.<P1] .= 1
    Y[Fs.>0.0 .&& P.>P1] .= 2
    # Correct Fs
    λ̇s   = zero(P)
    λ̇t   = zero(P)
    τiic = zero(P)
    Pc   = zero(P)
    λ̇s[Y.==2]  .= Fs[Y.==2] ./ (ηe + ξe*sin(ϕs)*sin(ψs))
    λ̇t[Y.==1]  .= Ft[Y.==1] ./ (ηe + ξe*tan(ϕt)*tan(ψt))
    τiic[Y.==1] = τii[Y.==1] .- λ̇t[Y.==1]*ηe
    Pc[Y.==1]   = P[Y.==1]   .+ λ̇t[Y.==1]*tan(ψt).*ξe
    τiic[Y.==2] = τii[Y.==2] .- λ̇s[Y.==2]*ηe
    Pc[Y.==2]   = P[Y.==2]   .+ λ̇s[Y.==2]*sin(ψs).*ξe




    p = plot()
    # p = plot!(P, yield)
    p = plot!(Ps, τiis, label="shear")
    p = plot!(Pt, τiit, label="tension")
    p = plot!(P, τii, label="trial")
    p = scatter!(Pc, τiic, label="corrected")
    p = plot!(ylim=(0, 3e8))
    display(p)
end

main_simple_yields()

