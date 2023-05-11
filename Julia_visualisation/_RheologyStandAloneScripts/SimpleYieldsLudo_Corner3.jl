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
    Pt = -0.5e7
    # Shear
    ϕs = 35.0*π/180
    ψs = 10*π/180
    # Tension
    ϕt = 85.0*π/180
    ψt = 85.0*π/180
    # Trial pressure and stress
    P    = collect(LinRange(-100e6, 300e6, 2000))
    (v, iP) = findmin(abs.(P.-Pt))
    Pshear  = P[iP:end]
    Ptens   = P[1:iP]
    τiis    = C*cos(ϕs) .+ Pshear.*sin(ϕs)
    α       = π - π/2 - ϕt 
    τs1     = C*cos(ϕs) .+ Pt.*sin(ϕs)
    d       = τs1*tan(α)
    Ct      = (d - Pt)/tan(α)
    τiit    = Ct .+ Ptens.*tan(ϕt)
    τii     = 1.1*C     .+ P.*sin(2*ϕs)
    # Evaluate trial yields 
    Fs = τii .- C*cos(ϕs) .- P.*sin(ϕs)
    Ft = τii .- Ct        .- P.*tan(ϕt)
    Fc = τii .- τs1
    # Quick and dirty mode selection
    Y = zeros(Int64, size(P))
    Y[Ft.>0.0 ]             .= 1
    Y[Fs.>0.0 ]             .= 2
    Y[Fs.>0.0 .&& Ft.>0.0 ] .= 3
    # @show Y
    # Correct Fs
    λ̇             = zero(P)
    τiic          = zero(P)
    Pc            = zero(P)
    sinψc         = zero(P)
    λ̇[Y.==2]     .= Fs[Y.==2] ./ (ηe + ξe*sin(ϕs)*sin(ψs))
    λ̇[Y.==1]     .= Ft[Y.==1] ./ (ηe + ξe*tan(ϕt)*tan(ψt))
    τiic[Y.==1]  .= τii[Y.==1] .- λ̇[Y.==1]*ηe
    Pc[Y.==1]    .= P[Y.==1]   .+ λ̇[Y.==1]*tan(ψt).*ξe
    τiic[Y.==2]  .= τii[Y.==2] .- λ̇[Y.==2]*ηe
    Pc[Y.==2]    .= P[Y.==2]   .+ λ̇[Y.==2]*sin(ψs).*ξe

    # CORNER: Brutal cowboys
    τiic[Y.==3]  .= τs1
    Pc[Y.==3]    .= Pt

    # # CORNER: Lambda style update
    λ̇[Y.==3]     .= Fc[Y.==3]./ηe
    sinψc[Y.==3] .= (Pt .- P[Y.==3]) ./ (ξe.*λ̇[Y.==3])
    τiic[Y.==3]  .= τii[Y.==3] .- λ̇[Y.==3]*ηe   
    Pc[Y.==3]    .= P[Y.==3]   .+ λ̇[Y.==3].*(sinψc[Y.==3]).*ξe

    # Viz
    p = plot()
    p = plot!(Pshear, τiis, label="shear")
    p = plot!(Ptens, τiit, label="tension")
    p = plot!(P,  τii,  label="trial")
    p = scatter!(Pc, τiic, label="corrected")
    p = plot!(ylim=(0, 3e8))
    display(p)
end

main_simple_yields()

