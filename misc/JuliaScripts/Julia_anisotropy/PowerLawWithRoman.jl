using Plots
plotlyjs()

function main()

    # Power-law creep diagrams
    logε̇ = LinRange(-17, -10, 100)
    ε̇    = 10.0.^logε̇

    #-----------------------------#
    η    = 1e20 # ---> C independent of n

    # Linear Case
    n     = 1
    C     = (2η)^-1
    τ_n1  = C^(-1/n) .* ε̇.^(1/n)

    # Power law
    n     = 5
    C     = (2η)^-1
    τ_n5  = C^(-1/n) .* ε̇.^(1/n)
    
    # Power law
    n     = 10
    C     = (2η)^-1
    τ_n10 = C^(-1/n) .* ε̇.^(1/n)

    #-----------------------------#
    ε̇ref  = 10^(-14) # ---> C defined such that there is an invariant point in the diagram
    τref  = 100e6

    n     = 1
    C     = ε̇ref * τref^(-n)
    τ_n1  = C^(-1/n) .* ε̇.^(1/n)

    n     = 5
    C     = ε̇ref * τref^(-n)
    τ_n5  = C^(-1/n) .* ε̇.^(1/n)

    n     = 10
    C     = ε̇ref * τref^(-n)
    τ_n10 = C^(-1/n) .* ε̇.^(1/n)

    #-----------------------------#
    plot(xlabel="log10 ε̇", ylabel="log10 τ" )
    plot!(log10.(ε̇), log10.(τ_n1),  label="n = 1")
    plot!(log10.(ε̇), log10.(τ_n5),  label="n = 5")
    plot!(log10.(ε̇), log10.(τ_n10), label="n = 10")

end

main()