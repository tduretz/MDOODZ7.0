using CairoMakie

function main()

    Ts = 600.
    Tl = 1200.
    α  = -0.1
    T  = LinRange(Ts, Tl, 100) 
    ϕ  = (exp.(α.*T) .- exp.(α*Ts)) ./ (exp.(α*Tl) .- exp.(α*Ts))

    f   = Figure(fontsize=25)
    ax1 = Axis(f[1, 1], title = L"$$Melt fraction", xlabel =L"$T$ [C]", ylabel = L"$ϕ$ [-]")
    lines!(ax1, T, ϕ )
    display(f)

end

main()