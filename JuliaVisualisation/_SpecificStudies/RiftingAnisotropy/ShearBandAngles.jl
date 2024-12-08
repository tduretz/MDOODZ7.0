using CairoMakie, Printf, StaticArrays
import Statistics:mean


function ShearBandAngle()

    ϕ = 35.0

    θ = LinRange(0, 90, 150)
    δ = LinRange(1, 10, 200)

    X = 1 .+ (δ .- 1).*sind.(2*θ') 
    
    ϕeff = asind.(sind.(ϕ)./X)

    α = 45 .- ϕeff./2 

    @show size(α )

    f = Figure(size = (500, 500), fontsize=18)
    ax1 = Axis(f[1, 1],  title = L"$$Shear band angle relative to $σ_1$", ylabel = L"$δ$ [-]", xlabel = L"$θ$ [$^\circ$]",  xticklabelcolor = :black, xaxisposition=:bottom)
    hm = heatmap!(ax1, θ, δ, α')
    Colorbar(f[1, 2], hm, label = L"$α$ [$^\circ$]", width = 20, labelsize = 18, ticklabelsize = 18 )
    colgap!(f.layout, 20)
    f

    Exx = 5
    Eyy = -Exx
    Exy = 1e-6

    # Exx = 0
    # Eyy = -Exx
    # Exy = 1

    x2  = 1
    delta = δ
    theta = θ' .*π/180
    x2 = 1
    x1 = -2 * Exy .* x2 ./ (Exx - Eyy - sqrt(Exx .^ 2 - 2 * Exx .* Eyy + 4 * Exy .^ 2 + Eyy .^ 2))
    
    x1 =  @. 195959179422654.0 * x2 .* (-0.125 * Exx .* delta .* (1.0 - cos(4.0 * theta)) + Exx .* delta .* sin(theta) .^ 4 + 0.125 * Exx .* (1.0 - cos(4.0 * theta)) - Exx .* sin(theta) .^ 4 + 0.5 * Exy .* delta .* (1.0 - cos(4.0 * theta)) - 0.5 * Exy .* (1.0 - cos(4.0 * theta)) + Exy + 0.125 * Eyy .* delta .* (1.0 - cos(4.0 * theta)) - Eyy .* delta .* sin(theta) .^ 4 - 0.125 * Eyy .* (1.0 - cos(4.0 * theta)) + Eyy .* sin(theta) .^ 4) ./ (48989794855663.6 * Exx .* delta .* (1.0 - cos(4.0 * theta)) - 97979589711327.1 * Exx .* delta - 48989794855663.6 * Exx .* (1.0 - cos(4.0 * theta)) + 48989794855663.6 * Exy .* delta .* (1.0 - cos(4.0 * theta)) - 391918358845308.0 * Exy .* delta .* sin(theta) .^ 4 - 48989794855663.6 * Exy .* (1.0 - cos(4.0 * theta)) + 391918358845308.0 * Exy .* sin(theta) .^ 4 - 48989794855663.6 * Eyy .* delta .* (1.0 - cos(4.0 * theta)) + 97979589711327.1 * Eyy .* delta + 48989794855663.6 * Eyy .* (1.0 - cos(4.0 * theta)) + 1.92e+15 * sqrt(0.0833333333333333 * Exx .^ 2 .* delta .^ 2 .* sin(theta) .^ 8 - 0.125 * Exx .^ 2 .* delta .^ 2 .* sin(theta) .^ 6 + 0.0729166666666667 * Exx .^ 2 .* delta .^ 2 .* sin(theta) .^ 4 - 0.0208333333333333 * Exx .^ 2 .* delta .^ 2 .* sin(theta) .^ 2 + 0.00260416666666667 * Exx .^ 2 .* delta .^ 2 - 0.166666666666667 * Exx .^ 2 .* delta .* sin(theta) .^ 8 + 0.25 * Exx .^ 2 .* delta .* sin(theta) .^ 6 - 0.125 * Exx .^ 2 .* delta .* sin(theta) .^ 4 + 0.0208333333333333 * Exx .^ 2 .* delta .* sin(theta) .^ 2 + 0.0833333333333333 * Exx .^ 2 .* sin(theta) .^ 8 - 0.125 * Exx .^ 2 .* sin(theta) .^ 6 + 0.0520833333333333 * Exx .^ 2 .* sin(theta) .^ 4 + 0.0416666666666667 * Exx .* Exy .* delta .^ 2 .* sin(theta) .^ 4 - 0.0208333333333333 * Exx .* Exy .* delta .^ 2 .* sin(theta) .^ 2 - 0.0416666666666667 * Exx .* Exy .* sin(theta) .^ 4 + 0.0208333333333333 * Exx .* Exy .* sin(theta) .^ 2 - 0.166666666666667 * Exx .* Eyy .* delta .^ 2 .* sin(theta) .^ 8 + 0.25 * Exx .* Eyy .* delta .^ 2 .* sin(theta) .^ 6 - 0.145833333333333 * Exx .* Eyy .* delta .^ 2 .* sin(theta) .^ 4 + 0.0416666666666667 * Exx .* Eyy .* delta .^ 2 .* sin(theta) .^ 2 - 0.00520833333333333 * Exx .* Eyy .* delta .^ 2 + 0.333333333333333 * Exx .* Eyy .* delta .* sin(theta) .^ 8 - 0.5 * Exx .* Eyy .* delta .* sin(theta) .^ 6 + 0.25 * Exx .* Eyy .* delta .* sin(theta) .^ 4 - 0.0416666666666667 * Exx .* Eyy .* delta .* sin(theta) .^ 2 - 0.166666666666667 * Exx .* Eyy .* sin(theta) .^ 8 + 0.25 * Exx .* Eyy .* sin(theta) .^ 6 - 0.104166666666667 * Exx .* Eyy .* sin(theta) .^ 4 + 0.333333333333333 * Exy .^ 2 .* delta .^ 2 .* sin(theta) .^ 8 - 0.5 * Exy .^ 2 .* delta .^ 2 .* sin(theta) .^ 6 + 0.208333333333333 * Exy .^ 2 .* delta .^ 2 .* sin(theta) .^ 4 - 0.666666666666667 * Exy .^ 2 .* delta .* sin(theta) .^ 8 + Exy .^ 2 .* delta .* sin(theta) .^ 6 - 0.5 * Exy .^ 2 .* delta .* sin(theta) .^ 4 + 0.0833333333333333 * Exy .^ 2 .* delta .* sin(theta) .^ 2 + 0.333333333333333 * Exy .^ 2 .* sin(theta) .^ 8 - 0.5 * Exy .^ 2 .* sin(theta) .^ 6 + 0.291666666666667 * Exy .^ 2 .* sin(theta) .^ 4 - 0.0833333333333333 * Exy .^ 2 .* sin(theta) .^ 2 + 0.0104166666666667 * Exy .^ 2 - 0.0416666666666667 * Exy .* Eyy .* delta .^ 2 .* sin(theta) .^ 4 + 0.0208333333333333 * Exy .* Eyy .* delta .^ 2 .* sin(theta) .^ 2 + 0.0416666666666667 * Exy .* Eyy .* sin(theta) .^ 4 - 0.0208333333333333 * Exy .* Eyy .* sin(theta) .^ 2 + 0.0833333333333333 * Eyy .^ 2 .* delta .^ 2 .* sin(theta) .^ 8 - 0.125 * Eyy .^ 2 .* delta .^ 2 .* sin(theta) .^ 6 + 0.0729166666666667 * Eyy .^ 2 .* delta .^ 2 .* sin(theta) .^ 4 - 0.0208333333333333 * Eyy .^ 2 .* delta .^ 2 .* sin(theta) .^ 2 + 0.00260416666666667 * Eyy .^ 2 .* delta .^ 2 - 0.166666666666667 * Eyy .^ 2 .* delta .* sin(theta) .^ 8 + 0.25 * Eyy .^ 2 .* delta .* sin(theta) .^ 6 - 0.125 * Eyy .^ 2 .* delta .* sin(theta) .^ 4 + 0.0208333333333333 * Eyy .^ 2 .* delta .* sin(theta) .^ 2 + 0.0833333333333333 * Eyy .^ 2 .* sin(theta) .^ 8 - 0.125 * Eyy .^ 2 .* sin(theta) .^ 6 + 0.0520833333333333 * Eyy .^ 2 .* sin(theta) .^ 4))
    s1_angle =  90 .+ atand.(x2./x1)


    s1_angle2 = zero(s1_angle)

    # s1_angle = zeros(length(Δ, Mak), length(θ))
    for i=1:length(Δ, Mak), j=1:length(θ)
        delta = δ[i]
        theta = θ[j]*π/180


        # Prediction of stress using viscosity tensor in Cartesian plane
        ηvep     = 1.
        two      = 2.
        n        = theta + π/2
        nx       = cos(n); ny = sin(n)
        N        = [nx; ny]
        d0       = 2*N[1]^2*N[2]^2
        d1       = N[1]*N[2] * (-N[1]^2 + N[2]^2)
        C_ANI    = [-d0 d0 -two*d1; d0 -d0 two*d1; -d1 d1 -two*(1/2-d0)]      ###### !!! - sign in D33 -a0
        C_ISO    = [1 0 0; 0 1 0; 0 0 two*1//2]
        ani      = 1. - 1.  ./ delta
        Dani     = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
        ε        = [Exx; Eyy; Exy]
        ε̇_eff    = ε
        τ        = 2*ηvep*Dani*ε̇_eff
        ev       = eigvecs([τ[1] τ[3]; τ[3] τ[2]])
        s1_angle[i,j] =  atand(ev[2,1]/ev[1,1]) 
    
        x2 = 1.
        x1 = x2 .* (-Exx .* delta .* sin(4.0 * theta) + Exx .* sin(4.0 * theta) + 4.0 * Exy .* delta .* cos(2.0 * theta) .^ 2 - 4.0 * Exy .* delta - 4.0 * Exy .* cos(2.0 * theta) .^ 2 + Eyy .* delta .* sin(4.0 * theta) - Eyy .* sin(4.0 * theta)) ./ (Exx .* delta .* cos(4.0 * theta) + Exx .* delta - Exx .* cos(4.0 * theta) + Exx - 8.0 * Exy .* delta .* sin(theta) .^ 3 .* cos(theta) + 8.0 * Exy .* delta .* sin(theta) .* cos(theta) .^ 3 + 8.0 * Exy .* sin(theta) .^ 3 .* cos(theta) - 8.0 * Exy .* sin(theta) .* cos(theta) .^ 3 - Eyy .* delta .* cos(4.0 * theta) - Eyy .* delta + Eyy .* cos(4.0 * theta) - Eyy + 4.0 * sqrt(0.25 * Exx .^ 2 .* delta .^ 2 .* cos(2.0 * theta) .^ 2 - 0.25 * Exx .^ 2 .* cos(2.0 * theta) .^ 2 + 0.25 * Exx .^ 2 + 0.5 * Exx .* Exy .* delta .^ 2 .* sin(4.0 * theta) - 0.5 * Exx .* Exy .* sin(4.0 * theta) - 0.5 * Exx .* Eyy .* delta .^ 2 .* cos(2.0 * theta) .^ 2 + 0.5 * Exx .* Eyy .* cos(2.0 * theta) .^ 2 - 0.5 * Exx .* Eyy - Exy .^ 2 .* delta .^ 2 .* cos(2.0 * theta) .^ 2 + Exy .^ 2 .* delta .^ 2 + Exy .^ 2 .* cos(2.0 * theta) .^ 2 - 0.5 * Exy .* Eyy .* delta .^ 2 .* sin(4.0 * theta) + 0.5 * Exy .* Eyy .* sin(4.0 * theta) + 0.25 * Eyy .^ 2 .* delta .^ 2 .* cos(2.0 * theta) .^ 2 - 0.25 * Eyy .^ 2 .* cos(2.0 * theta) .^ 2 + 0.25 * Eyy .^ 2))
        s1_angle2[i,j] =  atand.(x2./x1)
    end

    @show extrema(s1_angle)

    f = Figure(size = (500, 500), fontsize=18)
    ax1 = Axis(f[1, 1],  title = L"$$Angle of $σ_1$", ylabel = L"$δ$ [-]", xlabel = L"$θ$ [$^\circ$]",  xticklabelcolor = :black, xaxisposition=:bottom)
    hm = heatmap!(ax1, θ, δ, s1_angle')
    Colorbar(f[1, 2], hm, label = L"$α$ [$^\circ$]", width = 20, labelsize = 18, ticklabelsize = 18 )
    colgap!(f.layout, 20)

    ax2 = Axis(f[2, 1],  title = L"$$Angle of $σ_1$ ($δ = 3$)", ylabel = L"$δ$ [-]", xlabel = L"$θ$ [$^\circ$]",  xticklabelcolor = :black, xaxisposition=:bottom)
    ind = findmin(abs.(δ .- 3))[2]
    lines!(ax2, θ, s1_angle[ind,:])
    lines!(ax2, θ, s1_angle2[ind,:])

    # ind = findmin(abs.(δ .- 6))[2]
    # lines!(ax2, θ, s1_angle[ind,:])

    display(f)

    

end


ShearBandAngle()