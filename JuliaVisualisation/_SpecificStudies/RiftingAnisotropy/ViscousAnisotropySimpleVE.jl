using CairoMakie, Printf
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    δ          = 2.0 
    ηv         = 1 # normal viscosity
    G          = 1.0
    Δt         = 1e-4
    ηe         = G*Δt
    ηve        = (1/ηv + 1/ηe)^-1
    # Kinematics
    pure_shear = 0
    ε̇xxd       = pure_shear*5.0
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = ((1.0 - pure_shear)*5.0)
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    # Arrays
    θ          = LinRange( 0.0, π/2, 51 ) 
    τii_cart   = zero(θ)
    τii_cart_MD7 = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    τii_rot3   = zero(θ)
    τxx_cart   = zero(θ)
    τyy_cart   = zero(θ)
    τxy_cart   = zero(θ)
    τxx_rot    = zero(θ)
    τyy_rot    = zero(θ)
    τxy_rot    = zero(θ)

    τxx_cart  = zero(θ)
    τyy_cart  = zero(θ)
    τxy_cart  = zero(θ)

    nt = 200

    f = Figure(size = (500, 500), fontsize=18)
    ax1 = Axis(f[1, 1], title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")

    for it=1:nt
    # Loop over all orientations
    for i in eachindex(θ)
            τxx0 = τxx_cart[i] 
            τyy0 = τyy_cart[i] 
            τxy0 = τxy_cart[i] 
            # Transformation matrix: towards principal plane
            Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
            ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
            ε̇_rot       = Q*ε̇_tens*Q'  
            τ0_tens     = [τxx0 τxy0; τxy0 τyy0]
            τ0_rot      = Q*τ0_tens*Q' 
            # Viscous 
            ε_eff       = [ε̇_rot[1,1]+τ0_rot[1,1]/2/ηe; ε̇_rot[2,2]+τ0_rot[2,2]/2/ηe; ε̇_rot[1,2]+τ0_rot[1,2]/2/ηe*δ]
            D           = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/δ;]
            τ           = D*ε_eff
            Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
            τii_rot1[i] = sqrt(Y2)
            # Check strain rate components
            τxx         = τ[1]
            τyy         = τ[2]
            τxy         = τ[3]
            # ε̇xxd_chk    = τxx/2/ηv
            # ε̇xyd_chk    = τxy/2/(ηv/a^2)
            # @show (ε̇_rot[1,1], ε̇xxd_chk)
            # @show (ε̇_rot[1,2], ε̇xyd_chk)
            τxx_rot[i]  = τxx
            τyy_rot[i]  = τyy
            τxy_rot[i]  = τxy
            # Effective viscosity: compute stress invariant from strain rate invariant
            I2          = 0.5*(ε_eff[1]^2 + ε_eff[2]^2) + ε_eff[3]^2
            τii_rot2[i] = 2*ηve*sqrt(I2)
            # if one does the same by scaling the strain rate invariant by δ...
            I2          = 0.5*(ε_eff[1]^2 + ε_eff[2]^2) + ε_eff[3]^2/δ^2
            τii_rot3[i] = 2*ηve*sqrt(I2)
            # Rotate stress back
            τ_rot       = [τ[1] τ[3]; τ[3] τ[2]]
            τ           = Q'*τ_rot*Q
            J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
            τii_cart[i] = sqrt(J2)  
            τxx_cart[i] = τ[1,1]
            τyy_cart[i] = τ[2,2]
            τxy_cart[i] = τ[1,2]

            # Prediction of stress using viscosity tensor in Cartesian plane
            two      = 2.
            n        = θ[i] - π/2
            nx       = cos(n); ny = sin(n)
            N        = [nx; ny]
            d0       = 2*N[1]^2*N[2]^2
            d1       = N[1]*N[2] * (-N[1]^2 + N[2]^2)
            C_ANI    = [-d0 d0 -two*d1; d0 -d0 two*d1; -d1 d1 -two*(1/2-d0)]      ###### !!! - sign in D33 -a0
            C_ISO    = [1 0 0; 0 1 0; 0 0 two*1//2]
            ani      = 1. - 1.  ./ δ
            Dani     = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
            ε        = [ε̇xxd; ε̇yyd; ε̇xyd]
            τ0       = [τxx0; τyy0; τxy0]
            ε_eff    = ε .+ inv(Dani) * τ0 ./(2*ηe)
            τ_MD7_v1 = 2*ηve*Dani*ε_eff

            τii_cart_MD7[i]          = sqrt(0.5*(τ_MD7_v1[1]^2 + τ_MD7_v1[2]^2) + τ_MD7_v1[3]^2)

            # @show ' '
            # @show τxx_cart[i], τxy_cart[i]
            # @show τ_MD7_v1[1], τ_MD7_v1[3]

        end
    end


    lines!(ax1, θ*180/π, τii_cart )
    scatter!(ax1, θ*180/π, τii_cart_MD7 )
    lines!(ax1, θ*180/π, τii_rot1 )
    lines!(ax1, θ*180/π, τii_rot2 )
    lines!(ax1, θ*180/π, τii_rot3 )
    display(f)

end

main_simple_ani_vis()