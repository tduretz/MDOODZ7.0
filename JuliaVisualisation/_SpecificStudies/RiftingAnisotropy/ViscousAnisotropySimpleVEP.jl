using CairoMakie, Printf
import Statistics:mean


function main_simple_ani_vis()
    # Material parameters
    δ          = 3.0 
    ηv         = 1 # normal viscosity
    G          = 1.0 
    Δt         = 1e-1
    ηe         = G*Δt
    ηve        = (1/ηv + 1/ηe)^-1
    C          = 5 
    nt         = 100

    # Kinematics
    pure_shear = 1
    ε̇xxd       = pure_shear*5.0
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = ((1.0 - pure_shear)*5.0)
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    # Arrays
    θ          = LinRange( 0.0, π/2, 51 ) 
    τii_cart   = zero(θ)
    τii_cart2  = zero(θ)
    τii_cart_y = zero(θ)
    τii_min    = zero(θ)
    τii_cart_MD7 = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    τii_rot3   = zero(θ)
    plastic    = zero(θ)

    τxx_cart   = zero(θ)
    τyy_cart   = zero(θ)
    τxy_cart   = zero(θ)
    τxx_rot2  = zero(θ)
    τyy_rot2  = zero(θ)
    τxy_rot2  = zero(θ)

    τxx_rot    = zero(θ)
    τyy_rot    = zero(θ)
    τxy_rot    = zero(θ)
    σ1_angle   = zero(θ)


    for it=1:nt
    # Loop over all orientations
    for i in eachindex(θ)
            plastic[i] = 0
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
            ε̇_rot_eff   = [ε̇_rot[1,1]+τ0_rot[1,1]/2/ηe; ε̇_rot[2,2]+τ0_rot[2,2]/2/ηe; ε̇_rot[1,2]+τ0_rot[1,2]/2/ηe*δ]
            D           = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/δ;]
            τ           = D*ε̇_rot_eff
            Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
            τII         = sqrt(Y2)
            # Plasticty
            I2   = 0.5*(ε̇_rot_eff[1]^2 + ε̇_rot_eff[2]^2) + ε̇_rot_eff[3]^2
            F    = τII - C
            γ̇    = 0.0
            ηvep = ηve
            if F>0
                plastic[i] = 1
                γ̇    = F/ηve
                τII  = τII - γ̇*ηve
                ηvep = τII/2/sqrt(I2)
                τ[1]  = 2*ηvep*ε̇_rot_eff[1]
                τ[2]  = 2*ηvep*ε̇_rot_eff[2]
                τ[3]  = 2*ηvep*ε̇_rot_eff[3]/δ
                Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
                τIIc         = sqrt(Y2)
                @show τII, τIIc
            end
            D           = 2*ηvep*[1 0 0; 0 1 0; 0 0 1.0/δ;]
            τ           = D*ε̇_rot_eff
            @show  τ
            τii_rot1[i] = τII

            ########################################
            ########################################
            ########################################

            # Check strain rate components
            τxx         = τ[1]
            τyy         = τ[2]
            τxy         = τ[3]
            # ε̇xxd_chk    = τxx/2/ηv
            # ε̇xyd_chk    = τxy/2/(ηv/a^2)
            # @show (ε̇_rot[1,1], ε̇xxd_chk)
            # @show (ε̇_rot[1,2], ε̇xyd_chk)
            τxx_rot[i]   = τxx
            τyy_rot[i]   = τyy
            τxy_rot[i]   = τxy
            # Principal plane invariant
            I2           = 0.5*(ε̇_rot_eff[1]^2 + ε̇_rot_eff[2]^2) + ε̇_rot_eff[3]^2
            τii_rot2[i]  = 2*ηvep*sqrt(I2)
            # Cartesian invariant...
            I2           = 0.5*(ε̇_rot_eff[1]^2 + ε̇_rot_eff[2]^2) + ε̇_rot_eff[3]^2/δ^2
            τii_cart2[i] = 2*ηvep*sqrt(I2)
            # Rotate stress back
            τ_rot        = [τ[1] τ[3]; τ[3] τ[2]]
            τ            = Q'*τ_rot*Q
            J2           = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
            τii_cart[i]  = sqrt(J2)  
            τxx_cart[i]  = τ[1,1]
            τyy_cart[i]  = τ[2,2]
            τxy_cart[i]  = τ[1,2]

            ########################################
            ########################################
            ########################################

            # Prediction of stress using viscosity tensor in Cartesian plane
            two      = 2.
            n        = θ[i] + π/2
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
            ε̇_eff    = ε .+ inv(Dani) * τ0 ./(2*ηe)
            τ_MD7_v1 = 2*ηvep*Dani*ε̇_eff

            # Predicts stress in principal plane
            I2           = 0.5*(ε̇_eff[1]^2 + ε̇_eff[2]^2) + ε̇_eff[3]^2
            τii_rot3[i]  = 2*ηvep*sqrt(I2)

            # Predicts minimum stress
            I2           = 0.5*(ε̇_eff[1]^2 + ε̇_eff[2]^2) + ε̇_eff[3]^2/δ^2
            τii_min[i] = 2*ηvep*sqrt(I2)

            τii_cart_MD7[i] = sqrt(0.5*(τ_MD7_v1[1]^2 + τ_MD7_v1[2]^2) + τ_MD7_v1[3]^2)

            τii_cart_y[i]   = 1/(1+(δ-1)*sin(2*θ[i])) * (C)

            eta_ve = ηv
            delta  = δ
            theta  = θ[i]
            Exx = ε̇xxd
            Eyy = ε̇yyd
            Exy = ε̇xyd
            τii_cart_y[i] = (C^2/δ) / sqrt(C^2*cos(2*θ[i])^2 + C^2/δ^2*sin(2*θ[i])^2)
            #
            τii_cart_y[i] = C .* sqrt((0.25 * Exx .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Exx .^ 2 .* delta .^ 2 - 0.25 * Exx .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Exx .^ 2 + 0.5 * Exx .* Exy .* delta .^ 2 .* sin(4 * theta) - 0.5 * Exx .* Exy .* sin(4 * theta) - 0.5 * Exx .* Eyy .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.5 * Exx .* Eyy .* delta .^ 2 + 0.5 * Exx .* Eyy .* cos(2 * theta) .^ 2 - 0.5 * Exx .* Eyy - 1.0 * Exy .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + 1.0 * Exy .^ 2 .* delta .^ 2 + 1.0 * Exy .^ 2 .* cos(2 * theta) .^ 2 - 0.5 * Exy .* Eyy .* delta .^ 2 .* sin(4 * theta) + 0.5 * Exy .* Eyy .* sin(4 * theta) + 0.25 * Eyy .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Eyy .^ 2 .* delta .^ 2 - 0.25 * Eyy .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Eyy .^ 2) ./ (0.5 * Exx .^ 2 + 1.0 * Exy .^ 2 + 0.5 * Eyy .^ 2)) ./ delta

            τxx_rot2[i]  = C   * cos(2*θ[i]) 
            τyy_rot2[i]  =-C   * cos(2*θ[i])
            τxy_rot2[i]  = C/δ * sin(2*θ[i])
            # τii_cart_y[i]   = C/δ*sqrt(0.5* (-δ^4*cos(2*θ[i])^2 + δ^4 + δ^3*sin(4*θ[i]) + 2*δ^2*cos(2*θ[i])^2 - δ*sin(4*θ[i]) -1/2*cos(4*θ[i]) + 1/2 )) 

            τii_cart_y[i]  = C*sqrt((δ^2-1)*cos(2*θ[i])^2 + 1)/δ

            ev = eigvecs(τ)
            σ1_angle[i] = atand(ev[2,1]/ev[1,1])
            @show τ

            # @show ' '
            # @show abs((τxx_cart[i]-τ_MD7_v1[1])/τxx_cart[i]), abs((τxy_cart[i]- τ_MD7_v1[3])/τxy_cart[i])

        end
    end

    f = Figure(size = (500, 500), fontsize=18)
    ax1 = Axis(f[1, 1], title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    # Cartesian
    lines!(ax1, θ*180/π, τii_cart, label="τii_cart" )
    scatter!(ax1, θ[2:2:end]*180/π, τii_cart2[2:2:end], label="τii_cart2" )
    scatter!(ax1, θ[3:2:end]*180/π, τii_cart_MD7[3:2:end], label="τii_cart_MD7", marker=:dtriangle )
    lines!(ax1, θ*180/π, τii_cart_y, label="plastic yield", marker=:dtriangle )
    scatter!(ax1, θ[3:3:end]*180/π, τii_min[3:3:end], label="τii_min", marker=:xcross )
    # Principal plane
    lines!(ax1, θ*180/π, τii_rot2, label="τii_rot2" )
    scatter!(ax1, θ[2:2:end]*180/π, τii_rot1[2:2:end], label="τii_rot1" )
    scatter!(ax1, θ[3:3:end]*180/π, τii_rot3[3:3:end], label="τii_rot3", marker=:xcross )

    # axislegend()

    ax2 = Axis(f[1, 2], title="Plastic", xlabel="θ", ylabel="τᵢᵢ")
    scatter!(ax2, θ*180/π, plastic, label="plastic", marker=:xcross )

    ax3 = Axis(f[2, 1:2], title="Flow enveloppe", xlabel=L"$τ_{xx}$", ylabel=L"$τ_{xy}$", aspect=DataAspect())
    lines!(ax3, τxx_cart, τxy_cart, marker=:xcross )
    scatter!(ax3, τxx_rot, τxy_rot, marker=:xcross )
    scatter!(ax3, τxx_rot2, τxy_rot2, marker=:xcross )

    


    @info "Control points"
    @printf("| dir. angle | lay. angle | tauxx    | tauxy   | σ1 angle |\n")
    i = 1
    @printf("| %+1.4f   | %+1.4f   | %+1.4f | %+1.4f | %+1.4f |\n", (θ[i] - π/2)*180/π,  θ[i], τxx_cart[i], τxy_cart[i], σ1_angle[i])

    i = Int(round(length(θ)*1//3))
    @printf("| %+1.4f   | %+1.4f   | %+1.4f  | %+1.4f | %+1.4f |\n", (θ[i] - π/2)*180/π, θ[i], τxx_cart[i], τxy_cart[i], σ1_angle[i])

    i = Int(round(length(θ)*1//2))
    @printf("| %+1.4f   | %+1.4f   | %+1.4f  | %+1.4f | %+1.4f |\n", (θ[i] - π/2)*180/π, θ[i], τxx_cart[i], τxy_cart[i], σ1_angle[i])

    i = Int(round(length(θ)*2//3))
    @printf("| %+1.4f  | %+1.4f    | %+1.4f  | %+1.4f | %+1.4f |\n", (θ[i] - π/2)*180/π, θ[i], τxx_cart[i], τxy_cart[i], σ1_angle[i])

    display(f)

end

main_simple_ani_vis()