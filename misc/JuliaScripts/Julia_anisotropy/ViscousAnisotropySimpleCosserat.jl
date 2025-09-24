using Plots, Printf
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    ani_fac_v  = 3.0 
    a          = sqrt(ani_fac_v)
    ηv         = 1.0 # normal viscosity
    # Kinematics
    pure_shear = 0
    ε̇xxd       = pure_shear*5.0
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = (1.0 - pure_shear)*5.0 + 5.0/3
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    # Arrays
    θ          = LinRange( 0.0, π, 51 ) .- π/2
    τii_cart   = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    τii_rot3   = zero(θ)
    τxx_cart   = zero(θ)
    τxy_cart   = zero(θ)
    τxx_rot    = zero(θ)
    τxy_rot    = zero(θ)
    # Loop over all orientations
    for i in eachindex(θ)
        # Transformation matrix: towards principal plane
        Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
        ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
        ε̇_rot       = Q*ε̇_tens*Q'  
        # Viscous 
        ε           = [ε̇_rot[1,1]; ε̇_rot[2,2]; ε̇_rot[1,2]]
        D           = 2*ηv*[1 0 0; 0 1 0; 0 0 1.0/ani_fac_v;]
        τ           = D*ε
        Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac_v^2
        τii_rot1[i] = sqrt(Y2)
        # Check strain rate components
        τxx         = τ[1]
        τxy         = τ[3]
        ε̇xxd_chk    = τxx/2/ηv
        ε̇xyd_chk    = τxy/2/(ηv/a^2)
        @show (ε̇_rot[1,1], ε̇xxd_chk)
        @show (ε̇_rot[1,2], ε̇xyd_chk)
        τxx_rot[i]  = τxx
        τxy_rot[i]  = τxy
        # Effective viscosity: compute stress invariant from strain rate invariant
        I2          = 0.5*(ε[1]^2 + ε[2]^2) + ε[3]^2
        τii_rot2[i] = 2*ηv*sqrt(I2)
        # if one does the same by scaling the strain rate invariant by ani_fac_v...
        I2          = 0.5*(ε[1]^2 + ε[2]^2) + ε[3]^2/ani_fac_v^2
        τii_rot3[i] = 2*ηv*sqrt(I2)
        # Rotate stress back
        τ_rot       = [τ[1] τ[3]; τ[3] τ[2]]
        τ           = Q'*τ_rot*Q
        J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
        τii_cart[i] = sqrt(J2)  
        τxx_cart[i] = τ[1,1]
        τxy_cart[i] = τ[1,2]
    end
    p1 = plot(title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    p1 = plot!(θ*180/π, τii_cart, label="Cartesian")
    p1 = plot!(θ*180/π, τii_rot1, label="principal")
    p1 = plot!(θ*180/π, τii_rot2, label="τii_rot2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π, τii_rot3, label="τii_rot3", linewidth=0, marker =:cross)
    # p1 = plot(title="Stress components", xlabel="θ", ylabel="τᵢᵢ")
    # p1 = plot!(θ*180/π, τxx_cart, label="τxx_cart")
    # p1 = plot!(θ*180/π, τxy_cart, label="τxy_cart")
    # p1 = plot!(θ*180/π, τxx_rot, label="τxx_rot", linewidth=0, marker =:cross)
    # p1 = plot!(θ*180/π, τxy_rot, label="τxy_rot", linewidth=0, marker =:cross)
end

main_simple_ani_vis()