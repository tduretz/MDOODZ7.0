using Plots, Printf
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    ani_fac_v  = 3.0 
    ηv         = 1.0 # normal viscoisty
    # Kinematics
    εbg        = 1.0
    pure_shear = 1
    εxxd       = pure_shear*εbg
    εyyd       = -εxx
    εxy        = (1.0-pure_shear)*εbg
    # Arrays
    θ          = LinRange( 0.0, π, 51 ) .- π/2
    τii_cart   = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    τii_rot3   = zero(θ)
    # Loop over all orientations
    for i in eachindex(θ)
        # Transformation matrix: towards principal plane
        Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
        ε_tens      = [εxxd εxy; εxy εyyd]
        ε_rot       = Q*ε_tens*Q'  
        # Viscous 
        ε           = [ε_rot[1,1]; ε_rot[2,2]; ε_rot[1,2]]
        D           = 2*ηv*[1 0 0; 0 1 0; 0 0 1.0/ani_fac_v;]
        τ           = D*ε
        J2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac_v^2
        τii_rot1[i] = sqrt(J2)
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
    end
    p1 = plot(title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    p1 = plot!(θ*180/π, τii_cart, label="Cartesian")
    p1 = plot!(θ*180/π, τii_rot1, label="principal")
    p1 = plot!(θ*180/π, τii_rot2, label="τii_rot2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π, τii_rot3, label="τii_rot3", linewidth=0, marker =:cross)
end

main_simple_ani_vis()