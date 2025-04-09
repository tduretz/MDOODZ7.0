using Plots, Printf
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    ani_fac_v  = 4.0 
    # Kinematics
    pure_shear = 1
    ε̇xxd       = pure_shear*5.0
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = (1.0-pure_shear)*5.0 + 5.0/3
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    # Power law
    npwl      = 10
    τbg       = 2.0
    Bpwl      = 2^npwl*ε̇bg/τbg^(npwl)
    τ_chk     = 2*Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl)
    ηpwl      =   Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl-1)
    τ_chk     = 2*ηpwl*ε̇bg
    Cpwl      = (2*Bpwl^(-1/npwl))^(-npwl)
    if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
    # --------------------- TEST 1: Isotropic power law --------------------- #
    # TEST 1: Isotropic power law
    τxx = 2*ηpwl*ε̇xxd
    τxy = 2*ηpwl*ε̇xyd
    Y2  = τxx^2 + τxy^2
    # Check strain rates
    ε̇xxd_chk = Cpwl*Y2^((npwl-1)/2)*τxx
    ε̇xyd_chk = Cpwl*Y2^((npwl-1)/2)*τxy
    println("ε̇xx check: ", ε̇xxd ≈ ε̇xxd_chk) 
    println("ε̇xy check: ", ε̇xyd ≈ ε̇xyd_chk)
    # -------------------- TEST 2: Anisotropic power law -------------------- #
    # Arrays
    θ          = LinRange( 0.0, π, 51 ) .-0* π/2
    τii_cart1  = zero(θ) 
    τii_cart2  = zero(θ)
    ε̇ii_rot    = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    η_rot      = zero(θ)
    η_cart     = zero(θ)
    # Loop over all orientations
    for i in eachindex(θ)
        @show θ[i] 
        # Transformation matrix: towards principal plane
        Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
        ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
        ε̇_rot       = Q*ε̇_tens*Q'  
        # Viscous power law
        ε̇           = [ε̇_rot[1,1]; ε̇_rot[2,2]; ε̇_rot[1,2]]
        I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
        ηpwl        =   Bpwl^(-1.0/npwl)*sqrt(I2)^(1.0/npwl-1)
        D           = 2*ηpwl*[1 0 0; 0 1 0; 0 0 1.0/ani_fac_v;]
        τ           = D*ε̇
        # ---> Independent on orientation (objective) !!!!!!
        Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac_v^2
        τii_rot1[i] = sqrt(Y2)
        ε̇ii_rot[i]  = sqrt(I2)
        # Check strain rate components
        τxx         = τ[1]
        τxy         = τ[3]
        ε̇xxd_chk0   = Cpwl*Y2^((npwl-1)/2)*τxx
        ε̇xyd_chk0   = Cpwl*Y2^((npwl-1)/2)*(τxy*ani_fac_v)
        ε̇xxd_chk1   = τxx/2/ηpwl
        ε̇xyd_chk1   = τxy/2/ηpwl*ani_fac_v
        println( "ε̇xx = ", ε̇_rot[1,1], " ε̇xx check: ", ε̇xxd_chk0≈ε̇_rot[1,1], " ", ε̇xxd_chk1≈ε̇_rot[1,1])
        println( "ε̇xy = ", ε̇_rot[1,2], " ε̇xy check: ", ε̇xyd_chk0≈ε̇_rot[1,2], " ", ε̇xyd_chk1≈ε̇_rot[1,2])
        # Check strain rate invariant
        ε̇ii_pwl     = Cpwl * sqrt(Y2) ^(npwl)
        ε̇ii_pwl1    = sqrt(ε̇xxd_chk1^2 + ε̇xyd_chk1^2)
        println( "ε̇ii = ", sqrt(I2), " ε̇ii check: ", ε̇ii_pwl≈sqrt(I2), " ", ε̇ii_pwl1≈sqrt(I2))
        # Effective viscosity: compute stress invariant from strain rate invariant
        # ---> Predicts rotated stress invariant !!!!!!
        I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
        τii_rot2[i] = 2*ηpwl*sqrt(I2)
        # If one does the same by scaling the strain rate invariant by ani_fac_v^2...
        # ---> Predicts Cartesian stress invariant !!!!!!
        I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2/ani_fac_v^2
        τii_cart1[i] = 2*ηpwl*sqrt(I2)
        # Rotate stress back 
        τ_rot       = [τ[1] τ[3]; τ[3] τ[2]]
        τ           = Q'*τ_rot*Q
        # ---> Dependent on orientation (non-objective) !!!!!!
        J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
        τii_cart2[i] = sqrt(J2) 
        # Viscosities
        η_rot[i]  = ηpwl
        η_cart[i] = τii_cart1[i]/2/I2
    end
    p1 = plot(title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    p1 = plot!(θ*180/π, τii_cart1, label="τii_cart1")
    p1 = plot!(θ*180/π, τii_cart2, label="τii_cart2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π,  τii_rot1, label="τii_rot1")
    p1 = plot!(θ*180/π,  τii_rot2, label="τii_rot2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π,   ε̇ii_rot, label="ε̇ii_rot",  linewidth=0, marker =:circle)
    # p1 = plot!(θ*180/π,   η_rot, label="η_rot")
    # p1 = plot!(θ*180/π,  η_cart, label="η_cart")
end

main_simple_ani_vis()