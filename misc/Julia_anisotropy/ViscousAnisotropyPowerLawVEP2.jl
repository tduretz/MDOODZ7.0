using Plots, Printf
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    ani_fac    = 2.0 
    # Kinematics
    pure_shear = 1
    ε̇xxd       = pure_shear*5.0
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = (1.0-pure_shear)*5.0 + 5.0/3
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    P          = 1.0
    τxx        = 0.
    τyy        = 0.
    τxy        = 0.
    # Elasticity
    nt         = 1
    Δt         = 1e0
    G          = 1.0
    K          = 3.0 
    ηe         = G*Δt
    # Plasticity
    C          = 10.25
    fric       = 0*π/180
    dil        = 0*π/180
    a1         = 1.0
    a2         = 1.0
    a3         = 1.
    ηvp        = 0.0
    am         = 1//3*(a1+a2+a3)
    # Power law viscosity
    npwl       = 1
    τbg        = 2.0
    Bpwl       = 2^npwl*ε̇bg/τbg^(npwl)
    τ_chk      = 2*Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl)
    ηpwl       =   Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl-1)
    τ_chk      = 2*ηpwl*ε̇bg
    Cpwl       = (2*Bpwl^(-1/npwl))^(-npwl)
    if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
    # -------------------- TEST: Anisotropic power law -------------------- #
    # Arrays
    θ          = LinRange( 0.0, π, 26 ) .-0* π/2
    τii_cart1  = zero(θ) 
    τii_cart2  = zero(θ)
    ε̇ii_rot    = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    # Loop over all orientations
    for i in eachindex(θ)
        # τxx, τyy, τxy = -0.0, 0.0, 0.0
        τxx, τyy, τxy = -.2, 0.2, -0.7
        @show θ[i] 
        for it=1:nt
            τxx0, τyy0, τxy0 = τxx, τyy, τxy
            # Transformation matrix: towards principal plane
            Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
            ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
            τ0_tens     = [τxx0 τxy0; τxy0 τyy0]
            # ε̇_tens     += ε̇_tens .+ τ0_tens./2.0/ηe
            ε̇_rot       = Q*ε̇_tens*Q' 
            τ0_rot      = Q*τ0_tens*Q' 
            
            # VISCO ELASTICITY
            ε̇           = [ε̇_rot[1,1]; ε̇_rot[2,2]; ε̇_rot[1,2]]
    
            ε̇           = [ε̇_rot[1,1]+τ0_rot[1,1]/2/ηe; ε̇_rot[2,2]+τ0_rot[2,2]/2/ηe; ε̇_rot[1,2]+τ0_rot[1,2]/2/ηe*ani_fac]
            I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
            ε̇ii         = sqrt(I2)
            ηpwl        =  Bpwl^(-1.0/npwl)*sqrt(I2)^(1.0/npwl-1)
            ηve         = (1.0/ηpwl + 1.0/ηe).^(-1)
            ε̇pwl        = 0.
            τii         = 0.
            for iter=1:10
                τii        = 2*ηve*ε̇ii
                ε̇pwl       = Cpwl*τii^npwl
                r          = ε̇ii - τii/2/ηe - ε̇pwl
                ∂r∂ηve     = - ε̇ii/ηe - ε̇pwl*npwl/ηve
                ηve       -= r/∂r∂ηve
                if (abs(r)<1e-9) break; end
                @show (iter, r)
            end
            D           = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/ani_fac;]
            τ           = D*ε̇
            # ---> Independent on orientation (objective) !!!!!!
            Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac^2
            τii = sqrt(Y2); τxxt = τ[1]; τyyt = τ[2]; τzzt = -τxxt-τyyt; τxyt = τ[3]

            # # PLASTICITY
            γ̇ = 0.
            F = τii - cos(fric)*C - P*sin(fric)*am + sin(fric)*( a1*τxxt + a2*τyyt + a3*τzzt)/3
            if F>0
                γ̇    = F/(ηve + ηvp + K*Δt*sin(dil)*sin(fric)*am)
                τiic = τii - γ̇*ηve
                τxxc = τxxt*(1 - γ̇*ηve/τii)
                τyyc = τyyt*(1 - γ̇*ηve/τii)
                τzzc = τzzt*(1 - γ̇*ηve/τii)
                Pc   = P    + K*Δt*sin(dil)*γ̇ 
                F    = τiic - cos(fric)*C - Pc*sin(fric)*am - ηvp*γ̇ + sin(fric)*( a1*τxxt + a2*τyyt + a3*τzzt)/3
                @show (F)
                Y2   = τiic^2
                ηvep = τiic/2/ε̇ii
                ηve  = ηvep
                D    = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/ani_fac;]
                τ    = D*ε̇ 
                Y2   = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac^2
                τiic = sqrt(Y2)
                F    = τiic - cos(fric)*C - Pc*sin(fric)*am - ηvp*γ̇ + sin(fric)*( a1*τxxt + a2*τyyt + a3*τzzt)/3
                @show (F)
            end
            τii_rot1[i] = sqrt(Y2)
            ε̇ii_rot[i]  = sqrt(I2)

            # Check strain rate components
            τxx         = τ[1]
            τxy         = τ[3]
            # ε̇xxd_v = Cpwl*τii^((npwl-1))*τxx
            # ε̇xyd_v = Cpwl*τii^((npwl-1))*(τxy*ani_fac)
            ηpwl   =  Bpwl^(-1.0/npwl)*ε̇pwl^(1.0/npwl-1)
            ε̇xxd_v = τxx/2/ηpwl
            ε̇xyd_v = τxy/2/ηpwl*ani_fac
            ε̇xxd_e = (τxx - τ0_rot[1,1])/2/ηe
            ε̇xyd_e = (τxy - τ0_rot[1,2])/2/ηe*ani_fac
            ε̇xxd_p = γ̇*τxxt/τii/2
            ε̇xyd_p = γ̇*τxyt/τii/2*ani_fac
            fxx    = ε̇_rot[1,1] - ε̇xxd_v - ε̇xxd_e - ε̇xxd_p
            fxy    = ε̇_rot[1,2] - ε̇xyd_v - ε̇xyd_e - ε̇xyd_p
            @show (fxx, fxy)

            # Effective viscosity: compute stress invariant from strain rate invariant
            # ---> Predicts rotated stress invariant !!!!!!
            I2           = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
            τii_rot2[i]  = 2*ηve*sqrt(I2)
            # If one does the same by scaling the strain rate invariant by ani_fac^2...
            # ---> Predicts Cartesian stress invariant !!!!!!
            I2           = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2/ani_fac^2
            τii_cart1[i] = 2*ηve*sqrt(I2)
            # Rotate stress back 
            τ_rot        = [τ[1] τ[3]; τ[3] τ[2]]
            τ            = Q'*τ_rot*Q
            # ---> Dependent on orientation (non-objective) !!!!!!
            J2           = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
            τii_cart2[i] = sqrt(J2) 
        end 
    end
    p1 = plot(title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    p1 = plot!(θ*180/π, τii_cart1, label="τii_cart1")
    p1 = plot!(θ*180/π, τii_cart2, label="τii_cart2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π,  τii_rot1, label="τii_rot1")
    p1 = plot!(θ*180/π,  τii_rot2, label="τii_rot2", linewidth=0, marker =:cross)
    # p1 = plot!(θ*180/π,   ε̇ii_rot, label="ε̇ii_rot",  linewidth=0, marker =:circle)
end

main_simple_ani_vis()