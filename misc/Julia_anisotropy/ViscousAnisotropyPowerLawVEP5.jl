using Plots, Printf, HDF5
import Statistics:mean
gr()
# plotlyjs()

# Viusalise more things to check MD7 is correctly programmed

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function main_simple_ani_vis()
    # Material parameters
    ani_fac    = 2.0 
    # Kinematics
    pure_shear = 1
    ε̇xxd       = -pure_shear
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = (1.0-pure_shear)*1.0 + 0*5.0/3
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    τxx        = 0.
    τyy        = 0.
    τxy        = 0.
    P          = 1.0
    # Elasticity
    nt         = 10
    Δt         = 1e0
    G          = 1.0
    K          = 3.0  
    ηe         = G*Δt
    # Plasticity
    C          = 1.25
    fric       = 0*π/180
    dil        = 0*π/180
    ηvp        = 0.05
    # Power law
    npwl       = 3
    τbg        = 2.0
    Bpwl       = 2^npwl*ε̇bg/τbg^(npwl)
    τ_chk      = 2*Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl)
    ηpwl       =   Bpwl^(-1.0/npwl)*ε̇bg^(1.0/npwl-1)
    τ_chk      = 2*ηpwl*ε̇bg
    Cpwl       = (2*Bpwl^(-1/npwl))^(-npwl)
    if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
   
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # filename = string(path, @sprintf("Output%05d.gzip.h5", nt))        
    # model  = ExtractData( filename, "/Model/Params")
    # nvx, nvz = Int(model[4]), Int(model[5])
    # ncx, ncz = nvx-1, nvz-1
    # Ndx    = Float64.(reshape(ExtractData( filename, "/Centers/nx"),  ncx, ncz))
    # τxx  = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"),  ncx, ncz))
    # τxz  = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
    # τII  = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) );
    # τII_MD7 = mean(τII)
    # τxx_MD7 = mean(τxx)
    # τxy_MD7 = mean(τxz)
    # θ_MD7   = acosd(mean(Ndx)) + 90
    
    # -------------------- TEST: Anisotropic power law -------------------- #
    # Arrays 
    θ          = LinRange( 0.0, π, 51 )  # Orientation of layer (90 degree to director angle)
    τii_cart1  = zero(θ) 
    τii_cart2  = zero(θ)
    ε̇ii_rot    = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    η_rot      = zero(θ)
    τxx_cart   = zero(θ) 
    τxy_cart   = zero(θ) 
    τxx_rot    = zero(θ) 
    τxy_rot    = zero(θ) 
    # Loop over all layer orientations 
    for i in eachindex(θ)
        τxx, τyy, τxy = -.0, .0, 0.0
        @show θ[i] 
        for it=1:nt
            τxx0, τyy0, τxy0 = τxx, τyy, τxy
            # Transformation matrix: towards principal plane
            Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])] # transformation matrix
            ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
            τ0_tens     = [τxx0 τxy0; τxy0 τyy0]
            ε̇_rot       = Q*ε̇_tens*Q' # dev. strain rate --> from Cartesian to principal coordinates
            τ0_rot      = Q*τ0_tens*Q'# old dev. stress  --> from Cartesian to principal coordinates

            # VISCO ELASTICITY: Newton iteration because combined power-law and elasticity is non-linear
            #       Exx = εxx + τxx/2/G/Δt "effective strain rate" in principal coordinates
            #       Exy = εxy + τxy/2/(G/ani_fac)/Δt = εxy + ani_fac* τxy/2/G/Δt "effective strain rate" in principal coordinates
            ε̇           = [ε̇_rot[1,1]+τ0_rot[1,1]/2/ηe; ε̇_rot[2,2]+τ0_rot[2,2]/2/ηe; ε̇_rot[1,2]+τ0_rot[1,2]/2/ηe*ani_fac]
            I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
            ε̇ii         = sqrt(I2)
            ηpwl        =  Bpwl^(-1.0/npwl)*sqrt(I2)^(1.0/npwl-1)
            ηve         = (1.0/ηpwl + 1.0/ηe).^(-1)
            τii         = 0.
            r0          = 0.
            for iter=1:30 # Newton iteration
                τii        = 2*ηve*ε̇ii
                ε̇pwl       = Cpwl*τii^npwl
                r          = ε̇ii - τii/2/ηe - ε̇pwl
                if iter==1 r0 = r end
                ∂r∂ηve     = - ε̇ii/ηe - ε̇pwl*npwl/ηve
                ηve       -= r/∂r∂ηve
                @show (iter, abs(r)/abs(r0), r, r0)
                if (abs(r)/abs(r0)<1e-9 || r<1e-9) break; end
            end
            @show ηve, ηpwl
            # Effective visco-elastic tensor
            D           = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/ani_fac;]
            # Dev. stress tensor
            τ           = D*ε̇
            # ---> Anisotropic Dev. stress invariant in pricipal coordinates: Independent on orientation (objective) !!!!!!
            Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac^2
            τii = sqrt(Y2); τxxt = τ[1]; τyyt = τ[2]; τzzt = -τxxt-τyyt; τxyt = τ[3]

            # PLASTICITY
            γ̇ = 0.
            F = τii - cos(fric)*C - P*sin(fric)
            if F>0
                γ̇    = F/(ηve + ηvp + K*Δt*sin(dil)*sin(fric))
                τiic = τii  - γ̇*ηve
                Pc   = P    + K*Δt*sin(dil)*γ̇ 
                F    = τiic - cos(fric)*C - Pc*sin(fric) - ηvp*γ̇ 
                @show (F)
                Y2   = τiic^2
                ηvep = τiic/2/ε̇ii
                ηve  = ηvep
                D    = 2*ηve*[1 0 0; 0 1 0; 0 0 1.0/ani_fac;]
                τ    = D*ε̇ 
                Y2   = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*ani_fac^2
                τiic = sqrt(Y2)
                F    = τiic - cos(fric)*C - Pc*sin(fric) - ηvp*γ̇
                @show (F)
            end
            τxx_rot[i]  = τ[1]
            τxy_rot[i]  = τ[3]
            τii_rot1[i] = sqrt(Y2)
            ε̇ii_rot[i]  = sqrt(I2)
            η_rot[i]    = ηve

            # Check strain rate components
            τxx         = τ[1]
            τyy         = τ[2]
            τxy         = τ[3]
            ε̇xxd_v = Cpwl*τii^((npwl-1))*τxx
            ε̇xyd_v = Cpwl*τii^((npwl-1))*(τxy*ani_fac)
            ε̇xxd_e = (τxx-τ0_rot[1,1])/2/ηe
            ε̇xyd_e = (τxy-τ0_rot[1,2])/2/ηe*ani_fac
            ε̇xxd_p = γ̇*τxxt/τii/2
            ε̇xyd_p = γ̇*τxyt/τii/2*ani_fac
            fxx    = ε̇_rot[1,1] - ε̇xxd_v - ε̇xxd_e - ε̇xxd_p
            fxy    = ε̇_rot[1,2] - ε̇xyd_v - ε̇xyd_e - ε̇xyd_p
            @show (fxx, fxy)
            # Effective viscosity: compute stress invariant from strain rate invariant
            # ---> Predicts rotated stress invariant !!!!!!
            I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
            τii_rot2[i] = 2*ηve*sqrt(I2)
            # If one does the same by scaling the strain rate invariant by ani_fac^2...
            # ---> Predicts Cartesian stress invariant !!!!!!
            I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2/ani_fac^2
            τii_cart1[i] = 2*ηve*sqrt(I2)
            # Rotate stress back 
            τ_rot       = [τ[1] τ[3]; τ[3] τ[2]]
            τ           = Q'*τ_rot*Q
            # ---> Dependent on orientation (non-objective) !!!!!!
            J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
            τii_cart2[i] = sqrt(J2) 
            τxx_cart[i]  = τ[1,1]
            τxy_cart[i]  = τ[1,2]
            τxx = τ[1,1]
            τyy = τ[2,2]
            τxy = τ[1,2]
            # @show τII_MD7
            @show ηpwl
            # @show τxy_MD7
        end
    end
    p1 = plot(title="Stress invariant", xlabel="θ", ylabel="τᵢᵢ")
    p1 = plot!(θ*180/π, τii_cart1, label="τii_cart1")
    p1 = plot!(θ*180/π, τii_cart2, label="τii_cart2", linewidth=0, marker =:cross)
    p1 = plot!(θ*180/π,  τii_rot1, label="τii_rot1")
    p1 = plot!(θ*180/π,  τii_rot2, label="τii_rot2", linewidth=0, marker =:cross, xlim=(0,180), ylim=(0,4))
    p1 = plot!(θ*180/π,  τii_rot2, label="τii_rot2", linewidth=0, marker =:cross, xlim=(0,180), ylim=(0,4))
    # p2 = scatter!([θ_MD7],  [τII_MD7], label="MD7", marker =:circle)
    p2 = plot(title="Stress components", xlabel="θ", ylabel="τᵢⱼ")
    p2 = plot!(θ*180/π,  τxx_cart, label="τxx_cart",  linewidth=1)
    p2 = plot!(θ*180/π,  τxx_rot,  label="τxx_rot",   linewidth=1)
    p2 = plot!(θ*180/π,  τxy_cart, label="τxy_cart",  linewidth=1)
    p2 = plot!(θ*180/π,  τxy_rot,  label="τxy_rot",   linewidth=1)
    # p2 = scatter!([θ_MD7],  [τxx_MD7], label="MD7", marker =:circle, c=:blue)
    # p2 = scatter!([θ_MD7],  [τxy_MD7], label="MD7", marker =:circle, c=:green)

    display(plot(p1,p2))

end

main_simple_ani_vis()