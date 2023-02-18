using Plots, Printf
import Statistics:mean
include("./fcts_VE_trial.jl")
# MOTIVATED BY COMBINED RHEOLOGY:
# In previous version localiteration used one variable: ηpwl
# Adding creep mechanisms would result in adding new isolated viscosities (ηcst, ηlin...)
# To avoid this, I try to formulate the local iteration in terms of ηve_normal and ηve_shear
visu       = true
plasticity = false
showMD7    = true

function main() 

    # Parameters
    R         = 0*8.314
    εbg       = 1.0
    Δt        = 0.1
    nt        = 10
    My        = 1e6*365.25*24*3600
    PS        = 1
    
    npwl      = 3.0             # POWER LAW EXPONENT
    ani_fac_v = 3.0
    ani_fac_e = 2.0
    ani_fac_p = 1.0
    a1 = 1.0; a2 = 1.0; a3 = 1.0 # friction angle anisotropy
    a_v       = sqrt(ani_fac_v)  # VISCOUS anisotropy strength in Ray's formulation a_v = sqrt(η_n/η_s)
    a_e       = sqrt(ani_fac_e)  # ELASTIC anisotropy strength -------------------- a_e = sqrt(G_n/G_s)
    a_p       =     (ani_fac_p)  # PLASTIC anisotropy strength                      a_p = τxx/τxy
    # aniso_ang = LinRange( 0.0, π/2, 11 )
    aniso_ang = 135*π/180.0
    layer_ang = aniso_ang .- π/2 # in 2D the layers are orthogonal to the director
    
    # Numerics
    noisy   = 4
    tol     = 1.0e-12
    nitmax  = 20
    
    # Ambient conditions
    # Total strain rate
    εxx  = PS*εbg
    εyy  = -εxx*1.0 # Introduce a bit of compaction
    εzz  = -εxx-εyy
    εxy  = (1.0-PS)*εbg
    # Volumetric and deviatoric strain rates
    div  = εxx + εyy
    εxxd = εxx - 1.0/3.0*div
    εyyd = εyy - 1.0/3.0*div
    ezzd = -εxxd - εyyd
    P_ini = 50.0e6
    T     = 1.
    
    # Power law flow params
    τbg    = 2.0
    B      = 2^npwl*εbg/τbg^(npwl)
    τ_chk  = 2*B^(-1.0/npwl)*εbg^(1.0/npwl)
    ηpwl   = B^(-1.0/npwl)*εbg^(1.0/npwl-1)
    τ_chk  = 2*ηpwl*εbg
    if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
    C      = (2*B^(-1/npwl))^(-npwl)
    G      = 1e1
    K      = 7.0e10
    ηe     = G*Δt
    Coh    = 7.0e7  
    fric   = 1*30*π/180.0
    dil    = 1*10*π/180.0
    ηvp    = 1e18
    
    # Storage/visualisation
    τii_Ray = zeros(nt, length(aniso_ang))
    τxx_Ray = zeros(nt, length(aniso_ang))
    τyy_Ray = zeros(nt, length(aniso_ang))
    τxy_Ray = zeros(nt, length(aniso_ang))
    τii_MD7 = zeros(nt, length(aniso_ang))
    τxx_MD7 = zeros(nt, length(aniso_ang))
    τyy_MD7 = zeros(nt, length(aniso_ang))
    τxy_MD7 = zeros(nt, length(aniso_ang))
    p_MD7   = zeros(nt, length(aniso_ang))
    
    # Loop over all orientations
    for i=1:length(aniso_ang)
    
        # Loop over time and model V-E-VP stress build-up 
        τxx = 0.; τyy = 0.; τxy = 0.; P = P_ini;

        for it=1:nt
            if noisy>0 @printf("****** Time step %03d ******\n", it) end
            
            # Swap old deviatoric stresses
            τxx0 = τxx; τyy0 = τyy; τxy0 = τxy; P0 = P

            # Pressure update
            P = P0 - K*Δt*div
                
            ###############
    
            # Transformation matrix: towards principal plane
            Q       = [cos(layer_ang[i]) sin(layer_ang[i]); -sin(layer_ang[i]) cos(layer_ang[i])]
            @show ε_tens  = [εxxd εxy; εxy εyyd]
            τ0_tens = [τxx0 τxy0; τxy0 τyy0]
            @show ε_rot   = Q*ε_tens*Q'
            τ0_rot  = Q*τ0_tens*Q'
            @show cos(layer_ang[i])*sin(layer_ang[i])
            @show nx = cos(aniso_ang[i])
            @show ny = sin(aniso_ang[i])

            # ε        = [εxx; εyy; εxy]
            # τ0       = [τxx0; τyy0; τxy0]
            # ε_eff    = ε .+ inv(Dani_e) * τ0 ./(2*ηe)

            ############### VISCO-ELASTIC TRIAL
            τ_rot, ηve, a_ve, ηpwl = LocalViscoElasticTrialStress_Newton_ηve( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )
            τxx     = τ_rot[1,1]; τxy   = τ_rot[1,2]; τyy   =  τ_rot[2,2]; τzz   = -τxx-τyy
            J2 = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_ve^2
            εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
            εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
            εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2

            εii_eff_ve = sqrt( 1/2*(εxx_eff^2+εyy_eff^2) + εxy_eff^2/a_ve^2 )
            @printf("VE:  Check J2 definition: %2.2e\n", (J2 - (2*ηve)^2 * ( εii_eff_ve^2 ))/J2*100 )
            @printf("VE:  Check T2 definition: %2.2e\n", (sqrt(J2) - 2*ηve * εii_eff_ve  ) / sqrt(J2)*100 )
            @printf("VE:  sqrt(J2) = %2.2e\n", sqrt(J2))
            @printf("VE:  a_ve     = %2.2e\n", a_ve)
            @printf("VE:  ani_ve   = %2.2e\n", a_ve^2)
            ηvep    = ηve
            a_vep   = a_ve
            ############### PLASTIC CORRECTION
            if plasticity == true
                # Check yield condition
                J2      = 0.5*(τ_rot[1,1]^2 + τ_rot[2,2]^2) + τ_rot[1,2]^2*a_p^2
                J2t     = J2
                J2_corr = J2
                F       = sqrt(J2) - Coh - P*sin(fric)/3*(a1+a2+a3) +  sin(fric)/3*( a1*τxx + a2*τyy + a3*τzz)
                # Initial guess 
                γ̇  = 0.
                Fc    = 0.
                if F>0 
                    for iter=1:nitmax
                        # It would be nice to have this section derived with Symbolics.jl
                        eta_ve = ηve
                        Txx    = τxx
                        Tyy    = τyy
                        Tzz    = -τyy-τxx
                        Txy    = τxy
                        Tyx    = τxy
                        gdot   = γ̇
                        ap_n   = (1 - eta_ve .* gdot ./ sqrt(J2))
                        ap_s   = (1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(J2) .* a_ve .^ 2))
                        Jii    = J2
                        dt     = Δt
                        eta_vp = ηvp
                        J1     = a1*Txx + a2*Tyy + a3*Tzz
                        # J2_corr = 1.0 * Txx .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2 + 1.0 * Txy .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + 1.0 * Tyx .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + 1.0 * Tyy .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2;
                        # Fc = -Coh - P .* (a1 + a2 + a3) .* sin(fric) / 3 + (Txx .* a1 + Tyy .* a2 + a3 .* (-Txx - Tyy)) .* sin(fric) / 3 + 1.0 * sqrt(Txx .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2 + Txy .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + Tyx .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + Tyy .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2);
                        # dFdgdot = 1.0 * (-Txx .^ 2 .* eta_ve .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) - Txy .^ 2 .* a_p .^ 4 .* eta_ve .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) - Tyx .^ 2 .* a_p .^ 4 .* eta_ve .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) - Tyy .^ 2 .* eta_ve .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) ./ sqrt(Txx .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2 + Txy .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + Tyx .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2)) + 0.707106781186547) .^ 2 + Tyy .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 + Txy .^ 2 .* a_p .^ 2 + Tyx .^ 2 .* a_p .^ 2 + Tyy .^ 2) + 0.707106781186547) .^ 2);
                        J2_corr = Txx .^ 2 .* (1 - eta_ve .* gdot ./ sqrt(J2)) .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 .* (1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(J2) .* a_ve .^ 2)) .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 .* (1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(J2) .* a_ve .^ 2)) .^ 2 / 2 + Tyy .^ 2 .* (1 - eta_ve .* gdot ./ sqrt(J2)) .^ 2 / 2;
                        Fc = -Coh + sqrt(J2_corr) - P .* (a1 + a2 + a3) .* sin(fric) / 3 + (Txx .* a1 + Tyy .* a2 + a3 .* (-Txx - Tyy)) .* sin(fric) / 3;
                        dFdgdot = (-Txx .^ 2 .* eta_ve .* (1 - eta_ve .* gdot ./ sqrt(J2)) ./ (2 * sqrt(J2)) - Txy .^ 2 .* a_p .^ 4 .* eta_ve .* (1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(J2) .* a_ve .^ 2)) ./ (2 * sqrt(J2) .* a_ve .^ 2) - Tyx .^ 2 .* a_p .^ 4 .* eta_ve .* (1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(J2) .* a_ve .^ 2)) ./ (2 * sqrt(J2) .* a_ve .^ 2) - Tyy .^ 2 .* eta_ve .* (1 - eta_ve .* gdot ./ sqrt(J2)) ./ (2 * sqrt(J2))) ./ sqrt(J2_corr);
                        
                        J2_corr = Txx .^ 2 .* ap_n .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 .* ap_s .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 .* ap_s .^ 2 / 2 + Tyy .^ 2 .* ap_n .^ 2 / 2;
                        Fc = -Coh + sqrt(J2_corr) - P .* (a1 + a2 + a3) .* sin(fric) / 3 + (Txx .* a1 + Tyy .* a2 + a3 .* (-Txx - Tyy)) .* sin(fric) / 3;
                        dFdgdot = (-Txx .^ 2 .* ap_n .* eta_ve ./ (2 * sqrt(Jii)) - Txy .^ 2 .* a_p .^ 4 .* ap_s .* eta_ve ./ (2 * sqrt(Jii) .* a_ve .^ 2) - Tyx .^ 2 .* a_p .^ 4 .* ap_s .* eta_ve ./ (2 * sqrt(Jii) .* a_ve .^ 2) - Tyy .^ 2 .* ap_n .* eta_ve ./ (2 * sqrt(Jii))) ./ sqrt(J2_corr);
                        
                        # With out of plane
                        Pc   = P + K*Δt*γ̇*sin(dil)/3*(a1+a2+a3)

                        ap_n = 1 - eta_ve .* gdot ./ sqrt(Jii);

                        ap_s = 1 - a_p .^ 2 .* eta_ve .* gdot ./ (sqrt(Jii) .* a_ve .^ 2);
                        
                        J1_corr = Txx .* a1 .* ap_n + Tyy .* a2 .* ap_n + a3 .* (-Txx .* ap_n - Tyy .* ap_n);
                        
                        J2_corr = Txx .^ 2 .* ap_n .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 .* ap_s .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 .* ap_s .^ 2 / 2 + Tyy .^ 2 .* ap_n .^ 2 / 2 + Tzz .^ 2 .* ap_n .^ 2 / 2;
                        
                        Fc = -Coh + J1_corr .* sin(fric) / 3 + sqrt(J2_corr) - Pc .* (a1 + a2 + a3) .* sin(fric) / 3 - eta_vp .* gdot;
                        
                        dFdgdot = -K .* dt .* (a1 + a2 + a3) .^ 2 .* sin(dil) .* sin(fric) / 9 - eta_vp + (a3 .* (Txx .* eta_ve ./ sqrt(Jii) + Tyy .* eta_ve ./ sqrt(Jii)) - Txx .* a1 .* eta_ve ./ sqrt(Jii) - Tyy .* a2 .* eta_ve ./ sqrt(Jii)) .* sin(fric) / 3 + (-Txx .^ 2 .* ap_n .* eta_ve ./ (2 * sqrt(Jii)) - Txy .^ 2 .* a_p .^ 4 .* ap_s .* eta_ve ./ (2 * sqrt(Jii) .* a_ve .^ 2) - Tyx .^ 2 .* a_p .^ 4 .* ap_s .* eta_ve ./ (2 * sqrt(Jii) .* a_ve .^ 2) - Tyy .^ 2 .* ap_n .* eta_ve ./ (2 * sqrt(Jii)) - Tzz .^ 2 .* ap_n .* eta_ve ./ (2 * sqrt(Jii))) ./ sqrt(J2_corr);
                        @printf("J1=%2.4e J2=%2.4e\n", J1_corr, J2_corr)

                        εxx_p      = γ̇*( (a1/3-a3/3)*sin(fric) + τxx/sqrt(J2)/2)
                        txx1 = 2*ηve*(εxx_eff-εxx_p)
                        txx2 = sqrt(Txx .^ 2 .* ap_n .^ 2)

                        @printf("%2.4e %2.4e\n", txx1, txx2)
                        # @show (J2_corr1-J2_corr)/J2_corr
                        # @show (Fc-Fc1)/Fc
                        # @show (dFdgdot1-dFdgdot)/dFdgdot

                        dFdγ̇ = dFdgdot
                        # It would be nice to have this section derived with Symbolics.jl
                        if noisy>2 @printf("Newton VP It. %02d, r abs. = %2.2e r rel. = %2.2e\n", iter, Fc, Fc/F) end
                        if abs(Fc)/F < tol || abs(Fc)<tol
                            # error("sto")
                            break
                        end
                        γ̇    -= Fc/dFdγ̇
                    end
                    # Post-process plasticity
                    if abs(Fc)/F>1e-6 Fc, error("Yield failed 3")  end
                    εxx_p      = γ̇*( (a1/3-a3/3)*sin(fric) + τxx/sqrt(J2)/2)
                    εyy_p      = γ̇*( (a2/3-a3/3)*sin(fric) + τyy/sqrt(J2)/2)
                    εxy_p      = γ̇*τxy/sqrt(J2)/2*a_p^2
                    εyx_p      = γ̇*τxy/sqrt(J2)/2*a_p^2
                    D          = 2*ηve * [1 0 0 0; 0 1 0 0; 0 0 a_ve^-2 0; 0 0 0 a_ve^-2]
                    τ_vec      = D*([εxx_eff; εyy_eff; εxy_eff; εxy_eff] .- [εxx_p; εyy_p; εxy_p; εyx_p])
                    τ_rot      = [τ_vec[1] τ_vec[3]; τ_vec[4] τ_vec[2]]
                    J2_corr_p  = 0.5*(τ_rot[1,1]^2 + τ_rot[2,2]^2) + τ_rot[1,2]^2*a_p^2
                    @printf("VEP: sqrt(J2_corr_p ) = %2.2e\n", sqrt(J2_corr_p) )
                    @printf("VEP: sqrt(J2_corr   ) = %2.2e\n", sqrt(J2_corr)   )
                    # Compute a_vep 
                    axx       = (sqrt(J2t)        -       ηve*γ̇)/(sqrt(J2t))
                    axy       = (sqrt(J2t)*a_ve^2 - a_p^2*ηve*γ̇)/(sqrt(J2t)*a_ve^2)
                    η_mat     = ηve * [axx 0 0 0; 0 axx 0 0; 0 0 axy 0; 0 0 0 axy]
                    D         = 2*η_mat 
                    τ_vec     = D*([εxx_eff; εyy_eff; εxy_eff; εxy_eff]) 
                    # --- 
                    a_vep     = a_ve*sqrt(axx/axy)
                    ηvep      = ηve*axx
                    η_mat     = ηvep * [1 0 0 0; 0 1 0 0; 0 0 a_vep^-2 0; 0 0 0 a_vep^-2] 
                    D         = 2*η_mat
                    τ_vec     = D*([εxx_eff; εyy_eff; εxy_eff; εxy_eff])
                    # ---
                    τ_rot     = [τ_vec[1] τ_vec[3]; τ_vec[4] τ_vec[2]]
                    τxxc      = τ_rot[1,1]; τyyc = τ_rot[2,2]; τzzc = -τxxc-τyyc
                    J2_corr_p = 0.5*(τ_rot[1,1]^2 + τ_rot[2,2]^2) + τ_rot[1,2]^2*a_p^2
                    J1_corr_p = a1*τxxc + a2*τyyc + a3*τzzc
                    Pc        = P + K*Δt*γ̇*sin(dil)/3*(a1+a2+a3)
                    Fchk      = sqrt(J2_corr_p) - Coh - Pc*sin(fric)/3*(a1+a2+a3) +  sin(fric)/3*J1_corr_p - ηvp*γ̇
                end
            end
    
            ##############################
            # Rotate back
            # τ   = τ_rot
            τ   = Q'*τ_rot*Q
            τii = sqrt( 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2)
            τxx = τ[1,1]; τyy = τ[2,2]; τxy = τ[1,2]
            @printf("VEP: τii      = %2.6e\n", τii)
    
            ###############
            # Store
            τii_Ray[it,i] = τii
            τxx_Ray[it,i] = τxx
            τyy_Ray[it,i] = τyy
            τxy_Ray[it,i] = τxy
    
            # ###################################################################
            # THIS IS WHAT WILL HAPPEN IN MDOODZ: It still uses the C_ANI business we derived earlier...
            two      = 2.
            nx       = cos(aniso_ang[i]); ny = sin(aniso_ang[i])
            N        = [nx; ny]
            d0       = 2*N[1]^2*N[2]^2
            d1       = N[1]*N[2] * (-N[1]^2 + N[2]^2)
            C_ANI    = [-d0 d0 -two*d1; d0 -d0 two*d1; -d1 d1 -two*(1/2-d0)]      ###### !!! - sign in D33 -a0
            C_ISO    = [1 0 0; 0 1 0; 0 0 two*1//2]
            ani      = 1. - 1.  ./ (a_vep^2)
            Dani_vep = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
            ani      = 1. - 1.  ./ (a_e^2)
            Dani_e   = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
            ani      = 1. - 1.  ./ (a_v^2)
            Dani_v   = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
  
            # 1 ---------------------------------------------------------------
            ε        = [εxx; εyy; εxy]
            τ0       = [τxx0; τyy0; τxy0]
            ε_eff    = ε .+ inv(Dani_e) * τ0 ./(2*ηe)
            τ_MD7_v1 = 2*ηvep*Dani_vep*ε_eff
         
            # 3 ---------------------------------------------------------------
            # same as 2) but without mistake. Doesn't look very sexy... Can maybe be simplified...
            E     = ηe*[1 0 0; 0 1 0; 0 0 1] .+ ηpwl*Dani_v*inv(Dani_e)
            τ_MD7_v2 = 2*ηvep*(ηpwl+ηe)*inv(E)*(Dani_v*ε + Dani_v*inv(Dani_e) * (τ0 ./(2*ηe)))
            
            # @printf("max. rel. VE tangent : %2.2e\n ", maximum( (τ_MD7_v1.-τ_MD7_v2)./τ_MD7_v1)*100 ) 
                
            τii_MD7[it,i] = sqrt(0.5*(τ_MD7_v1[1]^2 + τ_MD7_v1[2]^2) + τ_MD7_v1[3]^2)
            @printf("VEP: τii      = %2.6e\n", τii_MD7[it,i])
            τxx_MD7[it,i] = τ_MD7_v1[1]
            τyy_MD7[it,i] = τ_MD7_v1[2]
            τxy_MD7[it,i] = τ_MD7_v1[3]
            p_MD7[it,i]   = P
    
        end
    end
    
    if visu
        if length(aniso_ang)==1
        p=plot([1:nt].*(Δt/My), τii_Ray[:,1], linewidth=1, marker =:circle, color=:green,  label=:none)
        p=plot!([1:nt].*(Δt/My), τii_MD7[:,1], linewidth=0, marker =:xcross, color=:green,  label=:none)
        else
        st1 = Int64(floor(  nt/3))
        st2 = Int64(floor(2*nt/3))
        st3 = nt
        p1=plot( aniso_ang.*(180/π), τii_Ray[st1,:], color=:green, label="1/3 loading")
        p1=plot!(aniso_ang.*(180/π), τii_Ray[st2,:], color=:blue,  label="2/3 loading")
        p1=plot!(aniso_ang.*(180/π), τii_Ray[st3,:], color=:red,   label="3/3 loading")
        if showMD7
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st1,:], marker =:circle, linewidth=0, color=:green, label=:none)
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st2,:], marker =:circle, linewidth=0, color=:blue,  label=:none)
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st3,:], marker =:circle, linewidth=0, color=:red,   label=:none)
        end
        p1=plot!(xlabel="n angle [°]", ylabel="τII [MPa]")
        
        i1 = 1  # angle 1
        i2 = Int64(floor(  length(aniso_ang)*1/2))  # angle 2
        i3 = Int64(floor(  length(aniso_ang)*2/3))  # angle 2
        p2=plot( [1:nt].*(Δt/My), τii_Ray[:,i1], color=:red,   label=@sprintf("n angle %2.2f°", aniso_ang[i1]*180/π) )
        p2=plot!([1:nt].*(Δt/My), τii_Ray[:,i2], color=:blue,  label=@sprintf("n angle %2.2f°", aniso_ang[i2]*180/π) )
        p2=plot!([1:nt].*(Δt/My), τii_Ray[:,i3], color=:green, label=@sprintf("n angle %2.2f°", aniso_ang[i3]*180/π) )
        if showMD7
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i1], linewidth=0, marker =:circle, color=:red,    label=:none)
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i2], linewidth=0, marker =:circle, color=:blue,   label=:none)
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i3], linewidth=0, marker =:circle, color=:green,  label=:none)
        end
        p2=plot!(xlabel="t [My]", ylabel="τII [MPa]")

        p4=plot(  aniso_ang.*(180/π), τxx_Ray[st3,:], color=:blue,  label="τxx Ray")
        p4=plot!( aniso_ang.*(180/π), τyy_Ray[st3,:], color=:green, label="τyy Ray")
        p4=plot!( aniso_ang.*(180/π), τxy_Ray[st3,:], color=:red,   label="τxy Ray")
        if showMD7
        p4=plot!( aniso_ang.*(180/π), τxx_MD7[st3,:], linewidth=0, marker =:circle, color=:blue,  label="τxx MD7")
        p4=plot!( aniso_ang.*(180/π), τyy_MD7[st3,:], linewidth=0, marker =:circle, color=:green, label="τyy MD7")
        p4=plot!( aniso_ang.*(180/π), τxy_MD7[st3,:], linewidth=0, marker =:circle, color=:red,   label="τxy MD7")
        end
        display(plot(p1,p2,p4, size=(1200,900), legend=:none))
        end
    end
    
end

main()  
