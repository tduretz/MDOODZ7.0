####################################################################

# This function only computes eta_pwl, it would require modification more mechanisms are added (Peierls)
function LocalViscoElasticTrialStress_Newton_ηve_ani_fac( ε_rot, τ0_rot, ηe, ani_fac_e, ani_fac_v, B, C, npwl, tol, nitmax, noisy )

    εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*ani_fac_e 
    εzz_eff    = -εxx_eff -εyy_eff
    # Construct initial guess for viscosity
    ε2         = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/ani_fac_v 
    ηpwl       = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)
    ηve_n      = 1/(1/ηpwl + 1/ηe) 
    ηve_s      = 1/(ani_fac_v/ηpwl + ani_fac_e/ηe) 
    fxx0       = 0.
    fxy0       = 0.
    for iter=1:nitmax
        
        τxx       = 2*ηve_n * εxx_eff
        τyy       = 2*ηve_n * εyy_eff
        τxy       = 2*ηve_s * εxy_eff
        τzz       = -τxx - τyy 

        Y2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*ani_fac_v 
        εii_pwl   = C * sqrt(Y2_v) ^(npwl)
        ηpwl      = B^(-1/npwl) * εii_pwl^((1-npwl)/npwl)

        # Residuals
        fxx       = 1 - ηve_n/ηe         - ηve_n/ηpwl               
        fxy       = 1 - ηve_s/(ηe/ani_fac_e) - ηve_s/(ηpwl/ani_fac_v) 
        Wxx       = 2*(εxx_eff*τxx + εyy_eff*τyy + εzz_eff*τzz)
        Wxy       = 2*εxy_eff*τxy
        ieta_pwl  = 1 ./ ηpwl

        if iter==1 fxx0 = fxx  end
        if iter==1 fxy0 = fxy  end
        @printf("Newton ani VE It. %02d: fxx = %2.2e --- fxy = %2.2e\n", iter,  fxx, fxy)
        if (abs(fxx/fxx0)<tol || abs(fxx)<tol) && (abs(fxy/fxy0)<tol || abs(fxy)<tol)
            break
        end

        # Jacobian, check: AnisotropicVE_William_v2.ipynb
        eta_ve_n  = ηve_n
        eta_ve_s  = ηve_s
        eta_e     = ηe
        dfxxde_n = Wxx .* eta_ve_n       .* ieta_pwl .* (1 - npwl) ./ (2 * Y2_v) - ieta_pwl - 1 ./ eta_e;
        dfxxde_s = Wxy .* ani_fac_v      .* eta_ve_n .* ieta_pwl .* (1 - npwl) ./ Y2_v;
        dfxyde_n = Wxx .* ani_fac_v      .* eta_ve_s .* ieta_pwl .* (1 - npwl) ./ (2 * Y2_v);
        dfxyde_s = Wxy .* ani_fac_v .^ 2 .* eta_ve_s .* ieta_pwl .* (1 - npwl) ./ Y2_v - ani_fac_e ./ eta_e - ani_fac_v .* ieta_pwl;

        f      = [fxx; fxy]
        J      = [dfxxde_n dfxxde_s; dfxyde_n dfxyde_s]
        dx     = -J\f
        ηve_n += dx[1]
        ηve_s += dx[2]
    end

    ηve  = ηve_n
    ani_fac_ve = ηve_n/ηve_s

    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 1.0/ani_fac_ve]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, ani_fac_ve, ηpwl

end


# This function only computes eta_pwl, it would require modification more mechanisms are added (Peierls)
function LocalViscoElasticTrialStress_Newton_ηve( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )

    εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2 
    εzz_eff    = -εxx_eff -εyy_eff
    # Construct initial guess for viscosity
    ε2         = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/a_v^2 
    ηpwl       = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)
    ηve_n      = 1/(1/ηpwl + 1/ηe) 
    ηve_s      = 1/(a_v^2/ηpwl + a_e^2/ηe) 
    fxx0       = 0.
    fxy0       = 0.
    for iter=1:nitmax
        
        τxx       = 2*ηve_n * εxx_eff
        τyy       = 2*ηve_n * εyy_eff
        τxy       = 2*ηve_s * εxy_eff
        τzz       = -τxx - τyy 

        Y2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_v^2 
        εii_pwl   = C * sqrt(Y2_v) ^(npwl)
        ηpwl      = B^(-1/npwl) * εii_pwl^((1-npwl)/npwl)

        # Residuals
        fxx       = 1 - ηve_n/ηe         - ηve_n/ηpwl               
        fxy       = 1 - ηve_s/(ηe/a_e^2) - ηve_s/(ηpwl/a_v^2) 
        Wxx       = 2*(εxx_eff*τxx + εyy_eff*τyy + εzz_eff*τzz)
        Wxy       = 2*εxy_eff*τxy
        ieta_pwl  = 1 ./ ηpwl

        if iter==1 fxx0 = fxx  end
        if iter==1 fxy0 = fxy  end
        @printf("Newton ani VE It. %02d: fxx = %2.2e --- fxy = %2.2e\n", iter,  fxx, fxy)
        if (abs(fxx/fxx0)<tol || abs(fxx)<tol) && (abs(fxy/fxy0)<tol || abs(fxy)<tol)
            break
        end

        # Jacobian, check: AnisotropicVE_William_v2.ipynb
        eta_ve_n  = ηve_n
        eta_ve_s  = ηve_s
        eta_e     = ηe
        dfxxde_n = Wxx .* eta_ve_n .* ieta_pwl .* (1 - npwl) ./ (2 * Y2_v) - ieta_pwl - 1 ./ eta_e;
        dfxxde_s = Wxy .* a_v .^ 2 .* eta_ve_n .* ieta_pwl .* (1 - npwl) ./ Y2_v;
        dfxyde_n = Wxx .* a_v .^ 2 .* eta_ve_s .* ieta_pwl .* (1 - npwl) ./ (2 * Y2_v);
        dfxyde_s = Wxy .* a_v .^ 4 .* eta_ve_s .* ieta_pwl .* (1 - npwl) ./ Y2_v - a_e .^ 2 ./ eta_e - a_v .^ 2 .* ieta_pwl;

        f      = [fxx; fxy]
        J      = [dfxxde_n dfxxde_s; dfxyde_n dfxyde_s]
        dx     = -J\f
        ηve_n += dx[1]
        ηve_s += dx[2]
    end

    ηve  = ηve_n
    a_ve = sqrt(ηve_n/ηve_s)

    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 a_ve^-2]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, a_ve, ηpwl

end

####################################################################

# This function only computes eta_pwl, it would require modification more mechanisms are added (Peierls)
function LocalViscoElasticTrialStress_NewtonDoppelBesser( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )

    εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2 
    # Construct initial guess for viscosity
    ε2         = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/a_v^2 
    ηpwl       = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)
    f0         = 0.
    a_ve       = 0.
    ηve        = 0.
    for iter=1:nitmax
        
        ηve       = (1/ηe + 1/ηpwl)^(-1)
        a_ve      = sqrt((ηe*a_v^2 + ηpwl*a_e^2)/(ηe + ηpwl))     
        τxx       = 2*ηve * εxx_eff
        τyy       = 2*ηve * εyy_eff
        τxy       = 2*ηve * εxy_eff / a_ve^2
        τzz       = -τxx - τyy 

        τ2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_v^2 
        εii_pwl   = C * sqrt(τ2_v) ^(npwl)
        ηpwl1     = B^(-1/npwl) * εii_pwl^((1-npwl)/npwl)

        f      = 1 - ηpwl1/ηpwl

        if iter==1 f0 = f;  end
        @printf("Newton ani VE It. %02d: fη_pwl = %2.2e\n", iter,  f)
        if (abs(f/f0)<tol || abs(f)<tol) 
            break
        end

        Eii_pwl   = εii_pwl  # DO YOU MANAGE TO AVOID THIS WITH SYMBOLICS.jl?
        Exx_eff   = εxx_eff
        Eyy_eff   = εyy_eff
        Exy_eff   = εxy_eff
        eta_ve    = ηve
        eta_e     = ηe
        eta_pwl   = ηpwl
        Txx       = τxx
        Tyy       = τyy
        Txy       = τxy
        dfde = -(1 - npwl) .* (0.5 * Exx_eff .^ 2 .* eta_ve .^ 3 ./ eta_pwl .^ 2 - Exy_eff .^ 2 .* a_e .^ 2 .* a_v .^ 2 .* eta_ve .^ 2 .* (eta_e + eta_pwl) .^ 2 ./ (a_e .^ 2 .* eta_pwl + a_v .^ 2 .* eta_e) .^ 3 + Exy_eff .^ 2 .* a_v .^ 2 .* eta_ve .^ 2 .* (eta_e + eta_pwl) ./ (a_e .^ 2 .* eta_pwl + a_v .^ 2 .* eta_e) .^ 2 + Exy_eff .^ 2 .* a_v .^ 2 .* eta_ve .^ 3 .* (eta_e + eta_pwl) .^ 2 ./ (eta_pwl .^ 2 .* (a_e .^ 2 .* eta_pwl + a_v .^ 2 .* eta_e) .^ 2) + 0.5 * Eyy_eff .^ 2 .* eta_ve .^ 3 ./ eta_pwl .^ 2 + 0.0625 * (-Txx - Tyy) .* (-2 * Txx .* eta_ve ./ eta_pwl .^ 2 - 2 * Tyy .* eta_ve ./ eta_pwl .^ 2)) ./ (0.5 * Exx_eff .^ 2 .* eta_ve .^ 2 + Exy_eff .^ 2 .* a_v .^ 2 .* eta_ve .^ 2 .* (eta_e + eta_pwl) .^ 2 ./ (a_e .^ 2 .* eta_pwl + a_v .^ 2 .* eta_e) .^ 2 + 0.5 * Eyy_eff .^ 2 .* eta_ve .^ 2 + 0.125 * (-Txx - Tyy) .^ 2) + 1 ./ eta_pwl;

        ηpwl -= f/dfde

    end


    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 a_ve^-2]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, a_ve, ηpwl

end

####################################################################

# As far as we are, we needed two variables (eta_ve, a_ve) for the Newton iteration
# DO YOU MANAGE TO REDUCE IT TO A ONE VARIABLE PROCEDURE?
function LocalViscoElasticTrialStress_NewtonBesser( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )

    εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2 
    # Construct initial guess for viscosity
    ε2         = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/a_v^2 
    ηpwl       = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)
    ηve        = (1/ηe + 1/ηpwl)^(-1)
    # Construct initial guess for a_ve
    a_ve         = a_e # 1/(1/a_v+1/a_e)
    feta0, fave0 = 0., 0.
    off = 1.0
    for iter=1:nitmax
        
        τxx       = 2*ηve * εxx_eff
        τyy       = 2*ηve * εyy_eff
        τxy       = 2*ηve * εxy_eff / a_ve^2
        τzz       = -τxx - τyy 

        τ2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_v^2 
        εii_pwl   = C * sqrt(τ2_v) ^(npwl)
        ηpwl      = B^(-1/npwl) * εii_pwl^((1-npwl)/npwl)

        feta      = ηve - (    1/ηe +     1/ηpwl)^(-1)
        fave      = a_ve - sqrt((ηe*a_v^2 + ηpwl*a_e^2)/(ηe + ηpwl)) 

        if iter==1 feta0 = feta; fave0 = fave; end
        @printf("Newton ani VE It. %02d: fη = %2.2e, fave = %2.2e\n", iter,  feta, fave)
        if (abs(feta/feta0)<tol || abs(feta)<tol) &&  (abs(fave/fave0)<tol || abs(fave)<tol)
            break
        end

        Eii_pwl   = εii_pwl  # DO YOU MANAGE TO AVOID THIS WITH SYMBOLICS.jl?
        Exx_eff   = εxx_eff
        Eyy_eff   = εyy_eff
        Exy_eff   = εxy_eff
        eta_ve    = ηve
        eta_e     = ηe
        eta_pwl   = ηpwl
        Txx       = τxx
        Tyy       = τyy
        Txy       = τxy

        dfetadeta = -B  .^(1   ./npwl)  .*Eii_pwl .^((npwl - 1)./npwl) .*(1 - npwl).*(0.5*Exx_eff.^2 .*eta_ve + Exy_eff.^2 .*a_v.^2 .*eta_ve./a_ve.^4 + 0.5*Eyy_eff.^2 .*eta_ve + 0.0625*(-4*Exx_eff - 4*Eyy_eff).*(-Txx - Tyy))./((B.^(1 ./npwl).*Eii_pwl.^((npwl - 1)./npwl) + 1 ./eta_e).^2 .*(0.5*Exx_eff.^2 .*eta_ve.^2 + Exy_eff.^2 .*a_v.^2 .*eta_ve.^2 ./a_ve.^4 + 0.5*Eyy_eff.^2 .*eta_ve.^2 + 0.125*(-Txx - Tyy).^2)) + 1;
        dfetadave = B.^(1 ./npwl).*Eii_pwl.^(-(1 - npwl)./npwl).*Exy_eff.*Txy.*a_v.^2 .*eta_ve.*(1 - npwl)./(a_ve.^3 .*(B.^(1 ./npwl).*Eii_pwl.^((npwl - 1)./npwl) + 1 ./eta_e).^2 .*(0.5*Exx_eff.^2 .*eta_ve.^2 + Exy_eff.^2 .*a_v.^2 .*eta_ve.^2 ./a_ve.^4 + 0.5*Eyy_eff.^2 .*eta_ve.^2 + 0.125*(-Txx - Tyy).^2));
        dfavedeta = -sqrt((a_e.^2 .*eta_pwl + a_v.^2 .*eta_e)./(eta_e + eta_pwl)).*(eta_e + eta_pwl).*(a_e.^2 .*eta_pwl.^2 ./(2*eta_ve.^2 .*(eta_e + eta_pwl)) + eta_pwl.^2 .*(-a_e.^2 .*eta_pwl - a_v.^2 .*eta_e)./(2*eta_ve.^2 .*(eta_e + eta_pwl).^2))./(a_e.^2 .*eta_pwl + a_v.^2 .*eta_e);
        dfavedave = 1;

        f     = [feta; fave]
        J     = [dfetadeta off*dfetadave; off*dfavedeta dfavedave]
        dx    = -J\f
        ηve  += dx[1]
        a_ve += dx[2]
    end

    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 a_ve^-2]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, a_ve, ηpwl

end


####################################################################

# Just works!
function LocalViscoElasticTrialStress_Picard( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )

    εxx_eff    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2 
    ε2         = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/a_v^2 # ACHTUNG: a^2 moves to denomina
    ηpwl       = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)
    ηve        = (1/ηe + 1/ηpwl)^(-1)
    ηve0       = ηve
    τ2         = 1.0

    @show ε_rot
    @show εxx_eff
    @show εyy_eff
    @show εxy_eff

    τxx, τyy, τxy, τzz = 0., 0., 0., 0.
    εii_pwl = 0.
    for iter=1:nitmax
        a_ve      = sqrt((ηe*a_v^2 + ηpwl*a_e^2)/(ηe + ηpwl)) 
        τxx       = 2*ηve * εxx_eff
        τyy       = 2*ηve * εyy_eff
        τxy       = 2*ηve * εxy_eff / a_ve^2
        τzz       = -τxx - τyy 

        τ2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_v^2 # ACHTUNG: a^2 moves to denomina
        εii_pwl   = C * sqrt(τ2_v) ^(npwl)
        ηpwl      = B^(-1/npwl) * εii_pwl^((1-npwl)/npwl)

        ηve       = (    1/ηe +     1/ηpwl)^(-1)
        err_norm  = abs(ηve-ηve0)/ηve0
        @printf("Picard VE It. %02d: Δηve = %2.2e\n", iter, err_norm)
        if err_norm<tol/100 
            break
        end
        ηve0 = ηve
    end

    a_ve    = sqrt((ηe*a_v^2 + ηpwl*a_e^2)/(ηe + ηpwl)) 
    @printf("ηve = %2.2e \n", ηve )
    @printf("a_ve PICARD = %2.2e\n", a_ve)

    ε2_eff  = 0.5*(εxx_eff^2    + εyy_eff^2   ) + εxy_eff^2   /a_ve^2
    τ2_ve   = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_ve^2
    τii     = 2*ηve*sqrt(ε2_eff)
    @printf("Δτii = %2.2e\n", (τii-sqrt(τ2_ve))/τii*100)
    # @printf("τii = %2.2e\n", τii)

    r_η = sqrt(ε2_eff) - τii/(2*ηe) - εii_pwl 
    @printf("%2.2e\n", r_η/sqrt(ε2_eff))

    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 a_ve^-2]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, a_ve, ηpwl

end

####################################################################

# Looked good but does not work when a_e != a_v
function LocalViscoElasticTrialStress_Newton( ε_rot, τ0_rot, ηe, a_e, a_v, B, C, npwl, tol, nitmax, noisy )

    a_ve    = a_v # trial value
    εxx_eff = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
    εyy_eff = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
    εxy_eff = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*a_e^2                 # Checked with python
    ε2      = 0.5*(ε_rot[1,1]^2 + ε_rot[2,2]^2) + ε_rot[1,2]^2/a_ve^2 # ACHTUNG: a^2 moves to denominator!!!
    ε2_eff  = 0.5*(εxx_eff^2    + εyy_eff^2   ) + εxy_eff^2   /a_ve^2 # ACHTUNG: a^2 moves to denominator!!!

    # Isolated viscosities
    ηpwl = B^(-1/npwl) * ε2^(0.5*(1-npwl)/npwl)

    # Initial guess
    ηve  = 1/(1/ηe + 1/ηpwl)
    res0 = 1.0

    for iter=1:nitmax

        # Stress vector
        τii     = 2*ηve*sqrt(ε2_eff)

        
        
        τxx       = 2*ηve * εxx_eff
        τyy       = 2*ηve * εyy_eff
        τxy       = 2*ηve * εxy_eff / a_ve^2
        τzz       = -τxx - τyy 
        τ2_v      = 0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2*a_v^2 # ACHTUNG: a^2 moves to denomina
        
        
        
        εii_pwl = C*sqrt(τ2_v)^(npwl)
        ηpwl    = 2^((1-npwl)/npwl)*B^(-1/npwl) * εii_pwl^(0.5*(1-npwl)/npwl)
        a_ve    = sqrt((ηe*a_v^2 + ηpwl*a_e^2)/(ηe + ηpwl))  
        ε2_eff  = 0.5*(εxx_eff^2    + εyy_eff^2   ) + εxy_eff^2   /a_ve^2 # ACHTUNG: a^2 moves to denominator!!!
    
        # Additive strain decomposition: Scalar
        r_η = sqrt(ε2_eff) - τii/(2*ηe) - εii_pwl 

        # Residual check
        res = max(max(abs(r_η)))/sqrt(ε2_eff)
        if iter == 1 res0 = res end
        if noisy>2 @printf("Newton VE It. %02d, r abs. = %2.2e r rel. = %2.2e\n", iter, res, res/res0) end
        if res/res0 < tol || res<tol
            break
        end

        # Jacobian (analytic partial derivatives)
        d11 = -sqrt(ε2_eff)/ηe - εii_pwl.*npwl./ηve

        # Newton updates
        ηve -= r_η/d11
    end
    @printf("a_ve NEWTON = %2.2e\n", a_ve)
    
    # A posteriori stress evaluation
    D       = 2*ηve * [1 0 0; 0 1 0; 0 0 a_ve^-2]
    τ_vec   = D*[εxx_eff; εyy_eff; εxy_eff]
    τ_rot   = [τ_vec[1] τ_vec[3]; τ_vec[3] τ_vec[2]]

    return τ_rot, ηve, a_ve, ηpwl
end