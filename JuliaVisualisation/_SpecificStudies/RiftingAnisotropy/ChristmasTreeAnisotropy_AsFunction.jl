using CairoMakie, Printf, StaticArrays
import Statistics:mean

# Try to implement anistropic friction but resulting invariants are not objective

@views function ChrismasTreeAniso()
    # Anisotropy parameters
    δ    = 3
    θ    = π/4 #LinRange( 0.0, π, 51 ) .-0* π/2
    # Domain
    ymin       = -150e3
    ymax       = 0e3
    ncy        = 100
    Δy         = (ymax - ymin)/ncy
    yc         = LinRange(ymin+Δy/2, ymax+Δy/2, ncy)
    nt         = 12
    Δt         = 1e12
    # Kinematics
    pure_shear = 1
    ε̇xxd       = pure_shear*1e-15
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = 0.
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    τxx        = zeros(ncy)
    τyy        = zeros(ncy)
    τxy        = zeros(ncy)
    # Gravity
    g             = -9.81
    ρ0            = 3300.
    ρ             = ρ0*ones(ncy)
    ρ[yc.>-35e3] .= 2800.
    Pl            = zeros(ncy)
    for i = ncy-1:-1:1
        Pl[i] += Pl[i+1] - g*Δy*ρ[i+1]
    end
    # Temperature 
    T0 = 273.15
    hr = 10e3
    qm = 30e-3
    κ  = 1e-6
    k  = 3.0
    h  = -35e3
    H0 = 3.0e-9/1.9
    T  = -(qm/k .- ρ0.*H0.*hr./k.*exp.(h./hr)).*yc .+ ρ0.*H0.*hr^2/k.*(1.0 .- exp.(yc./hr)) .+ T0
    T  .-= minimum(T)
    T  .+= T0
    T[T.>1370] .= 1370
    # T = 773*ones(ncy)
    # Elasticity
    G          = 3e10
    K          = 8e10  
    ηe         = G*Δt
    # Plasticity
    C          = 5e7
    fric       = 30*π/180
    dil        = 0*π/180
    ηvp        = 0*1e21
    # Pressure in extension
    X    = 1 / (sqrt((δ^2-1)*cos(2*θ)^2 + 1)/δ)
    ϕeff = asin(sin(fric) / X)
    @show ϕeff*180/π
    σH = @. -Pl*(1-sin(ϕeff))/(1+sin(ϕeff))
    P  = @. (Pl - σH)/2
    @show cumsum(ρ .* g .* Δy)
    # Power law
    R    = 8.314
    npwl = 3.5*ones(ncy)
    Qpwl = 530.0e3*ones(ncy)
    Apwl = 1.1000e-16*ones(ncy)
    npwl[yc.>-35e3] .= 3.3
    Qpwl[yc.>-35e3] .= 186.5e3
    Apwl[yc.>-35e3] .= 3.1623e-26

    Bpwl  = Apwl.^(-1.0./npwl) .* exp.( Qpwl./R./npwl./T )
    Cpwl       = (2*Bpwl).^(-npwl)
    # if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
    # # -------------------- TEST: Anisotropic power law -------------------- #
    # Arrays
    τii_cart1  = zeros(ncy) 
    τii_cart2  = zeros(ncy)
    τii_rot1   = zeros(ncy)
    τii_rot2   = zeros(ncy)
    #  Loop over all orientations
    for i in eachindex(θ)
        τxx .= .0
        τyy .= .0
        τxy .= .0
        @show θ[i] 
        for it=1:nt
            for ic=1:ncy
                τxx0 = τxx[ic]
                τyy0 = τyy[ic]
                τxy0 = τxy[ic]
                # Transformation matrix: towards principal plane
                Q           = @SMatrix([cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])])
                ε̇_tens      = @SMatrix([ε̇xxd ε̇xyd; ε̇xyd ε̇yyd])
                τ0_tens     = [τxx0 τxy0; τxy0 τyy0]
                ε̇_rot       = Q*ε̇_tens*Q' 
                τ0_rot      = Q*τ0_tens*Q'
                # VISCO ELASTICITY
                ε̇           = @SMatrix([ε̇_rot[1,1]+τ0_rot[1,1]/2/ηe; ε̇_rot[2,2]+τ0_rot[2,2]/2/ηe; ε̇_rot[1,2]+τ0_rot[1,2]/2/ηe*δ])
                I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
                ε̇ii         = sqrt(I2)
                ηpwl        = Bpwl[ic].*sqrt(I2).^(1.0./npwl[ic].-1)
                ηve         = (1.0./ηpwl .+ 1.0./ηe).^(-1)
                τii         = 0.
                r0          = 0.
                for iter=1:30
                    τii        = 2*ηve*ε̇ii
                    ε̇pwl       = Cpwl[ic]*τii.^npwl[ic]
                    r          = ε̇ii - τii/2/ηe - ε̇pwl
                    if iter==1 r0 = r end
                    ∂r∂ηve     = - ε̇ii/ηe - ε̇pwl*npwl[ic]/ηve
                    ηve       -= r/∂r∂ηve
                    # @show (iter, abs(r)/abs(r0))
                    if (abs(r)/abs(r0)<1e-6) break; end
                end
                D           = 2*ηve*@SMatrix([1 0 0; 0 1 0; 0 0 1.0/δ;])
                τ           = D*ε̇
                # ---> Independent on orientation (objective) !!!!!!
                Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
                τii = sqrt(Y2); τxxt = τ[1]; τyyt = τ[2]; τzzt = -τxxt-τyyt; τxyt = τ[3]

                # PLASTICITY
                γ̇ = 0.
                F  = τii - cos(fric)*C - P[ic]*sin(fric)
                if F>0
                    γ̇    = F/(ηve + ηvp + K*Δt*sin(dil)*sin(fric))
                    τiic = τii  - γ̇*ηve
                    Pc   = P[ic]    + K*Δt*sin(dil)*γ̇ 
                    τxxc = τxxt - 2*ηve*γ̇*(τxxt/τii/2)
                    τyyc = τyyt - 2*ηve*γ̇*(τyyt/τii/2)
                    τzzc = -τxxc-τyyc
                    F    = τiic - cos(fric)*C - Pc*sin(fric) - ηvp*γ̇ 
                    # @show (F)
                    Y2   = τiic^2
                    ηvep = τiic/2/ε̇ii
                    ηve  = ηvep
                    D    = 2*ηve*@SMatrix([1 0 0; 0 1 0; 0 0 1.0/δ;])
                    τ    = D*ε̇ 
                    Y2   = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
                    τiic = sqrt(Y2)

                    # @show τxxt, τyyt, τzzt
                    # @show τxxc, τyyc, τzzc
                    # @show τ[1], τ[2], -τ[1]-τ[2]
                    F    = τiic - cos(fric)*C - Pc*sin(fric) - ηvp*γ̇ 
                    # @show (F)  
                end
                τii_rot1[ic] = sqrt(Y2)


                # Check strain rate components
                τxx[ic]         = τ[1]
                τyy[ic]         = τ[2]
                τxy[ic]         = τ[3]
                ε̇xxd_v = Cpwl[ic] *τii^((npwl[ic] -1))*τxx[ic] 
                ε̇xyd_v = Cpwl[ic] *τii^((npwl[ic] -1))*(τxy[ic]/δ)
                ε̇xxd_e = (τxx[ic] - τ0_rot[1,1])/2/ηe
                ε̇xyd_e = (τxy[ic] - τ0_rot[1,2])/2/ηe*δ
                ε̇xxd_p = γ̇*τxxt/τii/2
                ε̇xyd_p = γ̇*τxyt/τii/2*δ
                fxx    = ε̇_rot[1,1] - ε̇xxd_v - ε̇xxd_e - ε̇xxd_p
                fxy    = ε̇_rot[1,2] - ε̇xyd_v - ε̇xyd_e - ε̇xyd_p
                # @show  (fxx, fxy)
                # Effective viscosity: compute stress invariant from strain rate invariant
                # ---> Predicts rotated stress invariant !!!!!!
                I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2
                τii_rot2[ic] = 2*ηve*sqrt(I2)
                # If one does the same by scaling the strain rate invariant by δ^2...
                # ---> Predicts Cartesian stress invariant !!!!!!
                I2          = 0.5*(ε̇[1]^2 + ε̇[2]^2) + ε̇[3]^2/δ^2
                τii_cart1[ic] = 2*ηve*sqrt(I2)
                # Rotate stress back 
                τ_rot       = @SMatrix([τ[1] τ[3]; τ[3] τ[2]])
                τ           = Q'*τ_rot*Q

                ε̇_rot       = @SMatrix([ε̇[1] ε̇[3]; ε̇[3] ε̇[2]])
                ε̇_eff       = Q'*ε̇_rot*Q

                ε̇_eff1 = [ε̇_eff[1,1]; ε̇_eff[2,2]; ε̇_eff[2,1]]
                # @show ε̇_eff1

                n     = θ[i] - π/2  
                nx    = cos(n)
                ny    = sin(n)
                d1    = 2*nx^2*ny^2
                d2    = nx*nx*(ny^2 - nx^2)
                ani   = 1.0 - 1.0 / δ
                ani   = 1.0 - 1.0 / δ
                A = [2.0-2.0*ani*d1 2.0*ani*d1 -2.0*ani*d2;
                     2.0*ani*d1 2.0-2.0*ani*d1 2.0*ani*d2;
                    -2.0*ani*d2     2.0*ani*d2 1.0  + 2.0*ani*(d1 - 0.5);]
            
                E = [ε̇xxd; ε̇yyd; ε̇xyd]
                τ0 = [τxx0/ηe; τyy0/ηe; τxy0/ηe]
                b = A\τ0
                Eeff = E + b./[1; 1; 2]

                Da11  = 2.0 - 2.0*ani*d1;
                Da12  = 2.0*ani*d1;
                Da13  =-2.0*ani*d2;
                Da22  = 2.0 - 2.0*ani*d1;
                Da23  = 2.0*ani*d2;
                Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
                a11   = Da33 * Da22 - Da23^2;
                a12   = Da13 * Da23 - Da33 * Da12;
                a13   = Da12 * Da23 - Da13 * Da22;
                a22   = Da33 * Da11 - Da13^2;
                a23   = Da12 * Da13 - Da11 * Da23;
                a33   = Da11 * Da22 - Da12^2;
                det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
                iDa11 = a11/det; 
                iDa12 = a12/det; 
                iDa13 = a13/det;
                iDa22 = a22/det; 
                iDa23 = a23/det;
                iDa33 = a33/det;
                Exx = ε̇xxd + (iDa11*τxx0 + iDa12*τyy0 + iDa13*τxy0)/ηe;
                Eyy = ε̇yyd + (iDa12*τxx0 + iDa22*τyy0 + iDa23*τxy0)/ηe;
                Exy = ε̇xyd + (iDa13*τxx0 + iDa23*τyy0 + iDa33*τxy0)/ηe/2.0;

                Eeff1 = [Exx;  Eyy ; Exy ]
                # @show  Eeff1
                # println(' ')

                
                # ---> Dependent on orientation (non-objective) !!!!!!
                J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
                τii_cart2[ic] = sqrt(J2) 
                τxx[ic] = τ[1,1]
                τyy[ic] = τ[2,2]
                τxy[ic] = τ[1,2]
            end
        end

        # f = Figure(size = (500, 500), fontsize=18)

        # #####################

        # ax2 = Axis(f[1, 1], xlabel = L"A) $T$ [$^\circ$C]", ylabel = L"$z$ [km]",  xticklabelcolor = :tomato, xaxisposition=:top) # , title = L"$$Pressure profile"
        # ax1 = Axis(f[1, 1], xlabel = L"$P$ [GPa]", xticklabelcolor = :blue,  xaxisposition=:bottom) # , title = L"$$Temperature profile"
        # hidespines!(ax2)
        # hideydecorations!(ax2)
        

        # lines!(ax2, T .- T0, -35*ones(size(T)), color=:gray, linestyle=:dash)

        # lines!(ax1, P./1e9, yc./1e3, color=:blue )
        # lines!(ax2, T .- T0, yc./1e3, color=:tomato)
        
        # #####################

        τii_yield     = P*sin(fric) .+ C*cos(fric)
        τii_yield_min = τii_yield*sqrt((δ^2-1)*cos(2*θ[i])^2 + 1)/δ


        # ax3 = Axis(f[1, 2], title = L"B) $$Stress profile ($\delta = 3$)", xlabel = L"$τ_\mathrm{II}$ [GPa]")
        # lines!(ax3, τii_rot1./1e9, yc./1e3, label=L"$\sqrt{Y_2^{'}}$", color=:blue )
        # xlims!(ax3, 0, .7)
        # # lines!(ax3, τii_rot2./1e9, yc./1e3 )
        # lines!(ax3, τii_cart1./1e9, yc./1e3, label=L"$\sqrt{J_2}$", color=:tomato )

        # lines!(ax3, τii_yield./1e9,     yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = 0.0$)", linestyle=:dash, color=:blue )
        # lines!(ax3, τii_yield_min./1e9, yc./1e3, label=L"$τ^\mathrm{fric}$ ($\theta = π/4$)", linestyle=:dash, color=:tomato )

        # # lines!(ax3, τii_cart2./1e9, yc./1e3 )
        # axislegend(position=:rb)
        # display(f)

        # save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/ChristmasTreeAnisotropy.png", f, px_per_unit = 4) 

        # @show τii_rot1./τii_cart1

        return τii_rot1, yc, τii_cart1, τii_yield, τii_yield_min, T, T0, P

    end
end

# ChrismasTreeAniso()
