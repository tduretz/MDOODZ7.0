
using Printf, Plots

function IsotropicDruckerPrager( γ̇, p )
    Tn2    = p.Tn2
    Ts2    = p.Ts2
    τii    = sqrt( 1/2*(Tn2 + Ts2))
    eta_ve = p.eta_ve;
    P      = p.P; 
    fric   = p.fric; 
    C      = p.C;
    eta_vp = p.eta_vp
    τiic   = τii  - γ̇*eta_ve
    F      = τiic - P*sin(fric) - C*cos(fric) - eta_vp*γ̇
    return F
end

function AnisotropicDruckerPrager( γ̇, p )
    eta_ve    = p.eta_ve;
    P         = p.P; 
    fric      = p.fric; 
    C         = p.C;
    eta_vp    = p.eta_vp
    Tn2       = p.Tn2
    Ts2       = p.Ts2
    ani_rat   = p.ani_rat
    am        = p.am
    Ji        = p.Ji
    a_p       = p.a_p
    τii       = sqrt( 1/2*(Tn2 + Ts2))
    J2_corr   = 0.
    J1_corr   = 0.
    ap_n      = 1.0 - eta_ve.*γ̇./τii;
    ap_s      = 1.0 - eta_ve.*γ̇./τii * ani_rat;
    ap_n = ap_s
    J2_corr   = 1/2*( Tn2*ap_n^2 + a_p^2*Ts2*ap_s^2)
    F         = sqrt(J2_corr) - C*cos(fric) - P.*am.*sin(fric) + J1_corr - eta_vp.*γ̇;
    # F    =  τii  - γ̇*eta_ve  - C*cos(fric) - P.*am.*sin(fric) + J1_corr - eta_vp.*γ̇;
    return F
end

function main()

    # Return mapping will fail if a_p != a_v (no root, F always larger than 0) 

    solve              = false
    isotropic_DP_ckeck = false

    eta = 1e23
    L   = 1e6
    V   = 1e-9
    T   = 100
    t   = L/V
    E   = 1/t
    S   = eta*E

    K = 0.00e+00; fric=0.523599; dil = 0.000000
    Ft = 6.6592e+00; Y2_v = 6.2088e+02; Y2_p = 1.5258e+02; eta_vp = 1.00e-02
    a_ve = 2.236068e+00; eta_ve = 2.721397e+01; G = 1.000000e+01;  dt = 1.262300e-04
    a_e = 2.236068; a_v = 2.236068; a_p = 1.000000
    Txx = -5.958186e+00; Tyy = 5.958186e+00; Tzz = -6.190604e-13; Txy = -1.082017e+01; P = 1.052000e+01
    gdot = 0.000000e+00
    C = 5.000000e-01; fric = 5.235988e-01
    a1 = 1.000000e+00; a2 = 1.000000e+00; a3 = 1.000000e+00
    # a_p = a_ve    # !!!!!!!!!!!! 
    
    dil = fric/2; K = 2G
    Tyx = Txy

    γ̇      = gdot
    F      = Ft
    Jii    = Y2_p
    J2_corr = Y2_p
    Δt     = 0.
    nitmax = 10
    noisy  = 3
    tol    = 1e-10

    println("--------------------")
    println(sqrt(Y2_v)*S)
    println(sqrt(Y2_p)*S)
    println(P*S)
    println(eta_ve*eta)
    println(dt*t)
    println(sqrt(Jii)*S)
    println("--------------------")

    # Linear isotropic yield
    @show Y2_p = 1/2*(Txx^2 + Tyy^2 + Tzz^2) + a_p^2*Txy^2
    F = sqrt(Y2_p) - C*cos(fric) - P*sin(fric)
    γ̇ = F/(eta_ve+eta_vp)
    @show sqrt(Y2_p) - γ̇*eta_ve
    @show F = sqrt(Y2_p) - γ̇*eta_ve - C*cos(fric) - P*sin(fric) - eta_vp .* gdot
    τ_ax = (0.0:1e7:6e9) ./S
    p_ax = (0.0:1e7:1e9) ./S
    Fdp  = zeros(length(p_ax), length(τ_ax))
    for j in eachindex(τ_ax),  i in eachindex(p_ax)
        Fdp[i,j]  = τ_ax[j] - C*cos(fric) - p_ax[i]*sin(fric) 
    end


    # ckecks J2
    am      = (a1 + a2 + a3)/3
    Tn2     = Txx^2+Tyy^2+Tzz^2
    Ts2     = Txy^2+Tyx^2
    Tii     = sqrt(Jii)
    τii     = sqrt(Jii)
    ani_rat = a_p^2/a_ve^2
    gdot    = γ̇
    Jii     = Y2_p
    J2_corr = Y2_p
    J2_corr = Txx .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1) .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2)) + 1) .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 .* (-a_p .^ 2 .* eta_ve .* gdot ./ (a_ve .^ 2 .* sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2)) + 1) .^ 2 / 2 + Tyy .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1) .^ 2 / 2 + Tzz .^ 2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1) .^ 2 / 2
    @show (1, sqrt(J2_corr))
    J2_corr = Y2_p*(1.0 - eta_ve * gdot/sqrt(Y2_p))^2
    @show (2, sqrt(J2_corr))
    @show (3, sqrt(Y2_p)*(1-eta_ve * gdot/sqrt(Y2_p)) )
    @show (4, sqrt(Y2_p) - eta_ve * gdot )
    @show (5, sqrt(Y2_p) - γ̇*eta_ve )
    ap_n      = 1.0 - eta_ve .*  γ̇./τii;
    ap_s      = 1.0 - eta_ve .*  γ̇./τii * ani_rat;
    J2_corr   = 1/2*( Tn2*ap_n^2 + a_p^2*Ts2*ap_s^2)
    @show (6, sqrt(J2_corr))
    @show ("J2", sqrt(Y2_p))
    @show F = sqrt(J2_corr) - C*cos(fric) - P*sin(fric) - eta_vp .* gdot
   
    # ckecks J1
    a1  = 1.0
    J1c = (Txx .* a1 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1) + Tyy .* a2 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1) + Tzz .* a3 .* (-eta_ve .* gdot ./ sqrt(Txx .^ 2 / 2 + Txy .^ 2 .* a_p .^ 2 / 2 + Tyx .^ 2 .* a_p .^ 2 / 2 + Tyy .^ 2 / 2 + Tzz .^ 2 / 2) + 1)) .* sin(fric) / 3
    @show ( 1, J1c*S )
    Ji   = sin(fric)/3*( a1*Txx + a2*Tyy + a3*Tzz)
    ap_n = 1. - eta_ve .* gdot ./ sqrt(Jii)
    J1c  = ap_n*Ji
    @show ( 2, J1c*S )

    if isotropic_DP_ckeck
        a1, a2, a3 = 1.0, 1.0, 1.0
        a_p = 1.
        a_ve = 1.0
        ani_rat = a_p^2/a_ve^2
    end
    
    # Pre compute
    am      = (a1 + a2 + a3)/3
    Tn2     = Txx^2+Tyy^2+Tzz^2
    Ts2     = Txy^2+Tyx^2
    Tii     = sqrt(Jii)
    ani_rat = a_p^2/a_ve^2

    # Make 2D maps of yield function
    γ̇v    = LinRange(0.0, 20, 500)
    Fiso  = zero(γ̇v)
    Fani  = zero(γ̇v)
    p     = (τii = sqrt(Y2_p), eta_ve = eta_ve, P = P, fric = fric, C = C, eta_vp = eta_vp, Tn2 = Tn2, Ts2 = Ts2, ani_rat = ani_rat, Ji = Ji, am = am, a_p=a_p)
    for i in eachindex(γ̇v)
        Fiso[i] =   IsotropicDruckerPrager( γ̇v[i], p )
        Fani[i] = AnisotropicDruckerPrager( γ̇v[i], p )
    end

    p=plot(γ̇v, Fiso, xlabel="gamma", ylabel="F")
    p=plot!(γ̇v, Fani, title=minimum(Fani))
    (v,imin) = findmin(Fani)
    gdot = γ̇v[imin]
     display(p)

    # gdot = a_ve .^ 2 .* (Jii .^ (5 // 2) .* (Txx .^ 2 .* a_ve .^ 2 + Txy .^ 2 .* a_p .^ 4 + Tyx .^ 2 .* a_p .^ 4 + Tyy .^ 2 .* a_ve .^ 2 + Tzz .^ 2 .* a_ve .^ 2) + sqrt(-Jii .^ 5 .* a_p .^ 2 .* (Txx .^ 2 .* Txy .^ 2 .* a_p .^ 4 - 2 * Txx .^ 2 .* Txy .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Txx .^ 2 .* Txy .^ 2 .* a_ve .^ 4 + Txx .^ 2 .* Tyx .^ 2 .* a_p .^ 4 - 2 * Txx .^ 2 .* Tyx .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Txx .^ 2 .* Tyx .^ 2 .* a_ve .^ 4 + Txy .^ 2 .* Tyy .^ 2 .* a_p .^ 4 - 2 * Txy .^ 2 .* Tyy .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Txy .^ 2 .* Tyy .^ 2 .* a_ve .^ 4 + Txy .^ 2 .* Tzz .^ 2 .* a_p .^ 4 - 2 * Txy .^ 2 .* Tzz .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Txy .^ 2 .* Tzz .^ 2 .* a_ve .^ 4 + Tyx .^ 2 .* Tyy .^ 2 .* a_p .^ 4 - 2 * Tyx .^ 2 .* Tyy .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Tyx .^ 2 .* Tyy .^ 2 .* a_ve .^ 4 + Tyx .^ 2 .* Tzz .^ 2 .* a_p .^ 4 - 2 * Tyx .^ 2 .* Tzz .^ 2 .* a_p .^ 2 .* a_ve .^ 2 + Tyx .^ 2 .* Tzz .^ 2 .* a_ve .^ 4))) ./ (Jii .^ 2 .* eta_ve .* (Txx .^ 2 .* a_ve .^ 4 + Txy .^ 2 .* a_p .^ 6 + Tyx .^ 2 .* a_p .^ 6 + Tyy .^ 2 .* a_ve .^ 4 + Tzz .^ 2 .* a_ve .^ 4))

    ap_n      = 1.0 - eta_ve .* gdot./Tii;
    ap_s      = 1.0 - eta_ve .* gdot./Tii * ani_rat;
    J2_corr   = 1/2*( Tn2*ap_n^2 + a_p^2*Ts2*ap_s^2)
    

    # p = heatmap(p_ax.*S./1e6, τ_ax.*S./1e6, Fdp'./1e6.*S)
    # p = contour!(p_ax.*S./1e6, τ_ax.*S./1e6, Fdp'./1e6.*S, levels=[0,0], c=:white)
    # p = plot!( [P].*S./1e6, [sqrt(Y2_p)].*S./1e6, markershape=:star, c=:white, label=:none)
    # p = plot!( [P].*S./1e6, [sqrt(Y2_p)-γ̇*eta_ve].*S./1e6, markershape=:star, c=:white, label=:none)

    # p = plot!( [P].*S./1e6, [sqrt(J2_corr)].*S./1e6, markershape=:star, c=:red, label=:none)


    # println(P.*S./1e6)
    # println(sqrt(Y2_p).*S./1e6)
    # display(p)

    # eta_ve_n = eta_ve
    # eta_ve_s = eta_ve_n/a_ve/a_ve
    # g_n = 0.0
    # g_s = 0.0
    # Pc  = P

    # for i=1:200

    #     Txxc    = Txx - g_n*eta_ve_n*Txx/2/Tii
    #     Tyyc    = Tyy - g_n*eta_ve_n*Tyy/2/Tii
    #     Tzzc    = Tzz - g_n*eta_ve_n*Tzz/2/Tii
    #     Txyc    = Txy - g_s*eta_ve_s*Txy/2/Tii
    #     Tyxc    = Txyc 
    #     Tn2     = Txxc^2 + Tyyc^2 + Tzzc^2
    #     Ts2     = Txyc^2 + Tyxc^2 
    #     J2_corr = 1/2*( Tn2 + a_p^2*Ts2)
    #     Fc      = sqrt(J2_corr) - C*cos(fric) - Pc.*am.*sin(fric) - eta_vp.*g_n;
    #     f_n     = Txxc/2/eta_ve_n + g_n*Txx/2/Tii
    #     f_s     = Txyc/2/eta_ve_s + g_s*Txy/2/Tii
    #     g_n    += f_n
    #     g_s    += f_s
    #     @show (f_n, f_s, Fc, g_n, g_s)
    # end


    if solve 
        γ̇ = 0;

        for iter=1:nitmax

            gdot   = γ̇
            dt     = Δt

            # With out-of-plane component
            Pc        = P + am*K*Δt*γ̇*sin(dil)

            ap_n      = 1.0 - eta_ve .* gdot./Tii;
            ap_s      = 1.0 - eta_ve .* gdot./Tii * ani_rat;
            
            J1_corr   = ap_n*Ji
            J2_corr   = 1/2*( Tn2*ap_n^2 + a_p^2*Ts2*ap_s^2)
            Fc        =  sqrt(J2_corr) - C*cos(fric) - Pc.*am.*sin(fric) + J1_corr - eta_vp.*gdot;
            @show  (sqrt(J2_corr),  C*cos(fric), Pc.*am.*sin(fric), J1_corr, eta_vp.*gdot)
    
            dapndgdot = -eta_ve/Tii
            dapsdgdot = dapndgdot*ani_rat
            dFdgdot   = -eta_vp + 1/2/sqrt(J2_corr) * ( a_p^2*Ts2*ap_s*dapsdgdot + Tn2*ap_n*dapndgdot ) - am^2*K*dt*sin(dil)*sin(fric) + dapndgdot*Ji
        
            dFdγ̇      = dFdgdot

            # It would be nice to have this section derived with Symbolics.jl
            if noisy>2 @printf("Newton VP It. %02d, r abs. = %2.2e r rel. = %2.2e, r ini. = %2.6e \n", iter, Fc, Fc/F, F) end
            if abs(Fc)/abs(F) < tol || abs(Fc)<tol
                break
            end
            γ̇    -= Fc/dFdγ̇
        end
    end

    # check incompressibility
    Exx =  γ̇*Txx/2/Tii
    Eyy =  γ̇*Tyy/2/Tii
    Ezz =  γ̇*Tzz/2/Tii
    @show Exx
    @show Eyy
    @show Ezz
    @show Exx+Eyy+Ezz

end

main()