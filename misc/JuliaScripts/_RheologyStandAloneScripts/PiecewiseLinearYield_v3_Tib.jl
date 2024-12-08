using Plots

line1_Dani(p,K,dt,η_ve,η_vp,p1,t1)       = (η_ve + η_vp)/(dt*K + 2.0/3.0*η_vp)*(-p + p1) + t1 
line2_Dani(p,K,dt,η_ve,η_vp,p1,t1)       = (η_ve + η_vp)/(dt*K + 2.0/3.0*η_vp)*(-p + p1) + t1 
line3_Dani(p,K,dt,η_ve,η_vp,p1,t1,φ,ψ)   = (η_ve + η_vp)/((dt*K + 2.0/3.0*η_vp)*sin(ψ))*(-p + p1) + t1  

function line(p, K, dt, η_ve, ψ, p1, t1)
    p2 = p1 + K*dt*sin(ψ)
    t2 = t1 - η_ve  
    a  = (t2-t1)/(p2-p1)
    b  = t2 - a*p2
    return a*p + b
end

function plastic_corrections!(τ, P, η_ve, η_vp, K, φ, ψ, C, σ_T, δσ_T, dt, pc1, τc1, pc2, τc2)
    ε_dev_pl  = 0.0
    ε_vol_pl  = 0.0
    domain_pl = 0.0

    l1    = line(P, K, dt, η_ve, π/2, pc1, τc1)
    l2    = line(P, K, dt, η_ve, π/2, pc2, τc2)
    l3    = line(P, K, dt, η_ve,   ψ, pc2, τc2)
    Pc = P
    τc = τ

    if max(τ - P*sin(φ) - C*cos(φ) , τ - P - σ_T , - P - (σ_T - δσ_T) ) > 0.0                                                         # check if F_tr > 0
        if τ <= τc1 
            # pressure limiter 
            dqdp = -1.0
            f    = - P - (σ_T - δσ_T) 
            λ̇    = f / (K*dt)                                                                                                                          # tensile pressure cutoff
            τc   = τ 
            Pc   = P - K*dt*λ̇*dqdp
            f    = - Pc - (σ_T - δσ_T) 
            domain_pl = 1.0
        elseif τc1 < τ <= l1    
            # corner 1 
            τc = τ - η_ve*(τ - τc1)/(η_ve)
            Pc = P - K*dt*(P - pc1)/(K*dt)
            domain_pl = 2.0
        elseif l1 < τ <= l2            # mode-1
            # tension
            dqdp = -1.0
            dqdτ =  1.0
            f    = τ - P - σ_T 
            λ̇    = f / (K*dt + η_ve) 
            τc   = τ - η_ve*λ̇*dqdτ
            Pc   = P - K*dt*λ̇*dqdp
            domain_pl = 3.0 
        elseif l2< τ <= l3 # 2nd corner
            # corner 2
            τc = τ - η_ve*(τ - τc2)/(η_ve)
            Pc = P - K*dt*(P - pc2)/(K*dt)
            domain_pl = 4.0
        elseif l3 < τ  
            # Drucker-Prager                                                              # Drucker Prager
            dqdp = -sin(ψ)
            dqdτ =  1.0
            f    = τ - P*sin(φ) - C*cos(φ) 
            λ̇    = f / (K*dt*sin(φ)*sin(ψ) + η_ve) 
            τc   = τ - η_ve*λ̇*dqdτ
            Pc   = P - K*dt*λ̇*dqdp
            @show domain_pl = 5.0 
        end
    end
    return Pc, τc
end

function main() 

    # Material properties
    G      = 2e10          # shear modulus
    K      = 1e10          # compressibility, 1/Pa
    φ      = 30.0*π/180    # friction angle, deg
    ψ      = 10.0*π/180    # dilation angle, deg
    σ_T    = 3e6           # tensile strength, always positive, Pa
    C      = 2.0*σ_T       # cohesion, typically 2x tensile strength, always positive, Pa
    δσ_T   = 1e6           # pressure limiter (if τII_tr<δσ_T) only pressure correction is applied, Pa
    # viscous regularization parameters
    η_vp   = 0.0           #           
    # Calculate corner coordinates 
    pc1    = -(σ_T - δσ_T)                           # p at the intersection of cutoff and Mode-1
    τc1    = δσ_T                                    # τII at the intersection of cutoff and Mode-1
    pc2    = -(σ_T - C*cos(φ))/(1.0-sin(φ))          # p at the intersection of Drucker-Prager and Mode-1
    τc2    = pc2 + σ_T                               # τII at the intersection of Drucker-Prager and Mode-1
    # Trial stresses
    τII_tr = [ 1e5 1.5e6 5e6   5e6 1e7  1e7 2e7]           # second invariant of trial stress, Pa
    P_tr   = [-3e6 -3e6 -1.5e6 0.0 5e6  4e6 2e7]           # trial pressure, Pa
    τII    = similar(τII_tr)
    P      = similar(P_tr)
    # numerical discretization
    dt     = 1e1          # time step, s 
    η_ve   = 1 / (1/(G*dt))

    p_tr1 = LinRange(-1e7, 0, 100)
    p_tr2 = LinRange(0, 1e7, 100)
    p_tr3 = LinRange(3e6, 6e6, 100)
    l1    = line.(p_tr1, K, dt, η_ve, π/2, pc1, τc1)
    l2    = line.(p_tr2, K, dt, η_ve, π/2, pc2, τc2)
    l3    = line.(p_tr3, K, dt, η_ve,   ψ, pc2, τc2)

    l1D    = line1_Dani.(p_tr1, K, dt, η_ve, η_vp, pc1, τc1)
    l2D    = line2_Dani.(p_tr2, K, dt, η_ve, η_vp, pc2, τc2)
    l3D    = line3_Dani.(p_tr3, K, dt, η_ve, η_vp, pc2, τc2,φ,ψ)

    P_end =  30e6
    p = plot(aspect_ratio=1)
    p = plot!([pc1, pc1, pc2, P_end],[0.0, τc1, τc2, P_end*sin(φ)+C*cos(φ)])
    p = plot!(p_tr1, l1, label="l1")
    p = plot!(p_tr1, l1D, label="l1D")
    p = plot!(p_tr2,  l2, label="l2")
    p = plot!(p_tr2, l2D, label="l2D")
    p = plot!(p_tr3,  l3, label="l3")
    p = plot!(p_tr3, l3D, label="l3D")

    for i in eachindex(P_tr)
        P[i], τII[i] = plastic_corrections!(τII_tr[i], P_tr[i], η_ve, η_vp, K, φ, ψ, C, σ_T, δσ_T, dt, pc1, τc1, pc2, τc2)
        p=plot!([P_tr[i], P[i]],[τII_tr[i], τII[i]])
    end
    display(p)

end

main()