using Plots

Inv2(Mij) = sqrt(1/2 * sum(Mij.^2))

function get_tau_ij(ϵij; Δt=1e10, η=1e20, G=1e10, Τij0=[0. 0.; 0. 0.])
    Eij = ϵij .+ Τij0 ./ (2 * G * Δt)
    ηVE = (1 / (G * Δt)+ 1 / η)^-1
    return 2 * ηVE .* Eij
end

n_steps = 20
Δt      = 1e10

function get_τ2_for_timesteps(ϵij; G=1e10, η=1e20)
    J2_values = zeros(n_steps+1)

    Τij = [0. 0.; 0. 0.]

    for i = 1:n_steps+1
        Τij0         = Τij
        Τij          = get_tau_ij(ϵij, Τij0=Τij0, G=G, η = η)
        J2_values[i] = Inv2(Τij)
    end

    return J2_values
end

pure_shear_ϵij = [1e-14 0; 0 -1e-14]
simple_shear_ϵij = [0 1e-14; 1e-14 0]

# PS - VE
τ2_pure_shear_ve  = get_τ2_for_timesteps(pure_shear_ϵij)

# SS - VE
τ2_simple_shear_ve  = get_τ2_for_timesteps(simple_shear_ϵij)

# PS - V
τ2_pure_shear_v  = get_τ2_for_timesteps(pure_shear_ϵij, G=1e20)

# SS - V
τ2_simple_shear_v  = get_τ2_for_timesteps(simple_shear_ϵij, G=1e20)

# PS - E
τ2_pure_shear_e  = get_τ2_for_timesteps(pure_shear_ϵij, η=1e40)

# SS - E
τ2_simple_shear_e  = get_τ2_for_timesteps(simple_shear_ϵij, η=1e40)

# Viz
plot(title="Stress build-up",   xlabel="Time [s]", ylabel="τ2 [Pa]")
plot!([0:n_steps].*Δt, τ2_pure_shear_ve,   label="PS - VE")
plot!([0:n_steps].*Δt, τ2_simple_shear_ve, label="SS - VE")
plot!([0:n_steps].*Δt, τ2_pure_shear_v,   label="PS - V")
plot!([0:n_steps].*Δt, τ2_simple_shear_v, label="SS - V")
plot!([0:n_steps].*Δt, τ2_pure_shear_e,   label="PS - E")
plot!([0:n_steps].*Δt, τ2_simple_shear_e, label="SS - E")