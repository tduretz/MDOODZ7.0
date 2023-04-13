using Plots

Inv2(Mij) = sqrt(1/2 * sum(Mij.^2))

function get_tau_ij(ϵij; Δt=1e10, η=1e20, G=1e10, Τij0=[0. 0.; 0. 0.])
    Eij = ϵij .+ Τij0 ./ (2 * G * Δt)
    ηVE = (1 / (G * Δt)+ 1 / η)^-1
    return 2 * ηVE .* Eij
end

n_steps = 20
Δt      = 1e10

function get_second_invariant(ϵij)
    J2_values = zeros(n_steps+1)

    Τij = [0. 0.; 0. 0.]

    for i = 1:n_steps+1
        Τij0         = Τij
        Τij          = get_tau_ij(ϵij, Τij0=Τij0)
        J2_values[i] = Inv2(Τij)
    end

    return J2_values
end

# PS
pure_shear_ϵij = [1e-14 0; 0 -1e-14]
τ2_pure_shear  = get_second_invariant(pure_shear_ϵij)

# SS
simple_shear_ϵij = [0 1e-14; 1e-14 0]
τ2_simple_shear  = get_second_invariant(simple_shear_ϵij)

# Viz
plot(title="Stress build-up",   xlabel="Time [s]", ylabel="τ2 [Pa]")
plot!([0:n_steps].*Δt, τ2_pure_shear,   label="PS - VE")
plot!([0:n_steps].*Δt, τ2_simple_shear, label="SS - VE")