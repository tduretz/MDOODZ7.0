function get_tau_ij(ϵij; Δt=1e10, η=1e20, G=1e10, Τij0=0)
    Eij = ϵij .+ Τij0 / (2 * G * Δt)
    ηVE = 1/G * Δt + 1 / η

    return 2 * ηVE .* Eij
end

pure_shear_ϵij = [1e-14 0; 0 -1e14]
get_tau_ij(pure_shear_ϵij)

simple_shear_ϵij = [0 1e-14; 2e-14 0]
get_tau_ij(simple_shear_ϵij)