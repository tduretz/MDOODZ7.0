using LinearAlgebra
using Plots
using MAT


# Anisotropy parameters
finite_strain = true
aniso_factor = 10.0

# Simple shear
dt = 1e-2

# Velocity gradients
dudx = 0.0
dudy = 1.0
dvdx = 0.0
dvdy = 0.0

μ_N   = 1.0 # normal viscosity
C_ISO = 2 * μ_N * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] # Visco tensor for isotropic flow
F_inc = [1.0 dt * dudy; 0.0 1.0] # deformation matrix

function main(Stefan)
    # Initialize arrays to store the values for plots
    iterations = []
    τxy_values = []
    Ny_values  = []
    Nx_values  = []
    strain     = []

    μ_S = μ_N / aniso_factor # shear viscosity

    # Initializing variables 
    ω̇12   = 0.5 * (dudy - dvdx) # rotation rate  
    ε̇12   = 0.5 * (dudy + dvdx) # strain rate component
    L     = [dudx dudy; dvdx dvdy] # deformation gradient tensor
    W     = [0.0 ω̇12; -ω̇12 0.0] # rotation rate matrix
    D     = [0.0 ε̇12; ε̇12 0.0] # strain rate matrix
    D_vec = [dudx, dvdy, ε̇12] # vector of strain rate
    τ_vec = C_ISO * D_vec # stress vector
    τxy   = τ_vec[3] # shear stress component
    γ     = 0.0 # measure of strain
    F     = [1.0 0.0; dt * dvdx 1.0] # shear deformation matrix

    # Director initialization
    Director  = [1.0, 1.0]
    Norm_dir  = sqrt(Director[1]^2 + Director[2]^2)
    Director  = Director ./ Norm_dir
    Director2 = copy(Director)
    Norm_dir  = sqrt(Director[1]^2 + Director[2]^2)
    Nx        = Director[1]
    Ny        = Director[2]
    Angle     = atan(Nx, Ny) * 180 / π
    Fol_x     = -Ny / Nx
    Fol_y     = 1.0

    # Time-stepping loop
    it = 1
    while Ny > -0.99
        it += 1

        # in those steps we rotate the director 
        Director_dot  = W * Director - D * Director * (transpose(Director) * Director) + (transpose(Director) * (D * Director)) * Director
        Director_dot2 = -transpose(L) * Director2
        Director     += Director_dot * dt
        Director2    += Director_dot2 * dt

        # we need to normalise the director every time it is updated
        Norm_dir   = norm(Director)
        Director ./= Norm_dir

        # once we know the director we compute anisotropy matrix
        a0 = 2 * Director[1]^2 * Director[2]^2
        a1 = Director[1] * Director[2] * (-Director[1]^2 + Director[2]^2)

        # build the matrix 
        C_ANI = [-a0 a0 2*a1; a0 -a0 -2*a1; a1 -a1 -1+2*a0]

        if finite_strain
            F      = F_inc * F # incrementation of deformation every step
            Aratio = get_A_ratio(F)
            μ_S    = μ_N / Aratio
        end

        S_C   = C_ISO + 2 * (μ_N - μ_S) * C_ANI
        τ_vec = S_C * D_vec
        τxy   = τ_vec[3]
        γ    += ε̇12 * dt
        Nx    = Director[1]
        Ny    = Director[2]
        Angle = atan(Nx, Ny) * 180 / π
        Fol_x = -Ny / Nx
        Fol_y = 1.0

        push!(iterations, it)
        push!(τxy_values, τxy)
        push!(Ny_values, Ny)
        push!(Nx_values, Nx)
        push!(strain, γ)

        #println("Iteration: $it, τxy: $τxy, Ny: $Ny, Nx: $Nx")
    end

    # Find the maximum length of your data arrays for setting the x-axis limits
    max_length = maximum(length.(values(Stefan)))


    # Plotting the reference lines (Stefan data)
    p1 = Plots.plot(iterations, τxy_values, label="τxy", xlabel="Iteration", ylabel="τxy", seriestype = :line)
    Plots.plot!(p1, iterations, Stefan["Tauxy"][1:length(iterations)], seriestype = :scatter, markersize = 1, label="Stefan")
        
    p2 = Plots.plot(iterations,  Nx_values, label= "Nx", xlabel="Iteration", ylabel= "Nx", seriestype = :line)
    Plots.plot!(p2, iterations, Stefan["Nx"][1:length(iterations)], seriestype = :scatter, markersize = 1, label="Stefan")
        
    p3 = Plots.plot(iterations,  Ny_values, label= "Ny", xlabel="Iteration", ylabel= "Ny", seriestype = :line)
    Plots.plot!(p3, iterations, Stefan["Ny"][1:length(iterations)], seriestype = :scatter, markersize = 1, label="Stefan")
        
    display(Plots.plot(p1, p2, p3, layout=(1, 3), legend=:topright))

    
    #Plots.savefig("output_plot.png")

end

function get_A_ratio(F)
    B      = F * transpose(F) # Left Cauchy-Green tensor
    BV, BE = eigen(B) # Spectral decomposition
    VE     = Diagonal(sqrt.(abs.(BV))) # Square roots of eigenvalues
    FS     = BE * VE # Scale eigenvectors with eigenvalues
    psa1   = FS[:, 2] # Finite strain axis 1
    psa2   = FS[:, 1] # Finite strain axis 2   
    return norm(psa1) / norm(psa2) # the ratio of the long axis of the strain ellipsoid divided by the small axis of the ellipsoid
end


file = MAT.matopen("/Users/romankulakov/Downloads/StefanConstantAnisotropyTest.mat")
Stefan = Dict()
Stefan["Tauxy"] = read(file, "Tauxy")
Stefan["Nx"] = read(file, "Nx")
Stefan["Ny"] = read(file, "Ny")
close(file)

# Insert the first value at the beginning of each array
for key in keys(Stefan)
    Stefan[key] = vcat(Stefan[key][1], Stefan[key])
end

main(Stefan)

