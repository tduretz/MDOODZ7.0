function strain_rate(velocity_field, compressible::Bool)
    nx, ny = size(velocity_field)
    strain_rate = zeros(Float64, (nx, ny, 2, 2))
    for i in 1:nx-1
        for j in 1:ny-1
            dx_u = velocity_field[i+1, j, 1] - velocity_field[i, j, 1]
            dy_u = velocity_field[i, j+1, 1] - velocity_field[i, j, 1]
            dx_v = velocity_field[i+1, j, 2] - velocity_field[i, j, 2]
            dy_v = velocity_field[i, j+1, 2] - velocity_field[i, j, 2]
            strain_rate[i, j, 1, 1] = 0.5*(dx_u + dx_u)
            strain_rate[i, j, 1, 2] = 0.5*(dy_u + dx_v)
            if compressible
                strain_rate[i, j, 2, 1] = 0.5*(dy_u + dx_v)
                strain_rate[i, j, 2, 2] = 0.5*(dy_v + dy_v)
            end
        end
    end
    return strain_rate
end
