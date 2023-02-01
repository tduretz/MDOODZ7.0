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


function strain_rate_2d_staggered(velocity_field_u, velocity_field_v, dx, dy)
    nx, ny = size(velocity_field_u)
    strain_rate = zeros(Float64, (nx, ny, 2, 2))
    for i in 1:nx-1
        for j in 1:ny-1
            du_dx = (velocity_field_u[i+1, j] - velocity_field_u[i, j]) / dx
            dv_dy = (velocity_field_v[i, j+1] - velocity_field_v[i, j]) / dy
            strain_rate[i, j, 1, 1] = du_dx
            strain_rate[i, j, 2, 2] = dv_dy
            strain_rate[i, j, 1, 2] = 0.5 * ((velocity_field_v[i+1, j] - velocity_field_v[i, j]) / dx + (velocity_field_u[i, j+1] - velocity_field_u[i, j]) / dy)
        end
    end
    return strain_rate
end


function calc_strain_rate(velocity_field_u, velocity_field_v, dx, dy, nt)
    nx, ny = size(velocity_field_u[:, :, 1])
    strain_rate = zeros(Float64, (nx, ny, 2, 2, nt))
    for t in 1:nt
        strain_rate[:, :, :, :, t] = strain_rate_2d_staggered(velocity_field_u[:, :, t], velocity_field_v[:, :, t], dx, dy)
    end
    return strain_rate
end


function solve_viscous_problem_explicit(velocity_field_u, velocity_field_v, dx, dy, dt, viscosity, nt)
    nx, ny = size(velocity_field_u[:, :, 1])
    for t in 1:nt
        strain_rate = strain_rate_2d_staggered(velocity_field_u[:, :, t], velocity_field_v[:, :, t], dx, dy)
        for i in 1:nx-1
            for j in 1:ny-1
                velocity_field_u[i, j, t+1] = velocity_field_u[i, j, t] + dt * viscosity * (strain_rate[i, j, 1, 1] + strain_rate[i+1, j, 1, 1]) / 2
                velocity_field_v[i, j, t+1] = velocity_field_v[i, j, t] + dt * viscosity * (strain_rate[i, j, 2, 2] + strain_rate[i, j+1, 2, 2]) / 2
            end
        end
    end
    return velocity_field_u, velocity_field_v
end

function solve_viscous_problem_implicit(velocity_field_u, velocity_field_v, dx, dy, dt, viscosity, nt)
    nx, ny = size(velocity_field_u[:, :, 1])
    for t in 1:nt
        strain_rate = strain_rate_2d_staggered(velocity_field_u[:, :, t], velocity_field_v[:, :, t], dx, dy)
        for i in 1:nx-1
            for j in 1:ny-1
                a_uu = viscosity * (strain_rate[i, j, 1, 1] + strain_rate[i+1, j, 1, 1]) / 2
                a_vv = viscosity * (strain_rate[i, j, 2, 2] + strain_rate[i, j+1, 2, 2]) / 2
                velocity_field_u[i, j, t+1] = (velocity_field_u[i, j, t] + dt * a_uu) / (1 + dt * viscosity)
                velocity_field_v[i, j, t+1] = (velocity_field_v[i, j, t] + dt * a_vv) / (1 + dt * viscosity)
            end
        end
    end
    return velocity_field_u, velocity_field_v
end


function strain_rate_3d(velocity_field, compressible::Bool)
    nx, ny, nz = size(velocity_field)
    strain_rate = zeros(Float64, (nx, ny, nz, 3, 3))
    for i in 1:nx-1
        for j in 1:ny-1
            for k in 1:nz-1
                dx_u = velocity_field[i+1, j, k, 1] - velocity_field[i, j, k, 1]
                dy_u = velocity_field[i, j+1, k, 1] - velocity_field[i, j, k, 1]
                dz_u = velocity_field[i, j, k+1, 1] - velocity_field[i, j, k, 1]
                dx_v = velocity_field[i+1, j, k, 2] - velocity_field[i, j, k, 2]
                dy_v = velocity_field[i, j+1, k, 2] - velocity_field[i, j, k, 2]
                dz_v = velocity_field[i, j, k+1, 2] - velocity_field[i, j, k, 2]
                dx_w = velocity_field[i+1, j, k, 3] - velocity_field[i, j, k, 3]
                dy_w = velocity_field[i, j+1, k, 3] - velocity_field[i, j, k, 3]
                dz_w = velocity_field[i, j, k+1, 3] - velocity_field[i, j, k, 3]
                strain_rate[i, j, k, 1, 1] = 0.5*(dx_u + dx_u)
                strain_rate[i, j, k, 1, 2] = 0.5*(dy_u + dx_v)
                strain_rate[i, j, k, 1, 3] = 0.5*(dz_u + dx_w)
                strain_rate[i, j, k, 2, 1] = 0.5*(dy_u + dx_v)
                strain_rate[i, j, k, 2, 2] = 0.5*(dy_v + dy_v)
                strain_rate[i, j, k, 2, 3] = 0.5*(dz_v + dy_w)
                if compressible
                    strain_rate[i, j, k, 3, 1] = 0.5*(dz_u + dx_w)
                    strain_rate[i, j, k, 3, 2] = 0.5*(dz_v + dy_w)
                    strain_rate[i, j, k, 3, 3] = 0.5*(dz_w + dz_w)
                end
            end
        end
    end
    return strain_rate
end
