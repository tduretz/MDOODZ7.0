function AddCountourQuivers!(PlotOnTop, ax1, coords, V, T, ϕ, σ1, ε̇1, PT, Fab, height, Lc, cm_y, group_phases, Δ, Mak)
    if PlotOnTop.quiver_origin
        xc2D = coords.c.x   .+ 0*coords.c.z'
        zc2D = 0*coords.c.x .+ coords.c.z'
        Mak.scatter!(ax1, xc2D[1:V.step:end,1:V.step:end][:]./Lc, zc2D[1:V.step:end,1:V.step:end][:]./Lc)
    end
    if PlotOnTop.ph_contours 
        Mak.contour!(ax1, coords.c_hr.x./Lc, coords.c_hr.z./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:black )  
    end
    if PlotOnTop.T_contours 
        Mak.contour!(ax1, coords.c.x./Lc, coords.c.z./Lc, T, levels=0:300:1250, linewidth = 4, color=:gray )  
    end  
    if PlotOnTop.PT_window 
        Mak.contour!(ax1, coords.c.x./Lc, coords.c.z./Lc, PT, linewidth = 4, color=:red )  
    end
    if PlotOnTop.fabric 
        Mak.arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, Fab.x[1:V.step:end,1:V.step:end], Fab.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5, linewidth=3)
    end 
    if PlotOnTop.σ1_axis
        Mak.arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, σ1.x[1:V.step:end,1:V.step:end], σ1.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5, color=:white, linewidth=5)
    end   
    if PlotOnTop.ε̇1_axis
        Mak.arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, ε̇1.x[1:V.step:end,1:V.step:end], ε̇1.z[1:V.step:end,1:V.step:end], arrowsize = 0, lengthscale=10Δ/1.5, color=:gray, linewidth=5)
    end
    if PlotOnTop.vel_vec
        Mak.arrows!(ax1, coords.c.x[1:V.step:end]./Lc, coords.c.z[1:V.step:end]./Lc, V.x[1:V.step:end,1:V.step:end]*cm_y, V.z[1:V.step:end,1:V.step:end]*cm_y, arrowsize = V.arrow, lengthscale = V.scale)
    end 
    if PlotOnTop.topo
        Mak.lines!(ax1, xv./Lc, height./Lc)
    end
    if PlotOnTop.ϕ_contours 
        Mak.contour!(ax1, coords.c.x./Lc, coords.c.z./Lc, ϕ, levels=0:1:0.1, linewidth = 6, color=:white, linestyle= :dashdot )  
    end 
end

function PrincipalStress(τxx, τzz, τxz, P)
    σ1   = (x=zeros(size(τxx)), z=zeros(size(τxx)) )
    τxzc = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end-0,1:end-1] .+ τxz[1:end-1,2:end-0] .+ τxz[2:end-0,2:end-0]) 
    for i in eachindex(τxzc)
        if P[i]>-1e-13
            σ = [-P[i]+τxx[i] τxzc[i]; τxzc[i] -P[i]+τzz[i]]
            v = eigvecs(σ)
            σ1.x[i] = v[1,1]
            σ1.z[i] = v[2,1]
        end
    end
    return σ1
end
