import CairoMakie

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    CairoMakie.save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

function ReadFile(path, step, scales, options, PlotOnTop)

    name     = @sprintf("Output%05d.gzip.h5", step)
    @info "Reading $(name)"
    filename = string(path, name)
    model    = ExtractData( filename, "/Model/Params")
    xc       = ExtractData( filename, "/Model/xc_coord")
    zc       = ExtractData( filename, "/Model/zc_coord")
    xv       = ExtractData( filename, "/Model/xg_coord")
    zv       = ExtractData( filename, "/Model/zg_coord")
    xvz      = ExtractData( filename, "/Model/xvz_coord")
    zvx      = ExtractData( filename, "/Model/zvx_coord")
    xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr")
    zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr")
    τzz_t    = ExtractData( filename, "TimeSeries/szzd_mean_time")
    P_t      = ExtractData( filename, "TimeSeries/P_mean_time")
    τxz_t    = ExtractData( filename, "TimeSeries/sxz_mean_time")
    t_t      = ExtractData( filename, "TimeSeries/Time_time")

    xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
    zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
    ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

    t      = model[1]
    tMy    = round(t/scales.tc, digits=6)
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1
    xmin, xmax = xv[1], xv[end]
    zmin, zmax = zv[1], zv[end]
    Lx, Lz    = (xmax-xmin)/scales.Lc, (zv[end]-zv[1])/scales.Lc
    Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
    (xmax-xmin)>1e3 ? length_unit="km" :  length_unit="m"

    # Overwrite coordinates array
    xc = LinRange(xmin+Δx/2, xmax-Δx/2, ncx)
    zc = LinRange(zmin+Δz/2, zmax-Δz/2, ncz)
    xv = LinRange(xmin, xmax, nvx)
    zv = LinRange(zmin, zmax, nvz)
    coords = ( c=(x=xc, z=zc), c_hr=(x=xc_hr, z=zc_hr), v=(x=xv, z=zv))

    @info "Model info"
    @show "Model apect ratio" Lx/Lz
    @show "Model time" t/My
    @show length_unit
    centroids, vertices = (ncx, ncz), (nvx, nvz)

    ph           = ExtractField(filename,  "/VizGrid/compo", centroids, false, 0)
    mask_air     = ph .== -1.00 
    ph_hr        = ExtractField(filename,  "/VizGrid/compo_hr", (ncx_hr, ncz_hr), false, 0)
    ph_dual_hr   = ExtractField(filename,  "/VizGrid/compo_dual_hr", (ncx_hr, ncz_hr), false, 0)
    group_phases = copy(ph_hr); ph_hr[ph_hr.==-1.00] .=   NaN; ph_dual_hr[ph_dual_hr.==-1.00] .=   NaN;
    ηc           = ExtractField(filename,  "/Centers/eta_n", centroids, true, mask_air)
    ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]   .= NaN
    P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]    .= NaN
    T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]    .= NaN
    d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]    .= NaN
    ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air]  .= NaN
    ε̇el   = Float64.(reshape(ExtractData( filename, "/Centers/eII_el"), ncx, ncz));         ε̇el[mask_air]  .= NaN
    ε̇pwl  = Float64.(reshape(ExtractData( filename, "/Centers/eII_pwl"), ncx, ncz));        ε̇pwl[mask_air] .= NaN
    ε̇lin  = Float64.(reshape(ExtractData( filename, "/Centers/eII_lin"), ncx, ncz));        ε̇lin[mask_air] .= NaN

    Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
    Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
    τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
    τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz))
    τyy   = -(τzz .+ τxx)
    τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
    ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
    ε̇zz   = Float64.(reshape(ExtractData( filename, "/Centers/ezzd"), ncx, ncz))
    ε̇yy   = -(ε̇xx .+ ε̇zz)
    ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
    τII   = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
    ε̇II   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
    εII_pl= Float64.(reshape(ExtractData( filename, "/Centers/strain_pl"), ncx, ncz))
    εII_pwl= Float64.(reshape(ExtractData( filename, "/Centers/strain_pwl"), ncx, ncz))
    εII   = εII_pl + εII_pwl
    τxzc  = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end,1:end-1] .+ τxz[1:end-1,2:end] .+ τxz[2:end,2:end]) 
    C     = Float64.(reshape(ExtractData( filename, "/Centers/cohesion"), ncx, ncz))
    ϕ     = ExtractField(filename, "/Centers/phi", centroids, false, 0)
    divu  = ExtractField(filename, "/Centers/divu", centroids, false, 0)
    εII_el   = Float64.(reshape(ExtractData( filename, "/Centers/strain_el"), ncx, ncz)); εII_el[mask_air]  .= NaN
    εII_lin  = Float64.(reshape(ExtractData( filename, "/Centers/strain_lin"), ncx, ncz)); εII_lin[mask_air]  .= NaN
    εII_pwl  = Float64.(reshape(ExtractData( filename, "/Centers/strain_pwl"), ncx, ncz)); εII_pwl[mask_air]  .= NaN
    εII_pl   = Float64.(reshape(ExtractData( filename, "/Centers/strain_pl"), ncx, ncz)); εII_pl[mask_air]  .= NaN
    εII = εII_pwl + εII_pl
    # εII   = Float64.(reshape(ExtractData( filename, "/Centers/strain"), ncx, ncz)); εII[mask_air]  .= NaN

    T_hr  = zeros(size(ph_hr)); T_hr[1:2:end-1,1:2:end-1] .= T; T_hr[2:2:end-0,2:2:end-0] .= T
    if options.LAB_color
        ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.<options.LAB_T] .= 2
        ph_hr[(ph_hr.==2 .|| ph_hr.==3) .&& T_hr.>options.LAB_T] .= 3
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.<options.LAB_T] .= 2
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==3) .&& T_hr.>options.LAB_T] .= 3
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.<options.LAB_T] .= 2
        ph_dual_hr[(ph_dual_hr.==2 .|| ph_dual_hr.==6) .&& T_hr.>options.LAB_T] .= 3
    end

    Fab = 0.
    if PlotOnTop.fabric
        δani    = ExtractField(filename, "/Centers/ani_fac", centroids, false, 0)
        Nx      = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
        Nz      = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
        Fab     = (x=-Nz./Nx, z=ones(size(Nz)))
        nrm     = sqrt.(Fab.x.^2 .+ Fab.z.^2)
        Fab.x ./= nrm
        Fab.z ./= nrm
    end
    height = 0.

    read_topo = true
    try Float64.(ExtractData( filename, "/Topo/z_grid"));
    catch 
        @warn "no topo data"
        read_topo = false
    end

    if read_topo
        height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
        Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
        Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
        Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
        Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
        x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
        z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
    end
    Vxc   = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
    Vzc   = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])
    V = (x=Vxc, z=Vzc, arrow=options.vel_arrow, step=options.vel_step, scale=options.vel_scale)
    σ1, ε̇1, PT = 0., 0., 0.
    if PlotOnTop.σ1_axis 
        σ1 = PrincipalStress(τxx, τzz, τxz, P) 
    end
    if PlotOnTop.ε̇1_axis 
        ε̇1 = PrincipalStress(ε̇xx, ε̇zz, ε̇xz, -1/3*divu) 
        @show extrema(ε̇1.x)
    end

    return (tMy=tMy,length_unit=length_unit, Lx=Lx, Lz=Lz, xc=xc, zc=zc, ε̇II=(ε̇II=ε̇II, ε̇el=ε̇el, ε̇lin=ε̇lin, ε̇pl=ε̇pl, ε̇pwl=ε̇pwl), τII=τII,
    coords=coords, V=V, T=T, ϕ=ϕ, σ1=σ1, ε̇1=ε̇1, PT=PT, Fab=Fab, height=height, all_phases=(ph=ph, ph_hr=ph_hr, ph_dual_hr=ph_dual_hr, group_phases=group_phases), εII=εII, C=C, Δ=Δ)
end
