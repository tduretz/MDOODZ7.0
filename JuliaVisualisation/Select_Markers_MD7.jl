import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
using CairoMakie#, GLMakie
using DataFramesMeta, CSV

const Mak = CairoMakie
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(true)
# Unit conversion
sec_year = 365*24*3600
Myear    = 1e6*sec_year
cm_year  = sec_year*100.

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 

function AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_year, group_phases, Δ)
    if PlotOnTop.T_contours 
        contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:black )  
    end 
    if PlotOnTop.ph_contours 
        contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
    end
    if PlotOnTop.fabric 
        arrows!(ax1, xc./Lc, zc./Lc, Fab.x, Fab.z, arrowsize = 0, lengthscale=Δ/1.5)
    end 
    if PlotOnTop.σ1_axis
        arrows!(ax1, xc./Lc, zc./Lc, σ1.x, σ1.z, arrowsize = 0, lengthscale=Δ/1.5, color=:white)
    end   
    if PlotOnTop.vel_vec
        arrows!(ax1, xc[1:V.step:end]./Lc, zc[1:V.step:end]./Lc, V.x[1:V.step:end,1:V.step:end]*cm_year, V.z[1:V.step:end,1:V.step:end]*cm_year, arrowsize = V.arrow, lengthscale = V.scale)
    end 
    if PlotOnTop.topo
        lines!(ax1, xv./Lc, height./Lc)
    end
end

function main()

    # Set the path to your files
    path = "/Users/lcandiot/Developer/MDOODZ7.0/cmake-exec/RiftingChenin/"
    # path ="/home/thanushika/software/Results/Phase_Diagram_Model/Marker/Oth5km/"

    # File numbers
    file_start = 00
    file_step  = 10
    file_end   = 200

    file_select_start = 0                             # First time step - Select marker coordinates from here
    file_select_end   = 200                             # Last time step - Select marker coordinates from here
    target_coords = (-50e3, -25e3, -15e3, -10e3)        # Coordinate bounds (min(x), max(x), min(y), max(y))

    # Select field to visualise
    field = :Phases
    # field = :Cohesion
    #  field = :Density
    # field = :Viscosity 
    # field = :PlasticStrainrate
    # field = :Stress
    # field = :StrainRate
    #  field = :Pressure
    # field = :Divergence
    # field = :Temperature
    # field = :Velocity_x
    # field = :Velocity_z
    # field = :Velocity
    # field = :GrainSize
    # field = :Topography
    # field = :TimeSeries
    # field = :AnisotropyFactor
    # field = :MeltFraction
    # field = :TimeSeries
    # field = :EffectiveFrictionTime


    # Switches
    printfig    = false  # print figures to disk
    printvid    = false
    framerate   = 3
    PlotOnTop = (
        ph_contours = false,  # add phase contours
        T_contours  = false,   # add temperature contours
        fabric      = false,  # add fabric quiver (normal to director)
        topo        = false,
        σ1_axis     = false,
        vel_vec     = false,
    )
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 
    vel_arrow   = 5
    vel_scale   = 1000000
    vel_step    = 5
    nap         = 0.1    # pause for animation 
    resol       = 1000
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    Lc = 1.0
    tc = Myear
    Vc = 1.0

    probe = (ϕeff = Float64.([]), t  = Float64.([]))
    cm_yearr = 100.0*3600.0*24.0*365.25

    # Select markers in last time step
    fmark    = string(path, @sprintf("Particles%05d.gzip.h5", file_select_end))
    xm       = Float64.(ExtractData( fmark, "/Particles/x"));  
    zm       = Float64.(ExtractData( fmark, "/Particles/z"));  
    phm      = Float64.(ExtractData( fmark, "/Particles/phase"));
    gen      = Float64.(ExtractData( fmark, "/Particles/generation"));        
    model    = ExtractData( fmark, "/Model/Params")
    idx_markers = []
    for idx in eachindex(xm)
        if xm[idx] >= target_coords[1] && xm[idx] <= target_coords[2] && zm[idx] >= target_coords[3] && zm[idx] <= target_coords[4]
            push!(idx_markers, idx)
        end
    end

    # Delete markers that have not been there from the start
    fmark      = string(path, @sprintf("Particles%05d.gzip.h5", file_select_start))
    xm         = Float64.(ExtractData( fmark, "/Particles/x"))
    idx_delete = findall(x -> x > length(xm), idx_markers)
    deleteat!(idx_markers, idx_delete)

    # Select markers and parameters for tracking
    IDs    = idx_markers[1:10]
    params = ("T [K]", "P [Pa]", "ph [Int]", "t [s]", "x [m]", "z [m]")
    store_types = Dict(
        "Time steps" => cld((file_end - file_start), file_step) + 1,
        "Markers"    => size(IDs, 1)
    )
    store_key = "Markers"

    # Initialize storage arrays
    data = if store_key == "Time steps"
        zeros(Float64, store_types[store_key], size(params, 1))
    elseif store_key == "Markers"
        (
            [Vector{Float64}(undef, store_types["Time steps"]) for _ in eachindex(IDs)],
            [Vector{Float64}(undef, store_types["Time steps"]) for _ in eachindex(IDs)],
            [0                                                 for _ in eachindex(IDs)],
            [Vector{Float64}(undef, store_types["Time steps"]) for _ in eachindex(IDs)],
            [Vector{Float64}(undef, store_types["Time steps"]) for _ in eachindex(IDs)],
            [Vector{Float64}(undef, store_types["Time steps"]) for _ in eachindex(IDs)]
        )
    end

    # Time loop
    ts_counter = 0
    for istep=file_start:file_step:file_end
        
        # Update time step counter
        ts_counter += 1

        # Figure
        f = Figure(size = (500, 500), fontsize=25)

        # Load Output file
        name     = @sprintf("Output%05d.gzip.h5", istep)
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
        P_t     = ExtractData( filename, "TimeSeries/P_mean_time")
        τxz_t    = ExtractData( filename, "TimeSeries/sxz_mean_time")
        t_t      = ExtractData( filename, "TimeSeries/Time_time")

        xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

        # Load marker file
        fmark    = string(path, @sprintf("Particles%05d.gzip.h5", istep))
        xm       = Float64.(ExtractData( fmark, "/Particles/x"));  
        zm       = Float64.(ExtractData( fmark, "/Particles/z"));  
        Tm       = Float64.(ExtractData( fmark, "/Particles/T"));  
        Pm       = Float64.(ExtractData( fmark, "/Particles/P"));  
        phm      = Float64.(ExtractData( fmark, "/Particles/phase"));
        gen      = Float64.(ExtractData( fmark, "/Particles/generation"));        
        model    = ExtractData( fmark, "/Model/Params")

        # Store marker data
        if store_key == "Time steps"
            data[ts_counter, 1] = Tm[IDs]
            data[ts_counter, 2] = Pm[IDs]
            data[ts_counter, 3] = phm[IDs]
            data[ts_counter, 4] = model[1]
        elseif store_key == "Markers"
            [data[1][Idx][ts_counter] = Tm[IDs[Idx]]  for Idx in eachindex(IDs)]
            [data[2][Idx][ts_counter] = Pm[IDs[Idx]]  for Idx in eachindex(IDs)]
            [data[3][Idx]             = phm[IDs[Idx]] for Idx in eachindex(IDs)]
            [data[4][Idx][ts_counter] = model[1]      for Idx in eachindex(IDs)]
            [data[5][Idx][ts_counter] = xm[IDs[Idx]]  for Idx in eachindex(IDs)]
            [data[6][Idx][ts_counter] = zm[IDs[Idx]]  for Idx in eachindex(IDs)]
        end

        if istep == file_end
            df = DataFrame()
            [df[!, param] = data[idx] for (idx, param) in enumerate(params)]
            @show df
            CSV.write("$(path)MarkerData.csv", df)
        end

        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        Lx>1e3 ? length_unit="km" :  length_unit="m"
        @info "Model info"
        @show "Model apect ratio" Lx/Lz
        @show "Model time" t/Myear
        @show length_unit
        centroids, vertices = (ncx, ncz), (nvx, nvz)

        ph           = ExtractField(filename,  "/VizGrid/compo", centroids, false, 0)
        mask_air     = ph .== -1.00 
        ph_hr        = ExtractField(filename,  "/VizGrid/compo_hr", (ncx_hr, ncz_hr), false, 0)
        ph_dual_hr   = ExtractField(filename,  "/VizGrid/compo_dual_hr", (ncx_hr, ncz_hr), false, 0)
        group_phases = copy(ph_hr); ph_hr[ph_hr.==-1.00] .= NaN
        ηc           = ExtractField(filename,  "/Centers/eta_n", centroids, true, mask_air)
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
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
        τxzc  = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end,1:end-1] .+ τxz[1:end-1,2:end] .+ τxz[2:end,2:end]) 
        C     = Float64.(reshape(ExtractData( filename, "/Centers/cohesion"), ncx, ncz))
        ϕ     = ExtractField(filename, "/Centers/phi", centroids, false, 0)
        divu  = ExtractField(filename, "/Centers/divu", centroids, false, 0)

        Fab = 0.
        if PlotOnTop.fabric
            δani  = ExtractField(filename, "/Centers/ani_fac", centroids, false, 0)
            Nx    = Float64.(reshape(ExtractData( filename, "/Centers/nx"), ncx, ncz))
            Nz    = Float64.(reshape(ExtractData( filename, "/Centers/nz"), ncx, ncz))
            Fab   = (x=-Nz./Nx, z=ones(size(Nz)))
            nrm   = sqrt.(Fab.x.^2 .+ Fab.z.^2)
            Fabx ./= nrm
            Fabz ./= nrm
        end
        height = 0.
        if PlotOnTop.topo
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
        V = (x=Vxc, z=Vzc, arrow=vel_arrow, step=vel_step, scale=vel_scale)
        σ1 = 0.
        if PlotOnTop.σ1_axis 
            σ1 = PrincipalStress(τxx, τzz, τxz, P) 
        end
        
        #####################################

        # Color palette for phase map
        cmap    = zeros(RGB{Float64}, 7)
        cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
        cmap[2] = RGBA{Float64}(207/255, 089/255, 087/255, 1.)  
        cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
        cmap[4] = RGBA{Float64}(213/255, 213/255, 209/255, 1.) 
        cmap[5] = RGBA{Float64}(117/255, 099/255, 097/255, 1.) 
        cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
        cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)

        # # Group phases for contouring
        # group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        # group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        # group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################
        empty!(f)
        f = Figure(size = (Lx/Lz*resol*1.2, resol), fontsize=25)

        if field==:Phases
            # Plot phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            # Plot markers on top
            sc = scatter!(ax1, xm[idx_markers[1:1:end]] ./ Lc, zm[idx_markers[1:1:end]] ./ Lc, color = :black, markersize = 4.0)
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_year, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end
        

        if field==:Viscosity
            ax1 = Axis(f[1, 1], title = L"$\eta$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ηc), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\eta$ [Pa.s]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Density
            ax1 = Axis(f[1, 1], title = L"$\rho$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ρc, colormap = (:turbo, α_heatmap))  
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\rho$ [kg.m$^{-3}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Stress
            ax1 = Axis(f[1, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII./1e6, colormap = (:turbo, α_heatmap)) 
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\tau_\textrm{II}$ [MPa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Pressure
            ax1 = Axis(f[1, 1], title = L"$P$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, P./1e9, colormap = (:turbo, α_heatmap)) #, colorrange=(1,1.2)1e4*365*24*3600
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label =  L"$P$ [GPa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Divergence
            ax1 = Axis(f[1, 1], title = L"∇⋅V at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, divu, colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:StrainRate
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [%$(length_unit)]", ylabel = L"$y$ [%$(length_unit)]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label =  L"$\dot{\varepsilon}_\textrm{II}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:PlasticStrainrate
            ε̇pl[ε̇pl.==0.0] .= 1e-30
            ax1 = Axis(f[1, 1], title = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇pl), colormap = (:turbo, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$\dot{\varepsilon}_\textrm{II}^\textrm{pl}$ [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity
            ax1 = Axis(f[1, 1], title = L"$V$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            Vx_BG = 0*xc .- 2*zc'  
            V     = sqrt.( (Vxc .- 0.0*Vx_BG).^2 + (Vzc).^2)
            hm = heatmap!(ax1, xc./Lc, zc./Lc, V, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            # xlims!(ax1, 0., 3.e-3)
            # ylims!(ax1, 0., 3.e-3)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$V$ [m.s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_x
            ax1 = Axis(f[1, 1], title = L"$Vx$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xv, zvx[2:end-1], Vx[2:end-1,:]*cm_yearr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vx$ [cm.yr$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Velocity_z
            ax1 = Axis(f[1, 1], title = L"$Vz$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xvz[2:end-1], zv, Vz[:,2:end-1]*cm_yearr, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = L"$Vz$ [cm.yr$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:GrainSize
            ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, α_heatmap), colorrange=(1, 3))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            Mak.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:AnisotropyFactor
            ax1 = Axis(f[1, 1], title = L"$δ_\textrm{ani}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, δani, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            # xminz, xmaxz = -0.4, 0.4
            # zminz, zmaxz = -0.17, 0.17
            # Lx = xmaxz - xminz
            # Lz = zmaxz - zminz
            # xlims!(ax1, -0.4, 0.4)
            # ylims!(ax1, -0.17, 0.17)
            Mak.Colorbar(f[1, 2], hm, label = L"$δ_\textrm{ani}$", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:MeltFraction
            ax1 = Axis(f[1, 1], title = L"ϕ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, ϕ, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            Mak.Colorbar(f[1, 2], hm, label = L"$ϕ$", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Cohesion
            ax1 = Axis(f[1, 1], title = L"$C$ [MPa] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, C./1e6, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            Mak.Colorbar(f[1, 2], hm, label = L"$C$ [MPa]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Temperature
            ax1 = Axis(f[1, 1], title = L"$T$ [C] at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, T, colormap = (:bilbao, α_heatmap))
            AddCountourQuivers!(PlotOnTop, ax1, xc, xv, zc, V, T, σ1, Fab, height, Lc, cm_yearear, group_phases, Δ)                
            Mak.Colorbar(f[1, 2], hm, label = L"$T$ [C]", width = 20, labelsize = 25, ticklabelsize = 14 )
            Mak.colgap!(f.layout, 20)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            if printfig Print2Disk( f, path, string(field), istep) end
        end

        if field==:Topography
            ax1 = Axis(f[1, 1], title = L"Topography at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$h$ [km]")
            @show mean(height)
            @show mean(z_mark)
            lines!(ax1, xv./Lc, height./Lc)
            scatter!(ax1, x_mark./Lc, z_mark./Lc)

            ax2 = Axis(f[2, 1], xlabel = L"$x$ [km]", ylabel = L"$Vx$ [km]")
            lines!(ax2, xv./Lc, Vx_grid./Vc)
            scatter!(ax2, x_mark./Lc, Vx_mark./Vc)

            ax3 = Axis(f[3, 1], xlabel = L"$x$ [km]", ylabel = L"$Vz$ [km]")
            lines!(ax3, xc./Lc, Vz_grid[2:end-1]./Vc)
            scatter!(ax3, x_mark/Lc, Vz_mark./Vc)
        end

        if field==:TimeSeries
            ax1 = Axis(f[1, 1], title = L"$τ_{xz}$", xlabel = L"$t$", ylabel = L"$\tau_{xz}$")
            τxz_t[1] = 0.
            @show .-τxz_t
            lines!(ax1, t_t, .-τxz_t./(.-P_t.+τzz_t))
        end

        if field==:EffectiveFrictionTime
            ϕ_eff = -mean(τxzc[:,end] ./ (τzz[:,end] .- P[:,end]) )
            if istep==0 ϕ_eff = 0. end
            push!(probe.ϕeff, ϕ_eff) 
            push!(probe.t, t) 
            if istep==file_end
                ax1 = Axis(f[1, 1], title = L"$ϕ_\mathrm{eff}$", xlabel = L"$t$", ylabel = L"$ϕ_\mathrm{eff}$")
                lines!(ax1, Float64.(probe.t), Float64.(probe.ϕeff))
            end

        end

        if field!=:EffectiveFrictionTime || istep!=file_end
            DataInspector(f)
            display(f)
            sleep(nap)
        end

    end

    yscale = Lz/Lx

    if printfig && printvid
        FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -pattern_type glob -i $(path)_$(field)/'*'.png -vf "scale=1080:1080*$(yscale)" -c:v libx264 -pix_fmt yuv420p -y "$(mov_name).mov"`)
    end

end

function PrincipalStress(τxx, τzz, τxz, P)
    σ1   = (x=zeros(size(τxx)), z=zeros(size(τxx)) )
    τxzc = 0.25*(τxz[1:end-1,1:end-1] .+ τxz[2:end-0,1:end-1] .+ τxz[1:end-1,2:end-0] .+ τxz[2:end-0,2:end-0]) 
    for i in eachindex(τxzc)
        if P[i]>1e-13
            σ = [-P[i]+τxx[i] τxzc[i]; τxzc[i] -P[i]+τzz[i]]
            v = eigvecs(σ)
            σ1.x[i] = v[1,1]
            σ1.z[i] = v[2,1]
        end
    end
    return σ1
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, path, field, istep; res=4)
     path1 = path*"/_$field/"
     mkpath(path1)
     save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()