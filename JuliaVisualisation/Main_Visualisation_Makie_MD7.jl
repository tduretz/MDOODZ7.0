import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

My = 1e6*365*24*3600

function main()

    # path ="/Users/imac/REPO_GIT/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD67_assembly_LR/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiverTom/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/NR03_ALE0/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/1_NR09/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/NR09_miroir/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x150/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x200/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt_inc/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt_inc_eyy0/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_LR_noPTmatrix/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_LR_noPTmatrix_el_inc_eyy0_1e-13/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_eta_drop/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_Ebg5e-14/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x130/"

    # File numbers
    file_start = 100
    file_step  = 20
    file_end   = 100

    Lc = 1

    field = :phases
    # field = :density
    # field = :viscosity
    # field = :pl_srate
    # field = :stress
    # field = :strain_rate
    # field = :pressure
    # field = :temperature
    # field = :velocity_x
    # field = :velocity_z
    field = :grain_size

    printfig    = true
    ph_contours = true
    T_contours  = true

    for istep=file_start:file_step:file_end
    
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model  = ExtractData( filename, "/Model/Params")
        xc     = ExtractData( filename, "/Model/xc_coord")
        zc     = ExtractData( filename, "/Model/zc_coord")
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")
        xv_hr  = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr  = ExtractData( filename, "/VizGrid/zviz_hr")
        xc_hr  = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr  = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr = length(xc_hr)
        ncz_hr = length(zc_hr)

        t   = model[1]
        tMy = round(t/My, digits=2)
        nvx = Int(model[4])
        nvz = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        @show "Model length" Lx
        @show "Model time" t/My

        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));          mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr)); 
        group_phases = copy( ph_hr); ph_hr[ph_hr.==-1.00] .= NaN
        ηc    = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz));          ηc[mask_air]  .= NaN
        ρc    = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz));          ρc[mask_air]  .= NaN
        P     = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz));              P[mask_air]   .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
        d     = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz));              d[mask_air]   .= NaN
        ε̇pl   = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz));         ε̇pl[mask_air] .= NaN
        Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(2*ε̇xx.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
        Vxc = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])


        #####################################

        # Color palette for phase map
        cmap    = zeros(RGB{Float64}, 7)
        cmap[1] = RGBA{Float64}(  244/255,  218/255, 205/255, 1.)  
        cmap[2] = RGBA{Float64}(  217/255,  99/255, 97/255, 1.)  
        cmap[3] = RGBA{Float64}(117/255,164/255, 148/255, 1.) 
        cmap[4] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cmap[5] = RGBA{Float64}(217/255,  99/255, 97/255, 1.) 
        cmap[6] = RGBA{Float64}(244/255,  218/255, 205/255, 1.) 
        cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        #####################################

        f = Figure(resolution = ( Lx/Lz*1000, 1000), fontsize=25, aspect = 2.0)

        if field==:phases
            ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = "x [km]", ylabel = "y [km]")
            hm = heatmap!(ax1, xc_hr./Lc, zc_hr./Lc, ph_hr, colormap = phase_colors)
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = "Phases", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig save(path*"Phase"*@sprintf("%05d", istep)*".png", f, px_per_unit = 4) end
        end

        if field==:stress
            ax1 = Axis(f[1, 1], title = "τII", xlabel = "x [km]", ylabel = "y [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo,0.85))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end            
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = "τII", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig save(path*"Stress"*@sprintf("%05d", istep)*".png", f, px_per_unit = 4) end
        end


        if field==:strain_rate
            ax1 = Axis(f[1, 1], title = "ε̇II", xlabel = "x [km]", ylabel = "y [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, 0.85))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = "ε̇II", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig save(path*"StrainRate"*@sprintf("%05d", istep)*".png", f, px_per_unit = 4) end
        end

        if field==:grain_size
            ax1 = Axis(f[1, 1], title = L"$d$ at $t$ = %$(tMy) Ma", xlabel = "x [km]", ylabel = "y [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(d.*1e6), colormap = (:turbo, 0.85), colorrange=(1, 3))
            if T_contours 
                contour!(ax1, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
            end
            if ph_contours 
                contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )  
            end
            xminz, xmaxz = -0.4, 0.4
            zminz, zmaxz = -0.17, 0.17
            Lx = xmaxz - xminz
            Lz = zmaxz - zminz
            xlims!(ax1, -0.4, 0.4)
            ylims!(ax1, -0.17, 0.17)
            colsize!(f.layout, 1, Aspect(1, Lx/Lz))
            GLMakie.Colorbar(f[1, 2], hm, label = "d", width = 20, labelsize = 25, ticklabelsize = 14 )
            GLMakie.colgap!(f.layout, 20)
            if printfig save(path*"GrainSize"*@sprintf("%05d", istep)*".png", f, px_per_unit = 4) end
        end

        DataInspector(f)
        display(f)
        
    end

end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

main()