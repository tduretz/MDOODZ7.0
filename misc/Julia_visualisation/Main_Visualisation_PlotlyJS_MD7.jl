import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, PlotlyJS, Printf, LazyGrids

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function main_PlotlyJS()

    cmy = 100*365.25*24*3600.0

    path ="/Users/imac/REPO_GIT/MDOODZ7.0/MDLIB/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD67_assembly_LR/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD7_assembly_LR/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiverTom/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x150/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x200/"

    # File numbers
    file_start = 405
    file_step  = 20
    file_end   = 405

    Lc = 1e3
    τc = 1e9
    Vc = 1.0/(cmy)

    zoom     = true
    fig_size = 2000
    field    = :phases
    # field = :density
    # field = :viscosity
    # field = :pl_srate
    # field = :pressure
    # field = :velocity_x
    # field = :velocity_z
    # field = :grain_size

    phase_contours = true

    for istep=file_start:file_step:file_end
    
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        # filename = string(path, @sprintf("Output_BeforeSolve%05d.gzip.h5", istep))
        model  = ExtractData( filename, "/Model/Params")
        xc     = ExtractData( filename, "/Model/xc_coord")
        zc     = ExtractData( filename, "/Model/zc_coord")
        xv     = ExtractData( filename, "/Model/xg_coord")
        zv     = ExtractData( filename, "/Model/zg_coord")
        xv_ph  = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_ph  = ExtractData( filename, "/VizGrid/zviz_hr")
        xc_ph  = 0.5.*(xv_ph[1:end-1] .+ xv_ph[2:end])
        zc_ph  = 0.5.*(zv_ph[1:end-1] .+ zv_ph[2:end])
        ncx_ph = length(xc_ph)
        ncz_ph = length(zc_ph)
        if zoom==false
            xmin, xmax = xc[1], xc[end]
            zmin, zmax = zc[1], zc[end]
        else
            xmin, xmax = -20e3, 20e3
            zmin, zmax = -20e3, 10e3
        end
        Lx, Lz = (xv[end]-xv[1])/Lc, (zv[end]-zv[1])/Lc
        Δx, Δz = xv[2]-xv[1],  zv[2]-zv[1]
        ###
        x_mark  = ExtractData( filename, "/Topo/x_mark");
        z_mark  = ExtractData( filename, "/Topo/z_mark");
        Vx_mark = ExtractData( filename, "/Topo/Vx_mark");
        Vz_mark = ExtractData( filename, "/Topo/Vz_mark");
        Vx_grid = ExtractData( filename, "/Topo/Vx_grid");
        Vz_grid = ExtractData( filename, "/Topo/Vz_grid");
        height  = ExtractData( filename, "/Topo/z_grid");
        xce      = [xc[1]-Δx; xc; xc[end]+Δx]; 
        zce      = [zc[1]-Δz; zc; zc[end]+Δz]; 
        xc2, zc2 = ndgrid(xc, zc)
        xvx2, zvx2 = ndgrid(xv, zce)
        xvz2, zvz2 = ndgrid(xce, zv)
        BCp   = ExtractData( filename, "/Centers/BCp");
        tag_s = ExtractData( filename, "/Flags/tag_s");
        tag_n = ExtractData( filename, "/Flags/tag_n");
        tag_u = ExtractData( filename, "/Flags/tag_u");
        tag_v = ExtractData( filename, "/Flags/tag_v");
        ###
        @show model
        t   = model[1]
        nvx = Int(model[4])
        nvz = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        ###
        ηc  = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz))
        ρc  = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz))
        P   = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz))
        d   = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz))
        ε̇pl = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz))
        ph  = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_ph, ncz_ph))
        Vx  = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz  = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))

        Vxc = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])

        if field==:phases      Xc = xc_ph; Zc = zc_ph; f = ph                end
        if field==:density     Xc = xc;    Zc = zc;    f = ρc                end
        if field==:viscosity   Xc = xc;    Zc = zc;    f = log10.(ηc)        end
        if field==:pl_srate    Xc = xc;    Zc = zc;    f = log10.(ε̇pl)       end
        if field==:pressure    Xc = xc;    Zc = zc;    f = P./τc             end
        if field==:grain_size  Xc = xc;    Zc = zc;    f = log10.(d)         end
        if field==:velocity_x  Xc = xv;    Zc = zc;    f = Vx[:,2:end-1]./Vc end
        if field==:velocity_z  Xc = xc;    Zc = zv;    f = Vz[2:end-1,:]./Vc end

        # if phase_contours  p = contour!(xc_ph./1e3, zc_ph./1e3, ph', c=:black) end
        # p = plot!(aspect_ratio=1)
        
        if @isdefined(cont)

        end

        @show minimum(Vx)
        @show maximum(Vx)
        @show maximum(Vx./Vc)

        p = Plot(
            [
            heatmap(x=Xc./Lc, y=Zc./Lc, z=f', colorscale=colors.turbo, colorbar_thickness=24), 
            # contour(x=xc_ph./Lc, y=zc_ph./Lc, z=ph', mode="lines", line_width=1, color="white", contours_coloring="none"),
            scatter(x=x_mark./Lc, y=z_mark./Lc, mode="markers"),
            scatter(x=xv./Lc, y=height./Lc, mode="markers"),
            # scatter(x=xc2[tag_n.==31]./Lc, y=zc2[tag_n.==31]./Lc, mode="markers"),
            # scatter(x=xc2[tag_n.==-1]./Lc, y=zc2[tag_n.==-1]./Lc, mode="markers"),
            scatter(x=xvx2[tag_u.==30]./Lc, y=zvx2[tag_u.==30]./Lc, mode="markers"),
            scatter(x=xvx2[tag_u.==-1]./Lc, y=zvx2[tag_u.==-1]./Lc, mode="markers"),
            heatmap(x=Xc./Lc, y=Zc./Lc, z=f', colorscale=colors.turbo, colorbar_thickness=24), 
            ],
            Layout(width=fig_size, height=fig_size*Lz/Lx,
            xaxis = attr(;
            range = [xmin, xmax]./Lc,
            constrain = "domain",
            ),
            yaxis = attr(;
                range = [zmin, zmax]./Lc,
                scaleanchor = "x",
                scaleratio = 1,
            )))
        
        display(p)
        sleep(0.5)
        @show P[1]

    end

end

main_PlotlyJS()