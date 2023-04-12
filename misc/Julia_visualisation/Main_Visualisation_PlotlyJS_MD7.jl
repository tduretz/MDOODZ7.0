using HDF5, PlotlyJS, Printf

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function main_PlotlyJS()

    path ="/Users/imac/REPO_GIT/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD67_assembly_LR/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD7_assembly_LR/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiverTom/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x150/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x200/"

    # File numbers
    file_start = 1
    file_step  = 10
    file_end   = 1

    Lc = 1

    # field = :phases
    field = :density
    field = :viscosity
    # field = :pl_srate
    field = :pressure
    # field = :velocity_x
    # field = :velocity_z
    # field = :grain_size

    phase_contours = true

    for istep=file_start:file_step:file_end
    
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
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
        xmin, xmax = xc[1], xc[end]
        zmin, zmax = zc[1], zc[end]

        @show model
        t   = model[1]
        nvx = Int(model[4])
        nvz = Int(model[5])
        ncx, ncz = nvx-1, nvz-1

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

        # p = plot()

        # if field==:phases      p = heatmap!(xc_ph./Lc, zc_ph./Lc, ph', c=:turbo)   end
        # if field==:density     p = heatmap!(xc./Lc, zc./Lc, (ρc'), c=:turbo) end
        # if field==:viscosity   p = heatmap!(xc./Lc, zc./Lc, log10.(ηc'), c=:turbo) end
        # if field==:pl_srate    p = heatmap!(xc./Lc, zc./Lc, log10.(ε̇pl'), c=:turbo) end
        # if field==:pressure    p = heatmap!(xc./Lc, zc./Lc, (P'), c=:turbo) end
        # if field==:grain_size  p = heatmap!(xc./Lc, zc./Lc, log10.(d'), c=:turbo) end
        # if field==:velocity_x  p = heatmap!(xc./Lc, zc./Lc, Vxc', c=:turbo) end
        # if field==:velocity_z  p = heatmap!(xc./Lc, zc./Lc, Vzc', c=:turbo) end

        # if phase_contours  p = contour!(xc_ph./1e3, zc_ph./1e3, ph', c=:black) end
        # p = plot!(aspect_ratio=1)

        p = Plot([PlotlyJS.heatmap(x=xc_ph./Lc, y=zc_ph./Lc, z=P', colorscale=colors.turbo, colorbar_thickness=24), 
        #PlotlyJS.contour(x=xc_ph./1e3, y=zc_ph./1e3, z=ph',
                #line_width=2, line_color="black")
                ],
        PlotlyJS.Layout(xaxis_range=[-3, 3], yaxis_range=[zmin/Lc, zmax/Lc], height=400, yaxis_scaleanchor="x", yaxis_scaleratio=1))
        
        display(p)
        sleep(0.1)
        @show P[1]

    end

end

main_PlotlyJS()