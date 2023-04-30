import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, Plots, Printf
gr()
plotlyjs()

function main()

    # path ="/Users/imac/REPO_GIT/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/MD67_assembly_LR/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiverTom/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x150/"
    # path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x200/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt_inc/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_el_gt_inc_eyy0/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_LR_noPTmatrix/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_LR_noPTmatrix_el_inc_eyy0_1e-13/"
    path ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/qcoe_x100_eta_drop/"

    # File numbers
    file_start = 960
    file_step  = 20
    file_end   = 960

    Lc = 1

    field = :phases
    # field = :density
    field = :viscosity
    # field = :pl_srate
    # field = :pressure
    # field = :temperature
    # field = :velocity_x
    # field = :velocity_z
    # field = :grain_size

    phase_contours = false

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

        @show model
        t   = model[1]
        nvx = Int(model[4])
        nvz = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        zmin, zmax = zv[1], zv[end]
        Lx, Lz = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc

        ηc  = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz))
        ρc  = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz))
        P   = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz))
        T   = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz))
        d   = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz))
        ε̇pl = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz))
        ph  = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_ph, ncz_ph))
        Vx  = Float64.(reshape(ExtractData( filename, "/VxNodes/Vx"), (ncx+1), (ncz+2)))
        Vz  = Float64.(reshape(ExtractData( filename, "/VzNodes/Vz"), (ncx+2), (ncz+1)))

        Vxc = 0.5 .* (Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
        Vzc = 0.5 .* (Vz[2:end-1,1:end-1] .+ Vz[2:end-1,2:end-0])

        p = plot()

        if field==:phases      p = heatmap!(xc_ph./Lc, zc_ph./Lc, ph', c=:turbo)   end
        if field==:density     p = heatmap!(xc./Lc, zc./Lc, (ρc'), c=:turbo, title="Density [kg/m3]") end
        if field==:viscosity   p = heatmap!(xc./Lc, zc./Lc, log10.(ηc'), c=:turbo) end
        if field==:pl_srate    p = heatmap!(xc./Lc, zc./Lc, log10.(ε̇pl'), c=:turbo) end
        if field==:pressure    p = heatmap!(xc./Lc, zc./Lc, (P'), c=:turbo) end
        if field==:temperature p = heatmap!(xc./Lc, zc./Lc, (T'.-273.15), c=:turbo) end
        if field==:grain_size  p = heatmap!(xc./Lc, zc./Lc, log10.(d'), c=:turbo) end
        if field==:velocity_x  p = heatmap!(xc./Lc, zc./Lc, Vxc', c=:turbo, clims=(-0.05, 0.05)) end
        if field==:velocity_z  p = heatmap!(xc./Lc, zc./Lc, Vzc', c=:turbo, clims=(-0.2, 0.1)) end

        if phase_contours  p = contour!(xc_ph./1e3, zc_ph./1e3, ph', c=:black) end
        p = plot!(aspect_ratio=Lx/Lz, xlims=(xmin, xmax), ylims=(zmin, zmax))
        display(p)
        sleep(0.1)
        @show P[1]
    end

end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

main()