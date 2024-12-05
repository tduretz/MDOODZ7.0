import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
using CairoMakie, GLMakie
Mak = CairoMakie
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

const y    = 365*24*3600
const My   = 1e6*y
const cm_y = y*100.

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 
`d`
function main()

    # Set the path to your files
    path ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"

    # File numbers
    istep = 0

    # Select field to visualise
    field = :Residual_x
    # field = :Divergence

    # Switches
    printfig    = false  # print figures to disk
    printvid    = false
    framerate   = 3
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 
    vel_arrow   = 5
    vel_scale   = 300000
    nap         = 0.1    # pause for animation 
    resol       = 500
    mov_name    = "$(path)/_$(field)/$(field)"  # Name of the movie
    Lx, Lz      = 1.0, 1.0

    # Scaling
    # Lc = 1000.
    # tc = My
    # Vc = 1e-9

    Lc = 1.0
    tc = My
    Vc = 1.0

    probe = (ϕeff = Float64.([]), t  = Float64.([]))

    # Time loop
    f = Figure(size = (Lx/Lz*resol*1.2, resol), fontsize=25)

    name     = @sprintf("Output%05d.gzip.h5", 1)
    filename = string(path, name)
    model    = ExtractData( filename, "/Model/Params")
    xc       = ExtractData( filename, "/Model/xc_coord")
    zc       = ExtractData( filename, "/Model/zc_coord")
    xv       = ExtractData( filename, "/Model/xg_coord")
    zv       = ExtractData( filename, "/Model/zg_coord")
    xvz      = ExtractData( filename, "/Model/xvz_coord")
    zvx      = ExtractData( filename, "/Model/zvx_coord")
    t      = model[1]
    tMy    = round(t/tc, digits=6)
    nvx    = Int(model[4])
    nvz    = Int(model[5])
    ncx, ncz = nvx-1, nvz-1
    centroids, vertices = (ncx, ncz), (nvx, nvz)

    name     = @sprintf("Residuals%05d.gzip.h5", istep)
    @info "Reading $(name)"
    filename = string(path, name)
     
    Vx    = Float64.(reshape(ExtractData( filename, "/VxNodes/ru"), (ncx+1), (ncz+2)))
    Vz    = Float64.(reshape(ExtractData( filename, "/VzNodes/rv"), (ncx+2), (ncz+1)))
    divu  = ExtractField(filename, "/Centers/rp", centroids, false, 0)

    # #####################################
    empty!(f)
    f = Figure(fontsize=25) #size = (Lx/Lz*resol*1.2, resol), 

    if field==:Residual_x
        ax1 = Axis(f[1, 1], title = L"$Vx$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
        hm = heatmap!(ax1, xv, zvx, Vx[:,:], colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label = L"$Vx$ [cm.yr$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
        # Mak.colgap!(f.layout, 20)
        if printfig Print2Disk( f, path, string(field), istep) end
    end

    if field==:Residual_z
        ax2 = Axis(f[1, 1], title = L"$Vz$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
        hm = heatmap!(ax2, xvz, zv, Vz, colormap = (:jet, α_heatmap))#, colorrange=(0., 0.6)
        # colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label = L"$Vz$ [cm.yr$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
        # Mak.colgap!(f.layout, 20)
        if printfig Print2Disk( f, path, string(field), istep) end
    end 

    if field==:Divergence
        ax3 = Axis(f[1, 1], title = L"∇⋅V at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
        hm = heatmap!(ax3, xc./Lc, zc./Lc, divu, colormap = (:turbo, α_heatmap))
        colsize!(f.layout, 1, Aspect(1, Lx/Lz))
        Mak.Colorbar(f[1, 2], hm, label =  L"∇⋅V [s$^{-1}$]", width = 20, labelsize = 25, ticklabelsize = 14 )
        Mak.colgap!(f.layout, 20)
        if printfig Print2Disk( f, path, string(field), istep) end
    end

    DataInspector(f)
    display(f)
    sleep(nap)

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