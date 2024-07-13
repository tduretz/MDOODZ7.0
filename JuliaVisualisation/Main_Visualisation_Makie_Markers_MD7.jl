import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)
My = 1e6*365*24*3600

function main()

    # Set the path to your files"
    path = "/Users/tduretz/REPO/MDOODZ7.0/RUNS/1_NR07/"
    path = "/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"

    # File numbers
    file_start = 1
    file_step  = 1
    file_end   = 1

    # Switches
    resolution  = 500

    # Scaling
    Lc = 1000.
    tc = My
    Vc = 1e-9

    # Time loop
    for istep=file_start:file_step:file_end

        # Read grid file
        filename = string(path, @sprintf("Output%05d.gzip.h5", istep))
        model    = ExtractData( filename, "/Model/Params")
        xc       = ExtractData( filename, "/Model/xc_coord")
        zc       = ExtractData( filename, "/Model/zc_coord")
        xv       = ExtractData( filename, "/Model/xg_coord")
        zv       = ExtractData( filename, "/Model/zg_coord")
        xv_hr    = ExtractData( filename, "/VizGrid/xviz_hr")
        zv_hr    = ExtractData( filename, "/VizGrid/zviz_hr")
        xc_hr    = 0.5.*(xv_hr[1:end-1] .+ xv_hr[2:end])
        zc_hr    = 0.5.*(zv_hr[1:end-1] .+ zv_hr[2:end])
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)
    
        # Read marker file
        fmark    = string(path, @sprintf("Particles%05d.gzip.h5", istep))
        xm       = Float64.(ExtractData( fmark, "/Particles/x"));  
        zm       = Float64.(ExtractData( fmark, "/Particles/z"));  
        phm      = Float64.(ExtractData( fmark, "/Particles/phase"));
        gen      = Float64.(ExtractData( fmark, "/Particles/generation"));        
        model    = ExtractData( fmark, "/Model/Params")

        # Model properties
        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        Lx, Lz    = model[2], model[3]
        Δx, Δz, Δ = Lx/ncx, Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        @show "Model apect ratio" Lx/Lz
        @show "Model time" t/My

        # Load unreasonable amount of data
        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));          mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr)); 
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
        height  = Float64.(ExtractData( filename, "/Topo/z_grid")); 
        height_finer_c  = Float64.(ExtractData( filename, "/Topo/z_finer_cent"));
        Vx_grid = Float64.(ExtractData( filename, "/Topo/Vx_grid"));
        Vz_grid = Float64.(ExtractData( filename, "/Topo/Vz_grid"));  
        Vx_mark = Float64.(ExtractData( filename, "/Topo/Vx_mark"));
        Vz_mark = Float64.(ExtractData( filename, "/Topo/Vz_mark"));
        x_mark  = Float64.(ExtractData( filename, "/Topo/x_mark"));
        z_mark  = Float64.(ExtractData( filename, "/Topo/z_mark"));
        tag_s   = Float64.(reshape(ExtractData( filename, "/Flags/tag_s"), nvx, nvz))
        tag_t   = Float64.(reshape(ExtractData( filename, "/Flags/tag_t_fine"), 2*nvx-2, 2*nvz-2))
        tag_n   = Float64.(reshape(ExtractData( filename, "/Centers/BCp"), ncx, ncz))
        tag_u   = Float64.(reshape(ExtractData( filename, "/Flags/tag_u"), nvx, nvz+1))

        # Figure
        f = Figure(size = (Lx/Lz*resolution, resolution), fontsize=25)
        ax1 = Axis(f[1, 1], title = L"Phases at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
        # heatmap!(ax1, xc./Lc, zc./Lc, T, linewidth = 4, color=:white )  
        
        nvx2 = 2*nvx - 1
        nvz2 = 2*nvz - 1
        xv2  = LinRange(xv[1], xv[end], nvx2) 
        zv2  = LinRange(zv[1], zv[end], nvz2) 
        xc2  = 0.5*(xv2[1:end-1] .+ xv2[2:end])
        zc2  = 0.5*(zv2[1:end-1] .+ zv2[2:end])
        x_finer_c = 0.5*(xv2[2:end] .+ xv2[1:end-1])

        for i in eachindex(xv2)
            lines!(ax1, xv2[i]./Lc.*ones(size(zv2)), zv2./Lc, color=:gray)
        end
        for i in eachindex(zv2)
            lines!(ax1, xv2./Lc, zv2[i]./Lc.*ones(size(xv2)), color=:gray)
        end

        for i in eachindex(xv)
            lines!(ax1, xv[i]./Lc.*ones(size(zv)), zv./Lc, color=:black)
        end
        for i in eachindex(zv)
            lines!(ax1, xv./Lc, zv[i]./Lc.*ones(size(xv)), color=:black)
        end
        
        # scatter!(ax1, xm[phm.==0]./Lc, zm[phm.==0]./Lc, color=:blue)  
        # scatter!(ax1, xm[phm.==1]./Lc, zm[phm.==1]./Lc, color=:red) 
        # scatter!(ax1, xm[phm.==2]./Lc, zm[phm.==2]./Lc, color=:red)  
        # scatter!(ax1, xm[phm.==3]./Lc, zm[phm.==3]./Lc, color=:green) 
        # scatter!(ax1, xm[phm.==4]./Lc, zm[phm.==4]./Lc, color=:orange)
        # scatter!(ax1, xm[phm.==5]./Lc, zm[phm.==5]./Lc, color=:cyan)
        # scatter!(ax1, xm[phm.==6]./Lc, zm[phm.==6]./Lc, color=:black)

        # scatter!(ax1, [-2.273688e+03]./Lc, [-2.279599e+03]./Lc, color=:black, marker=:cross )

        # xx = xv[1] + (0+1)*Δx/2
        # zz = zv[1] + (190+1)*Δz/2
        # scatter!(ax1, [xx]./Lc, [zz]./Lc, color=:black, marker=:cross )

        scatter!(ax1, [xc2[128]]./Lc, [zc2[186]]./Lc, color=:black, marker=:cross )

        @show sum(gen.==1)
        scatter!(ax1, xm[gen.==1]./Lc, zm[gen.==1]./Lc, color=:black) 

        # lines!(ax1, xv./Lc, height./Lc)
        lines!(ax1, x_finer_c./Lc, height_finer_c./Lc)

        # Add vertices flags
        xv2d = xv .+ 0.0.*zv'
        zv2d = 0.0.*xv .+ zv'
        xv2f = xv2d[:]
        zv2f = zv2d[:]
        tag_sf = tag_s[:]
        scatter!(ax1, xv2f[tag_sf.==-1]./Lc, zv2f[tag_sf.==-1]./Lc, color=:black)
        
        # Add centroid flags
        xc2d = xc .+ 0.0.*zc'
        zc2d = 0.0.*xc .+ zc'
        xc2f = xc2d[:]
        zc2f = zc2d[:]
        tag_nf = tag_n[:]
        scatter!(ax1, xc2f[tag_nf.==-1]./Lc, zc2f[tag_nf.==-1]./Lc, color=:blue)

        # Vx flags
        zce  = [zc[1]-Δz/2; zc; zc[end]+Δz/2]
        xc2d = xv .+ 0.0.*zv'
        zc2d = 0.0.*zce .+ zce'
        xc2f = xc2d[:]
        zc2f = zc2d[:]
        tag_uf = tag_u[:]
        scatter!(ax1, xc2f[tag_uf.==-1]./Lc, zc2f[tag_uf.==-1]./Lc, color=:pink)
        
        # Add centroid flags from the twice finer mesh used for remeshing
        xc2d = xc2 .+ 0.0.*zc2'
        zc2d = 0.0.*xc2 .+ zc2'
        xc2f = xc2d[:]
        zc2f = zc2d[:]
        tag_tf = tag_t[:]
        scatter!(ax1, xc2f[tag_tf.==-1]./Lc, zc2f[tag_tf.==-1]./Lc, color=:blue, marker=:xcross)

        hidedecorations!(ax1)
        # scatter!(ax1, x_mark./Lc, z_mark./Lc)

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

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

main()