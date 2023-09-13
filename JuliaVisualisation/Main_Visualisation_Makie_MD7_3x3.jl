import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

My = 1e6*365*24*3600

function main()
    path ="/Users/romankulakov/CLionProjects/MDOODZ70/cmake-exec/CollisionPolarCartesian/vis/"

    # File numbers
    # files_list = [280, 590, 960, 620, 1290, 2030, 950, 1960, 3070]  
    files_list = [10, 20, 30, 40, 50, 60, 220, 300, 380] 

    # Switches
    T_contours  = true  # add temperature contours
    α_heatmap   = 0.85   # transparency of heatmap 
    nap         = 0.3    # pause for animation 

    # Scaling
    Lc = 1000.
    tc = My

    # Figure
    f = Figure(fontsize=25)  # Adjust figure size for a 3x3 grid of plots

    # File loop
    for (plot_i, istep) in enumerate(files_list)
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
        ncx_hr, ncz_hr = length(xc_hr), length(zc_hr)

        t      = model[1]
        tMy    = round(t/tc, digits=6)
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xmin, xmax = xv[1], xv[end]
        Lx, Lz    = (xmax-xmin)/Lc, (zv[end]-zv[1])/Lc
        @show "Model apect ratio" Lx/Lz
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
        τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
        τzz   = Float64.(reshape(ExtractData( filename, "/Centers/szzd"), ncx, ncz))
        τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
        ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
        ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
        τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); τII[mask_air] .= NaN
        ε̇II   = sqrt.( 0.5*(2*ε̇xx.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN
        
        #####################################

        # Color palette for phase map
        cmap    = zeros(RGB{Float64}, 7)
        cmap[1] = RGBA{Float64}(210/255, 218/255, 205/255, 1.)  
        cmap[2] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
        cmap[3] = RGBA{Float64}(117/255, 164/255, 148/255, 1.) 
        cmap[4] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cmap[5] = RGBA{Float64}(217/255, 099/255, 097/255, 1.) 
        cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
        cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
        cgrad(cmap, length(cmap), categorical=true, rev=false)

        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 3

        # Plotting
        ax = Axis(f[div(plot_i-1, 3)+1, mod(plot_i-1, 3)+1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")  # Place the plot in the appropriate position in the grid
        hm = heatmap!(ax, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
        if T_contours 
            contour!(ax, xc./Lc, zc./Lc, T, levels=0:200:1400, linewidth = 4, color=:white )  
        end

        DataInspector(f)
        display(f)
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