import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)

My = 1e6*365*24*3600

function main()
    # Switches
     T_contours  = false  # add temperature contours
    tc = My
    field = :Phases
    field = :StrainRate
    # field = :Stress
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 


    Lc = 1000.0


    # 10Ma strong crust
    fileNames = [
        "/Volumes/T7/d/iso/Output00640.gzip.h5", # aniso 1
        "/Volumes/T7/d/aniso1_5/Output01140.gzip.h5", # aniso 1.5
        "/Volumes/T7/d/aniso2/Output01170.gzip.h5", # aniso 2
        "/Volumes/T7/d/aniso3/Output01170.gzip.h5", # aniso 3
        "/Volumes/T7/d/aniso6/Output01180.gzip.h5", # aniso 6
        "/Volumes/T7/d/aniso9/Output01450.gzip.h5", # aniso 9
    ]

    # 20Ma strong crust
    fileNames = [
        "/Volumes/T7/d/iso/Output01790.gzip.h5", # aniso 1
        "/Volumes/T7/d/aniso1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/d/aniso2/Output02960.gzip.h5", # aniso 2
        "/Volumes/T7/d/aniso3/Output02960.gzip.h5", # aniso 3
        "/Volumes/T7/d/aniso6/Output03240.gzip.h5", # aniso 6
        "/Volumes/T7/d/aniso9/Output03240.gzip.h5", # aniso 9
    ]

    filelabels = [
        "isotropic (δ = 1)",
        "anisotropic (δ = 1.5)",
        "anisotropic (δ = 2)",
        "anisotropic (δ = 3)",
        "anisotropic (δ = 6)",
        "anisotropic (δ = 9)",
    ]

    f = Figure(resolution = (1200, 1800), fontsize = 16)

    for (i, filename) in enumerate(fileNames)
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
        Lz/ncz, sqrt( (Lx/ncx)^2 + (Lz/ncz)^2)
        @show "filename" filename
        @show "Model time" t/My
    
        ph    = Float64.(reshape(ExtractData( filename, "/VizGrid/compo"), ncx, ncz));          mask_air = ph .== -1.00 
        ph_hr = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_hr, ncz_hr)); 
        group_phases = copy( ph_hr); ph_hr[ph_hr.==-1.00] .= NaN
        T     = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN
      
        #####################################
        
        # Group phases for contouring
        group_phases[ ph_hr.==4 .|| ph_hr.==0 .|| ph_hr.==5  .|| ph_hr.==1  ] .= 0
        group_phases[ ph_hr.==2 .|| ph_hr.==6   ]                             .= 1
        group_phases[ ph_hr.==3 ]                                             .= 2
        
        #####################################

        filelabel = filelabels[i]

        if field==:Phases
            # Color palette for phase map
            cmap    = zeros(RGB{Float64}, 7)
            cmap[1] = RGBA{Float64}(255/255, 255/255, 255/255, 1.)  
            cmap[2] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
            cmap[3] = RGBA{Float64}(217/255, 099/255, 097/255, 1.)  
            cmap[4] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
            cmap[5] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
            cmap[6] = RGBA{Float64}(244/255, 218/255, 205/255, 1.) 
            cmap[7] = RGBA{Float64}(223/255, 233/255, 219/255, 1.) 
            phase_colors = cgrad(cmap, length(cmap), categorical=true, rev=false)
        
            # Define the depth limit for the upper layer
            upper_layer_depth = -35.0 * Lc # Assuming Lc is the length scaling factor (1000.0 in your script)
        
            # Filter the coordinate arrays to only include the upper layer
            zc_hr_upper = zc_hr[zc_hr .>= upper_layer_depth]
            zc_upper = zc[zc .>= upper_layer_depth]
        
            # Similarly filter the data arrays for properties to include only the upper layer
            T_upper = T[:, zc .>= upper_layer_depth]
            group_phases_upper = group_phases[:, zc_hr .>= upper_layer_depth]

            # Create a subplot at the determined position
            ax = Axis(f[i, 1], title = "Phases at t = $(@sprintf("%.2f", tMy)) Ma, $filelabel",
                    xlabel = "x [km]", ylabel = "y [km]", aspect = DataAspect())
                
            # Use the filtered data for the heatmap and contours
            heatmap!(ax, xc_hr./Lc, zc_hr_upper./Lc, group_phases_upper, colormap = phase_colors)
            if T_contours 
                Makie.update_theme!(reset=true)
                contour!(ax, xc./Lc, zc_upper./Lc, T_upper, levels=0:200:1400, linewidth = 4 ) 
            end
        end  
        
        if field==:StrainRate
            ε̇xx   = Float64.(reshape(ExtractData( filename, "/Centers/exxd"), ncx, ncz))
            ε̇xz   = Float64.(reshape(ExtractData( filename, "/Vertices/exz"), nvx, nvz))
            ε̇II   = sqrt.( 0.5*(2*ε̇xx.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) ); ε̇II[mask_air] .= NaN

            ax1 = Axis(f[i, 1], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma, %$(filelabel)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )   
        end

        if field==:Stress
            τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
            τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
            τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); 
            τII[mask_air] .= NaN
            

            ax1 = Axis(f[i, 1], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma, %$(filelabel)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo, α_heatmap)) 
            contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )   
        end
    end

    DataInspector(f)
    display(f)
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

main()