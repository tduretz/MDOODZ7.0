import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(false)

My = 1e6*365*24*3600

function main()
    # Switches
     T_contours  = true  # add temperature contours
    tc = My
    field = :Phases
    #field = :StrainRate
    # field = :Stress
    α_heatmap   = 1.0 #0.85   # transparency of heatmap 


    Lc = 1000.0

     # 15Ma 
     fileNames = [
        "/Volumes/T7/dgh/transition1/Output01710.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dgh/transition2/Output01710.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output01880.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output02080.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition9/Output02430.gzip.h5", # aniso 9
        "/Volumes/T7/dgh/transition12/Output02450.gzip.h5", # aniso 12
        "/Volumes/T7/dgh/transition15/Output02680.gzip.h5", # aniso 15
    ]

    # 10Ma 
    fileNames = [
        "/Volumes/T7/dgh/transition1/Output01060.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dgh/transition2/Output01110.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output01190.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output01350.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition9/Output01370.gzip.h5", # aniso 9
        "/Volumes/T7/dgh/transition12/Output01370.gzip.h5", # aniso 12
        "/Volumes/T7/dgh/transition15/Output01440.gzip.h5", # aniso 15
    ]

    # 20Ma 
    fileNames = [
        "/Volumes/T7/dgh/transition1/Output02730.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dgh/transition2/Output02550.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output02920.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output03170.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition9/Output03610.gzip.h5", # aniso 9
        "/Volumes/T7/dgh/transition12/Output03600.gzip.h5", # aniso 12
        "/Volumes/T7/dgh/transition15/Output03900.gzip.h5", # aniso 15
    ]

    #
    #
    # aniso angle 85
    #
    # 10Ma 
    fileNames = [
        "/Volumes/T7/dghh/transition1/Output01060.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dghh/transition2/Output01090.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output01200.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output01290.gzip.h5", # aniso 6
        #"/Volumes/T7/dghh/transition9/Output01370.gzip.h5", # aniso 9
        "/Volumes/T7/dghh/transition12/Output01220.gzip.h5", # aniso 12
        "/Volumes/T7/dghh/transition15/Output01240.gzip.h5", # aniso 15
    ]

    # 15Ma 
    fileNames = [
        "/Volumes/T7/dghh/transition1/Output01710.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dghh/transition2/Output01670.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output01940.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output02230.gzip.h5", # aniso 6
        #"/Volumes/T7/dgh/transition9/Output02430.gzip.h5", # aniso 9
        "/Volumes/T7/dghh/transition12/Output02250.gzip.h5", # aniso 12
        #"/Volumes/T7/dghh/transition15/Output02680.gzip.h5", # aniso 15
    ]

     # 20Ma 
     fileNames = [
        "/Volumes/T7/dghh/transition1/Output02710.gzip.h5", # aniso 1
        # "/Volumes/T7/dgh/transition1_5/Output02660.gzip.h5", # aniso 1.5
        "/Volumes/T7/dghh/transition2/Output02560.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output02900.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output03220.gzip.h5", # aniso 6
        #"/Volumes/T7/dghh/transition9/Output03610.gzip.h5", # aniso 9
        "/Volumes/T7/dghh/transition12/Output03300.gzip.h5", # aniso 12
        #"/Volumes/T7/dghh/transition15/Output03900.gzip.h5", # aniso 15
    ]

    filelabels = [
        "isotropic (δ = 1, θ = 80)",
        # "anisotropic (δ = 1.5, θ = 80)",
        "anisotropic (δ = 2, θ = 80)",
        "anisotropic (δ = 3, θ = 80)",
        "anisotropic (δ = 6, θ = 80)",
        #"anisotropic (δ = 9, θ = 80)",
        "anisotropic (δ = 12, θ = 80)",
        #"anisotropic (δ = 15, θ = 80)",
    ]


    #
    #
    # compare angles
    #
    #

    # 10Ma 
    fileNames = [
        "/Volumes/T7/dgh/transition1/Output01060.gzip.h5", # aniso 1
        "/Volumes/T7/dgh/transition2/Output01110.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output01190.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output01350.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition12/Output01370.gzip.h5", # aniso 12

        "/Volumes/T7/dghh/transition1/Output01060.gzip.h5", # aniso 1
        "/Volumes/T7/dghh/transition2/Output01090.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output01200.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output01290.gzip.h5", # aniso 6
        "/Volumes/T7/dghh/transition12/Output01220.gzip.h5", # aniso 12
    ]    

     # 15Ma 
     fileNames = [
        "/Volumes/T7/dgh/transition1/Output01710.gzip.h5", # aniso 1
        "/Volumes/T7/dgh/transition2/Output01710.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output01880.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output02080.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition12/Output02450.gzip.h5", # aniso 12

        "/Volumes/T7/dghh/transition1/Output01710.gzip.h5", # aniso 1
        "/Volumes/T7/dghh/transition2/Output01670.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output01940.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output02230.gzip.h5", # aniso 6
        "/Volumes/T7/dghh/transition12/Output02250.gzip.h5", # aniso 12
    ]

    # 20Ma 
    fileNames = [
        "/Volumes/T7/dgh/transition1/Output02730.gzip.h5", # aniso 1
        "/Volumes/T7/dgh/transition2/Output02550.gzip.h5", # aniso 2
        "/Volumes/T7/dgh/transition3/Output02920.gzip.h5", # aniso 3
        "/Volumes/T7/dgh/transition6/Output03170.gzip.h5", # aniso 6
        "/Volumes/T7/dgh/transition12/Output03600.gzip.h5", # aniso 12

        "/Volumes/T7/dghh/transition1/Output02710.gzip.h5", # aniso 1
        "/Volumes/T7/dghh/transition2/Output02560.gzip.h5", # aniso 2
        "/Volumes/T7/dghh/transition3/Output02900.gzip.h5", # aniso 3
        "/Volumes/T7/dghh/transition6/Output03220.gzip.h5", # aniso 6
        "/Volumes/T7/dghh/transition12/Output03300.gzip.h5", # aniso 12
    ]

    filelabels = [
        "isotropic (δ = 1, θ = 80)",
        "anisotropic (δ = 2, θ = 80)",
        "anisotropic (δ = 3, θ = 80)",
        "anisotropic (δ = 6, θ = 80)",
        "anisotropic (δ = 12, θ = 80)",
        "isotropic (δ = 1, θ = 85)",
        "anisotropic (δ = 2, θ = 85)",
        "anisotropic (δ = 3, θ = 85)",
        "anisotropic (δ = 6, θ = 85)",
        "anisotropic (δ = 12, θ = 85)",
    ]

    #
    #
    # compare moho temperature
    #
    #


    fileNames = [
        "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3\\Output01160.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3\\Output01870.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3\\Output02710.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_50\\Output00970.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_50\\Output01690.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_50\\Output02670.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_20\\Output01170.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_20\\Output01830.gzip.h5", # aniso 1
    "C:\\Users\\rkulakov\\CLionProjects\\MDOODZ7.0\\cmake-exec\\NeckingReview\\aniso3_20\\Output02730.gzip.h5", # aniso 1
    ]

    filelabels = [
        "anisotropic (δ = 3, θ = 85, T moho = 500)",
        "anisotropic (δ = 3, θ = 85, T moho = 500)",
        "anisotropic (δ = 3, θ = 85, T moho = 500)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
        "anisotropic (δ = 3, θ = 85, T moho = 430)",
    ]


    f = Figure(resolution = (1800, 1800), fontsize = 16)

    if field==:Phases
        f = Figure(resolution = (1800, 800), fontsize = 20)
    end


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

        column = 1
        row = i

        if (i > 3)
            column = 2
            row = i - 3
        end


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
            ax = Axis(f[row, column], title = "Phases at t = $(@sprintf("%.2f", tMy)) Ma, $filelabel",
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

            ax1 = Axis(f[row, column], title = L"$\dot{\varepsilon}_\textrm{II}$ at $t$ = %$(tMy) Ma, %$(filelabel)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            hm = heatmap!(ax1, xc./Lc, zc./Lc, log10.(ε̇II), colormap = (:turbo, α_heatmap))
            contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )   
        end

        if field==:Stress
            τxx   = Float64.(reshape(ExtractData( filename, "/Centers/sxxd"), ncx, ncz))
            τxz   = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
            τII   = sqrt.( 0.5*(2*τxx.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); 
            τII[mask_air] .= NaN
            

            ax1 = Axis(f[row, column], title = L"$\tau_\textrm{II}$ at $t$ = %$(tMy) Ma, %$(filelabel)", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
            heatmap!(ax1, xc./Lc, zc./Lc, τII, colormap = (:turbo, α_heatmap)) 
            contour!(ax1, xc_hr./Lc, zc_hr./Lc, group_phases, levels=-1:1:maximum(group_phases), linewidth = 4, color=:white )   
        end
    end

    DataInspector(f)
    Print2Disk( f, string(field))
    display(f)
    
end

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function Print2Disk( f, field; res=4)
    save("/Volumes/T7/strain.png", f, px_per_unit = res) 
end

main()