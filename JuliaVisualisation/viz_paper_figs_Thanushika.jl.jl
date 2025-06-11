# Minimal frame work to make single nice visualisations in julia
using CairoMakie
using HDF5, Colors, ColorSchemes, MathTeXEngine, Printf

# Helper functions
function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function ExtractField(filename, field, size, mask_air, mask)
    field = try (Float64.(reshape(ExtractData( filename, field), size...)))
    catch 
        @warn "$field not found"
    end
    mask_air ? field[mask] .= NaN : nothing
    return field
end 

# Define main function
function run_main()

    # Set path
    path = "/home/thanushika/software/Results/Phase_Diagram_Model/New_Marker/WestPush/OceThick5/"

    # File numbers
    file_start = 0
    file_step  = 500
    file_end   = 1500

    # Initialize figure
    px_res   = 1000
    AR       = 3200.0 / 680.0       # Model aspect ratio
    fig_size = (px_res, cld(px_res, AR))
    fg1 = Figure(size = fig_size, figure_padding = 15)
    ax1 = Axis(
        fg1[1,1], 
        xlabel = L"$$x [km]", 
        ylabel = L"$$z [km]",
        aspect = AR
    )

    # File loop - This can very well be replace by just a single file
    for idx_file = file_start:file_step:file_end 
        
        # Construct file name and load data
        name  = @sprintf("Output%05d.gzip.h5", idx_file)
        filename = string(path, name)
        model = ExtractData( filename, "/Model/Params")
        t      = model[1]
        nvx    = Int(model[4])
        nvz    = Int(model[5])
        ncx, ncz = nvx-1, nvz-1
        xc    = ExtractData( filename, "/Model/xc_coord")
        zc    = ExtractData( filename, "/Model/zc_coord")

        # Mask the free air
        ph       = ExtractField(filename,  "/VizGrid/compo", (ncx, ncz), false, 0)
        mask_air = ph .== -1.00 
        T        = Float64.(reshape(ExtractData( filename, "/Centers/T"), ncx, ncz)) .- 273.15;    T[mask_air]   .= NaN

        # Fill Figure
        heatmap!(ax1, xc ./ 1e3, zc ./ 1e3, T)
        display(fg1)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# Run main
run_main();