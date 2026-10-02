import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using JuliaVisualisation
using HDF5, Printf, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG, Statistics
using CairoMakie, GLMakie
using CSV, Tables

Mak = GLMakie # or CairoMakie
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

const y = 365 * 24 * 3600
const My = 1e6 * y
const cm_y = y * 100.

@views function main()
    # Define paths to the two output files
    path1 = "/home/lilou/librairies/MDOODZ7.0/save_runs/SS_inclusion_n1_Pidam2_noEvol_20261001_182359/"
    path2 = "/home/lilou/librairies/MDOODZ7.0/save_runs/SS_inclusion_n1_Pidam2_evolaniso_20261001_181918/"

    # File numbers (assuming same time step for both files)
    file_start = 0
    file_step = 5
    file_end = 10

    # Scaling
    Lc = 1000
    tc = My
    Vc = 1.0
    τc = 1e6

    # Time loop
    for istep = file_start:file_step:file_end
        # Load data from the first file
        name1 = @sprintf("Output%05d.gzip.h5", istep)
        filename1 = string(path1, name1)
        model1 = ExtractData(filename1, "/Model/Params")
        xc1 = ExtractData(filename1, "/Model/xc_coord")
        zc1 = ExtractData(filename1, "/Model/zc_coord")
        nvx1    = Int(model1[4])
        nvz1    = Int(model1[5])
        ncx1, ncz1 = nvx1-1, nvz1-1


        τxx   = Float64.(reshape(ExtractData( filename1, "/Centers/sxxd"), ncx1, ncz1))
        τzz   = Float64.(reshape(ExtractData( filename1, "/Centers/szzd"), ncx1, ncz1))
        τxz   = Float64.(reshape(ExtractData( filename1, "/Vertices/sxz"), nvx1, nvz1))
        ε̇xx   = Float64.(reshape(ExtractData( filename1, "/Centers/exxd"), ncx1, ncz1))
        ε̇xz   = Float64.(reshape(ExtractData( filename1, "/Vertices/exz"), nvx1, nvz1))
        ε̇zz   = Float64.(reshape(ExtractData( filename1, "/Centers/ezzd"), ncx1, ncz1))
        τyy   = -(τzz .+ τxx)
        ε̇yy   = -(ε̇xx .+ ε̇zz)
        τII_1  = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); 
        ε̇II_1   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) );

        # Load data from the second file
        name2 = @sprintf("Output%05d.gzip.h5", istep)
        filename2 = string(path2, name2)
        model2 = ExtractData(filename2, "/Model/Params")
        xc2 = ExtractData(filename2, "/Model/xc_coord")
        zc2 = ExtractData(filename2, "/Model/zc_coord")
        nvx2    = Int(model2[4])
        nvz2    = Int(model2[5])
        ncx2, ncz2 = nvx2-1, nvz2-1


        τxx   = Float64.(reshape(ExtractData( filename2, "/Centers/sxxd"), ncx2, ncz2))
        τzz   = Float64.(reshape(ExtractData( filename2, "/Centers/szzd"), ncx2, ncz2))
        τxz   = Float64.(reshape(ExtractData( filename2, "/Vertices/sxz"), nvx2, nvz2))
        ε̇xx   = Float64.(reshape(ExtractData( filename2, "/Centers/exxd"), ncx2, ncz2))
        ε̇xz   = Float64.(reshape(ExtractData( filename2, "/Vertices/exz"), nvx2, nvz2))
        ε̇zz   = Float64.(reshape(ExtractData( filename2, "/Centers/ezzd"), ncx2, ncz2))
        τyy   = -(τzz .+ τxx)
        ε̇yy   = -(ε̇xx .+ ε̇zz)
        τII_2  = sqrt.( 0.5*(τxx.^2 .+ τyy.^2 .+ τzz.^2 .+ 0.5*(τxz[1:end-1,1:end-1].^2 .+ τxz[2:end,1:end-1].^2 .+ τxz[1:end-1,2:end].^2 .+ τxz[2:end,2:end].^2 ) ) ); 
        ε̇II_2   = sqrt.( 0.5*(ε̇xx.^2 .+ ε̇yy.^2 .+ ε̇zz.^2 .+ 0.5*(ε̇xz[1:end-1,1:end-1].^2 .+ ε̇xz[2:end,1:end-1].^2 .+ ε̇xz[1:end-1,2:end].^2 .+ ε̇xz[2:end,2:end].^2 ) ) );


        # Ensure the grids are the same
        if xc1 != xc2 || zc1 != zc2
            error("Grids from the two files do not match!")
        end

        # Compute differences
        ΔτII = τII_1 - τII_2
        Δε̇II = ε̇II_1 - ε̇II_2

       max_abs_ΔτII = maximum(abs.(ΔτII[.!isnan.(ΔτII)]))
       max_abs_Δε̇II = maximum(abs.(Δε̇II[.!isnan.(Δε̇II)]))

        # Create a figure with 2 rows and 3 columns
        f = Figure(size = (1500, 1000), fontsize=30)

        # Plot original stress field from file 1
        ax1 = Axis(f[1, 1],
                   title = L"$\tau_{II}$ reference",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm1 = heatmap!(ax1, xc1, zc1, τII_1, colormap = (:reds, 1.0), colorrange=(0.6,2))
        Colorbar(f[1, 2], hm1, label = L"$\tau_{II}$ [Pa]", height = Relative(0.5))

        # Plot original stress field from file 2
        ax2 = Axis(f[1, 3],
                   title = L"$\tau_{II}$ with anisoEvol",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm2 = heatmap!(ax2, xc1, zc1, τII_2, colormap = (:reds, 1.0), colorrange=(0.6,2))
        Colorbar(f[1, 4], hm2, label = L"$\tau_{II}$ [Pa]", height = Relative(0.5))

        # Plot stress difference
        ax3 = Axis(f[1, 5],
                   title = L"$\Delta \tau_{II}$",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm3 = heatmap!(ax3, xc1, zc1, ΔτII, colormap = (:bluesreds, 1.0), colorrange=(-max_abs_ΔτII,max_abs_ΔτII))
        Colorbar(f[1, 6], hm3, label = L"$\Delta \tau_{II}$ [Pa]", height = Relative(0.5))

        # Plot original strain-rate field from file 1
        ax4 = Axis(f[2, 1],
                   title = L"$\dot{\varepsilon}_{II}$ reference",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm4 = heatmap!(ax4, xc1, zc1, ε̇II_1, colormap = (:lapaz, 1.0), colorrange=(1.0,3))
        Colorbar(f[2, 2], hm4, label = L"$\dot{\varepsilon}_{II}$ [s$^{-1}$]", height = Relative(0.5))

        # Plot original strain-rate field from file 2
        ax5 = Axis(f[2, 3],
                   title = L"$\dot{\varepsilon}_{II}$ with anisoEvol",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm5 = heatmap!(ax5, xc1, zc1, ε̇II_2, colormap = (:lapaz, 1.0), colorrange=(1.0,3))
        Colorbar(f[2, 4], hm5, label = L"$\dot{\varepsilon}_{II}$ [s$^{-1}$]", height = Relative(0.5))

        # Plot strain-rate difference
        ax6 = Axis(f[2, 5],
                   title = L"$\Delta \dot{\varepsilon}_{II}$ ",
                   xlabel = L"$x$ [km]",
                   ylabel = L"$y$ [km]",
                   aspect = 1)
        hm6 = heatmap!(ax6, xc1, zc1, Δε̇II, colormap = (:bluesreds, 1.0), colorrange=(-max_abs_Δε̇II,max_abs_Δε̇II) )
        Colorbar(f[2, 6], hm6, label = L"$\Delta \dot{\varepsilon}_{II}$ [s$^{-1}$]", height = Relative(0.5))

        # Save the figure
      #   save("compare_stress_strainRate.png", f)
        save("/home/lilou/Desktop/presentations-rapports/pres_sete_2026/figures/compare_stress_strainRate.png", f) 

        display(f)
    end
end

main()