using CairoMakie

function main_func()
    fname = "/home/thanushika/software/MDOODZ7.0/IMPORT/Hydrated_Pdt.dat"
    #            = 793;                     // Resolution for temperature (MANTLE) []
    # model.PDMnP[pid]         = 793;                     // Resolution for pressure    (MANTLE) []
    # model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
    # model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
    # model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
    # model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
    # model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Hydrated_Pdt.dat", model.PDMnT[pid], model.PDMnP

    nT     = 793
    nP     = 793
    rho    = Array{Float64}(undef, nT*nP)
    rho_2D = read!(fname, rho) |> x -> reshape(x, nT, nP)

    T = collect(LinRange(373.0, 1373.0, nT))
    P = collect(LinRange(10.13e6, 5.5816e9, nP))

    fig = Figure(size = (500, 500), fontsize = 18)
    ax1 = Axis(
        fig[1, 1][1, 1],
        aspect = 1.0,
        xlabel = L"$T$ [°C]",
        ylabel = L"$P$ [kbar]"
    )
    hm = heatmap!(ax1, T .- 273.15, P ./ 1e8, rho_2D)
    Colorbar(fig[1,1][1,2], hm, label = L"$\rho$ [kg.m$^{-1}$]")
    rowsize!(fig.layout, 1, Fixed(300))
    display(fig)
end

main_func();
save("heatmap_output.png", fig)
