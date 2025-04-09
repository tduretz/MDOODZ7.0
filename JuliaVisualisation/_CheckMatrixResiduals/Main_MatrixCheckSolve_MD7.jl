using HDF5, Plots, SparseArrays, LinearAlgebra, Statistics
plotlyjs()

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function LoadSparseCompressedMatrix( filename, block )
    Ic   = ExtractData(filename, "/matrix/I$(block)")
    J    = ExtractData(filename, "/matrix/J$(block)")
    V    = ExtractData(filename, "/matrix/V$(block)")
    Ic .+= 1 
    J  .+= 1 
    return Ic, J, V
end

function AssembleSparseMartrix(Ic, J, V)
    # Decompress index I
    I    = zeros(size(J))
    m    = 0
    I[1] = 1
    p    = 1
    for k=2:length(Ic)        
        off = Ic[k] - Ic[k-1]
        for l=1:off
            I[m+l] = p
        end
        m = m + off
        p = p + 1
    end
    return sparse(I, J, V)
end


function main()

    solve = false
    γ     = 2.1e7

    # File
    path     ="/Users/tduretz/REPO/MDOODZ7.0/MDLIB/"
    filename = string(path, "Stokes_01cpu_step01_iter00.gzip.h5") 
    # path     ="/Users/tduretz/REPO/MDOODZ7.0/RUNS/RiverTom/"
    # filename = string(path, "Stokes_01cpu_step229_iter00.gzip.h5") 
    # Matrix block A
    IcA, JA, VA = LoadSparseCompressedMatrix( filename, "A" )
    MA          = AssembleSparseMartrix(IcA, JA, VA)
    # Matrix block B
    IcB, JB, VB = LoadSparseCompressedMatrix( filename, "B" )
    MB          = AssembleSparseMartrix(IcB, JB, VB)
    # Matrix block C
    IcC, JC, VC = LoadSparseCompressedMatrix( filename, "C" )
    MC          = AssembleSparseMartrix(IcC, JC, VC)
    # Matrix block D
    MD  = spdiagm( γ.*ones(size(MC,1)) )
    # RHS
    ru    = ExtractData(filename, "/matrix/rhs_mom")
    rp    = ExtractData(filename, "/matrix/rhs_cont")
    # Equation numbers
    equ   = ExtractData(filename,"/numbering/eqn_u");
    eqv   = ExtractData(filename,"/numbering/eqn_v");
    eqp   = ExtractData(filename,"/numbering/eqn_p");
    # Grid flags
    tag_s = ExtractData(filename,"/fields/tag_s");
    tag_n = ExtractData(filename,"/fields/tag_n");
    tag_u = ExtractData(filename,"/fields/tag_u");
    tag_v = ExtractData(filename,"/fields/tag_v");
    # Viscosity fields
    eta_s = ExtractData(filename,"/fields/eta_s");
    eta_n = ExtractData(filename,"/fields/eta_n");
    # Model parameters
    model = ExtractData(filename,"/model/params");
    nx, nz, Δx, Δz = model

    @show size(MA)
    @show size(MB)
    @show size(MC)
    @show size(MD)

    # Spies
    Msc = MA - MB*(MD*MC);
    @show norm(MA-MA')
    @show norm(MB+MC')
    @show norm(Msc-Msc')
    p = spy(MA-MA')
    display(p)


    # Solve
    if solve
        Mscf = cholesky(Msc)
        u    = zero(ru)
        p    = zero(rp)
        nit  = 5
        for it=1:nit
            ru1  = ru - MB*p - MB*MD*rp
            u    = Mscf\ru1
            dp   = MD*(rp-MC*u);
            p    = p + dp;
            meanDiv = mean(MC*u - rp)
            println("it.",  it, " Divergence = ", norm((MC*u - rp) .- meanDiv)/length(p))
        end
    end
    # model  = ExtractData(filename, "/Model/Params")
    # xc     = ExtractData(filename, "/Model/xc_coord")
    # zc     = ExtractData(filename, "/Model/zc_coord")
    # xv     = ExtractData(filename, "/Model/xg_coord")
    # zv     = ExtractData(filename, "/Model/zg_coord")
    # xv_ph  = ExtractData(filename, "/VizGrid/xviz_hr")
    # zv_ph  = ExtractData(filename, "/VizGrid/zviz_hr")
    # xc_ph  = 0.5.*(xv_ph[1:end-1] .+ xv_ph[2:end])
    # zc_ph  = 0.5.*(zv_ph[1:end-1] .+ zv_ph[2:end])
    # ncx_ph = length(xc_ph)
    # ncz_ph = length(zc_ph)

    # @show model
    # t   = model[1]
    # nvx = Int(model[4])
    # nvz = Int(model[5])
    # ncx, ncz = nvx-1, nvz-1

    # ηc   = Float64.(reshape(ExtractData( filename, "/Centers/eta_n"), ncx, ncz))
    # ρc   = Float64.(reshape(ExtractData( filename, "/Centers/rho_n"), ncx, ncz))
    # P    = Float64.(reshape(ExtractData( filename, "/Centers/P"), ncx, ncz))
    # d    = Float64.(reshape(ExtractData( filename, "/Centers/d"), ncx, ncz))
    # ph   = Float64.(reshape(ExtractData( filename, "/VizGrid/compo_hr"), ncx_ph, ncz_ph))
    # ε̇_pl = Float64.(reshape(ExtractData( filename, "/Centers/eII_pl"), ncx, ncz))
    # p = plot()
    # # p = heatmap!(xc./1e3, zc./1e3, log10.(ηc"), c=:turbo)
    # # p = heatmap!(xc./1e3, zc./1e3, (ηc"), c=:turbo)
    # p = heatmap!(xc_ph./1e3, zc_ph./1e3, ph", c=:turbo)
    # # p = heatmap!(xc./1e3, zc./1e3, log10.(ε̇_pl"), c=:turbo)
    # p = contour!(xc_ph./1e3, zc_ph./1e3, ph", levels=[-1], lw=10, c=:black)


    # # p = heatmap!(xc./1e3, zc./1e3, (P"), c=:turbo)
    # # p = heatmap!(xc./1e3, zc./1e3, (ρc"), c=:turbo)
    # # p = heatmap!(xc./1e3, zc./1e3, log10.(d"), c=:turbo)
    # p = plot!(aspect_ratio=:equal)
    # display(p)

end

main()