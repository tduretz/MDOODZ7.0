import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
using HDF5, GLMakie, Printf#, Colors, ColorSchemes, MathTeXEngine, LinearAlgebra, FFMPEG
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

y    = 365*24*3600
My   = 1e6*y
cm_y = y*100.

function main()
    x  = (min=-1e6, max=1e6)
    z  = (min=-6e5, max=0)
    L  = (x=x.max-x.min, z=z.max-z.min)
    nc = (x=500, z=300)
    Δ  = (x=L.x/nc.x, z=L.z/nc.z)
    xc = LinRange(x.min, x.max, nc.x) 
    zc = LinRange(z.min, z.max, nc.z) 
    ρc = 3300.0.*ones(nc...)
    # ρc[(0.0*xc .+ zc') .> -1.5e5] .= 2700
    hc = 400.0.*ones(nc.x)
    # ρc[abs.(xc .+ 0.0*zc') .< 2e5 .&& (0.0*xc .+ zc') .> -1.5e5] .= 2800

    for j=1:nc.z, i=1:nc.x
        # Introduce a circle
        z0 = -250e3 
        r  = 25e3
        x  = xc[i]
        z  = zc[j]
        if (x-0.0)^2 + (z-z0)^2 < r^2
            ρc[i,j] = 3450.0 
        end

        # Introduce a margin
        Lm  = 100e3 # margin center
        L0  = -2e5   # margin position
        HW  = 1e5   # lithosphere thickness
        HE  = 2e4   # lithosphere thickness
        # 2 points  
        p1  = (x=L0-Lm/2, z=-HW)
        p2  = (x=L0+Lm/2, z=-HE)
        a   = (p2.z - p1.z) / (p2.x - p1.x)
        b   = p2.z - a*p2.x
        zm  = a*x +b
        if (abs(x-L0)<Lm && z>zm && z<-HE)
            ρc[i,j] = 2700.0
        end

        # Introduce a margin
        Lm  = 100e3 # margin center
        L0  = 2e5   # margin position
        HW  = 1e5   # lithosphere thickness
        HE  = 2e4   # lithosphere thickness
        # 2 points  
        p1  = (x=L0-Lm/2, z=-HE)
        p2  = (x=L0+Lm/2, z=-HW)
        a   = (p2.z - p1.z) / (p2.x - p1.x)
        b   = p2.z - a*p2.x
        zm  = a*x +b
        if (abs(x-L0)<Lm && z>zm && z<-HE)
            ρc[i,j] = 2700.0
        end


    end

    # Compute west column load in units of P/g
    for i=2:nc.x
        PW     = sum(ρc[1,:]*Δ.z)
        PE     = sum(ρc[i,:]*Δ.z)
        Δh     = -(PE - PW) / ρc[i,1]
        hc[i] += Δh
    end

    ##################################################################
    resol       = 1000
    f = Figure(resolution = (L.x/L.z*resol, resol), fontsize=25)

    ax1 = Axis(f[1, 1], title = L"$h$", xlabel = L"$x$ [km]", ylabel = L"$h$ [m]")
    lines!(ax1, xc/1e3, hc)

    ax2 = Axis(f[2, 1], title = L"$ρ$", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    hm = heatmap!(ax2, xc/1e3, zc/1e3, ρc, colormap = (:turbo, 1.0))           
    colsize!(f.layout, 1, Aspect(1, L.x/L.z))
    GLMakie.Colorbar(f[2, 2], hm, label = L"$ρ$ [kg.m$^3$]", width = 20, labelsize = 25, ticklabelsize = 14 )
    GLMakie.colgap!(f.layout, 20)
    display(f)
    DataInspector(f)

end

main()