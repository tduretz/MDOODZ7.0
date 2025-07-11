using CairoMakie, Printf, LinearAlgebra, MathTeXEngine, GridGeometryUtils
import Statistics:mean

function main_simple_ani_vis()
    # Material parameters
    δ          = 2.0 
    ηv         = 1.0 # normal viscosity
    # Kinematics
    pure_shear = 1
    ε̇xxd       = pure_shear*.5
    ε̇yyd       = -ε̇xxd
    ε̇xyd       = ((1.0 - pure_shear)*5.0)
    ε̇bg        = sqrt(ε̇xxd^2 + ε̇xyd^2)
    # Arrays
    θ          = LinRange( 0.0, π, 51 ) 
    τii_cart   = zero(θ)
    τii_cart_MD7 = zero(θ)
    τii_rot1   = zero(θ)
    τii_rot2   = zero(θ)
    τii_rot3   = zero(θ)
    τxx_cart   = zero(θ)
    τxy_cart   = zero(θ)
    τxx_rot    = zero(θ)
    τxy_rot    = zero(θ)
    τ_cart     = zero(θ)
    τii_cart2  = zero(θ)
    τxx_cart2  = zero(θ)
    τyy_cart2  = zero(θ)
    τxy_cart2  = zero(θ)
    # Loop over all orientations
    for i in eachindex(θ)
        # Transformation matrix: towards principal plane
        Q           = [cos(θ[i]) sin(θ[i]); -sin(θ[i]) cos(θ[i])]
        ε̇_tens      = [ε̇xxd ε̇xyd; ε̇xyd ε̇yyd]
        ε̇_rot       = Q*ε̇_tens*Q'  
        # Viscous 
        ε           = [ε̇_rot[1,1]; ε̇_rot[2,2]; ε̇_rot[1,2]]
        D           = 2*ηv*[1 0 0; 0 1 0; 0 0 1.0/δ;]
        τ           = D*ε
        # a           = (δ^2)
        Y2          = 0.5*(τ[1]^2 + τ[2]^2) + τ[3]^2*δ^2
        τii_rot1[i] = sqrt(Y2)
        # Check strain rate components
        τxx         = τ[1]
        τxy         = τ[3]
        # ε̇xxd_chk    = τxx/2/ηv
        # ε̇xyd_chk    = τxy/2/(ηv/a^2)
        # @show (ε̇_rot[1,1], ε̇xxd_chk)
        # @show (ε̇_rot[1,2], ε̇xyd_chk)
        τxx_rot[i]  = τxx
        τxy_rot[i]  = τxy
        # Effective viscosity: compute stress invariant from strain rate invariant
        I2          = 0.5*(ε[1]^2 + ε[2]^2) + ε[3]^2
        τii_rot2[i] = 2*ηv*sqrt(I2)
        # if one does the same by scaling the strain rate invariant by δ...
        I2          = 0.5*(ε[1]^2 + ε[2]^2) + ε[3]^2/δ^2
        τii_cart2[i] = 2*ηv*sqrt(I2)
        # Rotate stress back
        τ_rot       = [τ[1] τ[3]; τ[3] τ[2]]
        τ           = Q'*τ_rot*Q
        J2          = 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2
        τii_cart[i] = sqrt(J2)  
        τxx_cart[i] = τ[1,1]
        τxy_cart[i] = τ[1,2]

        @show ' '
        @show τxx_cart[i], τxy_cart[i]

        # Prediction of stress using viscosity tensor in Cartesian plane
        two      = 2.
        n        = θ[i] - π/2
        nx       = cos(n); ny = sin(n)
        N        = [nx; ny]
        d0       = 2*N[1]^2*N[2]^2
        d1       = N[1]*N[2] * (-N[1]^2 + N[2]^2)
        C_ANI    = [-d0 d0 -two*d1; d0 -d0 two*d1; -d1 d1 -two*(1/2-d0)]      ###### !!! - sign in D33 -a0
        C_ISO    = [1 0 0; 0 1 0; 0 0 two*1//2]
        ani      = 1. - 1.  ./ δ
        Dani     = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
        ε        = [ε̇xxd; ε̇yyd; ε̇xyd]
        #τ0       = [τxx0; τyy0; τxy0]
        ε_eff    = ε #.+ inv(Dani) * τ0 ./(2*ηe)
        τ_MD7_v1 = 2*ηv*Dani*ε_eff

        τii_cart_MD7[i]          = sqrt(0.5*(τ_MD7_v1[1]^2 + τ_MD7_v1[2]^2) + τ_MD7_v1[3]^2)

        @show τ_MD7_v1[1], τ_MD7_v1[3]

        eta_ve = ηv
        delta  = δ
        theta  = θ[i]
        Exx = ε̇xxd
        Eyy = ε̇yyd
        Exy = ε̇xyd
        τ_cart[i] = 2.0 * eta_ve .* sqrt(0.25 * Exx .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Exx .^ 2 .* delta .^ 2 - 0.25 * Exx .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Exx .^ 2 + 0.5 * Exx .* Exy .* delta .^ 2 .* sin(4 * theta) - Exx .* Exy .* sin(4 * theta) / 2 - 0.5 * Exx .* Eyy .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.5 * Exx .* Eyy .* delta .^ 2 + 0.5 * Exx .* Eyy .* cos(2 * theta) .^ 2 - 0.5 * Exx .* Eyy - Exy .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + Exy .^ 2 .* delta .^ 2 + Exy .^ 2 .* cos(2 * theta) .^ 2 - 0.5 * Exy .* Eyy .* delta .^ 2 .* sin(4 * theta) + Exy .* Eyy .* sin(4 * theta) / 2 + 0.25 * Eyy .^ 2 .* delta .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Eyy .^ 2 .* delta .^ 2 - 0.25 * Eyy .^ 2 .* cos(2 * theta) .^ 2 + 0.25 * Eyy .^ 2) ./ delta

        τv           = τii_rot1[i]
        τxx_cart2[i] = τv   * cos(2*θ[i]) 
        τyy_cart2[i] =-τv   * cos(2*θ[i])
        τxy_cart2[i] = τv/δ * sin(2*θ[i])

        ev = eigvecs(τ)
        @show atand(ev[2,1]/ev[1,1])

    end

    f = Figure(size = (500, 400), fontsize=18)


    ncx, ncy = 1000, 1000    
    phase  = zeros(ncx, ncy)
    dx, dy = 1/ncx, 1/ncy
    xc     = LinRange(-1/2+dx/2, 1/2-dx/2, ncx)
    yc     = LinRange(-1/2+dy/2, 1/2-dy/2, ncy)

    ############################ Flat 
    layering = Layering(
            (0., 0.), 
            0.1, 
            0.5; 
            θ = 0,  
            perturb_amp=0*1.0, 
            perturb_width=1.0
        )

    # Set material geometry 
    phase .= 0.
    for I in CartesianIndices(phase)   # loop on centroids
        𝐱 = @SVector([xc[I[1]], yc[I[2]]])
        isin = inside(𝐱, layering)
        if isin 
            phase[I[1], I[2]] = 2
        end 
    end

    ax1 = Axis(f[3, 1], aspect=DataAspect(), title=L"$\theta=0$")
    heatmap!(ax1, xc, yc, phase, colormap=:grays)
    hidedecorations!(ax1)

    ############################ 45 
    
    layering = Layering(
            (0., 0.), 
            0.1, 
            0.5; 
            θ = π/4,  
            perturb_amp=0*1.0, 
            perturb_width=1.0
        )

    # Set material geometry 
    phase .= 0.
    for I in CartesianIndices(phase)   # loop on centroids
        𝐱 = @SVector([xc[I[1]], yc[I[2]]])
        isin = inside(𝐱, layering)
        if isin 
            phase[I[1], I[2]] = 2
        end 
    end

    ax1 = Axis(f[3, 2], aspect=DataAspect(), title=L"$\theta= 45^\circ$")
    heatmap!(ax1, xc, yc, phase, colormap=:grays)
    hidedecorations!(ax1)

    rowsize!(f.layout, 1, Relative(0.5))


    ############################ vertical 
    phase .= 0
    layering = Layering(
            (0., 0.), 
            0.1, 
            0.5; 
            θ = π/2,  
            perturb_amp=0*1.0, 
            perturb_width=1.0
        )

    # Set material geometry 
    phase .= 0.
    for I in CartesianIndices(phase)   # loop on centroids
        𝐱 = @SVector([xc[I[1]], yc[I[2]]])
        isin = inside(𝐱, layering)
        if isin 
            phase[I[1], I[2]] = 2
        end 
    end

    ax1 = Axis(f[3, 3], aspect=DataAspect(), title=L"$\theta= 90^\circ$")
    heatmap!(ax1, xc, yc, phase, colormap=:grays)
    hidedecorations!(ax1)

    ############################ 45 
    phase .= 0
    layering = Layering(
            (0., 0.), 
            0.1, 
            0.5; 
            θ = 3*π/4,  
            perturb_amp=0*1.0, 
            perturb_width=1.0
        )

    # Set material geometry 
    phase .= 0.
    for I in CartesianIndices(phase)   # loop on centroids
        𝐱 = @SVector([xc[I[1]], yc[I[2]]])
        isin = inside(𝐱, layering)
        if isin 
            phase[I[1], I[2]] = 2
        end 
    end

    ax1 = Axis(f[3, 4], aspect=DataAspect(), title=L"$\theta= 135^\circ$")
    heatmap!(ax1, xc, yc, phase, colormap=:grays)
    hidedecorations!(ax1)


   ############################ Flat 
   layering = Layering(
    (0., 0.), 
    0.1, 
    0.5; 
    θ = 0,  
    perturb_amp=0, 
    perturb_width=1.0
)

    # Set material geometry 
    phase .= 0.
    for I in CartesianIndices(phase)   # loop on centroids
        𝐱 = @SVector([xc[I[1]], yc[I[2]]])
        isin = inside(𝐱, layering)
        if isin 
            phase[I[1], I[2]] = 2
        end 
    end

    ax1 = Axis(f[3, 5], aspect=DataAspect(), title=L"$\theta= 0$")
    heatmap!(ax1, xc, yc, phase, colormap=:grays)
    hidedecorations!(ax1)
    rowsize!(f.layout, 1, Relative(0.5))

    ax1 = Axis(f[1, 1:5], title=L"$$Horizontal extension: variable $\tau_{II}$ and invariant $\tau_{II}'$", xlabel=L"$\theta$ [$^\circ$]", ylabel=L"$\tau_{II}$, $\tau_{II}'$  [-]")
    lines!(ax1, θ*180/π, τii_cart, label=L"$\tau_{II}$  ($\delta=2$, transformations)")
    scatter!(ax1, θ*180/π, τii_cart_MD7, label=L"$\tau_{II}$  ($\delta=2$, tensor)")

    # scatter!(ax1, θ*180/π, τii_cart_MD7 )
    lines!(ax1, θ*180/π, τii_rot1, label=L"$\tau_{II}'$ $(\delta=2)$")
    # text!(ax1, 50, 0.5; text=L"$\tau_{II} = \tau_{II}' \frac{\sqrt{(\delta^2-1)\cos{(2\theta)} + 1)}}{\delta} $", fontsize=12, rotation=π/2.7)

    tii = 1.0*sqrt.((δ^2-1)*cos.(2*θ).^2 .+ 1)/δ
    # scatter!(ax1, θ*180/π, tii, )

    scatter!(ax1, θ[1:5:end]*180/π, τii_rot2[1:5:end], label=L"$\tau_{II}$ $(\delta=1)$")#, label=L"$\tau_{II}$  $(\delta$=1)$" )
    # lines!(ax1, θ*180/π, τii_cart2 )
    # scatter!(ax1, θ*180/π, τii_cart, marker=:xcross, markersize=10 )
    # axislegend(position=:rb)

    f[2, 1:5] = Legend(f, ax1, framevisible = false, orientation = :horizontal, labelsize=20, nbanks=2)

    # ax3 = Axis(f[1, 1], title=L"$$A) Deviatoric stress components ($\delta = 2$)", xlabel=L"$\tau_{xx}'$ [-]", ylabel=L"$\tau_{xy}'$ [-]", aspect=DataAspect())
    # lines!(ax3, LinRange(0,1,10), zeros(10), color=:black)
    # text!(ax3, 0.3, 0.05; text=L"$\tau_{xx}' = \tau_{II}'$", fontsize=16)
    # lines!(ax3, zeros(10), LinRange(0,0.5,10), color=:black)
    # text!(ax3, -0.05, 0.0; text=L"$\tau_{xy}' = \frac{\tau_{II}'}{\delta}$", rotation=π/2, fontsize=16)

    # a = LinRange(0, 2π, 100)
    # txx = 1*cos.(a)
    # txy = 1*sin.(a)/δ

    # # lines!(ax3, τxx_cart, τxy_cart, marker=:xcross )
    # lines!(ax3, τxx_rot, τxy_rot, label=L"$\delta$=2$" )
    # # scatter!(ax3, txx, txy, marker=:xcross )

    save("/Users/tduretz/PowerFolders/_manuscripts/RiftingAnisotropy/Figures/ViscousAnisotropySimple_Revision1.png", f, px_per_unit = 4) 

    
    display(f)

end

main_simple_ani_vis()