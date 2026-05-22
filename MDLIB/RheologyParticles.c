// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2021  MDOODZ Developper team
//
// This file is part of MDOODZ.
//
// MDOODZ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MDOODZ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MDOODZ.  If not, see <http://www.gnu.org/licenses/>.
// =========================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "mdoodz-private.h"

#include "mdoodz-log.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseGrainSizeParticles( markers* particles, mat_prop *materials ){
    // Set initial particle grain size
#pragma omp parallel for shared( particles, materials )
    for( int k=0; k<particles->Nb_part; k++ ) {
        const int phase = particles->phase[k];
        if (phase != -1) particles->d[k] = materials->gs_ref[phase];
    }
}

// Per-marker init-from-finite-strain helper (design D1-D4, change
// aniso-init-from-finite-strain). For a single marker on an ani_fstrain==3
// phase with aniso_factor > 1, constructs F_init such that the per-marker
// step-0 δ derived from F_init via the existing FiniteStrainAspectRatio SVD
// + aniso_delta_fn pipeline equals exactly aniso_factor[phase]. Also sets
// aniso_delta_fs_prev = aniso_factor so the step-1 operator-split increment
// (δ_FS − aniso_delta_fs_prev) is the change in Hansen δ from step 0 to step
// 1 — no permanent +offset. Cold-limit telescope collapses to
//     aniso_delta_N → aniso_delta_fn(FS_AR_N)
// and Hansen-saturation δ ≤ 14.13 holds by construction.
//
// Math (design D1-D3):
//   γ_eff   = aniso_delta_fn_inv[phase](aniso_factor[phase])
//   FS_AR   = ((γ_eff + √(γ_eff² + 4)) / 2)²
//   γ_eq    = ln(FS_AR)
//   a, b    = exp(±γ_eq/2)   (principal stretches of pure-shear F)
//   θ_F     = aniso_angle[phase] − π/2     ← D3 corrigendum (foliation parallel)
//   c, s    = cos(θ_F), sin(θ_F)
//   F_init  = R(θ_F) · diag(a, b) · R(θ_F)^T (symmetric pure-shear)
//   Fxx = c²·a + s²·b,  Fzz = s²·a + c²·b,  Fxz = Fzx = c·s·(a − b)
//
// aniso_angle[phase] is in RADIANS at this point (InputOutput.c:1456 already
// did the deg → rad conversion). We subtract π/2 unconditionally to align
// the F principal extension with the easy-slip plane direction (foliation),
// not the director.
static inline void aniso_init_finite_strain_for_marker( int k, int phase,
                                                        mat_prop *materials,
                                                        markers *particles ) {
    const double aniso_factor = materials->aniso_factor[phase];
    const double aniso_angle  = materials->aniso_angle[phase];    // radians

    // Resolve γ_eff via per-phase analytic inverse. The caller has already
    // gated this on aniso_delta_fn_inv != NULL.
    const double gamma_eff = materials->aniso_delta_fn_inv[phase]( aniso_factor );

    // FS_AR(γ_eff) from the simple-shear-equivalent inversion
    //   γ_eff = √FS_AR − 1/√FS_AR     (forward in aniso_delta_fn)
    //   → √FS_AR = (γ_eff + √(γ_eff² + 4)) / 2
    const double sqrt_FS_AR = 0.5 * ( gamma_eff + sqrt( gamma_eff*gamma_eff + 4.0 ) );
    const double FS_AR      = sqrt_FS_AR * sqrt_FS_AR;
    const double gamma_eq   = log( FS_AR );

    // Pure-shear F_canonical principal stretches.
    const double a = exp(  0.5 * gamma_eq );
    const double b = exp( -0.5 * gamma_eq );

    // Rotation θ_F = aniso_angle − π/2 (foliation direction, not director).
    const double theta_F = aniso_angle - 0.5 * M_PI;
    const double c       = cos( theta_F );
    const double s       = sin( theta_F );

    // F_init = R(θ_F) · diag(a, b) · R(θ_F)^T   (symmetric, det = 1)
    particles->Fxx[k] = c*c*a + s*s*b;
    particles->Fzz[k] = s*s*a + c*c*b;
    particles->Fxz[k] = c*s*( a - b );
    particles->Fzx[k] = c*s*( a - b );

    // δ-state consistency: with F_init as built, the step-0 production δ from
    // the SVD-based FiniteStrainAspectRatio + aniso_delta_fn pipeline equals
    // aniso_factor by construction. Set aniso_delta_fs_prev to match — so the
    // step-1 operator-split increment is purely (δ_FS_1 − δ_FS_0), no offset.
    particles->aniso_delta_fs_prev[k] = aniso_factor;
    particles->aniso_delta[k]         = aniso_factor;
}

void InitialiseAnisoDeltaParticles( markers* particles, mat_prop *materials ){
    // Initialise per-marker viscous-anisotropy state on ani_fstrain == 3 phases
    // so the rifting "inherited cratonic fabric" interpretation is internally
    // consistent (design aniso-init-from-finite-strain). For each marker on an
    // ani_fstrain == 3 phase with aniso_factor > 1 and a non-NULL analytic
    // inverse aniso_delta_fn_inv (populated per aniso_db by FlowLaws.c::
    // ReadDataAnisotropy), the helper constructs:
    //   F_init = R(aniso_angle − π/2) · diag(a, b) · R(aniso_angle − π/2)^T
    // where (a, b) = exp(±γ_eq/2) and γ_eq = ln(FS_AR(γ_eff(aniso_factor))).
    // It then sets aniso_delta_fs_prev = aniso_factor and aniso_delta =
    // aniso_factor. With these three boundary conditions consistent, the
    // FiniteStrainAspectRatio operator split's cold-limit telescope collapses
    // to aniso_delta_N → aniso_delta_fn(FS_AR_N) exactly, the Hansen-saturation
    // ceiling holds by construction, and the inherited-fabric persistence in
    // cold quiescent rock is preserved by the strain-rate gate.
    //
    // For aniso_factor ≤ 1 (legacy isotropic init), aniso_db = 0 (no analytic
    // forward / inverse), or ani_fstrain != 3, the helper is skipped and the
    // legacy F = I, aniso_delta = 1, aniso_delta_fs_prev = 1 state from
    // PartInit is preserved.
    //
    // Guards (design D7):
    //  • aniso_factor < 1.0     → clamp to 1.0 with one-shot LOG_WARN per phase.
    //  • aniso_delta_fn_inv NULL (aniso_db = 0) → skip silently; legacy fallthrough.
    static int phase_logged[20] = {0};       // one-shot LOG_INFO per phase
    static int phase_warn_neg[20] = {0};     // one-shot LOG_WARN per phase

    for( int phase=0; phase<materials->Nb_phases && phase<20; phase++ ) {
        if ( materials->ani_fstrain[phase] != 3 ) continue;

        // guard aniso_factor < 1 (legacy invalid setting).
        if ( materials->aniso_factor[phase] < 1.0 ) {
            if ( !phase_warn_neg[phase] ) {
                LOG_WARN("InitialiseAnisoDelta: phase %d aniso_factor=%g < 1 (invalid); clamping to 1.0",
                         phase, materials->aniso_factor[phase]);
                phase_warn_neg[phase] = 1;
            }
            materials->aniso_factor[phase] = 1.0;
        }

        // aniso_db = 0 ⇒ aniso_delta_fn_inv NULL ⇒ skip helper.
        const int has_inverse = ( materials->aniso_delta_fn_inv[phase] != NULL );

        // per-phase one-shot LOG_INFO (suppressed for aniso_factor == 1).
        if ( !phase_logged[phase]
             && materials->aniso_factor[phase] > 1.0 + 1.0e-12
             && has_inverse ) {
            const double factor    = materials->aniso_factor[phase];
            const double angle_rad = materials->aniso_angle[phase];
            const double gamma_eff = materials->aniso_delta_fn_inv[phase]( factor );
            const double sqrt_FS_AR = 0.5 * ( gamma_eff + sqrt( gamma_eff*gamma_eff + 4.0 ) );
            const double FS_AR      = sqrt_FS_AR * sqrt_FS_AR;
            LOG_INFO("InitialiseAnisoDelta: phase %d ani_fstrain=3 aniso_factor=%g aniso_angle=%g deg -> gamma_eff=%g FS_AR=%g",
                     phase, factor, angle_rad * 180.0 / M_PI, gamma_eff, FS_AR);
            phase_logged[phase] = 1;
        }
    }

    // OpenMP review: each iteration touches only particles->{Fxx,Fxz,
    // Fzx,Fzz,aniso_delta,aniso_delta_fs_prev}[k] — fully per-marker, no shared
    // mutable state. materials and particles pointers are shared (read-only for
    // materials; per-k writes on disjoint indices for particles arrays). The
    // existing static schedule is preserved.
#pragma omp parallel for shared( particles, materials )
    for( int k=0; k<particles->Nb_part; k++ ) {
        const int phase = particles->phase[k];
        if ( phase == -1 )                                  continue;
        if ( materials->ani_fstrain[phase] != 3 )           continue;

        // Skip phases with no analytic inverse (aniso_db = 0 or unconfigured).
        if ( materials->aniso_delta_fn_inv[phase] == NULL ) {
            // Legacy fallthrough: keep PartInit defaults (F = I, aniso_delta = 1,
            // aniso_delta_fs_prev = 1). Preserve the historical aniso_factor copy
            // into aniso_delta so the few legacy SETs that relied on it still get
            // the per-marker initial δ they expect (degenerate without inverse).
            particles->aniso_delta[k] = materials->aniso_factor[phase];
            continue;
        }

        // aniso_factor ≈ 1: helper is a no-op (γ_eff=0, F=I). Skip so legacy
        // init1-style runs (aniso_factor = 1, F = I, δ = 1) are byte-identical.
        if ( materials->aniso_factor[phase] <= 1.0 + 1.0e-12 ) {
            particles->aniso_delta[k]         = materials->aniso_factor[phase];
            particles->aniso_delta_fs_prev[k] = materials->aniso_factor[phase];
            continue;
        }

        aniso_init_finite_strain_for_marker( k, phase, materials, particles );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// // Old deviatoric stress field and pressure --- DEPRECATED
// void OldDeviatoricStressesPressure( grid* mesh, markers* particles, scale scaling, params* model ) {

//     int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
//     double *sxx0,  *syy0, *szz0, *sxz0;
//     int    cent=1, vert=0, prop=1, interp=0;

//     Nx  = mesh->Nx;
//     Nz  = mesh->Nz;
//     Ncx = Nx-1;
//     Ncz = Nz-1;

//     sxx0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     syy0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     szz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     sxz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

//     P2Mastah( model, *particles, particles->sxxd,    mesh, sxx0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil, NULL);
//     P2Mastah( model, *particles, particles->szzd,    mesh, szz0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil, NULL);
//     P2Mastah( model, *particles, particles->syy,     mesh, syy0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil, NULL);
//     P2Mastah( model, *particles, particles->sxz,     mesh, mesh->sxz0,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil, NULL);

// #pragma omp parallel for shared( mesh, sxx0, syy0, szz0 ) private( k, k1, l, c0, c1, c2 ) firstprivate( Nx, Ncx, Ncz )
//     for ( k1=0; k1<Ncx*Ncz; k1++ ) {
//         k  = mesh->kp[k1];
//         l  = mesh->lp[k1];
//         c0 = k  + l*(Nx-1);
//         c1 = k  + l*(Nx);
//         c2 = k  + l*(Nx+1);

//         if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
//             mesh->p0_n[c0]       = -1.0/3.0*(sxx0[c0] + syy0[c0] + szz0[c0] );
//             mesh->sxxd0[c0]      =  mesh->p0_n[c0] + sxx0[c0];
//             mesh->szzd0[c0]      =  mesh->p0_n[c0] + szz0[c0];
//         }
//     }

//     InterpCentroidsToVerticesDouble( mesh->sxxd0,  mesh->sxxd0_s, mesh, model );
//     InterpCentroidsToVerticesDouble( mesh->szzd0,  mesh->szzd0_s, mesh, model );
//     InterpCentroidsToVerticesDouble( mesh->p0_n,   mesh->p0_s,    mesh, model );
//     InterpVerticesToCentroidsDouble( mesh->sxz0_n, mesh->sxz0,    mesh, model );

//     DoodzFree( sxx0 );
//     DoodzFree( szz0 );
//     DoodzFree( sxz0 );
// }

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// // Total stress field --- DEPRECATED
// void TotalStresses( grid* mesh, markers* particles, scale scaling, params* model ) {

//     int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
//     double *dsxx, *dsyy, *dszz, *dsxz, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0;
//     double *mdsxx, *mdszz, *mdsyy, *mdsxz;

//     Nx  = mesh->Nx;
//     Nz  = mesh->Nz;
//     Ncx = Nx-1;
//     Ncz = Nz-1;

//     dsxx = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dsyy = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dszz = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dsxz = DoodzCalloc(Nx *Nz , sizeof(DoodzFP));

//     mdsxx = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdsyy = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdszz = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdsxz = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

// #pragma omp parallel for shared( mesh, dsxx, dsyy, dszz ) private( k, k1, l, c0, c1, c2, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0  ) firstprivate( Nx, Ncx, Ncz )
//     for ( k1=0; k1<Ncx*Ncz; k1++ ) {
//         k  = mesh->kp[k1];
//         l  = mesh->lp[k1];
//         c0 = k  + l*(Nx-1);
//         c1 = k  + l*(Nx);
//         c2 = k  + l*(Nx+1);

//         if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
//             tyy      = -(mesh->sxxd[c0]  + mesh->szzd[c0]);
//             tyy0     = -(mesh->sxxd0[c0] + mesh->szzd0[c0]);
//             sxx0     = -mesh->p0_n[c0] + mesh->sxxd0[c0];
//             szz0     = -mesh->p0_n[c0] + mesh->szzd0[c0];
//             syy0     = -mesh->p0_n[c0] + tyy0;
//             sxx      = -mesh->p_in[c0] + mesh->sxxd[c0];
//             szz      = -mesh->p_in[c0] + mesh->szzd[c0];
//             syy      = -mesh->p_in[c0] + tyy;
//             dsxx[c0] = sxx - sxx0;
//             dsyy[c0] = syy - syy0;
//             dszz[c0] = szz - szz0;
//         }
//     }

//     // Vertex: shear stress change
//     for (k=0; k<Nx; k++) {
//         for (l=0; l<Nz; l++) {
//             c1 = k  + l*(Nx);
//             if (mesh->BCg.type[c1] !=30 ) dsxz[c1] =  mesh->sxz[c1] - mesh->sxz0[c1];
//         }
//     }

//     // Interpolate stress changes to markers
//     Interp_Grid2P_centroids2( *particles, mdsxx,  mesh, dsxx, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P_centroids2( *particles, mdszz,  mesh, dszz, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P_centroids2( *particles, mdsyy,  mesh, dsyy, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P( *particles, mdsxz,  mesh, dsxz, mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );

//     // Update marker stresses
//     ArrayPlusArray( particles->sxxd, mdsxx, particles->Nb_part );
//     ArrayPlusArray( particles->szzd, mdszz, particles->Nb_part );
//     ArrayPlusArray( particles->syy , mdsyy, particles->Nb_part );
//     ArrayPlusArray( particles->sxz,  mdsxz, particles->Nb_part );

//     // Freedom
//     DoodzFree( dsxx );
//     DoodzFree( dsyy );
//     DoodzFree( dszz );
//     DoodzFree( dsxz );
//     DoodzFree(mdsxx);
//     DoodzFree(mdszz);
//     DoodzFree(mdsyy);
//     DoodzFree(mdsxz);
// }

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrainII( grid* mesh, scale scaling, params model, markers* particles, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {

    double *strain_inc;
    int k, l, c1;
    DoodzFP *strain_inc_el, *strain_inc_pl, *strain_inc_pl_vol, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;

    //    printf("Accumulating strain\n");

    // Allocate
    strain_inc        = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_el     = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pl     = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pl_vol = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pwl    = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_exp    = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_lin    = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_gbs    = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));

    // Interpolate exz to cell centers
#pragma omp parallel for shared( mesh, model, strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs  ) private( c1 )
    for (c1=0; c1<(mesh->Nx-1)*(mesh->Nz-1); c1++) {

        strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);

        if (strain_inc[c1]<0) {
            LOG_ERR("negative strain increment at cell %d", c1);
            LOG_ERR("pl = %2.2e pwl = %2.2e exp = %2.2e lin = %2.2e gbs = %2.2e el = %2.2e" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
            exit(1);
        }

        strain_inc_el[c1]      = model.dt*mesh->eII_el[c1];
        strain_inc_pl[c1]      = model.dt*mesh->eII_pl[c1];
        strain_inc_pl_vol[c1]  = model.dt*mesh->div_u_pl[c1];
        strain_inc_pwl[c1]     = model.dt*mesh->eII_pwl[c1];
        strain_inc_exp[c1]     = model.dt*mesh->eII_exp[c1];
        strain_inc_lin[c1]     = model.dt*mesh->eII_lin[c1];
        strain_inc_gbs[c1]     = model.dt*mesh->eII_gbs[c1];
    }

    //---------------------------------------------------------------//

    double dE_tot, dE_el, dE_pl, dE_pl_vol, dE_pwl, dE_exp, dE_lin, dE_gbs;
    double dx, dz, dxm, dzm, dst, sumW;
    int    i_part, j_part, iSW, iNW, iSE, iNE;
    dx=mesh->dx;
    dz=mesh->dz;

#pragma omp parallel for shared( particles, strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pl_vol, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs, tag ) private( k, l, dE_tot, dE_el, dE_pl, dE_pl_vol, dE_pwl, dE_exp, dE_lin, dE_gbs, sumW, dst, dxm, dzm, iSW, iSE, iNW, iNE, i_part, j_part ) firstprivate( dx, dz, Nx, Nz )
    for (k=0;k<particles->Nb_part;k++) {

        dE_tot     = 0.0;
        dE_el      = 0.0;
        dE_pl      = 0.0;
        dE_pl_vol  = 0.0;
        dE_pwl     = 0.0;
        dE_exp     = 0.0;
        dE_lin     = 0.0;
        dE_gbs     = 0.0;


        if (particles->phase[k] != -1) {

            dst = (particles->x[k]-X_vect[0]);
            j_part   = ceil((dst/dx)) - 1;
            if (j_part<0) {
                j_part = 0;
            }
            if (j_part>Nx-2) {
                j_part = Nx-2;
            }

            // Get the line:
            dst = (particles->z[k]-Z_vect[0]);
            i_part   = ceil((dst/dz)) - 1;
            if (i_part<0) {
                i_part = 0;
            }
            if (i_part>Nz-2) {
                i_part = Nz-2;
            }

            dxm = (particles->x[k] - X_vect[j_part]);
            dzm = (particles->z[k] - Z_vect[i_part]);

            iSW = j_part+i_part*Nx;
            iSE = j_part+i_part*Nx+1;
            iNW = j_part+(i_part+1)*Nx;
            iNE = j_part+(i_part+1)*Nx+1;

            dE_tot     = 0.0;
            dE_el      = 0.0;
            dE_pl      = 0.0;
            dE_pl_vol  = 0.0;
            dE_pwl     = 0.0;
            dE_exp     = 0.0;
            dE_lin     = 0.0;
            dE_gbs     = 0.0;

            sumW         = 0.0;

            //            if (j_part>0 && j_part<Nx-2 && i_part>0 && i_part<Nz-2) {

            if (tag[iSW]!=30 && tag[iSW]!=31) {
                dE_tot    +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc[iSW];
                dE_el     +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSW];
                dE_pl     +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSW];
                dE_pl_vol +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pl_vol[iSW];
                dE_pwl    +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSW];
                dE_exp    +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSW];
                dE_lin    +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSW];
                dE_gbs    +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSW];
                sumW      +=  (1.0-dxm/dx) * (1.0-dzm/dz);
            }
            if (tag[iSE]!=30 && tag[iSE]!=31) {
                dE_tot    +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc[iSE];
                dE_el     +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSE];
                dE_pl     +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSE];
                dE_pl_vol +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pl_vol[iSE];
                dE_pwl    +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSE];
                dE_exp    +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSE];
                dE_lin    +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSE];
                dE_gbs    +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSE];
                sumW      +=  (dxm/dx) * (1.0-dzm/dz);
            }
            if (tag[iNW]!=30 && tag[iNW]!=31) {
                dE_tot    +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc[iNW];
                dE_el     +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_el[iNW];
                dE_pl     +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pl[iNW];
                dE_pl_vol +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pl_vol[iNW];
                dE_pwl    +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pwl[iNW];
                dE_exp    +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_exp[iNW];
                dE_lin    +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_lin[iNW];
                dE_gbs    +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_gbs[iNW];
                sumW      +=  (1.0-dxm/dx) * (dzm/dz);
            }
            if (tag[iNE]!=30 && tag[iNE]!=31) {
                dE_tot    +=  (dxm/dx) * (dzm/dz) * strain_inc[iNE];
                dE_el     +=  (dxm/dx) * (dzm/dz) * strain_inc_el[iNE];
                dE_pl     +=  (dxm/dx) * (dzm/dz) * strain_inc_pl[iNE];
                dE_pl_vol +=  (dxm/dx) * (dzm/dz) * strain_inc_pl_vol[iNE];
                dE_pwl    +=  (dxm/dx) * (dzm/dz) * strain_inc_pwl[iNE];
                dE_exp    +=  (dxm/dx) * (dzm/dz) * strain_inc_exp[iNE];
                dE_lin    +=  (dxm/dx) * (dzm/dz) * strain_inc_lin[iNE];
                dE_gbs    +=  (dxm/dx) * (dzm/dz) * strain_inc_gbs[iNE];
                sumW      += (dxm/dx)* (dzm/dz);
            }

            //            printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", sumW, dE_tot, strain_inc[iSW], strain_inc[iSE], strain_inc[iNW], strain_inc[iNE]);
            //
            //             if (dE_tot<0.0)exit(1);
            //            }

            //            if (j_part==0    && i_part==0   ) {
            //                dE_tot = strain_inc[iNE];
            //                dE_el  = strain_inc[iNE];
            //                dE_pl  = strain_inc[iNE];
            //                dE_pwl = strain_inc[iNE];
            //                dE_exp = strain_inc[iNE];
            //                dE_lin = strain_inc[iNE];
            //                dE_gbs = strain_inc[iNE];
            //            }
            //            if (j_part==Nx-2 && i_part==0   ) {
            //                dE_tot = strain_inc[iNW];
            //                dE_el  = strain_inc[iNW];
            //                dE_pl  = strain_inc[iNW];
            //                dE_pwl = strain_inc[iNW];
            //                dE_exp = strain_inc[iNW];
            //                dE_lin = strain_inc[iNW];
            //                dE_gbs = strain_inc[iNW];
            //            }
            //            if (j_part==0    && i_part==Nz-2) {
            //                dE_tot = strain_inc[iSE];
            //                dE_el  = strain_inc[iSE];
            //                dE_pl  = strain_inc[iSE];
            //                dE_pwl = strain_inc[iSE];
            //                dE_exp = strain_inc[iSE];
            //                dE_lin = strain_inc[iSE];
            //                dE_gbs = strain_inc[iSE];
            //            }
            //            if (j_part==Nx-2 && i_part==Nz-2) {
            //                dE_tot = strain_inc[iSW];
            //                dE_el  = strain_inc[iSW];
            //                dE_pl  = strain_inc[iSW];
            //                dE_pwl = strain_inc[iSW];
            //                dE_exp = strain_inc[iSW];
            //                dE_lin = strain_inc[iSW];
            //                dE_gbs = strain_inc[iSW];
            //            }

            //            if (j_part==0 && i_part>0 && i_part<Nz-2) {
            //
            //                if (tag[iSE]!=30 && tag[iSE]!=31) {
            //                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSE];
            //                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSE];
            //                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSE];
            //                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSE];
            //                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSE];
            //                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSE];
            //                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSE];
            //                    sumW +=  (1.0-dzm/dz);
            //                }
            //                if (tag[iNE]!=30 && tag[iNE]!=31) {
            //                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNE];
            //                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNE];
            //                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNE];
            //                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNE];
            //                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNE];
            //                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNE];
            //                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNE];
            //                    sumW +=  (dzm/dz);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            //
            //            if (j_part==Nx-2 && i_part>0 && i_part<Nz-2) {
            //
            //                if (tag[iSW]!=30 && tag[iSW]!=31) {
            //                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSW];
            //                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSW];
            //                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSW];
            //                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSW];
            //                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSW];
            //                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSW];
            //                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSW];
            //                    sumW +=  (1.0-dzm/dz);
            //                }
            //                if (tag[iNW]!=30 && tag[iNW]!=31) {
            //                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNW];
            //                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNW];
            //                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNW];
            //                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNW];
            //                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNW];
            //                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNW];
            //                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNW];
            //                    sumW +=  (dzm/dz);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }

            //            if (i_part==0 && j_part>0 && j_part<Nx-2) {
            //                if (tag[iNW]!=30 && tag[iNW]!=31) {
            //                    dE_tot += (1.0-dxm/dx)   * strain_inc[iNW];
            //                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iNW];
            //                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iNW];
            //                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iNW];
            //                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iNW];
            //                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iNW];
            //                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iNW];
            //                    sumW += (1.0-dxm/dx);
            //                }
            //                if (tag[iNE]!=30 && tag[iNE]!=31) {
            //                    dE_tot += (dxm/dx)   * strain_inc[iNE];
            //                    dE_el  += (dxm/dx)   * strain_inc_el[iNE];
            //                    dE_pl  += (dxm/dx)   * strain_inc_pl[iNE];
            //                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iNE];
            //                    dE_exp += (dxm/dx)   * strain_inc_exp[iNE];
            //                    dE_lin += (dxm/dx)   * strain_inc_lin[iNE];
            //                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iNE];
            //                    sumW += (dxm/dx);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            //            if (i_part==Nz-2 && j_part>0 && j_part<Nx-2) {
            //                if (tag[iSW]!=30 && tag[iSW]!=31) {
            //                    dE_tot += (1.0-dxm/dx)   * strain_inc[iSW];
            //                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iSW];
            //                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iSW];
            //                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iSW];
            //                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iSW];
            //                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iSW];
            //                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iSW];
            //                    sumW += (1.0-dxm/dx);
            //                }
            //                if ( tag[iSE] != 30 && tag[iSE] != 31 ) {
            //                    dE_tot += (dxm/dx)   * strain_inc[iSE];
            //                    dE_el  += (dxm/dx)   * strain_inc_el[iSE];
            //                    dE_pl  += (dxm/dx)   * strain_inc_pl[iSE];
            //                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iSE];
            //                    dE_exp += (dxm/dx)   * strain_inc_exp[iSE];
            //                    dE_lin += (dxm/dx)   * strain_inc_lin[iSE];
            //                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iSE];
            //                    sumW += (dxm/dx);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }

            if(sumW>1e-13) {
                dE_tot    /= sumW;
                dE_el     /= sumW;
                dE_pl     /= sumW;
                dE_pl_vol /= sumW;
                dE_pwl    /= sumW;
                dE_exp    /= sumW;
                dE_lin    /= sumW;
                dE_gbs    /= sumW;
            }


        }

        particles->strain[k]        += dE_tot;
        particles->strain_el[k]     += dE_el;
        particles->strain_pl[k]     += dE_pl;
        particles->strain_pl_vol[k] += dE_pl_vol;
        particles->strain_pwl[k]    += dE_pwl;
        particles->strain_exp[k]    += dE_exp;
        particles->strain_lin[k]    += dE_lin;
        particles->strain_gbs[k]    += dE_gbs;
    }

    //---------------------------------------------------------------//

    // Free
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pl_vol);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DeformationGradient ( grid mesh, scale scaling, params model, mat_prop materials, markers *particles ) {

    int k, l, cp, cu, cv, Nx, Nz;
    double dx, dz;
    double fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o;

    Nx = mesh.Nx;
    Nz = mesh.Nz;
    dx = mesh.dx;
    dz = mesh.dz;

    double *dudx, *dudz, *dvdx, *dvdz;
    DoodzFP *pdudx, *pdudz, *pdvdx, *pdvdz;

    dudx   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dvdz   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dudz   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    dvdx   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    pdudx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdudz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));

    // Compute dudx and dvdz (cell centers)
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx) + Nx;
            cv = k  + l*(Nx+1) + 1;
            dudx[cp] = 1.0/dx * ( mesh.u_in[cu+1]      - mesh.u_in[cu] );
            dvdz[cp] = 1.0/dz * ( mesh.v_in[cv+(Nx+1)] - mesh.v_in[cv] );
        }
    }

    // Compute dudx and dvdz (cell vertices)
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx);
            cv = k  + l*(Nx+1);
            dudz[cu] = 1.0/dz * ( mesh.u_in[cu+Nx] - mesh.u_in[cu] );
            dvdx[cu] = 1.0/dx * ( mesh.v_in[cv+1]  - mesh.v_in[cv] );
        }
    }

    // Interpolate from grid to particles
    Interp_Grid2P_centroids2( *(particles), pdudx, &mesh, dudx, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model  );
    Interp_Grid2P_centroids2( *(particles), pdvdz, &mesh, dvdz, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
    Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );

#pragma omp parallel for shared ( particles ) private ( k, fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o ) firstprivate( model, materials ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {

        // Deformation-gradient update: integrate dF/dt = L*F over one step.
        // The former first-order form F_new = (I + L*dt)*F_old is NOT
        // volume-consistent: det(I + L*dt) = 1 + tr(L)*dt + det(L)*dt^2, so
        // even for incompressible flow (tr(L)=0) it carries a spurious
        // det(L)*dt^2 term and det(F) drifts every step (geometric decay for
        // pure shear). We instead use the exact closed-form 2x2 matrix
        // exponential F_new = expm(L*dt)*F_old, for which
        // det(expm(L*dt)) = exp(tr(L)*dt) exactly — the physically correct
        // volume change. Write A = L*dt = m*I + B with m = tr(A)/2 and B
        // traceless (d = det(B)); then expm(A) = exp(m)*expm(B) with
        //   d < 0: q=sqrt(-d), expm(B) = cosh(q)*I + (sinh(q)/q)*B
        //   d > 0: p=sqrt(d),  expm(B) = cos(p)*I  + (sin(p)/p)*B
        //   d ~ 0:             expm(B) = I + B
        const double Axx = pdudx[k]*model.dt;
        const double Axz = pdudz[k]*model.dt;
        const double Azx = pdvdx[k]*model.dt;
        const double Azz = pdvdz[k]*model.dt;
        const double m   = 0.5*(Axx + Azz);          // A = m*I + B
        const double Bxx = Axx - m;                  // traceless part (Bzz = -Bxx)
        const double Bzz = Azz - m;
        const double dB  = Bxx*Bzz - Axz*Azx;        // det(B)
        double cI, sB;                               // expm(B) = cI*I + sB*B
        if ( dB < -1.0e-300 ) {
            const double q = sqrt(-dB);
            cI = cosh(q);  sB = sinh(q)/q;
        }
        else if ( dB > 1.0e-300 ) {
            const double p = sqrt(dB);
            cI = cos(p);   sB = sin(p)/p;
        }
        else {                                       // d ~ 0 (nilpotent B)
            cI = 1.0;      sB = 1.0;
        }
        const double em = exp(m);
        // fij = expm(L*dt) (off-diagonals share the traceless off-diag of A)
        fxx = em*(cI + sB*Bxx);
        fxz = em*(sB*Axz);
        fzx = em*(sB*Azx);
        fzz = em*(cI + sB*Bzz);
        fxx_o = particles->Fxx[k];
        fxz_o = particles->Fxz[k];
        fzx_o = particles->Fzx[k];
        fzz_o = particles->Fzz[k];
        particles->Fxx[k] = fxx*fxx_o + fxz*fzx_o;
        particles->Fxz[k] = fxx*fxz_o + fxz*fzz_o;
        particles->Fzx[k] = fzx*fxx_o + fzz*fzx_o;
        particles->Fzz[k] = fzx*fxz_o + fzz*fzz_o;
    }

    // Free
    DoodzFree( dudx );
    DoodzFree( dvdz );
    DoodzFree( dudz );
    DoodzFree( dvdx );
    DoodzFree( pdudx );
    DoodzFree( pdudz );
    DoodzFree( pdvdx );
    DoodzFree( pdvdz );

    LOG_INFO("-----> Deformation gradient tensor Updated");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FiniteStrainAspectRatio ( grid *mesh, scale scaling, params model, markers *particles, mat_prop *materials ) {

    int k, l, cp, cu, cv;
    double *FS_AR;
    int    cent=1, vert=0, prop=1, interp=0;

    // Largest finite aspect ratio we ever store. The downstream consumer
    // (AnisoFactorEvolv) clips with the per-phase ani_fac_max anyway; this is
    // only a last-resort guard so the function never emits a non-finite value
    // when the smaller singular value underflows for an extremely degenerate F.
    const double FS_AR_CAP = 1.0e12;

    FS_AR = DoodzCalloc (particles->Nb_part, sizeof(double));

    //#pragma omp parallel for shared ( particles ) private ( k ) schedule( static )
    for (k=0; k<particles->Nb_part; k++) {

        // Finite-strain aspect ratio = ratio of the eigenvalues of the right
        // stretch tensor U = sqrt(F^T F), i.e. sigma_max/sigma_min of F.
        //
        // The former implementation formed CG = F^T F, then
        // Det = CGxx*CGzz - CGxz*CGzx and s = sqrt(Det). Det equals det(F)^2
        // (>= 0 in exact arithmetic) but is built as a difference of large
        // nearly-equal products, so under finite precision it goes slightly
        // negative for ill-conditioned F and sqrt() returned a non-finite
        // value; the downstream e1/e2 divide could also divide by ~0.
        //
        // Instead use the stable closed-form 2x2 SVD of F directly:
        //   E = (Fxx+Fzz)/2,  Fd = (Fxx-Fzz)/2
        //   G = (Fzx+Fxz)/2,  H  = (Fzx-Fxz)/2
        //   Q = sqrt(E*E+H*H), R = sqrt(Fd*Fd+G*G)
        //   sigma_max = Q + R,  sigma_min = |Q - R|
        // Each sqrt has a non-negative argument by construction, so this never
        // yields a non-finite value for finite F.
        const double Fxx = particles->Fxx[k];
        const double Fxz = particles->Fxz[k];
        const double Fzx = particles->Fzx[k];
        const double Fzz = particles->Fzz[k];

        const double E  = 0.5*(Fxx + Fzz);
        const double Fd = 0.5*(Fxx - Fzz);
        const double G  = 0.5*(Fzx + Fxz);
        const double H  = 0.5*(Fzx - Fxz);
        const double Q  = sqrt(E*E + H*H);
        const double R  = sqrt(Fd*Fd + G*G);

        const double sigma_max = Q + R;
        const double sigma_min = fabs(Q - R);

        // Guard sigma_min: if it underflows, clamp the ratio to a large finite
        // value instead of returning a non-finite result.
        if ( sigma_min > sigma_max/FS_AR_CAP ) {
            FS_AR[k] = sigma_max/sigma_min;
            if ( FS_AR[k] > FS_AR_CAP ) FS_AR[k] = FS_AR_CAP;
        }
        else {
            FS_AR[k] = FS_AR_CAP;
        }

        // -------------------------------------------------------------------
        // ani_fstrain == 3 — per-marker δ-relaxation update (design D1, D6)
        // -------------------------------------------------------------------
        // Placed HERE, inside the FS_AR k-loop, because the per-marker FS_AR
        // exists only in this local array — there is no marker FS_AR field.
        // Two-field operator split:
        //   delta_FS   = aniso_delta_fn(FS_AR)                  (the ==2 value)
        //   delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev)
        //                                          (inherit the ==2 increment)
        //   delta_new  = 1 + (delta_prod - 1)·exp(-dt/tau_relax)  (relax → 1)
        //   clamp to [1, ani_fac_max]; store aniso_delta, aniso_delta_fs_prev.
        // Cold limit (tau → ∞): exp → 1 ⇒ delta_new = delta_prod ⇒ telescoping
        // collapses (since aniso_delta_0 = aniso_delta_fs_prev_0 = aniso_factor
        // by the init-from-finite-strain helper, design D4 / change aniso-init-
        // from-finite-strain) to
        //   aniso_delta_N → aniso_delta_fn(FS_AR_N)
        // exactly — the mineral-specific Hansen-saturation ceiling holds with
        // no permanent +offset. The "inherited fabric persists in cold static
        // rock" feature survives via the strain-rate gate below: in quiescent
        // cells τ_eff stays finite, so aniso_delta relaxes toward 1, but the
        // gate AND/OR the (zero-deformation) cold-limit means FS_AR stays at
        // its initial value, so δ stays at aniso_factor. Quiescent-hot limit
        // (delta_FS == aniso_delta_fs_prev): pure relaxation
        //   delta_new = 1 + (aniso_delta - 1)·exp(-dt/tau_eff).
        // With the strain-rate gate below, tau_eff → ∞ in actively-creeping cells
        // (gate suppresses DiSRX during deformation, per Boneh+21 static-anneal
        // protocol), so CPO build-up via FS_AR survives in shear bands.
        //
        // ONE-STEP LAG (note, not a defect): at the line-676 anisotropy pass
        // F/FS_AR and strain_pwl are end-of-previous-step values, while
        // particles->T is current. The relaxation is a slow process and this
        // is consistent with how ani_fstrain==2 already consumes FS_AR, so the
        // lag is acceptable.
        const int phase_d = particles->phase[k];
        if ( phase_d >= 0 && materials->ani_fstrain[phase_d] == 3 ) {
            // Production δ from the calibrated ani_fstrain==2 form. aniso_db>0
            // is auto-defaulted for ani_fstrain==3 (InputOutput.c / Main_DOODZ.c),
            // so aniso_delta_fn is populated; guard anyway.
            double delta_FS = 1.0;
            if ( materials->aniso_delta_fn[phase_d] != NULL ) {
                delta_FS = materials->aniso_delta_fn[phase_d]( FS_AR[k] );
            }
            // Resolve the effective L_relax for this phase.
            //   ani_relax_length >= 0  → explicit per-phase length
            //   ani_relax_length  < 0  → sentinel: inherit gs_ref[phase]
            // CAVEAT (constant-L_relax assumption, design D7a): with grain-size
            // evolution OFF, L_relax is a per-phase CONSTANT. The annealed
            // grain size that DiSRX actually grows to is itself T- and time-
            // dependent; holding it fixed introduces a bounded error that is
            // largest in the mid-temperature band and under episodic loading
            // (where the "caught partway" regime lives) — the regime analysis
            // bounds it at the ≤17 %-class on δ_final, never qualitative. With
            // the gs flow law ON, L_relax could instead track mesh->d_n; that
            // GSE-on path is deferred (design O3 — the -1 sentinel keeps it
            // non-breaking to add later).
            const double L_relax = ( materials->ani_relax_length[phase_d] >= 0.0 )
                                   ? materials->ani_relax_length[phase_d]
                                   : materials->gs_ref[phase_d];
            // Boneh-2021 DiSRX kinetics constants — mineral-specific, pulled
            // per-phase from the material database (set in ReadDataAnisotropy,
            // FlowLaws.c; olivine-calibrated there).
            const double tau_relax = DeltaRelaxationTau( particles->T[k], L_relax,
                                                         particles->strain_pwl[k],
                                                         materials->R, scaling,
                                                         materials->ani_relax_Q[phase_d],
                                                         materials->ani_relax_M0[phase_d],
                                                         materials->ani_relax_mu[phase_d],
                                                         materials->ani_relax_b[phase_d],
                                                         materials->ani_relax_drho_min[phase_d],
                                                         materials->ani_relax_drho_max[phase_d],
                                                         materials->ani_relax_eps_ref[phase_d] );
            // STRAIN-RATE GATE on DiSRX relaxation: Boneh-2021's M0/Q come from
            // STATIC-anneal experiments (Karato 89, Cooper-Kohlstedt 84, Speciale
            // 20); concurrent dislocation creep suppresses DiSRX (Hansen+12/14/16
            // build CPO at ε̇ ≈ 1e-5 with no observed DiSRX). The gate slows the
            // relaxation when this cell's dislocation strain rate exceeds the
            // per-phase threshold ani_relax_eps_max (default 1e-13 s^-1):
            //   f_gate = 1 / (1 + (eII_pwl / eps_max)^2),   tau_eff = tau_relax / f_gate
            // f_gate → 1 in quiescent rock (gate disabled, classic DiSRX),
            // f_gate → 0 in actively-deforming rock (tau_eff → ∞, CPO preserved).
            // A sentinel ani_relax_eps_max < 0 disables the gate (legacy behaviour).
            double tau_eff = tau_relax;
            const double eps_max = materials->ani_relax_eps_max[phase_d];
            if ( eps_max > 0.0 ) {
                // nearest-cell lookup of the centroid eII_pwl (scaled units)
                int icx = (int)floor( (particles->x[k] - model.xmin) / model.dx );
                int icz = (int)floor( (particles->z[k] - model.zmin) / model.dz );
                if (icx < 0) icx = 0;     if (icx > model.Nx-2) icx = model.Nx-2;
                if (icz < 0) icz = 0;     if (icz > model.Nz-2) icz = model.Nz-2;
                const int c2_m = icx + icz*(model.Nx-1);
                const double eII = mesh->eII_pwl[c2_m];
                if ( isfinite(eII) && eII >= 0.0 ) {
                    const double ratio  = eII / eps_max;
                    const double f_gate = 1.0 / ( 1.0 + ratio*ratio );
                    if ( f_gate > 1.0e-30 ) tau_eff = tau_relax / f_gate;
                    else                    tau_eff = 1.0e300;          // gate fully closed
                }
            }
            const double delta_prod = particles->aniso_delta[k]
                                      + ( delta_FS - particles->aniso_delta_fs_prev[k] );
            // Analytic exponential relaxation toward the isotropic limit δ_eq = 1.
            // Unconditionally stable for any dt; cannot overshoot below 1.
            double delta_new = 1.0 + ( delta_prod - 1.0 ) * exp( -model.dt / tau_eff );
            // Clamp to [1, ani_fac_max].
            if ( delta_new > materials->ani_fac_max[phase_d] ) delta_new = materials->ani_fac_max[phase_d];
            if ( delta_new < 1.0 )                             delta_new = 1.0;
            particles->aniso_delta[k]         = delta_new;
            particles->aniso_delta_fs_prev[k] = delta_FS;
        }
    }

    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_n,   mesh->BCp.type,  1, 0, interp, cent, model.interp_stencil, NULL);
    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_s,   mesh->BCg.type,  1, 0, interp, vert, model.interp_stencil, NULL);

    // ani_fstrain == 3: P2G the relaxed per-marker δ onto both the centroid
    // and vertex grids, mirroring exactly the two FS_AR P2Mastah calls above.
    // AnisoFactorEvolv's ani_fstrain==3 arm then consumes mesh->aniso_delta_n/_s
    // (passed in as its relaxed_delta argument) — so the relaxed δ flows
    // through the existing per-phase phase_perc weighting and harmonic/
    // geometric averaging in UpdateAnisoFactor.
    P2Mastah( &model, *particles, particles->aniso_delta, mesh, mesh->aniso_delta_n, mesh->BCp.type, 1, 0, interp, cent, model.interp_stencil, NULL);
    P2Mastah( &model, *particles, particles->aniso_delta, mesh, mesh->aniso_delta_s, mesh->BCg.type, 1, 0, interp, vert, model.interp_stencil, NULL);

    DoodzFree(FS_AR);

    LOG_INFO("-----> Finite strain aspect ratio updated");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrain( grid* mesh, scale scaling, params model, markers* particles ) {

    double *strain_inc;
    int Nx, Nz, k, l, c1;
    DoodzFP *strain_inc_mark, *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;

    Nx = mesh->Nx;
    Nz = mesh->Nz;

    // Allocate
    //    exz_n           = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc      = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_el   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pl   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pwl  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_exp  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_lin  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_gbs  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

    // Interpolate exz to cell centers
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c1 = k  + l*(Nx-1);

            strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);
            if (strain_inc[c1]<0) {
                LOG_INFO("negative strain increment");
                LOG_INFO("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
                exit(0);
            }

            strain_inc_el[c1]  = model.dt*mesh->eII_el[c1];
            strain_inc_pl[c1]  = model.dt*mesh->eII_pl[c1];
            strain_inc_pwl[c1] = model.dt*mesh->eII_pwl[c1];
            strain_inc_exp[c1] = model.dt*mesh->eII_exp[c1];
            strain_inc_lin[c1] = model.dt*mesh->eII_lin[c1];
            strain_inc_gbs[c1] = model.dt*mesh->eII_gbs[c1];
        }
    }

    // Interp increments to particles
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_el, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_el, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pl, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pwl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pwl, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_exp, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_exp, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_lin, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_lin, strain_inc_mark, particles->Nb_part );

    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_gbs, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_gbs, strain_inc_mark, particles->Nb_part );

    // Free
    //    DoodzFree(exz_n);
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
    DoodzFree(strain_inc_mark);

    LOG_INFO("-----> Accumulated strain updated");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateMaxPT ( scale scaling, params model, markers *particles ) {

    // Update maximum pressure and temperature on markers
    int k;

#pragma omp parallel for shared ( particles ) private ( k ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {

        if (particles->T[k] > particles->Tmax[k] ) particles->Tmax[k] = particles->T[k];
        if (particles->P[k] > particles->Pmax[k] ) particles->Pmax[k] = particles->P[k];

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleDensity( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *rho_inc_mark, *rho_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    rho_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    rho_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    //----------------------

    //    Interp_Grid2P_centroids2( *particles, particles->rho, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    //----------------------

    for (k=0;k<Ncx*Ncz;k++) {
        rho_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) rho_inc_grid[k] = mesh->rho_n[k] - mesh->rho0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, rho_inc_mark, mesh, rho_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->rho, rho_inc_mark, particles->Nb_part );

    DoodzFree(rho_inc_grid);
    DoodzFree(rho_inc_mark);

    //----------------------

    // for (int k=0; k<particles->Nb_part; k++ ) {
    //     particles->rho[k] = EvaluateDensity( particles->phase[k], particles->T[k], particles->P[k], particles->X[k],  particles->phi[k], &model, materials );
    // }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePhi( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *phi_inc_mark, *phi_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    phi_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    phi_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    for (k=0;k<Ncx*Ncz;k++) {
        phi_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) phi_inc_grid[k] = mesh->phi_n[k] - mesh->phi0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, phi_inc_mark, mesh, phi_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->phi, phi_inc_mark, particles->Nb_part );

    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->phi[k] <= 0.0) particles->phi[k] = 0.0;
        if (particles->phi[k] >= 1.0) particles->phi[k] = 1.0;
    }

    DoodzFree(phi_inc_grid);
    DoodzFree(phi_inc_mark);

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleX( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

//     DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr, *rho_part;
//     double *dT, *Tg0, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
//     int    cent=1, vert=0, prop=1, interp=0;
//     int Nx, Nz, Ncx, Ncz, k, c0, p;
//     Nx = mesh->Nx; Ncx = Nx-1;
//     Nz = mesh->Nz; Ncz = Nz-1;

//     // Allocations
//     Tm0        = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     Tg0        = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dT        = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

//     // Old temperature grid
// #pragma omp parallel for shared(mesh, Tg0) private(c0) firstprivate(Ncx,Ncz)
//     for ( c0=0; c0<Ncx*Ncz; c0++ ) {
//         if (mesh->BCt.type[c0] != 30) {
//             Tg0[c0] = mesh->X0_n[c0];
//             dT[c0]  = mesh->X_n[c0] - mesh->X0_n[c0];
//         }
//     }
//     Interp_Grid2P_centroids2( *particles, Tm0, mesh, Tg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );


//     /* CASE WITHOUT SUBGRID DIFFUSION: INTERPOLATE INCREMENT DIRECTLY */

//         // Interp increments to particles
//         Interp_Grid2P_centroids2( *particles, T_inc_mark, mesh, dT, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

//         // Increment temperature on particles
//         ArrayPlusArray( particles->X, T_inc_mark, particles->Nb_part );


//     // Freedom
//     DoodzFree(dT);
//     DoodzFree(Tg0);
//     DoodzFree(Tm0);
//     DoodzFree(T_inc_mark);

    //--------------------------------------------------
    DoodzFP *X_inc_mark, *X_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // // Total update (do not use!!!)
    // Interp_Grid2P_centroids2( *particles, particles->X, mesh, mesh->X_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

    // Incremental update
    X_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    X_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
#pragma omp parallel for shared ( X_inc_grid, mesh ) private( k ) firstprivate( Ncx, Ncz )
    for (k=0; k<Ncx*Ncz; k++) {
        X_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) X_inc_grid[k] = mesh->X_n[k] - mesh->X0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, X_inc_mark, mesh, X_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

    // Increment X on particles
    ArrayPlusArray( particles->X, X_inc_mark, particles->Nb_part );

    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->X[k] <= 0.0) particles->X[k] = 0.0;
        if (particles->X[k] >= 1.0) particles->X[k] = 1.0;
    }

    DoodzFree(X_inc_grid);
    DoodzFree(X_inc_mark);
}

//--------------------------------------------------
void UpdateParticleXpips( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *X_inc_mark, *X_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    //==================================================================
    // Diffuse Xreac array ---------------------------------------------
    //==================================================================

    //double diff = model.Reac_diff/(scaling->L*scaling->L/scaling->t);
    //double diff = 1.0e-15/(scaling.L*scaling.L/scaling.t);
    double DW, DE, DS, DN;
    double diff_time = model.dt;
    double dx = model.dx;
    double dz = model.dz;

    double X_wanted, valmin, valmax, deltaval, valt_test, Xtest;
    double d=1.0;
    int nb_it, nb_itmax;
    int i,c,l, c1, c2, c3;

    //    get min max kc
    double min_kc=1.0e30, max_kc=0.0;
    for ( k=0; k<materials->Nb_phases; k++) {
        if (materials->k_chem[k]<min_kc) min_kc = materials->k_chem[k];
        if (materials->k_chem[k]>max_kc) max_kc = materials->k_chem[k];
    }
    LOG_INFO("--> min kc = %2.2e m2.s-1 - max kc = %2.2e m2.s-1 ", min_kc*scaling.L*scaling.L/scaling.t, max_kc*scaling.L*scaling.L/scaling.t);


    double *knodes,*newkx,*newkz;

    knodes  = DoodzCalloc(Nx*Nz, sizeof(DoodzFP));
    newkx   = DoodzCalloc(Nx*(Nz+1), sizeof(DoodzFP));
    newkz   = DoodzCalloc((Nx+1)*Nz, sizeof(DoodzFP));

    int interp=0, vert=0;

    P2Mastah( &model, *particles, materials->k_chem, mesh, knodes, mesh->BCg.type,  0, 0, interp, vert, (&model)->interp_stencil, NULL);

    // Interp_P2N ( *particles, materials->k_chem,  mesh, knodes, mesh->xg_coord,  mesh->zg_coord, 0, 0, &model );
    MinMaxArray(knodes, scaling.L*scaling.L/scaling.t, Nx*Nz, "knodes" );

    // fill newkx
    for( l=0; l<Nz+1; l++) {
        for( k=0; k<Nx; k++) {

            c1 = k + l*Nx;     // current value
            c3 = k + (l-1)*Nx; // value below

            if (l==0)          newkx[c1] = knodes[c1];
            if (l==Nz)         newkx[c1] = knodes[c3];
            if (l>0 && l<Nz)   newkx[c1] = 0.5*(knodes[c1] + knodes[c3]);

        }
    }
    MinMaxArray(newkx, scaling.L*scaling.L/scaling.t, Nx*(Nz+1), "newkx" );

    // fill newkz
    for( l=0; l<Nz; l++) {
        for( k=0; k<Nx+1; k++) {

            c1 = k   + l*(Nx+1);     // current value
            c2 = k   + l*(Nx);       // value on nodes
            c3 = k   + l*(Nx) - 1 ;     // value on nodes left

            //            if (k==0)          newkz[c1] = knodes[c2];
            //            if (k==Nx)         newkz[c1] = knodes[c3];

            //          if (k==0 && model.isperiodic_x == 1)    newkz[c1] = knodes[c2 + Nx-1];
            if (k==0 && model.periodic_x == 1)    newkz[c1] = 0.5*(knodes[c2] + knodes[c2+(Nx-2)]);
            if (k==0 && model.periodic_x == 0)    newkz[c1] = knodes[c2];

            //          if (k==Nx && model.isperiodic_x == 1)   newkz[c1] = knodes[c3 - Nx-1];
            if (k==Nx && model.periodic_x == 1)   newkz[c1] = 0.5*(knodes[c3 - (Nx-2)] + knodes[c3]);
            if (k==Nx && model.periodic_x == 0)   newkz[c1] = knodes[c3];

            if (k>0 && k<Nz)   newkz[c1] = 0.5*(knodes[c2] + knodes[c3]);

        }
    }
    MinMaxArray(newkz, scaling.L*scaling.L/scaling.t, (Nx+1)*Nz, "newkz" );



    double dt   = 0.249*model.dx*model.dz/max_kc, time=0.0;
    double dt_special = 0.0;


    LOG_INFO("\n--------");
    LOG_INFO("--> model time  = %f ", model.dt);
    LOG_INFO("--> dt_explicit = %f ", dt);

    if (dt>=diff_time) dt = diff_time;

    int nstep   = floor(diff_time/dt);

    double vN, vS, vE, vW, vC;

    double *temp;

    temp  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    //    printf("\n--------\n");
    LOG_INFO("--> dt_chosen   = %f ", dt);
    //    printf("--> K           = %f \n", diff );
    LOG_INFO("--> nstep       = %d ", nstep);
    //    printf("--> K           = %2.2e m2.s-1 \n", diff*scaling.L*scaling.L/scaling.t );


    if (diff_time - (nstep*dt) > 0.0){
        dt_special = diff_time - (nstep*dt);
    }

    for (i=0; i<=nstep; i++) {

        ArrayEqualArray(temp, mesh->X_n,Ncx * Ncz);

        if (i==nstep){
            time += dt_special;
            dt = dt_special;
        }
        else{
            time += dt;
        }

        for( l=0; l<Ncz; l++) {
            for( k=0; k<Ncx; k++) {

                c = k + l*Ncx;

                c1  = k   + (l+1)*Nx;
                c3  = k   + l*(Nx+1)+1;

                //                DW    = mesh->kc_x[c1]         ;
                //                DE    = mesh->kc_x[c1 + 1]     ;
                //                DS    = mesh->kc_z[c3]         ;
                //                DN    = mesh->kc_z[c3 + Nx+1]  ;

                DW    = newkx[c1]         ;
                DE    = newkx[c1 + 1]     ;
                DS    = newkz[c3]         ;
                DN    = newkz[c3 + Nx+1]  ;


                vC = temp[c];

                if (k==0 && model.periodic_x == 1)     vW = temp[c + Ncx-1];
                if (k==0 && model.periodic_x == 0)     vW = temp[c]        ;
                if (k>0)                                 vW = temp[c - 1]    ;

                if (k==Ncx-1 && model.periodic_x == 1) vE = temp[c-Ncx+1]  ;
                if (k==Ncx-1 && model.periodic_x == 0) vE = temp[c]        ;
                if (k<Ncx-1)                             vE = temp[c+1]      ;

                if (l==0) vS = temp[c];
                if (l>0)  vS = temp[c - Ncx];

                if (l==Ncz-1) vN = temp[c];
                if (l<Ncz-1)  vN = temp[c + Ncx];

                //mesh->X_n[c] += dt*diff*((vE-2.0*vC+vW)/(model.dx*model.dx)+(vN-2.0*vC+vS)/(model.dz*model.dz));
                //mesh->X_n[c] += dt*1.0e-15/(scaling.L*scaling.L/scaling.t)*((vE-2.0*vC+vW)/(model.dx*model.dx)+(vN-2.0*vC+vS)/(model.dz*model.dz));

                mesh->X_n[c] += dt*( DE*(vE-vC)/(dx*dx) - DW*(vC-vW)/(dx*dx) + DN*(vN-vC)/(dz*dz) - DS*(vC-vS)/(dz*dz));

            }
        }
    }

    LOG_INFO("--> end of loops, time = %2.2e s ", time*scaling.t);
    LOG_INFO("--> model time  = %2.2e s ", model.dt*scaling.t);
    LOG_INFO("\n--------");
    MinMaxArrayTag( mesh->X_n,             1.0, Ncx*Ncz, "Xreac_n ", mesh->BCp.type );

    free(temp);
    free(knodes);
    free(newkx);
    free(newkz);
    //==================================================================
    //==================================================================

    // // Total update (do not use!!!)
    //    Interp_Grid2P( *particles, particles->X, mesh, mesh->X_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Incremental update
    X_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    X_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
#pragma omp parallel for shared ( X_inc_grid, mesh ) private( k ) firstprivate( Ncx, Ncz )
    for (k=0; k<Ncx*Ncz; k++) {
        X_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) X_inc_grid[k] = mesh->X_n[k] - mesh->X0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, X_inc_mark, mesh, X_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->X, X_inc_mark, particles->Nb_part );

    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->X[k] <= 0.0) particles->X[k] = 0.0;
        if (particles->X[k] >= 1.0) particles->X[k] = 1.0;
    }

    DoodzFree(X_inc_grid);
    DoodzFree(X_inc_mark);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleGrainSize( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

//     DoodzFP *d_inc_mark, *d_inc_grid;
//     int Nx, Nz, Ncx, Ncz, k;
    const int Nx = mesh->Nx; //Ncx = Nx-1;
    const int Nz = mesh->Nz; //Ncz = Nz-1;

//     d_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     d_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

//     for (k=0;k<Ncx*Ncz;k++) {
//         d_inc_grid[k] =0.0;
//         if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) d_inc_grid[k] = mesh->d_n[k] - mesh->d0_n[k];
//     }

//     // Interp increments to particles
//     Interp_Grid2P_centroids2( *particles, d_inc_mark, mesh, d_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

//     MinMaxArrayPart( particles->d, scaling.L, particles->Nb_part, "d on markers", particles->phase ) ;

//     // Increment grain size on particles
//     ArrayPlusArray( particles->d, d_inc_mark, particles->Nb_part );

//     MinMaxArrayPart( particles->d, scaling.L, particles->Nb_part, "d on markers", particles->phase ) ;


// #pragma omp parallel for shared ( particles )
//     for (k=0;k<particles->Nb_part;k++) {
//         if (particles->d[k]<1.0e-16) particles->d[k] = 1.0e-16;
//     }

//     DoodzFree(d_inc_mark);
//     DoodzFree(d_inc_grid);

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, particles->d, mesh, mesh->d_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleEnergy( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr, *rho_part;
    double *dT, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
    int    cent=1, vert=0, prop=1, interp=0;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // Allocations
    Tm0        = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    dT         = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

    // Old temperature grid
#pragma omp parallel for shared(mesh, dT) private(c0) firstprivate(Ncx,Ncz)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        if (mesh->BCt.type[c0] != 30) {
            dT[c0]  = mesh->T[c0] - mesh->T0_n[c0];
        }
    }
    Interp_Grid2P_centroids2( *particles, Tm0, mesh, mesh->T0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // SUBGRID
    if ( model.subgrid_diffusion >= 1 ) { /* CASE WITH SUBGRID DIFFUSION */

        LOG_INFO("Subgrid diffusion for temperature update");
        dTgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dTmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        rho_part = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        // Get density
        Interp_Grid2P_centroids2( *particles, rho_part, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Tm0,dTms) private(k,p,dtm) firstprivate(materials,dx,dz,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p = particles->phase[k];
                dtm     = materials->Cp[p] * rho_part[k]/ (materials->k[p] * (1.0/dx/dx + 1.0/dz/dz));
                dTms[k] = -( particles->T[k] - Tm0[k] ) * (1.0-exp(-d*model.dt/dtm));
            }
        }

        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dTms,     mesh, dTgs,   mesh->BCp.type,  1, 0, interp, cent, 1, NULL);

        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dTgs, dTgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) dTgr[c0] = dT[c0] - dTgs[c0];
        }

        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dTmr, mesh, dTgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Final temperature update on markers
#pragma omp parallel for shared(particles,dTms,dTmr,T_inc_mark) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) T_inc_mark[k]    = dTms[k] + dTmr[k];
            if (particles->phase[k] != -1) particles->T[k]  = particles->T[k] + T_inc_mark[k];
        }
        DoodzFree(dTms);
        DoodzFree(dTmr);
        DoodzFree(dTgs);
        DoodzFree(dTgr);
        DoodzFree(rho_part);
    }
    else {  /* CASE WITHOUT SUBGRID DIFFUSION: INTERPOLATE INCREMENT DIRECTLY */

        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, T_inc_mark, mesh, dT, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Increment temperature on particles
        ArrayPlusArray( particles->T, T_inc_mark, particles->Nb_part );
    }

    // Freedom
    DoodzFree(dT);
    DoodzFree(Tm0);
    DoodzFree(T_inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePressure( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *P_inc_mark;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    double d=1.0, dtm;
    int    cent=1, vert=0, prop=1, interp=0;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // Compute increment
#pragma omp parallel for shared(mesh) private(c0)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        mesh->dp[c0] = 0.0;
        if (mesh->BCp.type[c0] != 30 ) {
            mesh->dp[c0] = (mesh->p_in[c0]-mesh->p0_n[c0]);
        }
    }

    if ( model.subgrid_diffusion >= 2 ) {

        LOG_INFO("Subgrid diffusion for pressure update");
        double *Pg0  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *Pm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *etam = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        Interp_Grid2P_centroids2( *particles, etam,  mesh, mesh->eta_phys_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1  , Nz-1 , mesh->BCp.type, &model);

        /* -------------- */
        // Old Pressure grid
#pragma omp parallel for shared(mesh, Pg0) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) Pg0[c0] = mesh->p0_n[c0];
        }
        Interp_Grid2P_centroids2( *particles, Pm0, mesh, Pg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        /* -------------- */

        //        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );



        // Old temperature grid --> markers
        //        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Pm0,dPms) private(k,p,dtm) firstprivate(materials,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p         = particles->phase[k];
                dtm       = etam[k]/ (materials->bet[p]);
                dPms[k]   = -( particles->P[k] - Pm0[k]) * (1.0 - exp(-d*model.dt/dtm));
            }
        }

        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dPms,     mesh, dPgs,   mesh->BCp.type,  1, 0, interp, cent, 1, NULL);

        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dPgs, dPgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCp.type[c0] != 30 ) dPgr[c0] = mesh->dp[c0] - dPgs[c0];
        }

        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dPmr, mesh, dPgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Final temperature update on markers
#pragma omp parallel for shared(particles,dPms,dPmr) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) particles->P[k] = particles->P[k] + dPms[k] + dPmr[k];
        }

        DoodzFree(Pg0);
        DoodzFree(Pm0);
        DoodzFree(dPms);
        DoodzFree(dPmr);
        DoodzFree(dPgs);
        DoodzFree(dPgr);
        DoodzFree(etam);
    }
    else {

        P_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, P_inc_mark, mesh, mesh->dp, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Increment pressure on particles
        ArrayPlusArray( particles->P, P_inc_mark, particles->Nb_part );

        DoodzFree(P_inc_mark);

    }

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleStress( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {

    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, p, k1, c3;
    DoodzFP *mdsxxd, *mdszzd, *mdsxz, *dsxxd, *dszzd, *dsxz, d=1.0, dtaum;
    DoodzFP *dtxxgs, *dtzzgs, *dtxzgs, *dtxxgr, *dtzzgr, *dtxzgr, *txxm0, *tzzm0, *txzm0, *dtxxms, *dtzzms, *dtxzms, *dtxxmr, *dtzzmr, *dtxzmr, *etam;
    double *dudx_n, *dvdz_n, *dudz_s, *dvdx_s, *om_s, *om_n, *dudz_n, *dvdx_n, *dudx_s, *dvdz_s;
    double angle, tzz, txx, txz, dx, dz, dt;
    double *txz_n, *txx_s, *tzz_s, *dtxxg0, *dtzzg0, *dtxzg0;
    int    cent=1, vert=0, prop=1, interp=0;

    Nx = model->Nx;
    Nz = model->Nz;
    dx = model->dx;
    dz = model->dz;
    dt = model->dt;

    dudx_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dvdz_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dudz_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dvdx_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dudz_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dvdx_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dudx_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dvdz_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));

#pragma omp parallel for shared ( mesh, om_s, dudz_s, dvdx_s ) \
private ( k, l, k1, c1, c3 )                                   \
firstprivate( model )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        c3 = k + l*(Nx+1);
        if ( mesh->BCg.type[c1] != 30 ) {
            dudz_s[c1] = (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
            dvdx_s[c1] = (mesh->v_in[c3+1 ] - mesh->v_in[c3])/dx;
        }
    }

#pragma omp parallel for shared ( mesh, dudx_n, dvdz_n ) \
private ( k, l, k1, c0, c1, c2 )                         \
firstprivate( model )
    for ( k1=0; k1<(Nx-1)*(Nz-1); k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            dudx_n[c0]  = (mesh->u_in[c1+1+Nx]     - mesh->u_in[c1+Nx] )/dx;
            dvdz_n[c0]  = (mesh->v_in[c2+1+(Nx+1)] - mesh->v_in[c2+1]  )/dz;
        }
    }

    InterpCentroidsToVerticesDouble( dudx_n, dudx_s, mesh, model );
    InterpCentroidsToVerticesDouble( dvdz_n, dvdz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dudz_n, dudz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dvdx_n, dvdx_s, mesh, model );

// Rotate director directly on particles
    if ( model->anisotropy == 1 && model->advection==1) {

#pragma omp parallel for shared( particles, mesh ) firstprivate( dt, model ) private( k )
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                double nx = particles->nx[k];
                double nz = particles->nz[k];
                double mdudx =  Centers2Particle( particles, dudx_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, k, model->periodic_x );
                double mdudz = Vertices2Particle( particles, dudz_s,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, k );
                double mdvdz =  Centers2Particle( particles, dvdz_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, k, model->periodic_x );
                double mdvdx = Vertices2Particle( particles, dvdx_s,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, k );
                particles->nx[k] += dt*(-(mdudx - mdvdz)*nx*nz - mdvdx*nz*nz + mdudz*nx*nx)*nz;
                particles->nz[k] += dt*( (mdudx - mdvdz)*nx*nz + mdvdx*nz*nz - mdudz*nx*nx)*nx;
                double norm = sqrt( pow(particles->nx[k],2) + pow(particles->nz[k],2));
                particles->nx[k] /= norm;
                particles->nz[k] /= norm;
            }
        }
    }

    // Marker stress update
    if ( model->elastic==1 ) {

        Nx = mesh->Nx; Ncx = Nx-1;
        Nz = mesh->Nz; Ncz = Nz-1;

        if ( model->subgrid_diffusion == 2 ) {

            // Alloc
            dtxxgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtzzgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtxzgs = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
            dtxxgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtzzgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtxzgr = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
            txxm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            tzzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            txzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxxms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtzzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxxmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtzzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            etam   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

            // Old stresses grid --> markers
            Interp_Grid2P_centroids2( *particles, txxm0, mesh, mesh->sxxd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P_centroids2( *particles, tzzm0, mesh, mesh->szzd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P( *particles, txzm0, mesh, mesh->sxz0,       mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);
            Interp_Grid2P( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);


            MinMaxArray(mesh->eta_phys_s, scaling->eta, Nx*Nz, "eta grid  ");
            MinMaxArrayTag( mesh->eta_phys_s, scaling->eta, Nx*Nz,   "eta_s", mesh->BCg.type );
            MinMaxArray(etam, scaling->eta, particles->Nb_part, "eta phys part  ");

            // Interp_Grid2P_eta( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);


            if ( model->subgrid_diffusion == 2 ) {

                LOG_INFO("Subgrid diffusion for stress tensor component update");

                // Compute subgrid stress increments on markers
    #pragma omp parallel for shared(particles,txxm0,tzzm0,txzm0,dtxxms,dtzzms,dtxzms,etam) private(k,p,dtaum) firstprivate(materials,model,d)
                for ( k=0; k<particles->Nb_part; k++ ) {
                    if (particles->phase[k] != -1) {
                        p         = particles->phase[k];
                        dtaum     = etam[k] / materials->G[p];
                        dtxxms[k] = -( particles->sxxd[k] - txxm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        dtzzms[k] = -( particles->szzd[k] - tzzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        dtxzms[k] = -( particles->sxz[k]  - txzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        if (isinf(dtxxms[k])) {
                            LOG_INFO("Inf dtxxms[k]: sxxd=%2.2e txxm0=%2.2e exp(-d*dt/dtaum)=%2.2e d=%2.2e dt=%2.2e dtaum=%2.2e etam=%2.2e G=%2.2e", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum), d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S);
                            LOG_INFO("%2.2e %2.2e %2.2e %2.2e %2.2e", d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S);
                            exit(1);
                        }
                        if (isnan(dtxxms[k])) {
                            LOG_INFO("NaN dtxxms[k]: sxxd=%2.2e txxm0=%2.2e exp(-d*dt/dtaum)=%2.2e d=%2.2e dt=%2.2e dtaum=%2.2e etam=%2.2e G=%2.2e", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum), d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S);
                            exit(1);
                        }
                    }
                }

                // Subgrid stress increments markers --> grid
                P2Mastah( model, *particles, dtxxms,     mesh, dtxxgs,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil, NULL);
                P2Mastah( model, *particles, dtzzms,     mesh, dtzzgs,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil, NULL);
                P2Mastah( model, *particles, dtxzms,     mesh, dtxzgs,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil, NULL);

                // Remaining stress increments on the grid
    #pragma omp parallel for shared(mesh,dtxxgs,dtxxgr,dtzzgs,dtzzgr) private(c0) firstprivate(Ncx,Ncz)
                for ( c0=0; c0<Ncx*Ncz; c0++ ) {
                    if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtxxgr[c0] = (mesh->sxxd[c0] - mesh->sxxd0[c0]) - dtxxgs[c0];
                    if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtzzgr[c0] = (mesh->szzd[c0] - mesh->szzd0[c0]) - dtzzgs[c0];
                }
    #pragma omp parallel for shared(mesh,dtxzgs,dtxzgr) private(c0) firstprivate(Nx,Nz)
                for ( c0=0; c0<Nx*Nz; c0++ ) {
                    if (mesh->BCg.type[c0]!=30) dtxzgr[c0] = (mesh->sxz[c0]-mesh->sxz0[c0]) - dtxzgs[c0];
                }

                // Remaining stress increments grid --> markers
                Interp_Grid2P_centroids2( *particles, dtxxmr, mesh, dtxxgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
                Interp_Grid2P_centroids2( *particles, dtzzmr, mesh, dtzzgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
                Interp_Grid2P(           *particles, dtxzmr, mesh, dtxzgr, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type         );

                // Final stresses update on markers
    #pragma omp parallel for shared(particles,dtxxms,dtzzms,dtxzms,dtxxmr,dtzzmr,dtxzmr) private(k)
                for ( k=0; k<particles->Nb_part; k++ ) {
                    if (particles->phase[k] != -1) particles->sxxd[k]  = particles->sxxd[k] + dtxxms[k] + dtxxmr[k];
                    if (particles->phase[k] != -1) particles->szzd[k]  = particles->szzd[k] + dtzzms[k] + dtzzmr[k];
                    if (particles->phase[k] != -1) particles->sxz[k]   = particles->sxz[k]  + dtxzms[k] + dtxzmr[k];
                }
            }

            // Free
            DoodzFree( dtxxgs );
            DoodzFree( dtzzgs );
            DoodzFree( dtxzgs );
            DoodzFree( dtxxgr );
            DoodzFree( dtzzgr );
            DoodzFree( dtxzgr );
            DoodzFree( txxm0  );
            DoodzFree( tzzm0  );
            DoodzFree( txzm0  );
            DoodzFree( dtxxms );
            DoodzFree( dtzzms );
            DoodzFree( dtxzms );
            DoodzFree( dtxxmr );
            DoodzFree( dtzzmr );
            DoodzFree( dtxzmr );
            DoodzFree( etam   );
        }

        if (model->subgrid_diffusion==0 || model->subgrid_diffusion==1){

            LOG_INFO("No subgrid diffusion for stress tensor component update");

            // Alloc
            dsxxd  = DoodzCalloc((Nx-1)*(Nz-1), sizeof(double));
            dszzd  = DoodzCalloc((Nx-1)*(Nz-1), sizeof(double));
            dsxz   = DoodzCalloc((Nx)*(Nz), sizeof(double));
            mdsxxd = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            mdszzd = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            mdsxz  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

            // Cell: normal stress change
            for (k=0; k<Nx-1; k++) {
                for (l=0; l<Nz-1; l++) {
                    c0 = k  + l*(Nx-1);
                    if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dsxxd[c0] = mesh->sxxd[c0] - mesh->sxxd0[c0];
                    if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dszzd[c0] = mesh->szzd[c0] - mesh->szzd0[c0];
                }
            }

            // Vertex: shear stress change
            for (k=0; k<Nx; k++) {
                for (l=0; l<Nz; l++) {
                    c1 = k  + l*(Nx);
                    if (mesh->BCg.type[c1] !=30 ) dsxz[c1] =  mesh->sxz[c1] - mesh->sxz0[c1];
                }
            }

            // Interpolate stress changes to markers
            Interp_Grid2P_centroids2( *particles, mdsxxd, mesh, dsxxd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P_centroids2( *particles, mdszzd, mesh, dszzd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P( *particles, mdsxz,  mesh, dsxz,  mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );

            // Update marker stresses
            ArrayPlusArray(  particles->sxxd, mdsxxd, particles->Nb_part );
            ArrayPlusArray(  particles->szzd, mdszzd, particles->Nb_part );
            ArrayPlusArray(  particles->sxz,  mdsxz,  particles->Nb_part );

            // Free
            DoodzFree(dsxxd);
            DoodzFree(dszzd);
            DoodzFree(dsxz);
            DoodzFree(mdsxxd);
            DoodzFree(mdszzd);
            DoodzFree(mdsxz);
        }

        // Large strain: stress rotation/deformation
        if (model->elastic==1) {

            double *wxzm   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            P2Mastah( model, *particles, wxzm,     mesh, mesh->wxz,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil, NULL);

#pragma omp parallel for shared( particles, wxzm ) firstprivate( dt, model ) private( k, angle, txx, tzz, txz )
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) {
                    double txx = particles->sxxd[k];
                    double tzz = particles->szzd[k];
                    double txz = particles->sxz[k];
                    if (model->stress_rotation==1) { // Jaumann rate
                        particles->sxxd[k] -=  dt*2.0*txz*wxzm[k];
                        particles->szzd[k] -= -dt*2.0*txz*wxzm[k];
                        particles->sxz[k]  -=  dt*(tzz-txx)*wxzm[k];
                    }
                    if (model->stress_rotation==2) { // Analytical rotation
                        double angle = dt*wxzm[k];
                        particles->sxxd[k] =  ( txx*cos(angle) + txz*sin(angle))*cos(angle) + ( txz*cos(angle) + tzz*sin(angle))*sin(angle);
                        particles->szzd[k] = -(-txx*sin(angle) + txz*cos(angle))*sin(angle) + (-txz*sin(angle) + tzz*cos(angle))*cos(angle);
                        particles->sxz[k]  = -( txx*cos(angle) + txz*sin(angle))*sin(angle) + ( txz*cos(angle) + tzz*sin(angle))*cos(angle);
                    }
                }
            }
            DoodzFree(wxzm);
        }
    }

    DoodzFree(dudx_n);
    DoodzFree(dvdz_n);
    DoodzFree(dvdx_s);
    DoodzFree(dudz_s);
    DoodzFree(dudz_n);
    DoodzFree(dvdx_n);
    DoodzFree(dvdz_s);
    DoodzFree(dudx_s);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleDivThermal( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *inc_mark, *inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    for (k=0;k<Ncx*Ncz;k++) {
        inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) inc_grid[k] = mesh->divth_n[k] - mesh->divth0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, inc_mark, mesh, inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->divth, inc_mark, particles->Nb_part );

    DoodzFree(inc_grid);
    DoodzFree(inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateAdvectionMode( scale scaling, params* model ) {

    double real_time = model->time;
    double t_wanted = model->time_advection_start;

    LOG_INFO("real_time = %2.2e and t_wanted = %2.2e", real_time * scaling.t, t_wanted * scaling.t);

    if (real_time>t_wanted){
        LOG_INFO("real_time > t_wanted => model advection set to 1!");
        model->advection = 1;
        LOG_INFO("model.advection = %d", model->advection);
    }
    else {
        LOG_INFO("real_time < t_wanted => model advection set to 0!");
        model->advection = 0;
        LOG_INFO("model.advection = %d", model->advection);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePhase( grid* mesh, scale scaling, params* model, markers* particles, mat_prop* materials ) {

    int k;

    double real_time = model->time;
    double t_wanted = model->time_switch_phase;

    LOG_INFO("real_time = %2.2e and t_wanted = %2.2e", real_time * scaling.t, t_wanted * scaling.t);

    if (real_time>t_wanted){
    LOG_INFO("real_time > t_wanted => Phase switch!");

    // Change phase if t > t_wanted
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->phase[k] == model->phase_switch_init) particles->phase[k] = model->phase_switch_final;
        if (particles->dual[k] == model->phase_switch_init) particles->dual[k] = model->phase_switch_final;
    }

    }

}