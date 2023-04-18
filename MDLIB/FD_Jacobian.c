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
#include "stdbool.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "cholmod.h"
#include "mdoodz-private.h"

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

void PhaseRheologyLoop_v1( int is_centroid, double sign, double denom, double Exx, double Ezz, double Exz, double P, double ani_fstrain, double ani_fac_e, double d0, double d1, double angle, double lx2, double lxlz, int c, double** vol,
               double* G, double* T, double* P0, double T0, double* gs0, double* phi0, double* X0, double* txx0, double* tzz0, double* txz0, double* beta, double* div,
               double* strain, double* dil, double* fric, double* C,
               params* model, mat_prop* materials, scale* scaling,
               double* detadE, double* ddivpdE, double* danidE, double* drhodP, double**eta_phase ) {
        
    double txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, detadp, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Wtot, Wel, Wdiss;
    double ani_vep=1.0;
    // Loop on phases
    for (int p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(vol[p][c])>min_fraction ) is_phase_active = true;

        if ( is_phase_active==true ) {

            if (model->anisotropy==0) ViscosityConcise(      p,                                           G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 0, is_centroid, 0 );
            if (model->anisotropy==1) ViscosityConciseAniso( p, lxlz, lx2, angle, ani_fstrain, G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 0, is_centroid, 0 );
            if ( model->eta_average == ARITHMETIC) {
                *detadE      += sign*vol[p][c] * eta_vep/(denom); 
            }
            
            if ( model->eta_average == HARMONIC) {
                *detadE      += sign*vol[p][c] * eta_vep/(denom) / pow(eta_phase[p][c],2.0);
            }
            
            if ( model->eta_average == GEOMETRIC) {
                *detadE      += sign*vol[p][c] * eta_vep/(denom) / eta_phase[p][c];
            }

            // Anisotropic factor
            *danidE     += sign*vol[p][c] * ani_vep/(denom);

            // Plastic divergence (default arithmetic average)
            if (is_centroid > 0) { // if centroid, compute all derivatives of divp
                *ddivpdE     += sign*vol[p][c] * div_pl/(denom);
            }

            // Density (default arithmetic average)
            if (is_centroid == 2) { // if centroid and pressure is varied compute drhodp
                *drhodP      += sign*vol[p][c] * rho/(denom);
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DerivativesOnTheFly_n( double* detadexx, double* detadezz, double* detadgxz, double* detadp, double* ddivpdexx, double* ddivpdezz, double* ddivpdgxz, double* ddivpdp, double* danidexx, double* danidezz, double* danidgxz, double* danidp, double* drhodp, int c0, double Exx_ref, double Ezz_ref, double Exz_ref, double P_ref, double ani, double ani_e, double d1, double d2, double angle, double lx2, double lxlz, grid *mesh, mat_prop *materials, params *model, scale *scaling ) {
    const double tol = 1e-5;
    double pert_xx   = tol*fabs(Exx_ref) + tol/1e10;
    double pert_zz   = tol*fabs(Ezz_ref) + tol/1e10;
    double pert_xz   = tol*fabs(Exz_ref) + tol/1e10;
    double pert_p    = tol*fabs(P_ref)   + tol/1e10;

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Exx
    PhaseRheologyLoop_v1(  1, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadexx, ddivpdexx, danidexx, NULL, mesh->phase_eta_n );

    // 2) Negative perturbation in Exx
    PhaseRheologyLoop_v1( 1, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadexx, ddivpdexx, danidexx, NULL, mesh->phase_eta_n );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Ezz
    PhaseRheologyLoop_v1( 1, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadezz, ddivpdezz, danidezz, NULL, mesh->phase_eta_n );

    // 2) Negative perturbation in Ezz
    PhaseRheologyLoop_v1( 1, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadezz, ddivpdezz, danidezz, NULL, mesh->phase_eta_n );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
    PhaseRheologyLoop_v1(  1, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadgxz, ddivpdgxz, danidgxz, NULL, mesh->phase_eta_n );

    // 2) Negative perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
    PhaseRheologyLoop_v1( 1, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadgxz, ddivpdgxz, danidgxz, NULL, mesh->phase_eta_n );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in P
    PhaseRheologyLoop_v1(  2, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadp, ddivpdp, danidp, drhodp, mesh->phase_eta_n );
    // 2) Negative perturbation in P
    PhaseRheologyLoop_v1( 2, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, ani_e, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                        mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                        mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                        mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                        mesh->bet_n, mesh->div_u,
                        mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                        model, materials, scaling,
                        detadp, ddivpdp, danidp, drhodp, mesh->phase_eta_n );

    // ----------------------------------------------------------------------------------------------------------------------------------------------------//

    if ( model->eta_average == HARMONIC ) {
        *detadexx *= pow(mesh->eta_n[c0] , 2);
        *detadezz *= pow(mesh->eta_n[c0] , 2);
        *detadgxz *= pow(mesh->eta_n[c0] , 2);
        *detadp   *= pow(mesh->eta_n[c0] , 2);
    }
    if ( model->eta_average == GEOMETRIC ) {
        *detadexx *= mesh->eta_n[c0];
        *detadezz *= mesh->eta_n[c0];
        *detadgxz *= mesh->eta_n[c0];
        *detadp   *= mesh->eta_n[c0];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DerivativesOnTheFly_s( double* detadexx, double* detadezz, double* detadgxz, double* detadp, double* danidexx, double* danidezz, double* danidgxz, double* danidp, int c1, double Exx_ref, double Ezz_ref, double Exz_ref, double P_ref, double ani, double ani_e, double d1, double d2, double angle, double lx2, double lxlz, grid *mesh, mat_prop *materials, params *model, scale *scaling ) {
    const double tol = 1e-7;
    const double pert_xx   = tol*fabs(Exx_ref) + tol/1e10;
    const double pert_zz   = tol*fabs(Ezz_ref) + tol/1e10;
    const double pert_xz   = tol*fabs(Exz_ref) + tol/1e10;
    const double pert_p    = tol*fabs(P_ref)   + tol/1e10;

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Exx
    PhaseRheologyLoop_v1(  0, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadexx, NULL, danidexx, NULL, mesh->phase_eta_s );

    // 2) Negative perturbation in Exx
    PhaseRheologyLoop_v1( 0, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadexx, NULL, danidexx, NULL, mesh->phase_eta_s );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Ezz
    PhaseRheologyLoop_v1(  0, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadezz, NULL, danidezz, NULL, mesh->phase_eta_s );

    // 2) Negative perturbation in Ezz
    PhaseRheologyLoop_v1( 0, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadezz, NULL, danidezz, NULL, mesh->phase_eta_s );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//

    // 1) Positive perturbation in Exz
    PhaseRheologyLoop_v1(  0, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadgxz, NULL, danidgxz, NULL, mesh->phase_eta_s );

    // 2) Negative perturbation in Exz
    PhaseRheologyLoop_v1( 0, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadgxz, NULL, danidgxz, NULL, mesh->phase_eta_s );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//
    // 1) Positive perturbation in P
    PhaseRheologyLoop_v1(  0, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadp, NULL, danidp, NULL, mesh->phase_eta_s );

    // 2) Negative perturbation in P
    PhaseRheologyLoop_v1( 0, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, ani_e, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                        mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                        mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                        mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                        mesh->bet_s, mesh->div_u_s,
                        mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                        model, materials, scaling,
                        detadp, NULL, danidp, NULL, mesh->phase_eta_s );

    //----------------------------------------------------------------------------------------------------------------------------------------------------//
    if ( model->eta_average == HARMONIC ) {
        *detadexx *= pow(mesh->eta_s[c1] , 2);
        *detadezz *= pow(mesh->eta_s[c1] , 2);
        *detadgxz *= pow(mesh->eta_s[c1] , 2);
        *detadp   *= pow(mesh->eta_s[c1] , 2);
    }
    if ( model->eta_average == GEOMETRIC ) {
        *detadexx *= mesh->eta_s[c1];
        *detadezz *= mesh->eta_s[c1];
        *detadgxz *= mesh->eta_s[c1];
        *detadp   *= mesh->eta_s[c1];
    }
}
/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/