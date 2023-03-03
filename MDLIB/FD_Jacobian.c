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

            if (model->anisotropy==0) ViscosityConcise(      p,                                           G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 0, is_centroid );
            if (model->anisotropy==1) ViscosityConciseAniso( p, lxlz, lx2, angle, ani_fstrain, ani_fac_e, G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &ani_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 0, is_centroid );
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
    const double tol = 1e-7;
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

#if 0 // To be deleted

void PhaseRheologyLoop( int is_centroid, double sign, double denom, double Exx, double Ezz, double Exz, double P, double ani_fstrain, double d0, double d1, double angle, double lx2, double lxlz, int c, double** vol,
               double* G, double* T, double* P0, double T0, double* gs0, double* phi0, double* X0, double* txx0, double* tzz0, double* txz0, double* beta, double* div,
               double* strain, double* dil, double* fric, double* C,
               params* model, mat_prop* materials, scale* scaling,
               double* detadE, double* ddivpdE, double* drhodP, double**eta_phase ) {
        
    double txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, dnew, detadp, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r;

    // Loop on phases
    for (int p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(vol[p][c])>min_fraction ) is_phase_active = true;

        if ( is_phase_active==true ) {

            if (model->anisotropy==0) ViscosityConcise(      p,                                      G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep,   &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, 0, is_centroid );
            // if (model->anisotropy==1) ViscosityConciseAniso( p, lxlz, lx2, angle, ani_fstrain, G[c], T[c], P, gs0[c], phi0[c], X0[c], Exx, Ezz, Exz, txx0[c], tzz0[c], txz0[c], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &ani_vep, &ani_e, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &dnew, strain[c], dil[c], fric[c], C[c], P0[c], T0, &Xreac, &OverS, &Pcorr, &rho, beta[c], div[c], &div_el, &div_pl, &div_r, 0, is_centroid );
            if ( model->eta_average == ARITHMETIC) {
                detadE[c]      += sign*vol[p][c] * eta_vep/(denom); 
            }
            
            if ( model->eta_average == HARMONIC) {
                detadE[c]      += sign*vol[p][c] * eta_vep/(denom) / pow(eta_phase[p][c],2.0);
            }
            
            if ( model->eta_average == GEOMETRIC) {
                detadE[c]      += sign*vol[p][c] * eta_vep/(denom) / eta_phase[p][c];
            }

            // Plastic divergence (default arithmetic average)
            if (is_centroid > 0) { // if centroid, compute all derivatives of divp
                ddivpdE[c]     += sign*vol[p][c] * div_pl/(denom);
            }

            // Density (default arithmetic average)
            if (is_centroid == 2) { // if centroid and pressure is varied compute drhodp
                drhodP[c]      += sign*vol[p][c] * rho/(denom);
    //            if ( drhodP[c]<0.0 ) { //&&
    //                printf("cell %d, perturbed rho = %2.10e, P = %2.10e phase = %d, drhodP = %2.2e, P = %2.10e --- denom = %2.2e \n", c, rho*scaling->rho, P*scaling->S, p, drhodP[c], P, denom);
    //            }
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ViscosityDerivatives0( grid *mesh, mat_prop *materials, params *model, scale *scaling ) {

    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, eta_vep, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss;
    int average = model->eta_average, unsplit_diff_reac = model->unsplit_diff_reac;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double Pcorr, rho;
    double Exx, Ezz, Exz, el, etae, ani, d1, d2, angle, lxlz, lx2;
    double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
    double tol = 1e-7;

    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    InterpCentroidsToVerticesDouble( mesh->div_u,   mesh->div_u_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
    InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, el, etae, ani, d1, d2, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2 ) firstprivate( unsplit_diff_reac, materials, scaling, average, model, Ncx, Ncz, tol )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {

        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);

        // Initialise arrays to 0
        mesh->detadexx_n[c0]       = 0.0;
        mesh->ddivpdexx_n[c0]      = 0.0;
        mesh->detadezz_n[c0]       = 0.0;
        mesh->ddivpdezz_n[c0]      = 0.0;
        mesh->detadgxz_n[c0]       = 0.0;
        mesh->ddivpdgxz_n[c0]      = 0.0;
        mesh->detadp_n[c0]         = 0.0;
        mesh->ddivpdp_n[c0]        = 0.0;
        mesh->drhodp_n[c0]         = 0.0;

        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {
            
            //----------------------------------------------------------//
            if ( model->elastic==1 ) etae      = model->dt*mesh->mu_n[c0];
            else                       etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->anisotropy == 0 ) {
                ani = 0.0; d1   = 0.0; d2   = 0.0; angle = 0.0; lxlz = 0.0; lx2 = 0.0;
            }
            else {
                // Anisotropy
                // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
                // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_n[c0];
                d1      = mesh->d1_n[c0];
                d2      = mesh->d2_n[c0];
                angle   = mesh->angle_n[c0];
                lxlz    = 0.5*d1;
                lx2     = pow( cos(angle), 2);
            }
            //----------------------------------------------------------//
//            Exx = mesh->exxd[c0]  + mesh->sxxd0[c0] /etae/2.0;
//            Ezz = mesh->ezzd[c0]  + mesh->szzd0[c0] /etae/2.0;
//            Exz = mesh->exz_n[c0] + mesh->sxz0_n[c0]/etae/2.0;
//            gxz = 2.0*Exz;
            
            Da11  = 2.0 - 2.0*ani*d1;
            Da12  = 2.0*ani*d1;
            Da13  = 2.0*ani*d2;
            Da22  = 2.0 - 2.0*ani*d1;
            Da23  =-2.0*ani*d2;
            Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;
            
            double Exx_ref   = mesh->exxd[c0]      + (iDa11*mesh->sxxd0[c0] + iDa12*mesh->szzd0[c0] + iDa13*mesh->sxz0_n[c0])/etae;
            double Ezz_ref   = mesh->ezzd[c0]      + (iDa12*mesh->sxxd0[c0] + iDa22*mesh->szzd0[c0] + iDa23*mesh->sxz0_n[c0])/etae;
            double Exz_ref   = mesh->exz_n[c0]     + (iDa13*mesh->sxxd0[c0] + iDa23*mesh->szzd0[c0] + iDa33*mesh->sxz0_n[c0])/2.0/etae;
            double P_ref     = mesh->p_in[c0];
            double gxz_ref   = 2.0*Exz_ref;
            double pert_xx   = tol*fabs(Exx_ref) + tol/1e10;
            double pert_zz   = tol*fabs(Ezz_ref) + tol/1e10;
            double pert_xz   = tol*fabs(Exz_ref) + tol/1e10;
            double pert_p    = tol*fabs(P_ref)   + tol/1e10;
            
//            printf("exxd %2.2e ezzd %2.2e eyyd %2.2e \n", Exx_ref, Ezz_ref, - Exx_ref - Ezz_ref);
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Exx
            PhaseRheologyLoop(  1, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadexx_n, mesh->ddivpdexx_n, NULL, mesh->phase_eta_n );
            
            // 2) Negative perturbation in Exx
            PhaseRheologyLoop( 1, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadexx_n, mesh->ddivpdexx_n, NULL, mesh->phase_eta_n );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Ezz
            PhaseRheologyLoop( 1, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadezz_n, mesh->ddivpdezz_n, NULL, mesh->phase_eta_n );

            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 1, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadezz_n, mesh->ddivpdezz_n, NULL, mesh->phase_eta_n );

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
            PhaseRheologyLoop(  1, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadgxz_n, mesh->ddivpdgxz_n, NULL, mesh->phase_eta_n );

            // 2) Negative perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
            PhaseRheologyLoop( 1, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadgxz_n, mesh->ddivpdgxz_n, NULL, mesh->phase_eta_n );

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in P
            PhaseRheologyLoop(  2, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadp_n, mesh->ddivpdp_n, mesh->drhodp_n, mesh->phase_eta_n );

            // 2) Negative perturbation in P
            PhaseRheologyLoop( 2, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, d1, d2, angle, lx2, lxlz, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n, mesh->T0_n[c0],
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadp_n, mesh->ddivpdp_n, mesh->drhodp_n, mesh->phase_eta_n );
            // if ( mesh->drhodp_n[c0]<0.0 ) {
            //     printf("%2.2e\n", mesh->drhodp_n[c0]);
            //     exit(1);
            // }

            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            if ( model->eta_average == HARMONIC ) {
                mesh->detadexx_n[c0] *= pow(mesh->eta_n[c0] , 2);
                mesh->detadezz_n[c0] *= pow(mesh->eta_n[c0] , 2);
                mesh->detadgxz_n[c0] *= pow(mesh->eta_n[c0] , 2);
                mesh->detadp_n[c0]   *= pow(mesh->eta_n[c0] , 2);
            }
            if ( model->eta_average == GEOMETRIC ) {
                mesh->detadexx_n[c0] *= mesh->eta_n[c0];
                mesh->detadezz_n[c0] *= mesh->eta_n[c0];
                mesh->detadgxz_n[c0] *= mesh->eta_n[c0];
                mesh->detadp_n[c0]   *= mesh->eta_n[c0];
            }
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
        }
    }

    #pragma omp parallel for shared( mesh ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, el, etae, ani, d1, d2, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2 ) firstprivate( unsplit_diff_reac, materials, scaling, average, model, Nx, Nz, tol )
    for ( k1=0; k1<Nx*Nz; k1++ ) {

        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;

        // Initialise arrays to 0
        mesh->detadexx_s[c1]       = 0.0;
        mesh->detadezz_s[c1]       = 0.0;
        mesh->detadgxz_s[c1]       = 0.0;
        mesh->detadp_s[c1]         = 0.0;

        if ( mesh->BCg.type[c1] != 30 ) {

            //----------------------------------------------------------//
            if ( model->elastic==1   ) etae      = model->dt*mesh->mu_s[c1];
            else           etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->anisotropy == 0 ) {
                ani = 0.0; d1   = 0.0; d2   = 0.0; angle = 0.0; lxlz = 0.0; lx2 = 0.0;
            }
            else {
                // Anisotropy
                // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
                // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_s[c1];
                d1      = mesh->d1_s[c1];
                d2      = mesh->d2_s[c1];
                angle   = mesh->angle_s[c1];
                lxlz    = 0.5*d1;
                lx2     = pow( cos(angle), 2);
            }
            //----------------------------------------------------------//
            //            Exx = mesh->exxd_s[c1] + mesh->sxxd0_s[c1]/etae/2.0;
            //            Ezz = mesh->ezzd_s[c1] + mesh->szzd0_s[c1]/etae/2.0;
            //            Exz = mesh->exz[c1]    + mesh->sxz0[c1]   /etae/2.0;
            //            gxz = 2.0*Exz;

            Da11  = 2.0 - 2.0*ani*d1;
            Da12  = 2.0*ani*d1;
            Da13  = 2.0*ani*d2;
            Da22  = 2.0 - 2.0*ani*d1;
            Da23  =-2.0*ani*d2;
            Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;
            double Exx_ref = mesh->exxd_s[c1]  + (iDa11*mesh->sxxd0_s[c1] + iDa12*mesh->szzd0_s[c1] + iDa13*mesh->sxz0[c1])/etae;
            double Ezz_ref = mesh->ezzd_s[c1]  + (iDa12*mesh->sxxd0_s[c1] + iDa22*mesh->szzd0_s[c1] + iDa23*mesh->sxz0[c1])/etae;
            double Exz_ref = mesh->exz[c1]     + (iDa13*mesh->sxxd0_s[c1] + iDa23*mesh->szzd0_s[c1] + iDa33*mesh->sxz0[c1])/2.0/etae;
            double P_ref     = mesh->P_s[c1];
            double gxz_ref   = 2.0*Exz_ref;
            double pert_xx   = tol*fabs(Exx_ref) + tol/1e10;
            double pert_zz   = tol*fabs(Ezz_ref) + tol/1e10;
            double pert_xz   = tol*fabs(Exz_ref) + tol/1e10;
            double pert_p    = tol*fabs(P_ref)   + tol/1e10;

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Exx
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadexx_s, NULL, NULL, mesh->phase_eta_s );

            // 2) Negative perturbation in Exx
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadexx_s, NULL, NULL, mesh->phase_eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Ezz
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadezz_s, NULL, NULL, mesh->phase_eta_s );
            
            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadezz_s, NULL, NULL, mesh->phase_eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Exz
            PhaseRheologyLoop(  0, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadgxz_s, NULL, NULL, mesh->phase_eta_s );
        
            // 2) Negative perturbation in Exz
            PhaseRheologyLoop( 0, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadgxz_s, NULL, NULL, mesh->phase_eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            // 1) Positive perturbation in P
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadp_s, NULL, NULL, mesh->phase_eta_s );
            
            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, d1, d2, angle, lx2, lxlz, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s, 0.0,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadp_s, NULL, NULL, mesh->phase_eta_s );

            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            if ( model->eta_average == HARMONIC ) {
                mesh->detadexx_s[c1] *= pow(mesh->eta_s[c1] , 2);
                mesh->detadezz_s[c1] *= pow(mesh->eta_s[c1] , 2);
                mesh->detadgxz_s[c1] *= pow(mesh->eta_s[c1] , 2);
                mesh->detadp_s[c1]   *= pow(mesh->eta_s[c1] , 2);
            }
            if ( model->eta_average == GEOMETRIC ) {
                mesh->detadexx_s[c1] *= mesh->eta_s[c1];
                mesh->detadezz_s[c1] *= mesh->eta_s[c1];
                mesh->detadgxz_s[c1] *= mesh->eta_s[c1];
                mesh->detadp_s[c1]   *= mesh->eta_s[c1];
            }
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
        }

    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ViscosityDerivatives( grid *mesh, mat_prop *materials, params *model, scale *scaling ) {

    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, eta_vep, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss;
    int average = model->eta_average, unsplit_diff_reac = model->unsplit_diff_reac;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double Pcorr, rho;
    double Exx, Ezz, Exz, el, etae, ani, d1, d2, angle, lxlz, lx2;
    double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
    double tol = 1e-7;

    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    InterpCentroidsToVerticesDouble( mesh->div_u,   mesh->div_u_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
    InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, el, etae, ani, d1, d2, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2 ) firstprivate( unsplit_diff_reac, materials, scaling, average, model, Ncx, Ncz, tol )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {

        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);

        // Initialise arrays to 0
        mesh->detadexx_n[c0]       = 0.0;
        mesh->ddivpdexx_n[c0]      = 0.0;
        mesh->detadezz_n[c0]       = 0.0;
        mesh->ddivpdezz_n[c0]      = 0.0;
        mesh->detadgxz_n[c0]       = 0.0;
        mesh->ddivpdgxz_n[c0]      = 0.0;
        mesh->detadp_n[c0]         = 0.0;
        mesh->ddivpdp_n[c0]        = 0.0;
        mesh->drhodp_n[c0]         = 0.0;

        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {
            
            //----------------------------------------------------------//
            if ( model->elastic==1 ) etae      = model->dt*mesh->mu_n[c0];
            else                       etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->anisotropy == 0 ) {
                ani = 0.0; d1   = 0.0; d2   = 0.0; angle = 0.0; lxlz = 0.0; lx2 = 0.0;
            }
            else {
                // Anisotropy
                // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
                // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_n[c0];
                d1      = mesh->d1_n[c0];
                d2      = mesh->d2_n[c0];
                angle   = mesh->angle_n[c0];
                lxlz    = 0.5*d1;
                lx2     = pow( cos(angle), 2);
            }
            //----------------------------------------------------------//
//            Exx = mesh->exxd[c0]  + mesh->sxxd0[c0] /etae/2.0;
//            Ezz = mesh->ezzd[c0]  + mesh->szzd0[c0] /etae/2.0;
//            Exz = mesh->exz_n[c0] + mesh->sxz0_n[c0]/etae/2.0;
//            gxz = 2.0*Exz;
            
            Da11  = 2.0 - 2.0*ani*d1;
            Da12  = 2.0*ani*d1;
            Da13  = 2.0*ani*d2;
            Da22  = 2.0 - 2.0*ani*d1;
            Da23  =-2.0*ani*d2;
            Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;
            double Exx_ref   = mesh->exxd[c0]      + (iDa11*mesh->sxxd0[c0] + iDa12*mesh->szzd0[c0] + iDa13*mesh->sxz0_n[c0])/etae;
            double Ezz_ref   = mesh->ezzd[c0]      + (iDa12*mesh->sxxd0[c0] + iDa22*mesh->szzd0[c0] + iDa23*mesh->sxz0_n[c0])/etae;
            double Exz_ref   = mesh->exz_n[c0]     + (iDa13*mesh->sxxd0[c0] + iDa23*mesh->szzd0[c0] + iDa33*mesh->sxz0_n[c0])/2.0/etae;
            double P_ref     = mesh->p_in[c0];

            // Compute derivatives on the fly
            double detadexx =0., detadezz =0., detadgxz =0., detadp =0.;
            double ddivpdexx=0., ddivpdezz=0., ddivpdgxz=0., ddivpdp=0.;
            double danidexx =0., danidezz =0., danidgxz =0., danidp =0.;
            double drhodp=0.;
            DerivativesOnTheFly_n( &detadexx, &detadezz, &detadgxz, &detadp, &ddivpdexx, &ddivpdezz, &ddivpdgxz, &ddivpdp, &danidexx, &danidezz, &danidgxz, &danidp, &drhodp, c0, Exx_ref, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, mesh, materials, model, scaling );
            mesh->detadexx_n[c0]  = detadexx;
            mesh->detadezz_n[c0]  = detadezz;
            mesh->detadgxz_n[c0]  = detadgxz;
            mesh->detadp_n[c0]    = detadp;
            mesh->ddivpdexx_n[c0] = ddivpdexx;
            mesh->ddivpdezz_n[c0] = ddivpdezz;
            mesh->ddivpdgxz_n[c0] = ddivpdgxz;
            mesh->ddivpdp_n[c0]   = ddivpdp;
            mesh->drhodp_n[c0]    = drhodp;
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
        }
    }

    #pragma omp parallel for shared( mesh ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, el, etae, ani, d1, d2, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2 ) firstprivate( unsplit_diff_reac, materials, scaling, average, model, Nx, Nz, tol )
    for ( k1=0; k1<Nx*Nz; k1++ ) {

        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;

        // Initialise arrays to 0
        mesh->detadexx_s[c1]       = 0.0;
        mesh->detadezz_s[c1]       = 0.0;
        mesh->detadgxz_s[c1]       = 0.0;
        mesh->detadp_s[c1]         = 0.0;

        if ( mesh->BCg.type[c1] != 30 ) {

            //----------------------------------------------------------//
            if ( model->elastic==1   ) etae      = model->dt*mesh->mu_s[c1];
            else           etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->anisotropy == 0 ) {
                ani = 0.0; d1   = 0.0; d2   = 0.0; angle = 0.0; lxlz = 0.0; lx2 = 0.0;
            }
            else {
                // Anisotropy
                // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
                // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_s[c1];
                d1      = mesh->d1_s[c1];
                d2      = mesh->d2_s[c1];
                angle   = mesh->angle_s[c1];
                lxlz    = 0.5*d1;
                lx2     = pow( cos(angle), 2);
            }
            //----------------------------------------------------------//

            Da11  = 2.0 - 2.0*ani*d1;
            Da12  = 2.0*ani*d1;
            Da13  = 2.0*ani*d2;
            Da22  = 2.0 - 2.0*ani*d1;
            Da23  =-2.0*ani*d2;
            Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;

            double Exx_ref = mesh->exxd_s[c1]  + (iDa11*mesh->sxxd0_s[c1] + iDa12*mesh->szzd0_s[c1] + iDa13*mesh->sxz0[c1])/etae;
            double Ezz_ref = mesh->ezzd_s[c1]  + (iDa12*mesh->sxxd0_s[c1] + iDa22*mesh->szzd0_s[c1] + iDa23*mesh->sxz0[c1])/etae;
            double Exz_ref = mesh->exz[c1]     + (iDa13*mesh->sxxd0_s[c1] + iDa23*mesh->szzd0_s[c1] + iDa33*mesh->sxz0[c1])/2.0/etae;
            double P_ref     = mesh->P_s[c1];
            double gxz_ref   = 2.0*Exz_ref;
            
            // Compute derivatives on the fly
            double detadexx=0., detadezz=0., detadgxz=0., detadp=0.;
            double danidexx =0., danidezz =0., danidgxz =0., danidp =0.;
            DerivativesOnTheFly_s( &detadexx, &detadezz, &detadgxz, &detadp, &danidexx, &danidezz, &danidgxz, &danidp, c1, Exx_ref, Ezz_ref, Exz_ref, P_ref, ani, d1, d2, angle, lx2, lxlz, mesh, materials, model, scaling );
            mesh->detadexx_s[c1] = detadexx;
            mesh->detadezz_s[c1] = detadezz;
            mesh->detadgxz_s[c1] = detadgxz;
            mesh->detadp_s[c1]   = detadp;
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
        }

    }

}
#endif