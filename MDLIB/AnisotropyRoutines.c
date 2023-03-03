// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2018  MDOODZ Developper team
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
#include "mdoodz-private.h"
#include "include/mdoodz.h"
#include "RheologyDensity.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

// // Second invariant
// void Invii( Tensor2D *T ) {
//     T->ii = sqrt( 0.5*( pow(T->xx,2) + pow(T->zz,2) + pow(T->yy,2)) + pow(T->xz,2) );
// }

// // Squared second invariant
// void Invii2( Tensor2D *T ) {
//     T->ii2 =  0.5*( pow(T->xx,2) + pow(T->zz,2) + pow(T->yy,2)) + pow(T->xz,2);
// }

// Squared second invariant of deviatoric stress
double Y2( Tensor2D *T, double ani_fac ) {
    return 0.5*( pow(T->xx,2) + pow(T->zz,2) + pow(T->yy,2)) + ani_fac*pow(T->xz,2);
}

// Squared second invariant of deviatoric strain rate
double I2( Tensor2D *E, double ani_fac ) {
    return 0.5*( pow(E->xx,2) + pow(E->zz,2) + pow(E->yy,2)) + pow(E->xz,2)/ani_fac;
}

// Solution of 2x2 system: TODO: move to MiscFunction.c, make x and f vectors and A matrix
void Solve2x2(double *x1, double *x2, double f1, double f2, double a11, double a12, double a21, double a22) {
  const double det = a11 * a22 - a12 * a21; // determinant
  const double ai  =  1.0 / det * a22;      // inverse matrix components
  const double bi  = -1.0 / det * a12;
  const double ci  = -1.0 / det * a21;
  const double di  =  1.0 / det * a11;
  *x1 = ai * f1 + bi * f2;                  // inverse times rhs
  *x2 = ci * f1 + di * f2;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ViscosityConciseAniso( int phase, double lxlz, double lx2, double angle, double ani_fstrain, double ani_fac_e, double G, double T, double P, double d0, double phi, double X0, double Exx, double Ezz, double Exz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *Txx, double *Tzz, double *Txz, double* eta_vep, double* ani_vep, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp, double* Eii_lin, double* Eii_gbs, double* Eii_cst, double *d, double strain_acc, double dil, double fric, double C, double P0, double T0,  double *X1, double *OverS, double *Pcorr, double *rho, double beta, double div, double *div_el, double *div_pl, double *div_r, double *Wtot, double *Wel, double *Wdiss, int post_process, int centroid ) {
    // !!!!!!!!!!!!!!!!!!!!!!!!
    // ACTUNG FOR GSE:: d is now d0 and d1 is now d
    // !!!!!!!!!!!!!!!!!!!!!!!!
    // General paramaters
    const double tol    = 1.0e-11, R = materials->R, dt = model->dt, minEta = model->mineta, maxEta = model->maxeta;
    const int    nitmax = 20, noisy = 0;
    double eta = 0.0, eta_el, eta_cst;
    int    plastic = 0, constant = 0, dislocation = 0, peierls = 0, diffusion = 0, gbs = 0, elastic = model->iselastic, kinetics = 0, is_pl = 0;
    double eta_pwl = 0.0, B_pwl = 0.0, C_pwl = 0.0, Q_pwl = materials->Qpwl[phase], V_pwl = materials->Vpwl[phase], n_pwl = materials->npwl[phase], m_pwl = materials->mpwl[phase], r_pwl = materials->rpwl[phase], A_pwl = materials->Apwl[phase], f_pwl = materials->fpwl[phase], a_pwl = materials->apwl[phase], F_pwl = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase], t_pwl = materials->tpwl[phase];
    const double E_lin = materials->Qlin[phase], V_lin = materials->Vlin[phase], n_lin = materials->nlin[phase], m_lin = materials->mlin[phase], r_lin = materials->rlin[phase], A_lin = materials->Alin[phase], f_lin = materials->flin[phase], a_lin = materials->alin[phase], F_lin = materials->Flin[phase];
    double eta_ve_n, eta_ve_s, ani_fac_v, ani_fac_p, a_e, a_v, a_p, Tii;
    Tensor2D E_rot, T_rot, T0_rot;
    // Clean up
    *eta_vep   = 0.0, *ani_vep = 0.0;
    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0, *Eii_lin = 0.0, *Eii_pl = 0.0, *Eii_pwl = 0.0, *Eii_el = 0.0, *Eii_gbs = 0, *Eii_cst = 0.0;
    *Txx     = 0.0, *Tzz     = 0.0, *Txz    = 0.0;
    *Wtot    = 0.0, *Wdiss   = 0.0, *Wel    = 0.0, *d       = 0.0;
    *div_el  = 0.0, *div_pl  = 0.0, *div_r  = 0.0;
    *X1      = 0.0, *rho     = 0.0, *OverS  = 0.0, *Pcorr   = 0.0;
    // P corr will be corrected if plasticity feedbacks on pressure (dilation)
    *Pcorr = P;
    // Constant grain size
    *d = d0;
    // Activate deformation mechanisms
    if ( materials->cstv[phase] !=0                  ) constant    = 1;
    if ( materials->pwlv[phase] !=0                  ) dislocation = 1;
     // Turn of elasticity for the initialisation step  (viscous flow stress)
    if ( model->step    == 0                         ) elastic     = 0;
    // Visco-plastic incompressible limit
    if ( elastic==0                 ) { G = 1e1; dil = 0.0;}; //K = 1e1;
    // Zero C limit
    if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;
     // Precomputations
    if ( dislocation == 1 ) {
      B_pwl = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Q_pwl + P*V_pwl)/R/n_pwl/T ) * pow(d0, m_pwl/n_pwl) * pow(f_pwl, -r_pwl/n_pwl) * exp(-a_pwl*phi/n_pwl);
      C_pwl   = pow(2.0*B_pwl, -n_pwl);
    }
    // Isolated anisotropy factors
    if ( materials->ani_fstrain[phase] == 0 ) {
        ani_fac_v = materials->ani_fac_v[phase];
        ani_fac_p = materials->ani_fac_p[phase];
    }
    else {
        ani_fac_v = ani_fstrain;
        ani_fac_p = ani_fstrain;
    }
    a_e = sqrt(ani_fac_e), a_v = sqrt(ani_fac_v), a_p = sqrt(ani_fac_p);
    // Transform strain rate: Q*E*Qt
    const double nz2 = 1.0-lx2;
    E_rot.xx  =   lx2*Exx +  nz2*Ezz +   2.*lxlz*Exz;
    E_rot.zz  =   nz2*Exx +  lx2*Ezz -   2.*lxlz*Exz;
    E_rot.xz  = -lxlz*Exx + lxlz*Ezz + (lx2-nz2)*Exz;
    E_rot.yy  = -E_rot.xx - E_rot.zz;                 // !!! out-of-plane component !!!
    // Local iterations to determine eta_vep and ani_vep
    double f0_n, f0_s, f_n, f_s, dfnden, dfndes, dfsden, dfsdes;
    double I2_v, Y2_v, Eii_vis, W_n, W_s, ieta_pwl, deta_n, deta_s;
    // Construct initial guess for viscosity: Isolated viscosities 
    I2_v     = I2( &E_rot, ani_fac_v ); 
    if (constant)    eta_cst  = materials->eta0[phase];
    if (dislocation) eta_pwl  = B_pwl * pow( I2_v, 0.5*(1-n_pwl)/n_pwl );
    if (elastic)     eta_el   = G*dt; 

    // Cases where no local iterations are needed
    // if (dislocation) eta_ve_n = eta_pwl;
    // if (dislocation) eta_ve_s = eta_pwl/ani_fac_v;
    // if (constant) eta_ve_n = eta_cst;
    // if (constant) eta_ve_s = eta_cst/ani_fac_v;
    // T_rot.xx = 2.0 * eta_ve_n * E_rot.xx;
    // T_rot.zz = 2.0 * eta_ve_n * E_rot.zz;
    // T_rot.xz = 2.0 * eta_ve_s * E_rot.xz;
    // T_rot.yy = -T_rot.xx - T_rot.xx; 

    // Define viscosity bounds
    double eta_up_n  = 1.0e100 / scaling->eta, eta_up_s  = 1.0e100 / scaling->eta;
    if (constant)    eta_up_n = MINV(eta_up_n, eta_cst);
    if (dislocation) eta_up_n = MINV(eta_up_n, eta_pwl);
    if (elastic)     eta_up_n = MINV(eta_up_n, eta_el);
    //-----------
    if (constant)    eta_up_s = MINV(eta_up_s, eta_cst/ani_fac_v);
    if (dislocation) eta_up_s = MINV(eta_up_s, eta_pwl/ani_fac_v);
    if (elastic)     eta_up_s = MINV(eta_up_s, eta_el /ani_fac_e );
    // Initial guess
    eta_ve_n = eta_up_n;
    eta_ve_s = eta_up_s;
    // Iterations for visco-elastic trial
    if (noisy) printf("Start local VE iterations for phase %d (cst: %d --- el: %d --- pwl: %d)\n", phase, constant, elastic, dislocation);
    // if (noisy) printf("ani_fac_v = %1.1f  ani_fac_p = %1.1f\n", ani_fac_v, ani_fac_p );    // DELETE
    // if (noisy) printf("e_n=%2.2e e_s=%2.2e eta_pwl=%2.2e\n", eta_ve_n, eta_ve_s, eta_pwl); // DELETE
    for (int it = 0; it < nitmax; it++) {
      // Compute stress in principal plane
      T_rot.xx = 2.0 * eta_ve_n * E_rot.xx;
      T_rot.zz = 2.0 * eta_ve_n * E_rot.zz;
      T_rot.xz = 2.0 * eta_ve_s * E_rot.xz;
      T_rot.yy = -T_rot.xx - T_rot.zz; 
      Y2_v     = Y2( &T_rot, ani_fac_v);
      // ---------- Strain rates
      if (constant)    *Eii_cst = sqrt(Y2_v)  / 2.0 / eta_cst;
      if (dislocation) *Eii_pwl = C_pwl * pow(Y2_v, 0.5*n_pwl);
      Eii_vis = *Eii_pwl + *Eii_cst;
      // ---------- Viscosities
      // eta_pwl = pow( B_pwl, -1.0/n_pwl ) * pow( *Eii_pwl, (1-n_pwl)/n_pwl );
      eta_pwl = B_pwl * pow( *Eii_pwl, (1-n_pwl)/n_pwl );
      // ---------- Residuals
      f_n = 1.0; f_s = 1.0;
      if (constant)    { f_n -= eta_ve_n/eta_cst; f_s -= eta_ve_s/(eta_cst/ani_fac_v); }
      if (elastic)     { f_n -= eta_ve_n/eta_el;  f_s -= eta_ve_s/(eta_el /ani_fac_e); }
      if (dislocation) { f_n -= eta_ve_n/eta_pwl; f_s -= eta_ve_s/(eta_pwl/ani_fac_v); }
      W_n       = 2.0*(E_rot.xx*T_rot.xx + E_rot.zz*T_rot.zz + E_rot.yy*T_rot.yy);
      W_s       = 2.0*E_rot.xz*T_rot.xz;
      ieta_pwl  = 1.0 / eta_pwl;
      // Make some noise!!!!!
      if (noisy) printf("It. %02d: f_n = %2.2e --- f_s = %2.2e\n", it, f_n, f_s);
      // Exit criteria
      if ( fabs(f_n) < tol / 100 && fabs(f_s) < tol / 100 ) {
        if (it > 10) printf("V-E L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
        break;
      } else if (it == nitmax - 1 && (fabs(f_n) > tol || fabs(f_s) > tol) ) {
        printf("Visco-Elastic iterations failed!\n");
        exit(0);
      }
      // Compute Jacobian
      dfnden = 0.0; dfndes = 0.0; dfsden = 0.0; dfsdes = 0.0;
      if (constant) {
        dfnden -= 1.0       / eta_cst;
        dfsdes -= ani_fac_v / eta_cst;
      }
      if (elastic) {
        dfnden -= 1.0       / eta_el;
        dfsdes -= ani_fac_e / eta_el;
      }
      if (dislocation) {
        dfnden += W_n *                     eta_ve_n * ieta_pwl * (1.0 - n_pwl) / (2.0 * Y2_v)  - ieta_pwl;
        dfndes += W_s * ani_fac_v         * eta_ve_n * ieta_pwl * (1.0 - n_pwl) / Y2_v;
        dfsden += W_n * ani_fac_v         * eta_ve_s * ieta_pwl * (1.0 - n_pwl) / (2.0 * Y2_v);
        dfsdes += W_s * pow(ani_fac_v, 2) * eta_ve_s * ieta_pwl * (1.0 - n_pwl) / Y2_v - ani_fac_v * ieta_pwl;
      }
      // Inverse of the Jacobian
      Solve2x2(&deta_n, &deta_s, f_n, f_s, dfnden, dfndes, dfsden, dfsdes);
      // Coupled update of normal and shear viscosities
      eta_ve_n -= deta_n;
      eta_ve_s -= deta_s;
    }
    Tii      = sqrt(Y2_v);
    *eta_vep = eta_ve_n;
    *ani_vep = eta_ve_n/eta_ve_s;

    // Iterations for visco-elasto-viscoplastic correction
    double eta_vp0 = materials->eta_vp[phase], n_vp = materials->n_vp[phase], eta_vp = materials->eta_vp[phase];
    double a1=materials->axx[phase], a2=materials->azz[phase], a3=materials->ayy[phase];
    double Coh = C, K = 1.0 / beta;
    double Y2_p = Y2( &T_rot, ani_fac_p);
    double Ft   = sqrt(Y2_p) - Coh - P*sin(fric)/3*(a1+a2+a3) +  sin(fric)/3.0*( a1*T_rot.xx + a2*T_rot.zz + a3*T_rot.yy);
    double Fc = Ft;
    double ap_n, ap_s, eta_ve=eta_ve_n, gdot=0.0, a_ve=sqrt(eta_ve_n/eta_ve_s), eta_pl, divp;
    double Y1_p_c, Y2_p_c, dFdgdot;
    double txx = T_rot.xx, tzz=T_rot.zz, tyy=T_rot.yy, txz=T_rot.xz;
    double Pc = P;
    if (Ft>1e-17) {
      is_pl = 1;
      if (noisy) printf("Start local VP iterations (cst: %d --- el: %d --- pwl: %d)\n", constant, elastic, dislocation);
      // printf("Start local VP iterations (cst: %d --- el: %d --- pwl: %d)\n", constant, elastic, dislocation);
      // printf("Plastic!! Ft = %2.4e --- Y2_v = %2.4e Y2_p = %2.4e\n", Ft, Y2_v, Y2_p);  
      for (int it = 0; it < nitmax; it++) {

        ap_n   = 1 - eta_ve*gdot/sqrt(Y2_p);
        ap_s   = 1 - pow(a_p, 2)*eta_ve*gdot/(sqrt(Y2_p)*pow(a_ve, 2));
        Y1_p_c = a1*ap_n*txx + a2*ap_n*tzz + a3*(-ap_n*txx - ap_n*tzz);
        Y2_p_c = pow(a_p, 2)*pow(ap_s, 2)*pow(txz, 2) + (1.0/2.0)*pow(ap_n, 2)*pow(txx, 2) + (1.0/2.0)*pow(ap_n, 2)*pow(tyy, 2) + (1.0/2.0)*pow(ap_n, 2)*pow(tzz, 2);
        Pc     = P + K*dt*gdot*sin(dil)/3.*(a1+a2+a3);
        Fc     = -Coh - 1.0/3.0*Pc*(a1 + a2 + a3)*sin(fric) + (1.0/3.0)*Y1_p_c*sin(fric) + sqrt(Y2_p_c) - eta_vp*gdot;

        // double exxp      = gdot*( (a1/3-a3/3)*sin(fric) + txx/sqrt(Y2_p)/2);
        // double ezzp      = gdot*( (a2/3-a3/3)*sin(fric) + tzz/sqrt(Y2_p)/2);
        // double exzp      = gdot*txz/sqrt(Y2_p)*a_p*a_p;
        // T_rot.xx = 2.0 * eta_ve_n * (E_rot.xx - exxp);
        // T_rot.zz = 2.0 * eta_ve_n * (E_rot.zz - ezzp);
        // T_rot.xz = 2.0 * eta_ve_s * (E_rot.xz - 0.5*exzp);
        // T_rot.yy = -T_rot.xx-T_rot.zz;
        // double Y2_p_c1   = Y2(&T_rot, ani_fac_p);
        // double txx1 = sqrt(pow(ap_n, 2)*pow(txx, 2));
        // double tzz1 = sqrt(pow(ap_n, 2)*pow(tzz, 2));
        // double txz1 = sqrt(pow(ap_s, 2)*pow(txz, 2));

        // Make some noise!!!!!
        if (noisy) printf("It. %02d: f_n = %2.2e\n", it, fabs(Fc/Ft));
        // printf("It. %02d: f_n = %2.2e\n", it, fabs(Fc/Ft));

        if ( fabs(Fc/Ft) < 1e-8 ) {
        if (it > 10) printf("VP L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
        // if (model->step>0)exit(0);
        break;
        } else if (it == nitmax - 1 && fabs(Fc/Ft) > 1e-8 ) {
          printf("Viscoplastic iterations failed! tol=%2.2e\n", tol);
          exit(0);
        }
        // Newton update
        dFdgdot = -1.0/9.0*K*dt*pow(a1 + a2 + a3, 2)*sin(dil)*sin(fric) - eta_vp + (1.0/3.0)*(a3*(eta_ve*txx/sqrt(Y2_p) + eta_ve*tzz/sqrt(Y2_p)) - a1*eta_ve*txx/sqrt(Y2_p) - a2*eta_ve*tzz/sqrt(Y2_p))*sin(fric) + (-pow(a_p, 4)*ap_s*eta_ve*pow(txz, 2)/(sqrt(Y2_p)*pow(a_ve, 2)) - 1.0/2.0*ap_n*eta_ve*pow(txx, 2)/sqrt(Y2_p) - 1.0/2.0*ap_n*eta_ve*pow(tyy, 2)/sqrt(Y2_p) - 1.0/2.0*ap_n*eta_ve*pow(tzz, 2)/sqrt(Y2_p))/sqrt(Y2_p_c);
        gdot    -= Fc/dFdgdot;
      }
  
      // Y2_p_c   = Y2(&T_rot, ani_fac_p);
      // Y1_p_c   = a1*T_rot.xx + a2*T_rot.zz + a3*T_rot.yy;
      // double Fchk0 = sqrt(Y2_p_c) - Coh - Pc*sin(fric)/3.0*(a1+a2+a3) +  sin(fric)/3.0*Y1_p_c - eta_vp*gdot;

      double axx   = (sqrt(Y2_p)             -            eta_ve*gdot)/(sqrt(Y2_p));
      double axz   = (sqrt(Y2_p)*pow(a_ve,2) - pow(a_p,2)*eta_ve*gdot)/(sqrt(Y2_p)*pow(a_ve,2));
      // T_rot.xx = 2.0 * eta_ve_n*axx * E_rot.xx;
      // T_rot.zz = 2.0 * eta_ve_n*axx * E_rot.zz;
      // T_rot.xz = 2.0 * eta_ve_s*axz * E_rot.xz;
      // T_rot.yy = -T_rot.xx -T_rot.zz;
      // Y2_p_c   = Y2(&T_rot, ani_fac_p);
      // Y1_p_c   = a1*T_rot.xx + a2*T_rot.zz + a3*T_rot.yy;
      // double Fchk1 = sqrt(Y2_p_c) - Coh - Pc*sin(fric)/3.0*(a1+a2+a3) +  sin(fric)/3.0*Y1_p_c - eta_vp*gdot;

      *eta_vep = eta_ve*axx;
      *ani_vep = a_ve*a_ve*(axx/axz);
      eta_ve_n = *eta_vep;
      eta_ve_s = *eta_vep/(*ani_vep);
      T_rot.xx = 2.0 * eta_ve_n * E_rot.xx;
      T_rot.zz = 2.0 * eta_ve_n * E_rot.zz;
      T_rot.xz = 2.0 * eta_ve_s * E_rot.xz;
      T_rot.yy = -T_rot.xx - T_rot.zz; 
      *div_pl  = gdot*sin(dil)/3.*(a1+a2+a3);
      *Eii_pl  = gdot/2.0;
      Tii      = sqrt(Y2_p);
      // Y2_p_c   = Y2(&T_rot, ani_fac_p);
      // Y1_p_c   = a1*T_rot.xx + a2*T_rot.zz + a3*T_rot.yy;
      // double Fchk2 = sqrt(Y2_p_c) - Coh - Pc*sin(fric)/3*(a1+a2+a3) +  sin(fric)/3.0*Y1_p_c - eta_vp*gdot;
      // printf("Ft = %2.6e --- Fc = %2.2e --- Fchk0 = %2.6e --- Fchk1 = %2.6e --- Fchk2 = %2.6e\n", Ft, Fc, Fchk0, Fchk1, Fchk2);
      // exit(1);
    }
    // Back-transform stress: Qt*E*Q
    *Txx =   lx2*T_rot.xx +  nz2*T_rot.zz -   2.*lxlz*T_rot.xz;
    *Tzz =   nz2*T_rot.xx +  lx2*T_rot.zz +   2.*lxlz*T_rot.xz;
    *Txz =  lxlz*T_rot.xx - lxlz*T_rot.zz + (lx2-nz2)*T_rot.xz;
    // Update effective viscosity and anisotropy factor
    *Pcorr   = Pc;
    eta      = *eta_vep;
    eta_pl   = Tii/gdot;
    divp     = gdot*sin(dil);
    if ( post_process == 1) {

      //-------- Post-Processing
  if ( post_process == 1) {

    // Strain rates: VEP partitioning
    if ( dislocation == 1 )              eta_pwl  = B_pwl * pow( *Eii_pwl, 1.0/n_pwl - 1.0 ); 

    // Compute dissipative strain rate components
    if (elastic) *Eii_el   = Tii*Tii/(4*eta_el*eta_el);
    if (is_pl)   *Eii_pl   = gdot/2.0;

    // Partitioning of volumetric strain
    *div_pl  = divp;
    *div_el  = - (Pc - P0) / (K*dt);
    // *div_r   = divr;

    double inv_eta_diss = 0.0;
    if (dislocation== 1)  inv_eta_diss += (1.0/eta_pwl);
    if (constant   == 1)  inv_eta_diss += (1.0/eta_cst);
    if (is_pl      == 1)  inv_eta_diss += (1.0/eta_pl ); 
    eta        = 1.0/(inv_eta_diss);

    // Viscoplastic overstress
    *OverS = eta_vp*gdot;

    // Work(s)
    if (centroid==1) {
      *Wtot  = Tii*Tii/eta_ve;
      *Wel   = Tii*Tii/eta_el;
      *Wdiss = Tii*Tii/eta;
      // printf("Wtot = %2.6e Wdiss = %2.6e Wel = %2.6e\n", *Wtot, *Wdiss, *Wel);
      // printf("Wtot - (Wdiss + Wel) = %2.6e\n", *Wtot - (*Wdiss + *Wel));
      // printf("eta = %2.6e --- eta_ve = %2.6e --- eta_ve1 = %2.6e\n", eta, eta_ve, 1.0/(1./eta+1./eta_el));
    }
  }
      // *div_el  = - (Pc - P0) / (K*dt);
      // double Tii = sqrt( *Txx* *Txx + *Tzz * *Tzz + (-*Txx - *Tzz)*(-*Txx - *Tzz));
      // // double Eii = sqrt( Exx*Exx+Ezz*Ezz+(-Exx-Ezz)*(-Exx-Ezz));
      // if (elastic==1) *Eii_el  = Tii*Tii/(4*eta_el*eta_el); //fabs(Eii-*Eii_pwl-*Eii_pl); // this is a big shortcut, ideally one should compute E^e = 1/(2G)*Dani_e*(T-T0)/dt
    }
    return eta;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGridAniso( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {

  // printf("NonNewtonianViscosityGridAniso!\n");

  int    p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1;
  double eta, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
  double el = 0.0, Wtot, Wel, Wdiss;
  int    average = model->eta_avg, unsplit_diff_reac = model->unsplit_diff_reac;
  double Xreac;
  double OverS;
  double Pcorr, rho;
  double Exx, Ezz, Exz, eta_e, ani, d1, d2, angle, lxlz, lx2, lay_ang;
  double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
  if (model->iselastic == 1) el = 1.0;

  double etaVE, VEcoeff, d0; // to delete

  Nx  = mesh->Nx;
  Ncx = Nx - 1;
  Nz  = mesh->Nz;  Ncz = Nz-1;

  // Stuff to be interpolated to vertices
  InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
  InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

  // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Wtot, Wel, Wdiss, Xreac,OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e, ani, d1, d2, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2, lay_ang ) firstprivate( el, unsplit_diff_reac, materials, scaling, average, model, Ncx, Ncz )
  for ( k1=0; k1<Ncx*Ncz; k1++ ) {

    //    for ( l=0; l<Ncz; l++ ) {
    //        for ( k=0; k<Ncx; k++ ) {

    // printf("average = %d\n", average);
    k      = mesh->kp[k1];
    l      = mesh->lp[k1];
    c0     = k  + l*(Ncx);

    mesh->eta_n[c0]       = 0.0;
    mesh->eta_phys_n[c0]  = 0.0;
    mesh->VE_n[c0]        = 0.0;
    mesh->sxxd[c0]        = 0.0;
    mesh->szzd[c0]        = 0.0;
    mesh->sxz_n[c0]       = 0.0;
    mesh->eII_el[c0]      = 0.0;
    mesh->eII_pl[c0]      = 0.0;
    mesh->eII_pwl[c0]     = 0.0;
    mesh->eII_exp[c0]     = 0.0;
    mesh->eII_lin[c0]     = 0.0;
    mesh->eII_gbs[c0]     = 0.0;
    mesh->eII_cst[c0]     = 0.0;
    mesh->d_n[c0]         = 0.0;
    mesh->p_corr[c0]      = 0.0;                // Achtung baby
    mesh->div_u_el[c0]    = 0.0;
    mesh->div_u_pl[c0]    = 0.0;
    mesh->div_u_r[c0]     = 0.0;
    mesh->Wtot[c0]        = 0.0;
    mesh->Wel[c0]         = 0.0;
    mesh->Wdiss[c0]       = 0.0;
    //        X                     =  mesh->Xreac_n[c0]; // Save X first
    //        if (model->progress_transform==1) mesh->Xreac_n[c0]    = 0.0;
    if ( unsplit_diff_reac == 0 ) mesh->X_n[c0]        = 0.0;
    mesh->OverS_n[c0]    = 0.0;
    mesh->aniso_factor_n[c0]   = 0.0;
    // mesh->aniso_factor_e_n[c0] = 0.0;

    if ( model->density_change == 1 ) {
      mesh->rho_n[c0]  = 0.0;
    }

    // Loop on grid nodes
    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {

      //----------------------------------------------------------//
      if ( model->iselastic==1 ) eta_e      = model->dt*mesh->mu_n[c0];
      else                       eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
        // Anisotropy
        // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
        // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_n[c0];
        d1      = mesh->d1_n[c0];  // d1 = 2*lx^2*lz^2 
        d2      = mesh->d2_n[c0];  // d2 = lx*lz*(pow(-lx, 2.0) + pow(lz, 2.0));
        angle   = mesh->angle_n[c0];
        lay_ang = angle - M_PI/2.0;
        lxlz    = cos(lay_ang)*sin(lay_ang);
        lx2     = cos(lay_ang)*cos(lay_ang);
        // lxlz    = 0.5*d1;
        // lx2     = pow( cos(angle), 2);
      //----------------------------------------------------------//
      // if ( model->iselastic==1 ) {
      //   Exx = mesh->exxd[c0]  + mesh->sxxd0[c0] /eta_e/2.0;
      //   Ezz = mesh->ezzd[c0]  + mesh->szzd0[c0] /eta_e/2.0;
      //   Exz = mesh->exz_n[c0] + mesh->sxz0_n[c0]/eta_e/2.0;
      // }
      // else {
      //   Exx = mesh->exxd[c0] ;
      //   Ezz = mesh->ezzd[c0] ;
      //   Exz = mesh->exz_n[c0];
      // }

      double aniS_e, aniS_vep;
      aniS_e = 1.0 - 1.0 / mesh->aniso_factor_e_n[c0];
      EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], d1, d2, aniS_e, eta_e, model->iselastic ); 
      
      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_n[p][c0])>min_fraction ) is_phase_active = true;

        if ( is_phase_active ) {
          eta =  ViscosityConciseAniso( p, lxlz, lx2, angle, mesh->FS_AR_n[c0], mesh->aniso_factor_e_n[c0], mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], Exx, Ezz, Exz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &ani_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], mesh->p0_n[c0], mesh->T0_n[c0], &Xreac, &OverS, &Pcorr, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 1 );
          mesh->phase_eta_n[p][c0] = eta_vep;

          switch ( average ) {
            case 0 :
              // ARITHMETIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * eta_vep;
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * eta;
              mesh->aniso_factor_n[c0]   += mesh->phase_perc_n[p][c0] * ani_vep;
              // mesh->aniso_factor_e_n[c0] += mesh->phase_perc_n[p][c0] * ani_e;
              break;
            case 1 :
              // HARMONIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * 1.0/eta_vep;
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * 1.0/eta;
              mesh->aniso_factor_n[c0]   += mesh->phase_perc_n[p][c0] * 1.0/ani_vep;
              // mesh->aniso_factor_e_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/ani_e;
              break;
            case 2 :
              // GEOMETRIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * log(eta_vep);
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * log(eta);
              mesh->aniso_factor_n[c0]   += mesh->phase_perc_n[p][c0] * log(ani_vep);
              // mesh->aniso_factor_e_n[c0] += mesh->phase_perc_n[p][c0] * log(ani_e);
              break;
          }

        //   mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
          mesh->eII_el[c0]     += mesh->phase_perc_n[p][c0] * eII_el;
          mesh->eII_pl[c0]     += mesh->phase_perc_n[p][c0] * eII_pl;
          mesh->eII_pwl[c0]    += mesh->phase_perc_n[p][c0] * eII_pwl;
          mesh->eII_exp[c0]    += mesh->phase_perc_n[p][c0] * eII_exp;
          mesh->eII_lin[c0]    += mesh->phase_perc_n[p][c0] * eII_lin;
          mesh->eII_gbs[c0]    += mesh->phase_perc_n[p][c0] * eII_gbs;
          mesh->eII_cst[c0]    += mesh->phase_perc_n[p][c0] * eII_cst;
          mesh->d_n[c0]        += mesh->phase_perc_n[p][c0] * 1.0/dnew;

          mesh->Wtot[c0]       += mesh->phase_perc_n[p][c0] * Wtot;
          mesh->Wdiss[c0]      += mesh->phase_perc_n[p][c0] * Wdiss;
          mesh->Wel[c0]        += mesh->phase_perc_n[p][c0] * Wel;

          if (mesh->Wdiss[c0]<0.0) {printf("negative dissipation: you crazy! --> Wdiss = %2.2e\n", mesh->Wdiss[c0]*scaling->S*scaling->E); }

          mesh->p_corr[c0]      += mesh->phase_perc_n[p][c0] * Pcorr;
          mesh->div_u_el[c0]    += mesh->phase_perc_n[p][c0] * div_el;
          mesh->div_u_pl[c0]    += mesh->phase_perc_n[p][c0] * div_pl;
          mesh->div_u_r[c0]     += mesh->phase_perc_n[p][c0] * div_r;

          if ( unsplit_diff_reac == 0 ) mesh->X_n[c0]         += mesh->phase_perc_n[p][c0] * Xreac;
          mesh->OverS_n[c0]     += mesh->phase_perc_n[p][c0] * OverS;

          // Volume changes
          if ( model->density_change == 1 ) {
            mesh->rho_n[c0]     += mesh->phase_perc_n[p][c0] * (rho);
          }
        }
      }

      mesh->d_n[c0]          = 1.0/mesh->d_n[c0];

      // HARMONIC AVERAGE
      if ( average == 1 ) {
        mesh->eta_n[c0]            = 1.0/mesh->eta_n[c0];
        mesh->eta_phys_n[c0]       = 1.0/mesh->eta_phys_n[c0];
        mesh->aniso_factor_n[c0]   = 1.0/mesh->aniso_factor_n[c0];
        // mesh->aniso_factor_e_n[c0] = 1.0/mesh->aniso_factor_e_n[c0];

        if (isinf (mesh->eta_phys_n[c0]) ) {
          printf("Inf: Problem on cell centers:\n");
          for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
          printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", eta, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
          printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
          printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
          exit(1);
        }
        if (isnan (mesh->eta_phys_n[c0]) ) {
          printf("NaN: Problem on cell centers:\n");
          printf("progress_transform %d\n", model->progress_transform);
          for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
          printf("eta=%2.2e G=%2.2e T=%2.2e P=%2.2e d=%2.2e phi=%2.2e %2.2e %2.2e %2.2e %2.2e\n", eta*scaling->eta, mesh->mu_n[c0]*scaling->S, mesh->T[c0]*scaling->T, mesh->p_in[c0]*scaling->S, mesh->d0_n[c0]*scaling->L, mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
          printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
          printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
          exit(1);
        }
      }
      // GEOMETRIC AVERAGE
      if ( average == 2 ) {
        mesh->eta_n[c0]            = exp(mesh->eta_n[c0]);
        mesh->eta_phys_n[c0]       = exp(mesh->eta_phys_n[c0]);
        mesh->aniso_factor_n[c0]   = exp(mesh->aniso_factor_n[c0]);
        // mesh->aniso_factor_e_n[c0] = exp(mesh->aniso_factor_e_n[c0]);
      }

      // Effective strain rate
      // double aniS_e, aniS_vep;
      // if ( model->aniso_fstrain  == 0 ) aniS_e = 1.0 - 1.0 / mesh->aniso_factor_e_n[c0];
      // if ( model->aniso_fstrain  == 1 ) aniS_e = 1.0 - 1.0 / mesh->FS_AR_n[c0];
      // EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], d1, d2, aniS_e, eta_e, model->iselastic ); 
      // Final stress update
      aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
      // if ( model->aniso_fstrain  == 0 ) aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
      // if ( model->aniso_fstrain  == 1 ) aniS_vep = 1.0 - 1.0 / mesh->FS_AR_n[c0];
      // Da11  = 2.0 - 2.0*aniS_vep*d1;
      // Da12  = 2.0*aniS_vep*d1;
      // Da13  =-2.0*aniS_vep*d2;
      // Da22  = 2.0 - 2.0*aniS_vep*d1;
      // Da23  = 2.0*aniS_vep*d2;
      // Da33  = 1.0  + 2.0*aniS_vep*(d1 - 0.5);
      // Normal stress
      const double Dxx = Exx*(1.0 - aniS_vep*d1) + Ezz*aniS_vep*d1 - 2.0*Exz*aniS_vep*d2;
      const double Dzz = Ezz*(1.0 - aniS_vep*d1) + Exx*aniS_vep*d1 + 2.0*Exz*aniS_vep*d2;
      mesh->sxxd[c0] = 2.0 * mesh->eta_n[c0] * Dxx;
      mesh->szzd[c0] = 2.0 * mesh->eta_n[c0] * Dzz;
    }
  }

// Calculate vertices viscosity
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Wtot, Wel, Wdiss, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e, ani, d1, d2, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2, lay_ang ) firstprivate( el, unsplit_diff_reac, materials, scaling, average, model, Nx, Nz )
  for ( k1=0; k1<Nx*Nz; k1++ ) {

    k  = mesh->kn[k1];
    l  = mesh->ln[k1];
    c1 = k + l*Nx;

    mesh->sxz[c1]              = 0.0;
    mesh->eta_phys_s[c1]       = 0.0;
    mesh->eta_s[c1]            = 0.0;
    mesh->aniso_factor_s[c1]   = 0.0;
    // mesh->aniso_factor_e_s[c1] = 0.0;

    if (unsplit_diff_reac == 0) mesh->X_s[c1]        = 0.0;

    if ( mesh->BCg.type[c1] != 30 ) {

      if ( model->iselastic==1   ) eta_e      = model->dt*mesh->mu_s[c1];
      else           eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
        // // Anisotropy
        // if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
        // if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_s[c1];
        d1      = mesh->d1_s[c1];
        d2      = mesh->d2_s[c1];
        angle   = mesh->angle_s[c1];
        lay_ang = angle - M_PI/2.0;
        lxlz    = cos(lay_ang)*sin(lay_ang);
        lx2     = cos(lay_ang)*cos(lay_ang);
        // lxlz    = 0.5*d1;
        // lx2     = pow( cos(angle), 2);
      //----------------------------------------------------------//
      // Effective strain rate
      double aniS_e, aniS_vep;
      aniS_e = 1.0 - 1.0 / mesh->aniso_factor_e_s[c1];
      // if ( model->aniso_fstrain  == 0 ) aniS_e = 1.0 - 1.0 / mesh->aniso_factor_e_s[c1];
      // if ( model->aniso_fstrain  == 1 ) aniS_e = 1.0 - 1.0 / mesh->FS_AR_s[c1];
      EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], d1, d2, aniS_e, eta_e, model->iselastic ); 

      // if ( model->iselastic==1   ) {
      //   Exx = mesh->exxd_s[c1] + mesh->sxxd0_s[c1]/eta_e/2.0;
      //   Ezz = mesh->ezzd_s[c1] + mesh->szzd0_s[c1]/eta_e/2.0;
      //   Exz = mesh->exz[c1]    + mesh->sxz0[c1]   /eta_e/2.0;
      // }
      // else {
      //   Exx = mesh->exxd_s[c1];
      //   Ezz = mesh->ezzd_s[c1];
      //   Exz = mesh->exz[c1]   ;
      // }
      
      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_s[p][c1])>min_fraction ) is_phase_active = true;

        if ( is_phase_active ) {

          eta =  ViscosityConciseAniso( p, lxlz, lx2, angle, mesh->FS_AR_s[c1], mesh->aniso_factor_e_s[c1], mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], Exx, Ezz, Exz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &ani_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], mesh->p0_s[c1], 0.0, &Xreac, &OverS, &Pcorr, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 0 );
          mesh->phase_eta_s[p][c1] = eta_vep;

          switch ( average ) {
            case 0 :
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * eta_vep;
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * eta;
              mesh->aniso_factor_s[c1]   += mesh->phase_perc_s[p][c1] * ani_vep;
              // mesh->aniso_factor_e_s[c1] += mesh->phase_perc_s[p][c1] * ani_e;
              break;
            case 1:
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * 1.0/eta_vep;
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * 1.0/eta;
              mesh->aniso_factor_s[c1]   += mesh->phase_perc_s[p][c1] * 1.0/ani_vep;
              // mesh->aniso_factor_e_s[c1] += mesh->phase_perc_s[p][c1] * 1.0/ani_e;
              break;
            case 2:
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * log(eta_vep);
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * log(eta);
              mesh->aniso_factor_s[c1]   += mesh->phase_perc_s[p][c1] * log(ani_vep);
              // mesh->aniso_factor_e_s[c1] += mesh->phase_perc_s[p][c1] * log(ani_e);
              break;
          }
          if (unsplit_diff_reac == 0) mesh->X_s[c1]        += mesh->phase_perc_s[p][c1] * Xreac;
        }
      }
      // HARMONIC AVERAGE
      if (average == 1) {
        mesh->eta_s[c1]            = 1.0/mesh->eta_s[c1];
        mesh->eta_phys_s[c1]       = 1.0/mesh->eta_phys_s[c1];
        mesh->aniso_factor_s[c1]   = 1.0/mesh->aniso_factor_s[c1];
        // mesh->aniso_factor_e_s[c1] = 1.0/mesh->aniso_factor_e_s[c1];
        if (isinf (mesh->eta_phys_s[c1]) ) {
          printf("Inf: Problem on cell vertices:\n");
          for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
          printf("%2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
          printf("x=%2.2e z=%2.2e\n", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
          exit(1);
        }
        if (isnan (mesh->eta_phys_s[c1]) ) {
          printf("Nan: Problem on cell vertices:\n");
          for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
          printf("%2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1],  mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
          printf("x=%2.2e z=%2.2e\n", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
          exit(1);
        }
      }
      // GEOMETRIC AVERAGE
      if (average == 2) {
        mesh->eta_s[c1]            = exp(mesh->eta_s[c1]);
        mesh->eta_phys_s[c1]       = exp(mesh->eta_phys_s[c1]);
        mesh->aniso_factor_s[c1]   = exp(mesh->aniso_factor_s[c1]);
        // mesh->aniso_factor_e_s[c1] = exp(mesh->aniso_factor_e_s[c1]);
      }

      // // Effective strain rate
      // double aniS_e, aniS_vep;
      // if ( model->aniso_fstrain  == 0 ) aniS_e = 1.0 - 1.0 / mesh->aniso_factor_e_s[c1];
      // if ( model->aniso_fstrain  == 1 ) aniS_e = 1.0 - 1.0 / mesh->FS_AR_s[c1];
      // EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], d1, d2, aniS_e, eta_e, model->iselastic ); 

      // Final stress update
      aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
      // if ( model->aniso_fstrain  == 0 ) aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
      // if ( model->aniso_fstrain  == 1 ) aniS_vep = 1.0 - 1.0 / mesh->FS_AR_s[c1];
      const double Dxz  = -Exx*aniS_vep*d2 + Ezz*aniS_vep*d2 + 2*Exz*(aniS_vep*(d1 - 0.5) + 0.5); 
      // Shear stress
      mesh->sxz[c1]  = 2.0 * mesh->eta_s[c1] * Dxz;

      // if (l==0) printf("Exx = %2.2f, Ezz = %2.2f, Exz = %2.2f\n", Exx, Ezz, Exz);
        


    }
  }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseDirectorVector (grid* mesh, markers* particles, params* model, mat_prop* materials, double scaling_L) {

  int    cent=1, vert=0, prop=1, interp=0;
  int k;
  double angle, norm;
  const double Earth_radius = 6370e3/scaling_L;

#pragma omp parallel for shared( particles ) private( angle, norm )
  for (k=0; k<particles->Nb_part; k++) {

    if ( particles->phase[k] != -1 ) {

      // Set up director vector
      if (model->particle_aniso_angle) {
        angle           = particles->aniso_angle[k];
      } else {
        angle           = materials->aniso_angle[particles->phase[k]];
      }
      if (model->polar == 1 ) {
        angle = angle - atan(particles->x[k] / (particles->z[k] ));
      }
      particles->nx[k]  = cos(angle);
      particles->nz[k]  = sin(angle);
      norm              = sqrt(particles->nx[k]*particles->nx[k] + particles->nz[k]*particles->nz[k]);
      particles->nx[k] /= norm;
      particles->nz[k] /= norm;
    }
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NormalizeDirector( grid* mesh, DoodzFP* nx_n, DoodzFP* nz_n, DoodzFP* nx_s, DoodzFP* nz_s, params *model ) {

  int Nx, Nz, k;
  double nx, nz, norm;
  Nx = model->Nx;
  Nz = model->Nz;

#pragma omp parallel for shared ( mesh, nx_n, nz_n ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
  for ( k=0; k<(Nx-1)*(Nz-1); k++ ) {
    if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
      nx            = nx_n[k];
      nz            = nz_n[k];
      norm          = sqrt(nx*nx + nz*nz);
      nx            = nx/norm;
      nz            = nz/norm;
      nx_n[k]       = nx;
      nz_n[k]       = nz;
    }
  }

#pragma omp parallel for shared ( mesh, nx_s, nz_s ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
  for ( k=0; k<Nx*Nz; k++ ) {
    if (mesh->BCg.type[k] != 30) {
      nx            = nx_s[k];
      nz            = nz_s[k];
      norm          = sqrt(nx*nx + nz*nz);
      nx            = nx/norm;
      nz            = nz/norm;
      nx_s[k]       = nx;
      nz_s[k]       = nz;
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
