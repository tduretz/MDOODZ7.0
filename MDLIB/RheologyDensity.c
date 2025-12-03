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
#include "RheologyDensity.h"

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

void EffectiveStrainRate( double* Exx, double* Ezz, double* Exz, double exx, double ezz, double exz, double Txx, double Tzz, double Txz, double d1, double d2, double ani, double eta_e, int elastic ) {
  if (elastic) {
    const double Da11  = 2.0 - 2.0*ani*d1;
    const double Da12  = 2.0*ani*d1;
    const double Da13  =-2.0*ani*d2;
    const double Da22  = 2.0 - 2.0*ani*d1;
    const double Da23  = 2.0*ani*d2;
    const double Da33  = 1.0  + 2.0*ani*(d1 - 0.5);
    const double a11   = Da33 * Da22 - pow(Da23,2);
    const double a12   = Da13 * Da23 - Da33 * Da12;
    const double a13   = Da12 * Da23 - Da13 * Da22;
    const double a22   = Da33 * Da11 - pow(Da13,2);
    const double a23   = Da12 * Da13 - Da11 * Da23;
    const double a33   = Da11 * Da22 - pow(Da12,2);
    const double det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
    const double iDa11 = a11/det; 
    const double iDa12 = a12/det; 
    const double iDa13 = a13/det;
    const double iDa22 = a22/det; 
    const double iDa23 = a23/det;
    const double iDa33 = a33/det;
    *Exx = exx + (iDa11*Txx + iDa12*Tzz + iDa13*Txz)/eta_e;
    *Ezz = ezz + (iDa12*Txx + iDa22*Tzz + iDa23*Txz)/eta_e;
    *Exz = exz + (iDa13*Txx + iDa23*Tzz + iDa33*Txz)/eta_e/2.0;
  }
  else {
    *Exx = exx;
    *Ezz = ezz; 
    *Exz = exz;
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void LocalIterationViscoElastic(LocalIterationMutables mutables, LocalIterationParams params) {
  const int    nitmax = 20;
  const double tol    = 1.0e-11;
  double eta_ve  = *mutables.eta;
  double Tii     = *mutables.Tii;
  double Eii_cst = 0., Eii_pwl = 0., Eii_gbs = 0., Eii_exp = 0., Eii_lin = 0.;
  if (params.noisy==1) printf("Start local iteration cycle VE for phase %d:\n", params.phase); 
  for (int it = 0; it < nitmax; it++) {
    // Function evaluation at current effective viscosity
    Tii     = 2.0 * eta_ve * params.Eii;
    if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
    if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
    if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
    if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                     // Peierls - power law
    if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(params.d, -params.m_lin); // !!! gs - dependence !!!
    const double            Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
    // Residual check
    const double r_eta_ve = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;
    const double res_eta  = fabs(r_eta_ve / params.Eii);
    if (params.noisy==1) printf("%d %2.4e %2.4e %2.4e\n", params.phase, r_eta_ve, Eii_pwl, Tii); // REMOVE
    if (res_eta < tol / 100) {
      if (it > 10) printf("L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
      break;
    } else if (it == nitmax - 1 && res_eta > tol) {
      printf("Visco-Elastic iterations failed!\n");
      exit(0);
    }
    // Analytical derivative of function
    double dfdeta = 0.0;
    if (params.elastic)     dfdeta += -params.Eii / params.eta_el;
    if (params.peierls)     dfdeta += -(Eii_exp) * (params.ST + params.n_exp) / eta_ve;
    if (params.diffusion)   dfdeta += -(Eii_lin) *params.n_lin / eta_ve;
    if (params.dislocation) dfdeta += -(Eii_pwl) *params.n_pwl / eta_ve;
    if (params.constant)    dfdeta += -params.Eii / params.eta_cst;
    // Update viscosity
    eta_ve -= r_eta_ve / dfdeta;
  }

  *mutables.eta     = eta_ve;
  *mutables.Tii     = Tii;
  *mutables.Eii_cst = Eii_cst;
  *mutables.Eii_exp = Eii_exp;
  *mutables.Eii_gbs = Eii_gbs;
  *mutables.Eii_pwl = Eii_pwl;
  *mutables.Eii_lin = Eii_lin;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void LocalIterationViscoElasticGrainSize(LocalIterationMutables mutables, LocalIterationParams params) {
  // Local iterations
  const int   nitmax = 20;
  const double tol   = 1.0e-11;
  double d_ve        = *mutables.d1;
  double eta_ve      = *mutables.eta;
  double Tii         = *mutables.Tii;
  const double p     = params.pg;
  const double Eii   = params.Eii;
  const double n_pwl = params.n_pwl;
  const double Ag    = params.Ag;
  const double Bg    = params.lam / (params.cg*params.gam); 
  const double C_pwl = params.C_pwl;
  double eta_up      = params.eta_up;
  double eta_lo      = params.eta_lo;
  double Eii_cst = 0., Eii_pwl = 0., Eii_gbs = 0., Eii_exp = 0., Eii_lin = 0.; 
  double Eii_vis, dddeta;

  // Function evaluation at current effective viscosity
  Tii     = 2.0 * eta_lo * params.Eii;
  if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
  if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
  d_ve    = pow( 2*Bg*params.Eii*Eii_pwl*eta_ve*p/params.Ag, (-1.0/(p+1)) );
  if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
  if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
  if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
  dddeta  = -1.0/2.0*Ag*d_ve*pow(2*Eii*eta_ve, -n_pwl)*(2*Bg*C_pwl*Eii*n_pwl*p*pow(2*Eii*eta_ve, n_pwl)/Ag + 2*Bg*C_pwl*Eii*p*pow(2*Eii*eta_ve, n_pwl)/Ag)/(Bg*C_pwl*Eii*eta_ve*p*(p + 1));
  Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
  // Residual check
  const double r_eta_lo = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;


  // Function evaluation at current effective viscosity
  Tii     = 2.0 * eta_up * params.Eii;
  if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
  if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
  d_ve    = pow( 2*Bg*params.Eii*Eii_pwl*eta_ve*p/params.Ag, (-1.0/(p+1)) );
  if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
  if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
  if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
  dddeta  = -1.0/2.0*Ag*d_ve*pow(2*Eii*eta_ve, -n_pwl)*(2*Bg*C_pwl*Eii*n_pwl*p*pow(2*Eii*eta_ve, n_pwl)/Ag + 2*Bg*C_pwl*Eii*p*pow(2*Eii*eta_ve, n_pwl)/Ag)/(Bg*C_pwl*Eii*eta_ve*p*(p + 1));
  Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
  // Residual check
  const double r_eta_up = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;

  // printf("%2.2e %2.2e\n", r_eta_lo, r_eta_up);
  // if ( (r_eta_up>0. && r_eta_lo>0.) || (r_eta_up<0. && r_eta_lo<0.) ) {printf("ouch!\n"); exit(113);}

  // 1 residual formulations --- Alternate bisection and Newton
  // if (params.noisy==1) printf("Start local iteration cycle GSE for phase %d:\n", params.phase);
  // for (int it = 0; it < nitmax; it++) {

  //   if (it<6) {
  //     eta_ve = 0.5*(eta_up + eta_lo);
  //     // Function evaluation at current effective viscosity
  //     Tii     = 2.0 * eta_ve * params.Eii;
  //     if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
  //     if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
  //     d_ve    = pow( 2*Bg*params.Eii*Eii_pwl*eta_ve*p/params.Ag, (-1.0/(p+1)) );
  //     if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
  //     if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
  //     if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
  //     const double            dddeta  = -1.0/2.0*Ag*d_ve*pow(2*Eii*eta_ve, -n_pwl)*(2*Bg*C_pwl*Eii*n_pwl*p*pow(2*Eii*eta_ve, n_pwl)/Ag + 2*Bg*C_pwl*Eii*p*pow(2*Eii*eta_ve, n_pwl)/Ag)/(Bg*C_pwl*Eii*eta_ve*p*(p + 1));
  //     const double            Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
  //     // Residual check
  //     const double r_eta_ve = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;
  //     const double res_eta  = fabs(r_eta_ve / params.Eii);
  //     if (params.noisy==1) printf("Bisect. %d %d %2.4e %2.4e %2.4e\n", it, params.phase, r_eta_ve, Eii_pwl, Tii); // REMOVE
  //     if (res_eta < tol / 100) {
  //       if (it > 10) printf("L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
  //       if (it > 10) exit(1);
  //       break;
  //     } else if (it == nitmax - 1 && res_eta > tol) {
  //       printf("Visco-Elastic iterations failed!\n");
  //       exit(0);
  //     }
  //     if ( (r_eta_ve>0. && r_eta_lo>0.) || (r_eta_ve<0. && r_eta_lo<0.) ) {
  //       eta_lo = eta_ve; 
  //     }
  //     else {
  //       eta_up = eta_ve;
  //     }
  //   }
  //   else {

  //     // Function evaluation at current effective viscosity
  //   Tii     = 2.0 * eta_ve * params.Eii;
  //   if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
  //   if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
  //   d_ve    = pow( 2*Bg*params.Eii*Eii_pwl*eta_ve*p/params.Ag, (-1.0/(p+1)) );
  //   if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
  //   if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
  //   if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
  //   const double            dddeta  = -1.0/2.0*Ag*d_ve*pow(2*Eii*eta_ve, -n_pwl)*(2*Bg*C_pwl*Eii*n_pwl*p*pow(2*Eii*eta_ve, n_pwl)/Ag + 2*Bg*C_pwl*Eii*p*pow(2*Eii*eta_ve, n_pwl)/Ag)/(Bg*C_pwl*Eii*eta_ve*p*(p + 1));
  //   const double            Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
  //   // Residual check
  //   const double r_eta_ve = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;
  //   const double res_eta  = fabs(r_eta_ve / params.Eii);
  //   if (params.noisy==1) printf("Newton %d %d %2.4e %2.4e %2.4e\n", it, params.phase, r_eta_ve, Eii_pwl, Tii); // REMOVE
  //   if (res_eta < tol / 100) {
  //     if (it > 10) printf("L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
  //     if (it > 10) exit(1);
  //     break;
  //   } else if (it == nitmax - 1 && res_eta > tol) {
  //     printf("Visco-Elastic iterations failed!\n");
  //     exit(0);
  //   }
  //   // Analytical derivative of function
  //   double dfdeta = 0.0;
  //   if (params.elastic)     dfdeta += -params.Eii / params.eta_el;
  //   if (params.constant)    dfdeta += -params.Eii / params.eta_cst;
  //   if (params.peierls)     dfdeta += -(Eii_exp) * (params.ST + params.n_exp) / eta_ve;
  //   if (params.diffusion)   dfdeta += -(Eii_lin) *params.n_lin / eta_ve;
  //   if (params.dislocation) dfdeta += -(Eii_pwl) *params.n_pwl / eta_ve;
  //   dfdeta += params.m_lin*Eii_lin*dddeta/d_ve;
  //   // Update viscosity
  //   // printf("dfdeta = %2.2e\n");
  //   eta_ve -= r_eta_ve / dfdeta;

  //   }

  // }

  // 1 residual formulations --- full Newton (see LocalIterationGrainSize.ipynb)
  if (params.noisy==1) printf("Start local iteration cycle GSE for phase %d:\n", params.phase);
  for (int it = 0; it < nitmax; it++) {
    // Function evaluation at current effective viscosity
    Tii     = 2.0 * eta_ve * params.Eii;
    if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
    if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
    d_ve    = pow( 2*Bg*params.Eii*Eii_pwl*eta_ve*p/params.Ag, (-1.0/(p+1)) );
    if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
    if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
    if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
    const double            dddeta  = -1.0/2.0*Ag*d_ve*pow(2*Eii*eta_ve, -n_pwl)*(2*Bg*C_pwl*Eii*n_pwl*p*pow(2*Eii*eta_ve, n_pwl)/Ag + 2*Bg*C_pwl*Eii*p*pow(2*Eii*eta_ve, n_pwl)/Ag)/(Bg*C_pwl*Eii*eta_ve*p*(p + 1));
    const double            Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
    // Residual check
    const double r_eta_ve = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;
    const double res_eta  = fabs(r_eta_ve / params.Eii);
    if (params.noisy==1) printf("%d %2.4e %2.4e %2.4e\n", params.phase, r_eta_ve, Eii_pwl, Tii); // REMOVE
    if (res_eta < tol / 100) {
      if (it > 15) printf("L.I. Warnung: more that 10 local iterations, there might be a problem...\n");
      // if (it > 15) exit(1);
      break;
    } else if (it == nitmax - 1 && res_eta > tol) {
      printf("Visco-Elastic iterations failed!\n");
      exit(0);
    }
    // Analytical derivative of function
    double dfdeta = 0.0;
    if (params.elastic)     dfdeta += -params.Eii / params.eta_el;
    if (params.constant)    dfdeta += -params.Eii / params.eta_cst;
    if (params.peierls)     dfdeta += -(Eii_exp) * (params.ST + params.n_exp) / eta_ve;
    if (params.diffusion)   dfdeta += -(Eii_lin) *params.n_lin / eta_ve;
    if (params.dislocation) dfdeta += -(Eii_pwl) *params.n_pwl / eta_ve;
    // only activate Newton after a couple of steps
    if (it>1) dfdeta += params.m_lin*Eii_lin*dddeta/d_ve;
    // Update viscosity
    eta_ve -= r_eta_ve / dfdeta;
  }

  // 2 residuals formulations
  // if (params.noisy==1) printf("Start local iteration cycle GSE for phase %d:\n", params.phase); 
  // for (int it = 0; it < nitmax; it++) {
  //   // Function evaluation at current effective viscosity
  //   Tii     = 2.0 * eta_ve * params.Eii;
  //   if (params.constant)    Eii_cst = Tii  / 2.0 / params.eta_cst;
  //   if (params.dislocation) Eii_pwl = params.C_pwl * pow(Tii, params.n_pwl);
  //   if (params.gbs)         Eii_gbs = params.C_gbs * pow(Tii, params.n_gbs);
  //   if (params.peierls)     Eii_exp = params.C_exp * pow(Tii, params.ST + params.n_exp);                 // Peierls - power law
  //   if (params.diffusion)   Eii_lin = params.C_lin * pow(Tii, params.n_lin) * pow(d_ve, -params.m_lin);// !!! gs - dependence !!!
  //   const double            Eii_vis = Eii_pwl + Eii_exp + Eii_lin + Eii_gbs + Eii_cst;
  //   // Residual check
  //   const double r_eta_ve = params.Eii - params.elastic * Tii  / (2.0 * params.eta_el) - Eii_vis;
  //   const double d_it     = exp(log(params.Ag * params.gam / (params.lam * (1.0 / params.cg) * Tii * Eii_pwl *params.pg)) / (1.0 + params.pg));
  //   const double r_d      = d_ve - d_it;
  //   const double res_eta  = fabs(r_eta_ve / params.Eii);
  //   const double res_d    = fabs(r_d);
  //   if (params.noisy==1) printf("%d %d %2.4e %2.4e %2.4e %2.4e\n", it, params.phase, res_eta, res_d, Eii_pwl, Tii); // REMOVE
  //   if (res_eta < tol / 100 && res_d < tol) {
  //     if (it > 20) printf("L.I. GSE Warnung: more that 20 local iterations, there might be a problem...\n");
  //     break;
  //   } else if (it == nitmax - 1 && res_eta > tol) {
  //     printf("Visco-Elastic iterations failed!\n");
  //     exit(0);
  //   }
  //   // Analytical derivative of function
  //   double dr_eta_deta = 0.0;
  //   if (params.elastic)     dr_eta_deta += -params.Eii / params.eta_el;
  //   if (params.peierls)     dr_eta_deta += -(Eii_exp) * (params.ST + params.n_exp) / eta_ve;
  //   if (params.diffusion)   dr_eta_deta += -(Eii_lin) *params.n_lin / eta_ve;
  //   if (params.dislocation) dr_eta_deta += -(Eii_pwl) *params.n_pwl / eta_ve;
  //   if (params.constant)    dr_eta_deta += -params.Eii / params.eta_cst;
  //   const double dr_eta_dd = Eii_lin * params.m_lin / (d_ve);
  //   const double dr_d_dd   = 1.0;
  //   const double dr_d_deta = -2.0 * params.Eii * (Eii_pwl) *d_it * eta_ve * params.lam * params.pg * (-0.5 * params.Ag * params.cg * params.gam * params.n_pwl / (params.Eii * (Eii_pwl) *pow(eta_ve, 2) * params.lam * params.pg) - 0.5 * params.Ag * params.cg * params.gam / (params.Eii * (Eii_pwl) *pow(eta_ve, 2) * params.lam * params.pg)) / (params.Ag * params.cg * params.gam * (params.pg + 1.0));
  //   // Inverse of the Jacobian: TODO: call function Solve2x2()
  //   const double det       = dr_eta_deta * dr_d_dd - dr_eta_dd * dr_d_deta; // determinant
  //   const double ai        =  1.0 / det * dr_d_dd;                          // inverse matrix components
  //   const double bi        = -1.0 / det * dr_eta_dd;
  //   const double ci        = -1.0 / det * dr_d_deta;
  //   const double di        =  1.0 / det * dr_eta_deta;
  //   const double deta      = (ai * r_eta_ve + bi * r_d);// inverse times rhs
  //   const double dd        = (ci * r_eta_ve + di * r_d);
  //   // Coupled update of viscosity and grain size
  //   eta_ve -= deta;
  //   d_ve    -= dd;
  // }

  if (isnan(eta_ve)) {
    printf("Local iterations from grain size and viscoisity exploded -  guess is probably incorrect");
  }

  *mutables.eta     = eta_ve;
  *mutables.d1      = d_ve;
  *mutables.Tii     = Tii;
  *mutables.Eii_cst = Eii_cst;
  *mutables.Eii_exp = Eii_exp;
  *mutables.Eii_gbs = Eii_gbs;
  *mutables.Eii_pwl = Eii_pwl;
  *mutables.Eii_lin = Eii_lin;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ViscosityConcise( int phase, double G, double T, double P, double d, double phi, double X0, double Exx, double Ezz, double Exz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *Txx, double *Tzz, double *Txz, double* etaVE, double* VEcoeff, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp , double* Eii_lin, double* Eii_gbs, double* Eii_cst, double *d1, double strain_acc, double dil, double fric, double C, double P0, double T0,  double *X1, double *OverS, double *Pcorr, double *rho, double beta, double div, double *div_el, double *div_pl, double *div_r, double *Wtot, double *Wel, double *Wdiss, int post_process, int centroid, int final_update, int index ) {

  // printf("P0 in 1 = %.4e\n", P0);
  // General paramaters
  int    it, nitmax = 20;
  int noisy = 0;
  double eta = 0.0, R = materials->R;
  double min_eta = model->min_eta, max_eta = model->max_eta;
  double TmaxPeierls = (1200.0 + zeroC) / scaling->T;// max. T for Peierls

  double dt;
  if (model->step==0) dt = 1e100;
  else dt = model->dt;

  // Parameters for deformation map calculations
  int    plastic = 0, tensile = 0, constant = 0, dislocation = 0, peierls = 0, diffusion = 0, gbs = 0, elastic = model->elastic, kinetics = 0;
  double tol = model->eta_tol, res_pl = 0.0, Tii = 0.0, Tii0 = sqrt(Txx0 * Txx0 + Txz0 * Txz0);
  double eta_lo = 0.0, eta_up = 0.0, eta_ve = 0.0;
  double eta_pwl = 0.0, eta_exp = 0.0, eta_vep = 0.0, eta_lin = 0.0, eta_el = 0.0, eta_gbs = 0.0, eta_cst = 0.0, eta_pl = 0.0;
  double Eii = 0.0;
  double X   = X0;

  // Flow law parameters from input file
  double Tyield = 0.0, F_trial = 0.0, gdot = 0.0, dQdtxx = 0.0, dQdtzz = 0.0, dQdtxz = 0.0;
  int    is_pl  = 0;
  double Ea_pwl = materials->Qpwl[phase], Va_pwl = materials->Vpwl[phase], n_pwl = materials->npwl[phase], m_pwl = materials->mpwl[phase], r_pwl = materials->rpwl[phase], A_pwl = materials->Apwl[phase], f_pwl = materials->fpwl[phase], a_pwl = materials->apwl[phase], F_pwl = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase], t_pwl = materials->tpwl[phase];
  double Ea_lin = materials->Qlin[phase], Va_lin = materials->Vlin[phase], n_lin = materials->nlin[phase], m_lin = materials->mlin[phase], r_lin = materials->rlin[phase], A_lin = materials->Alin[phase], f_lin = materials->flin[phase], a_lin = materials->alin[phase], F_lin = materials->Flin[phase];
  double Ea_gbs = materials->Qgbs[phase], Va_gbs = materials->Vgbs[phase], n_gbs = materials->ngbs[phase], m_gbs = materials->mgbs[phase], r_gbs = materials->rgbs[phase], A_gbs = materials->Agbs[phase], f_gbs = materials->fgbs[phase], a_gbs = materials->agbs[phase], F_gbs = materials->Fgbs[phase];
  double B_pwl = 0.0, B_lin = 0.0, B_exp = 0.0, B_gbs = 0.0, C_pwl = 0.0, C_lin = 0.0, C_exp = 0.0;
  double Ea_exp = materials->Qexp[phase], S_exp = materials->Sexp[phase], E_exp = materials->Eexp[phase], F_exp = 0.0;
  double gamma = materials->Gexp[phase], ST=0.0, n_exp = materials->nexp[phase];
  double Exx_lin = 0.0, Ezz_lin = 0.0, Exz_lin = 0.0, Exx_exp = 0.0, Ezz_exp = 0.0, Exz_exp = 0.0, Exx_gbs = 0.0, Ezz_gbs = 0.0, Exz_gbs = 0.0, Exx_cst = 0.0, Ezz_cst = 0.0, Exz_cst = 0.0;
  double Exx_pl = 0.0, Exx_pwl = 0.0, Ezz_pl = 0.0, Ezz_pwl = 0.0, Exz_pl = 0.0, Exz_pwl = 0.0;
  int    gs = materials->gs[phase];
  double Ag;
  double Qkin = materials->Qkin[phase], Skin = materials->Skin[phase], kkin = materials->kkin[phase], Vkin, dG = -1.0, rho_eq, rho0;

  double eta_vp0 = materials->eta_vp[phase], n_vp = materials->n_vp[phase], eta_vp = materials->eta_vp[phase], T_st = materials->T_st[phase];
  double dQdTii = 1.0,  dQdP = 0.0, K = 1.0 / beta;
  double F_trial0 = F_trial;
  double dFdgdot, divp = 0.0, Pc = P;
  double Pmin = 1e5/scaling->S, Tiimin = 1e5/scaling->S;

  // Mix of constant viscosity
  int    phase_two           = materials->phase_two[phase];
  int    constant_mix        = materials->phase_mix[phase];
  int    mix_avg             = model->diffuse_avg;
  double rho1                = materials->rho[phase];
  double rho2                = materials->rho[phase];
  double tau_kin = materials->tau_kin[phase], Pr = materials->Pr[phase], dPr = materials->dPr[phase];
  bool   phase_chemical_production = model->chemical_production && materials->reac_soft[phase];
  double divr = 0.0;
  if (model->diffuse_X == 0) constant_mix = 0;

  // Partial melting  
  int    Melting = model->melting;
  double phi0 = phi, phi1 = phi, Ql1 = 0.;

  //------------------------------------------------------------------------//

  // Initialise strain rate invariants to 0
  *etaVE   = 0.0, *VEcoeff = 0.0;
  *Eii_exp = 0.0, *Eii_lin = 0.0, *Eii_pl = 0.0, *Eii_pwl = 0.0, *Eii_el = 0.0, *Eii_gbs = 0, *Eii_cst = 0.0;
  *Txx     = 0.0, *Tzz     = 0.0, *Txz    = 0.0;
  *Wtot    = 0.0, *Wel     = 0.0, *Wdiss  = 0.0, *d1      = 0.0;
  *div_el  = 0.0, *div_pl  = 0.0, *div_r  = 0.0;
  *X1      = 0.0, *rho     = 0.0, *OverS  = 0.0, *Pcorr   = 0.0;

  //------------------------------------------------------------------------//

  // Invariants
  double Eyy = -(Exx + Ezz);// definition of deviatoric tensor
  Eii   = sqrt(1.0 / 2.0 * (Exx * Exx + Ezz * Ezz + Eyy * Eyy) + Exz * Exz);
  if (Eii * scaling->E < 1e-30) Eii = 1e-30 / scaling->E;

  //    printf("E --> %2.2e %2.2e %2.2e\n", Exx, Ezz, Eyy); // if (fabs(Eyy)>1e-6)
  //    if (fabs(Eyy)>1e-11) exit(1);

  // P corr will be corrected if plasticity feedbacks on pressure (dilation)
  *Pcorr = P;

  //------------------------------------------------------------------------//
// printf("%d\n", materials->cstv[phase]);

  // Activate deformation mechanisms
  if ( materials->cstv[phase]   !=0                  ) constant    = 1;
  // printf("%d\n", materials->cstv[phase]);
  if ( materials->pwlv[phase]   !=0                  ) dislocation = 1;
  if ( materials->expv[phase]   !=0 && T<TmaxPeierls ) peierls     = 1;
  if ( materials->linv[phase]   !=0                  ) diffusion   = 1;
  if ( materials->gbsv[phase]   !=0                  ) gbs         = 1;
  if ( materials->gs[phase]     !=0                  ) gs          = 1;
  if ( materials->kin[phase]    !=0                  ) kinetics    = 1;
  if ( materials->plast[phase]  !=0                  ) plastic     = 1;
  if ( materials->tensile[phase]!=0                  ) tensile     = 1;

  // Turn of elasticity for the initialisation step  (viscous flow stress)
  if ( model->step    == 0                         ) elastic     = 0;

  // Constant grain size initially
  *d1 = materials->gs_ref[phase];

  // Tensional cut-off
  // if ( model->gz<0.0 && P<0.0     ) { P = 0.0; } // printf("Aie aie aie P < 0 !!!\n"); exit(122);
  // if ( model->gz<0.0 && P<0.0     ) fric = 0.0;

  // Visco-plastic limit
  if ( elastic==0                 ) { G = 1e1; dil = 0.0; K = 1e10;}; //K = 1e1;

  // Zero C limit
  if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;

  // Avoid redundant plasticity formulation
  if (plastic == 1 && tensile == 1){ printf("You activitate both mode-II and combined mode-I/mode-II plasticity. Choose one and deactivate the other\n"); exit(0); }
  
  //------------------------------------------------------------------------//

  // Precomputations
  if ( dislocation == 1 ) {
    B_pwl = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl + P*Va_pwl)/R/n_pwl/T ) * pow(d, m_pwl/n_pwl) * pow(f_pwl, -r_pwl/n_pwl) * exp(-a_pwl*phi/n_pwl);
    C_pwl   = pow(2.0*B_pwl, -n_pwl);
  }
  if ( diffusion == 1 ) {
    if (m_lin>0.0 && d<1e-13/scaling->L) {
      printf("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!\n", d*scaling->L);
      exit(1);
    }
    B_lin = F_lin * pow(A_lin,-1.0/n_lin) * exp( (Ea_lin + P*Va_lin)/R/n_lin/T ) * pow(f_lin, -r_lin/n_lin) * exp(-a_lin*phi/n_lin); // * pow(d, m_lin/n_lin) !!!!!!!!!!!!!!!!!!!!!!!!
    C_lin = pow(2.0*B_lin, -n_lin);
  }
  if ( gbs == 1 ) {
    B_gbs = F_gbs * pow(A_gbs,-1.0/n_gbs) * exp( (Ea_gbs + P*Va_gbs)/R/n_gbs/T ) * pow(d, m_gbs/n_gbs) * pow(f_gbs, -r_gbs/n_gbs) * exp(-a_gbs*phi/n_gbs);
  }
  if ( peierls   == 1 ) {
    ST                           = Ea_exp/R/T * 2.0*gamma*(1-gamma);
    // new
    double Arr_exp = exp(-Ea_exp/R/T*pow(1.0-gamma,2.0));
    F_exp = pow( pow(2.0,1.0-ST-n_exp) / pow(sqrt(3.0), ST+n_exp+1.0), 1.0/(ST+n_exp));
    B_exp = F_exp * ( pow(gamma*S_exp, ST/(ST+n_exp)) / pow( E_exp*Arr_exp, 1.0/(ST+n_exp)) );
    C_exp = pow(2.0*B_exp, -(ST+n_exp));
  }

  double cos_fric = cos(fric);
  double sin_fric = sin(fric);
  double sin_dil  = sin(dil);

      // printf("P0 in 2 = %.4e\n", P0);
  if ( plastic==1 ) {

    Tyield     = C*cos_fric + P*sin_fric;

    // Von-Mises cut-off
    if (materials->Slim[phase] < Tyield) {
      C          = materials->Slim[phase];
      Tyield     = materials->Slim[phase];
    }
  }

  //------------------------------------------------------------------------//
  // Reaction stuff: 1. Update reaction progress
  X = X0;

  // Reaction stuff: 2. Mixture rheology (Huet et al., 2014)
  if (phase_chemical_production) {
      X     = (X0*tau_kin - 0.5*dt*erfc((P - Pr)/dPr) + dt)/(dt + tau_kin);
      HuetAveragingModel( &B_pwl, &C_pwl, &n_pwl, phase, R, T, t_pwl, X, pre_factor, materials );
  }

  if ( Melting == 1 ) PartialMelting( &phi1, &Ql1, P, T, phi0, tau_kin, dt, materials->melt[phase], scaling );
  
  // Set pointer value
  *X1 = X;

  //------------------------------------------------------------------------//

  // Isolated viscosities
  eta_el   = G*dt; 
  if ( constant    == 1 )              eta_cst  = materials->eta0[phase]; 
  if ( constant_mix== 1 && mix_avg==0) eta_cst  =     X0*materials->eta0[phase] + (1-X0)*materials->eta0[phase_two];
  if ( constant_mix== 1 && mix_avg==1) eta_cst  = pow(X0/materials->eta0[phase] + (1-X0)/materials->eta0[phase_two], -1.0);
  if ( constant_mix== 1 && mix_avg==2) eta_cst  = exp(X0*log(materials->eta0[phase]) + (1-X0)*log(materials->eta0[phase_two]));
  if ( dislocation == 1 )              eta_pwl  = B_pwl * pow( Eii, 1.0/n_pwl - 1.0 ); 

  if (gs == 1) {
    double Kg = materials->Kpzm[phase], Qg = materials->Qpzm[phase];
    Ag = Kg * exp(-Qg / R / T);
    const double gam         = materials->Gpzm[phase];
    const double lam         = materials->Lpzm[phase];
    const double cg          = materials->cpzm[phase];
    const double pg          = materials->ppzm[phase];
    Tii          = 2.*eta_el*Eii;
    *Eii_pwl     = Tii/2.0/eta_pwl;
    d        = exp(log(Ag * gam / (lam * (1.0 / cg) * Tii * (*Eii_pwl) *pg)) / (1.0 + pg));
  }

  if ( diffusion   == 1 )              eta_lin  = B_lin * pow( Eii, 1.0/n_lin - 1.0 ) * pow(d, m_lin/n_lin);
  if ( gbs         == 1 )              eta_gbs  = B_gbs * pow( Eii, 1.0/n_gbs - 1.0 ); 
  if ( peierls     == 1 )              eta_exp  = B_exp * pow( Eii, 1.0 / (ST + n_exp) - 1.0); 
  //------------------------------------------------------------------------//

  // Viscoelasticity
  *Eii_pl = 0.0;
   eta_up =  1.0e100 / scaling->eta;
   eta_lo = -1.0e100 / scaling->eta;
  if (constant == 1)    {eta_up = MINV(eta_up, eta_cst); eta_lo = MAXV(eta_lo, eta_cst);}
  if (dislocation == 1) {eta_up = MINV(eta_up, eta_pwl); eta_lo = MAXV(eta_lo, eta_pwl);}
  if (elastic == 1)     {eta_up = MINV(eta_up, eta_el) ; eta_lo = MAXV(eta_lo, eta_el) ;}
  if (peierls == 1)     {eta_up = MINV(eta_up, eta_exp); eta_lo = MAXV(eta_lo, eta_exp);}
  if (diffusion == 1)   {eta_up = MINV(eta_up, eta_lin); eta_lo = MAXV(eta_lo, eta_lin);}
  if (gbs == 1)         {eta_up = MINV(eta_up, eta_gbs); eta_lo = MAXV(eta_lo, eta_gbs);}
  
  //------------------------------------------------------------------------//

  // Initial guess 
  eta_ve = eta_up; // better
  Tii    = 2.0 * eta_ve * Eii;

  if (gs == 1) {
    double Kg = materials->Kpzm[phase], Qg = materials->Qpzm[phase];
    Ag = Kg * exp(-Qg / R / T);
    const double gam         = materials->Gpzm[phase];
    const double lam         = materials->Lpzm[phase];
    const double cg          = materials->cpzm[phase];
    const double pg          = materials->ppzm[phase];
    const double Eii_pwl     = Tii/2.0/eta_pwl;
    const double d_it        = exp(log(Ag * gam / (lam * (1.0 / cg) * Tii * (Eii_pwl) *pg)) / (1.0 + pg));
    *d1 = d;
    if (noisy) printf("d guess = %2.2e, Tii = %2.2e Eiipwl=%2.2e\n", d_it*scaling->L, Tii*scaling->S, Eii_pwl*scaling->E);
  }
  
  LocalIterationParams params = {
          .noisy       = noisy,
          .phase       = phase,
          .Eii         = Eii,
          .gam         = materials->Gpzm[phase],
          .lam         = materials->Lpzm[phase],
          .cg          = materials->cpzm[phase],
          .pg          = materials->ppzm[phase],
          .Ag          = Ag,
          .C_pwl       = C_pwl,
          .n_pwl       = n_pwl,
          .C_lin       = C_lin,
          .n_lin       = n_lin,
          .m_lin       = m_lin,
          .ST          = ST,
          .n_exp       = n_exp,
          .elastic     = elastic,
          .eta_el      = eta_el,
          .peierls     = peierls,
          .C_exp       = C_exp,
          .dislocation = dislocation,
          .diffusion   = diffusion,
          .constant    = constant,
          .eta_cst     = eta_cst,
          .gbs         = gbs,
          .d           = d,
          .n_gbs       = n_gbs,
          .eta_up      = eta_up,
          .eta_lo      = eta_lo};
  if (gs) {
    LocalIterationViscoElasticGrainSize((LocalIterationMutables){.eta = &eta_ve, .Tii = &Tii, .Eii_cst=Eii_cst, .Eii_lin=Eii_lin, .Eii_gbs=Eii_gbs, .Eii_pwl=Eii_pwl, .Eii_exp=Eii_exp, .d1 = d1}, params);
  } else {
    LocalIterationViscoElastic((LocalIterationMutables){.eta = &eta_ve, .Tii = &Tii, .Eii_cst=Eii_cst, .Eii_lin=Eii_lin, .Eii_gbs=Eii_gbs, .Eii_pwl=Eii_pwl, .Eii_exp=Eii_exp}, params);
  }
  if (noisy) printf("d  = %2.2e, Tii = %2.2e \n", *d1*scaling->L, Tii*scaling->S);  //------------------------------------------------------------------------//

  // if (centroid==0 && index==65) printf("eta_ve = %2.4e Eii_pwl = %2.4e Eii=%2.4e C_pwl=%2.4e\n", eta_ve, *Eii_pwl, Eii, C_pwl);

  is_pl = 0;
  dQdP = 0.0;
  //------------------------------------------------------------------------//
  if (plastic == 1) {
    // Check yield stress
    F_trial = Tii - Tyield;

    // printf("Tii=%2.2e F=%2.2e C=%2.2e cos_fric=%2.2e P=%2.2e sin_fric=%2.2e\n", Tii*scaling->S, F_trial*scaling->S, C*scaling->S, cos_fric, P*scaling->S, sin_fric);

    double Tiic;

    if (F_trial > 1e-17) {

      // Initial guess - eta_vp = 0
      is_pl    = 1;
      eta_vp   = eta_vp0 * pow(Eii, 1.0/n_vp - 1);
      gdot     = F_trial / ( eta_ve + eta_vp + K*dt*sin_fric*sin_dil);
      dQdP     = -sin_dil; 
      F_trial0 = F_trial;

      // Return mapping --> find plastic multiplier rate (gdot)
      for (it=0; it<nitmax; it++) {

        dQdP    = -sin_dil;
        divp    = -gdot*dQdP;
        Pc      = P + K*dt*divp; // P0 - k*dt*(div-divp) = P + k*dt*divp
        eta_vp  = eta_vp0 * pow(fabs(gdot), 1.0/n_vp - 1.0);
        Tyield  = Tyield + gdot*eta_vp;
        Tiic    = Tii - eta_ve*gdot;
        F_trial = Tiic - Tyield;

        // Residual check
        res_pl = fabs(F_trial);
        if ( noisy>0 ) printf("%02d Viscoplastic iterations It., F = %2.2e Frel = %2.2e --- n_vp = %2.2e, eta_vp = %2.2e\n", it, res_pl, res_pl/F_trial0, n_vp, eta_vp*scaling->eta);
        if ( res_pl < tol || res_pl/F_trial0 < tol ) break;
        dFdgdot  = - eta_ve - eta_vp/n_vp - K*dt*sin_fric*sin_dil;
        gdot    -= F_trial / dFdgdot;

      }
      if ( noisy>0 && it==nitmax-1 && (res_pl > tol || res_pl/F_trial0 > tol)  ) { printf("Visco-Plastic iterations failed!\n"); exit(0);}

      // In case return mapping has failed (because of tension), return to a von Mises minimum stress 
      if (Tiic<0.0) {
        //printf("Aie, tension!\n"); exit(1);
        F_trial = Tii - Tiimin;
        gdot    = F_trial /  (eta_ve);
        Tiic    = Tii - eta_ve*gdot;
        Tiic    = Tiimin;
        Pc      = P;
      }

      *Pcorr  = Pc;
      eta_vep = Tiic / (2.0*Eii);
      Tii     = Tiic;
    }
  }

      double lam_dot = 0.0;


  // ----------------- Combined mode-I mode-II formulation of Popov et al. 2025, GMD
  if (tensile == 1)
  {
    // Step 1: Get tensile yield parameters
    double k        = sin_fric;
    double kq       = sin_dil;
    double c        = C * cos_fric;
    double a        = sqrt(1.0 + k*k);
    double b        = sqrt(1.0 + kq*kq);
    double py       = (T_st + c/a) / (1.0 - k/a);     // Center of yield cap circle
    double Ry       = py - T_st;                      // Radius of yield cap circle
    double pd       = py - Ry * k/a;                  // Pressure coordinate of the delimiter point (yield)
    double tau_d    = k * pd + c;                     // Stress coordinate of the delimiter point (yield)
    double pq       = pd + kq * tau_d;                // Center of potential flow cap circle
    // double Rq       = pq - T_st;                      // Radius of potential flow cap circle
    double R_haty   = sqrt(Tii*Tii + (P - py)* (P - py));
    double R_hatq   = sqrt(Tii*Tii + (P - pq)* (P - pq));
    int    on_yield = 0, on_flow = 0;
    gdot = 0.0, divp = 0.0;
    on_yield = (Tii *(py - pd)  >= tau_d * (py - P));

    is_pl = 0;

    // Step 2: Evaluate yield function
    if (on_yield)     // On Drucker-Prager yield
    {
      F_trial = Tii - k * P - c;
    } else                                       // On Cap
    {
      F_trial = a * (R_haty - Ry);
    }

    // Step 3: Start correction if above the yield
    if (F_trial > 1e-17)
    {
      // Initialize
      double Tiic;
      double d2QdTii2, d2QdP2, d2QdTP, d2QdPT;
      double dFdTii, dFdP;
      double R1, R2, R3, J11, J12, J13, J21, J22, J23, J31, J32, J33, dTii, dP, dlam_dot, R_norm, R0_norm, atol = 1e-16, rtol = 1e-14;

      // Newton iterations
      is_pl = 1, Pc = P, Tiic = Tii, lam_dot = 0.0;
      for (int it = 0; it < nitmax; it++)
      {

        // Yield and flow conditions
        on_yield = (Tiic *(py - pd) * 1e40) >= (tau_d * (py - Pc));
        on_flow  = (Tiic *(pq - pd) * 1e40) >= (tau_d * (pq - Pc));

        // Regularization viscosity
        double abs_lam_dot = fabs(lam_dot);
        if (abs_lam_dot < 1e-30) abs_lam_dot = 1e-30;
        eta_vp  = eta_vp0; // * pow(abs_lam_dot, 1.0/n_vp - 1.0);

        // Evaluate yield function
        R_haty   = sqrt(Tiic*Tiic + (Pc - py)* (Pc - py));
        R_hatq   = sqrt(Tiic*Tiic + (Pc - pq)* (Pc - pq));
        if (on_yield)     // On Drucker-Prager yield
        {
          F_trial = Tiic - k * Pc - c;
          dFdTii = 1.0;
          dFdP   = -k;
        }else  // On Cap
        {
          F_trial = a * (R_haty - Ry);
          dFdTii  = a * Tiic      / R_haty;
          dFdP    = a * (Pc - py) / R_haty;
        }

        // Macauly brackets
        double F = F_trial, Fpos, dFdTii_eff, dFdP_eff;
        if (F > 1e-17)
        {
          Fpos        = F;
          dFdTii_eff  = dFdTii;
          dFdP_eff    = dFdP;
        } else
        {
          Fpos       = 0.0;
          dFdTii_eff = 0.0;
          dFdP_eff   = 0.0;
        }

        
        
        // Compute flow potential derivatives
        if (on_flow)
        {
          dQdTii = 0.5       * 2; 
          dQdP   = -kq; // <- This might be -kq depending on the sign convention for pressure
          d2QdTii2 = 0.0;
          d2QdP2   = 0.0;
          d2QdTP   = 0.0;
          d2QdPT   = 0.0;
        }else
        {
          dQdTii   = 0.5 * b * Tiic / R_hatq          *       2;
          dQdP     = -b * (Pc - pq) / R_hatq;
          d2QdTii2 = 0.5 * b * (R_hatq*R_hatq - Tiic*Tiic) / (R_hatq*R_hatq*R_hatq);
          d2QdP2   = -b * (R_hatq*R_hatq - (Pc - pq)*(Pc - pq)) / (R_hatq*R_hatq*R_hatq);
          d2QdTP   = -0.5 * b * Tiic * (Pc - pq) / (R_hatq*R_hatq*R_hatq);
          d2QdPT   =  b * Tiic * (Pc - pq) / (R_hatq*R_hatq*R_hatq);
        }

        // Local linear system of equations - Residual form
        R1 = (Tii - Tiic) / 2.0 / eta_ve - lam_dot * dQdTii;
        R2 = (Pc - P) / (K * dt)         - lam_dot * dQdP;
        R3 = Fpos                        - lam_dot * eta_vp;

        // Stop if Newton residual did not converge
        if (it == 0)
        {
          R0_norm = sqrt(R1*R1 + R2*R2 + R3*R3);
        }
        R_norm = sqrt(R1*R1 + R2*R2 + R3*R3);

         if (centroid) {
        if (on_yield)
        {
          printf("No. iter = %05d; mode = 2; R1 = %.4e; R2 = %.4e; R3 = %.4e; R_norm = %.4e; R/R0 = %.4e\n", it, R1, R2, R3, R_norm, R_norm / R0_norm);
        }else
        {
          printf("No. iter = %05d; mode = 1; R1 = %.4e; R2 = %.4e; R3 = %.4e; R_norm = %.4e; R/R0 = %.4e\n", it, R1, R2, R3, R_norm, R_norm / R0_norm);
        }
      }

        
        // printf("R_norm = %.4e; R0_norm = %.4e; R/R0 = %.4e; nitmax = %d; atol = %.4e; rtol = %.4e;\n", R_norm, R0_norm, R_norm / R0_norm, nitmax, atol, rtol);
        if (R_norm < atol || R_norm / R0_norm < rtol) {
            if (noisy > 0) {printf("Visco-Plastic Newton converged in %d iters\n", it+1);};
            break;
        }
        
        // Assemble Jacobian
        J11 = -1.0 / 2.0 / eta_ve - lam_dot * d2QdTii2;
        J12 =                     - lam_dot * d2QdTP;
        J13 =                                -dQdTii;
        
        J21 =                  -lam_dot * d2QdPT;
        J22 =  1.0 / (K * dt) - lam_dot * d2QdP2;
        J23 =                            -dQdP;
        
        J31 = dFdTii_eff;
        J32 = dFdP_eff;
        J33 = -eta_vp;

        // Solve the linear system J * dx = R
        int ok = solve_local_3x3(
            J11, J12, J13,
            J21, J22, J23,
            J31, J32, J33,
            R1, R2, R3,
            &dTii, &dP, &dlam_dot
        );
        // Sanity check
        if (!ok) {
            fprintf(stderr, "ERROR: Singular Jacobian in Newton local solve.\n");
            exit(1);
        }

        // // Apply line search to find suitable Newton steps
        // double Tiic_ls, Pc_ls, lam_dot_ls, alpha, R0_norm_ls, R_norm_ls, F_trial_ls;
        // double c_ls = 1e-4;
        double alpha = 1.0;
        // R0_norm_ls = sqrt(R1*R1 + R2*R2 + R3*R3);

        // // Find good fraction of a step
        // for (int it_ls = 0; it_ls < 100; it_ls++)
        // {
        //   // Compute trial updates
        //   Tiic_ls    = Tiic    + dTii     * alpha;
        //   Pc_ls      = Pc      + dP       * alpha;
        //   lam_dot_ls = lam_dot + dlam_dot * alpha;
        //   lam_dot_ls = fmax(lam_dot_ls, 0.0);

        //   // Yield and flow conditions
        //   on_yield = (Tiic_ls *(py - pd) >= tau_d * (py - Pc_ls));
        //   on_flow  = (Tiic_ls *(pq - pd) >= tau_d * (pq - Pc_ls));

        //   // Regularization viscosity
        //   double abs_lam_dot_ls = fabs(lam_dot_ls);
        //   if (abs_lam_dot_ls < 1e-30) abs_lam_dot_ls = 1e-30;
        //   eta_vp  = eta_vp0; // * pow(abs_lam_dot_ls, 1.0/n_vp - 1.0);

        //   // Evaluate yield function
        //   R_haty   = sqrt(Tiic_ls*Tiic_ls + (Pc_ls - py)* (Pc_ls - py));
        //   R_hatq   = sqrt(Tiic_ls*Tiic_ls + (Pc_ls - pq)* (Pc_ls - pq));
        //   if (on_yield)     // On Drucker-Prager yield
        //   {
        //     F_trial_ls = Tiic_ls - k * Pc_ls - c;
        //     dFdTii     = 1.0;
        //     dFdP       = -k;
        //   }else  // On Cap
        //   {
        //     F_trial_ls = a * (R_haty - Ry);
        //     dFdTii     = a * Tiic_ls      / R_haty;
        //     dFdP       = a * (Pc_ls - py) / R_haty;
        //   }

        //   // Macauly brackets
        //   F = F_trial_ls;
        //   if (F > 1e-17)
        //   {
        //     Fpos        = F;
        //     dFdTii_eff  = dFdTii;
        //     dFdP_eff    = dFdP;
        //   } else
        //   {
        //     Fpos       = 0.0;
        //     dFdTii_eff = 0.0;
        //     dFdP_eff   = 0.0;
        //   }

        //   // Compute flow potential derivatives
        //   if (on_flow)
        //   {
        //     dQdTii = 0.5;
        //     dQdP   = -kq; // <- This might be -kq depending on the sign convention for pressure
        //     d2QdTii2 = 0.0;
        //     d2QdP2   = 0.0;
        //     d2QdTP   = 0.0;
        //     d2QdPT   = 0.0;
        //   }else
        //   {
        //     dQdTii   = 0.5 * b * Tiic_ls / R_hatq;
        //     dQdP     = -b * (Pc_ls - pq) / R_hatq;
        //     d2QdTii2 = 0.5 * b * (R_hatq*R_hatq - Tiic_ls*Tiic_ls) / (R_hatq*R_hatq*R_hatq);
        //     d2QdP2   = -b * (R_hatq*R_hatq - (Pc_ls - pq)*(Pc_ls - pq)) / (R_hatq*R_hatq*R_hatq);
        //     d2QdTP   = -0.5 * b * Tiic_ls * (Pc_ls - pq) / (R_hatq*R_hatq*R_hatq);
        //     d2QdPT   =  b * Tiic_ls * (Pc_ls - pq) / (R_hatq*R_hatq*R_hatq);
        //   }

        //   // Compute residuals
        //   R1 = (Tii - Tiic_ls) / 2.0 / eta_ve - lam_dot_ls * dQdTii;
        //   R2 = (Pc_ls - P) / (K * dt)         - lam_dot_ls * dQdP;
        //   R3 = Fpos                           - lam_dot_ls * eta_vp;
        //   R_norm_ls = sqrt(R1*R1 + R2*R2 + R3*R3);

        //   // printf("Line search it. = %02d: R0 = %.8e; R = %.8e; alpha = %.2e\n", it_ls, R0_norm_ls, R_norm_ls, alpha);
        //   // Accept alpha, if the residual is small enough
        //   if (R_norm_ls < (1.0 - c_ls * alpha) * R0_norm_ls)
        //   {
        //     // printf("Line search converged in %03d. Accepted alpha = %.2e\n", it_ls, alpha);
        //     break;
        //   }else
        //   {
        //     alpha *= 0.5;
        //   }
        //   // Restrict alpha to a minimum value
        //   if (alpha < 0.1)
        //   {
        //     alpha = 0.1;
        //     break;
        //   }
        // }
        
        // Correct stress, pressure and multiplier
        Tiic    += dTii     * alpha;
        Pc      += dP       * alpha;
        lam_dot += dlam_dot * alpha;
        lam_dot  = fmax(lam_dot, 0.0);

        // Failure statement
        if ( noisy>0 && it==nitmax-1 && (R_norm > atol || R_norm/R0_norm > rtol)  ) { printf("Visco-Plastic iterations failed!\n"); exit(0);}
      }

      // Otherwise store values
      *Pcorr  = Pc;
      eta_vep = Tiic / (2.0*Eii);
      Tii     = Tiic;
      gdot    = 2*lam_dot * dQdTii;
      divp    = 1*lam_dot * dQdP;

      if (centroid)    printf("Pc = %.4e, Pc = %.4e\n", Pc, (P  - K * dt * lam_dot * dQdP));


    }
  }


  // ----------------- Reaction volume changes, computation of updated density only on centroid nodes
  if ( centroid > 0 ) {
    if (final_update==1) P = Pc;
    // if (model->chemical_production == 1 || Melting == 1) {
    //   *rho = EvaluateDensity( phase, T,  P,  X,  phi1, model, materials );                     
    // }
    if (model->density_variations==1) {
      *rho = EvaluateDensity( phase, T,  P,  X,  phi1, model, materials );  
    }
    else {
      rho_eq = EvaluateDensity( phase, T, P, X, phi1, model, materials );
      *rho   = rho_eq;
      if ( kinetics == 1 ) {
        rho0          = EvaluateDensity( phase, T0, P0, X, phi, model, materials );
        if ( materials->kin[phase]==9 ) {
          tau_kin = 0.0;
          dG      = Interpolate2Ddata( T0, P0, model->kin_Tmin, model->kin_Tmax, model->kin_Pmin, model->kin_Pmax, model->kin_nT, model->kin_nP, model->kin_dG );
          if ( dG>0.0 ) {
            Vkin          = kkin*T*( exp(-Qkin/R/T) * (1.0 - exp(-dG/R/T)) ); // growth rate [m/s]
            tau_kin       = log(1-0.6666666) * pow(-2*Skin*Vkin, -1);
          }
        }
        if ( materials->kin[phase]==9 || materials->kin[phase]==1 ) {
          *rho          = 1.0/(tau_kin + dt) * (tau_kin*rho0 + dt*rho_eq);
        }
      }
    }
  }

  if (is_pl == 0) {
    (*etaVE)    = eta_ve;
    (*div_pl)   = 0.0;
  }
  else {
    (*etaVE)    = eta_vep;
    (*div_pl)   = divp;
    eta_pl      = Tii/gdot;
  }

  /*----------------------------------------------------*/
  /*----------------------------------------------------*/
  /*----------------------------------------------------*/

  // Deviatoric stress
  *Txx       = 2.0*(*etaVE)*Exx;
  *Tzz       = 2.0*(*etaVE)*Ezz;
  *Txz       = 2.0*(*etaVE)*Exz;

  //-------- Post-Processing
  if ( post_process == 1) {

    // Strain rates: VEP partitioning
    if ( dislocation == 1 )              eta_pwl  = B_pwl * pow( *Eii_pwl, 1.0/n_pwl - 1.0 ); 
    if ( diffusion   == 1 )              eta_lin  = B_lin * pow( *Eii_lin, 1.0/n_lin - 1.0 ) * pow(*d1, m_lin/n_lin);
    if ( gbs         == 1 )              eta_gbs  = B_gbs * pow( *Eii_gbs, 1.0/n_gbs - 1.0 ); 
    if ( peierls     == 1 )              eta_exp  = B_exp * pow( *Eii_exp, 1.0 / (ST + n_exp) - 1.0); 

    // Compute strain rate invariants
    if (elastic) {
      const double Tyy    = -(*Txx  + *Tzz  );  
      const double Tyy0   = -( Txx0 +  Tzz0 );  
      const double Exx_el = (*Txx - Txx0)/2/eta_el;
      const double Ezz_el = (*Tzz - Tzz0)/2/eta_el;
      const double Eyy_el = ( Tyy - Tyy0)/2/eta_el;
      const double Exz_el = (*Txz - Txz0)/2/eta_el;
      *Eii_el   = sqrt(0.5*(Exx_el*Exx_el + Ezz_el*Ezz_el + Eyy_el*Eyy_el) + Exz_el*Exz_el);
      double Exx_pl = gdot * (*Txx)/2.0/Tii;
      double Exz_pl = gdot * (*Txz)/2.0/Tii;
      double div_pl = gdot * dQdP;
      double Exx_real = Exx - Txx0 / 2.0 / eta_el;
      double Exz_real = Exz - Txz0 / 2.0 / eta_el;
      if (centroid == 1 && gdot>0.0)
      {
        printf("-------------\n");
        printf("dQdTii = %.4e\n", dQdTii);
        printf("Exxel = %.4e, Exxpl = %.4e; Exx = %.4e; Exx_diff = %.4e\n", Exx_el, Exx_pl, Exx_real, Exx_real - Exx_el - Exx_pl);
        printf("Exzel = %.4e, Exzpl = %.4e; Exz = %.4e; Exz_diff = %.4e\n", Exz_el, Exz_pl, Exz_real, Exz_real - Exz_el - Exz_pl);
        printf("div_e + div_p = %.4e; div = %.4e; divp = %.4e; dive = %.4e; div_diff = %.4e\n", -(Pc - P0) / (K*dt) + gdot * sin_dil, div, gdot * sin_dil, -(Pc - P0) / (K*dt), div + (Pc - P0) / (K*dt) - gdot * sin_dil);
        printf("gdot = %.4e; sin_dil = %.4e\n", gdot, sin_dil);
        printf("Pc = %.4e; Pc1 = %.4e; P = %.4e; P0 = %.4e\n", Pc, P  - K * dt * lam_dot * dQdP, P , P0);
        printf("Model.step = %05d\n", model->step);
      }
    }
    if (is_pl && plastic) {
      *Eii_pl   = gdot/2.0;
    } else if (is_pl && tensile)
    {
      *Eii_pl   = gdot;
    } else
    {
      *Eii_pl   = 0.0;
    }

    // Partitioning of volumetric strain
    *div_pl  = divp;
    *div_el  = - (Pc - P0) / (K*dt);
    *div_r   = divr;

    double inv_eta_diss = 0.0;
    if (peierls    == 1)  inv_eta_diss += (1.0/eta_exp);
    if (dislocation== 1)  inv_eta_diss += (1.0/eta_pwl);
    if (diffusion  == 1)  inv_eta_diss += (1.0/eta_lin);
    if (constant   == 1)  inv_eta_diss += (1.0/eta_cst);
    if (is_pl      == 1)  inv_eta_diss += (1.0/eta_pl ); 
    eta  = 1.0/(inv_eta_diss);
    // if (T*scaling->T<673.0) {
      //  printf("constant = %d eta = %2.2e eta_pl = %2.2e %2.2e %2.2e %2.2e\n", constant, eta*scaling->eta, eta_pl*scaling->eta, eta_pwl*scaling->eta, eta_pl*scaling->eta, eta_cst*scaling->eta);
    // }

    // Viscoplastic overstress
    *OverS = eta_vp*gdot;

    // Work(s)
    if (centroid==1) {
      *Wtot  = Tii*Tii/eta_ve;
      *Wel   = Tii*Tii/eta_el;
      *Wdiss = Tii*Tii/eta;
      if (*Wdiss<0.0) {
       printf("Wtot = %2.6e Wdiss = %2.6e Wel = %2.6e\n", *Wtot, *Wdiss, *Wel);
       printf("Wtot - (Wdiss + Wel) = %2.6e\n", *Wtot - (*Wdiss + *Wel));
       printf("eta = %2.6e --- eta_ve = %2.6e --- eta_ve1 = %2.6e\n", eta, eta_ve, 1.0/(1./eta+1./eta_el));
       printf("constant = %d %2.2e %2.2e %2.2e %2.2e, Tii=%2.2e P=%2.2e gdot=%2.2e\n", constant, eta*scaling->eta, eta_pwl*scaling->eta, eta_pl*scaling->eta, eta_cst*scaling->eta, Tii*scaling->S, Pc*scaling->S, gdot*scaling->E);
       exit(1);
      }
    }
  }

  //-------- Post-Processing

  // Viscosity limiter
  if( *etaVE > max_eta ) {
    *etaVE = max_eta;
  }

  if( *etaVE < min_eta ) {
    *etaVE = min_eta;
  }

  // Viscosity limiter
  if( eta > max_eta ) {
    eta = max_eta;
  }

  if( eta < min_eta ) {
    eta = min_eta;
  }

  //    printf("e --> %2.2e %2.2e %2.2e\n", Exx, Ezz, Eyy); // if (fabs(Eyy)>1e-6)
  //    printf("T --> %2.2e %2.2e %2.2e\n", *Txx, *Tzz, -(*Txx)-(*Tzz)); // if (fabs( -(*Txx)-(*Tzz))>1e-6)

  return eta;

}

int solve_local_3x3(
    double J11, double J12, double J13,
    double J21, double J22, double J23,
    double J31, double J32, double J33,
    double R1,  double R2,  double R3,
    double *dTii, double *dP, double *dlam)
{
    // Compute determinant of J
    double detJ =
        J11 * (J22 * J33 - J23 * J32)
      - J12 * (J21 * J33 - J23 * J31)
      + J13 * (J21 * J32 - J22 * J31);

    // --- Singularity / robustness check
    const double det_tol = 1e-200;  // you can tune this based on your scaling
    if (fabs(detJ) < det_tol) {
        fprintf(stderr,
                "Warning: local 3x3 Jacobian nearly singular in plastic corrector, detJ = %.3e\n",
                detJ);
        // Option 1: return failure so caller can handle it (e.g. cut timestep, bail out, etc.)
        return 0;
    }

    double invdet = 1.0 / detJ;

    // dx = -inv(J) * R, using adjugate(J)

    *dTii = -invdet * (
          (J22 * J33 - J23 * J32) * R1 +
          (J13 * J32 - J12 * J33) * R2 +
          (J12 * J23 - J13 * J22) * R3
    );

    *dP   = -invdet * (
          (J23 * J31 - J21 * J33) * R1 +
          (J11 * J33 - J13 * J31) * R2 +
          (J13 * J21 - J11 * J23) * R3
    );

    *dlam = -invdet * (
          (J21 * J32 - J22 * J31) * R1 +
          (J12 * J31 - J11 * J32) * R2 +
          (J11 * J22 - J12 * J21) * R3
    );

    return 1;  // success
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGrid( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling, int final_update ) {

  int    p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1;
  double eta, txx1, tzz1, txz1, etaVE, VEcoeff = 0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
  double Wtot, Wel, Wdiss;
  int    average = model->eta_average;
  double Xreac;
  double OverS;
  double Pcorr, rho;
  double Exx, Ezz, Exz, eta_e;
  double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
  double el = 0.0;
  if (model->elastic == 1) el = 1.0;

  Nx  = mesh->Nx;
  Ncx = Nx - 1;
  Nz  = mesh->Nz;  Ncz = Nz-1;

  // Stuff to be interpolated to vertices
  InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
  InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

  // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Xreac,OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e, Wtot, Wel, Wdiss ) firstprivate( el, materials, scaling, average, model, Ncx, Ncz )
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
    mesh->X_n[c0]         = 0.0;
    mesh->OverS_n[c0]     = 0.0;

    if ( model->density_variations == 1 ) {
      mesh->rho_n[c0]  = 0.0;
    }

    // Loop on grid nodes
    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {

      //----------------------------------------------------------//
      if ( model->elastic==1 ) eta_e      = model->dt*mesh->mu_n[c0];
      else                       eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
      Exx = mesh->exxd[c0]  + el*mesh->sxxd0[c0] /eta_e/2.0;
      Ezz = mesh->ezzd[c0]  + el*mesh->szzd0[c0] /eta_e/2.0;
      Exz = mesh->exz_n[c0] + el*mesh->sxz0_n[c0]/eta_e/2.0;

      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_n[p][c0])>min_fraction ) is_phase_active = true;

        if ( is_phase_active==true ) {
          // printf("P0 before centers = %.4e\n", mesh->p0_n[c0]);
          eta =  ViscosityConcise( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], Exx, Ezz, Exz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], mesh->p0_n[c0], mesh->T0_n[c0], &Xreac, &OverS, &Pcorr, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 1, final_update, c0 );
          mesh->phase_eta_n[p][c0] = etaVE;

          switch ( average ) {
            case 0 :
              // ARITHMETIC AVERAGE
              mesh->eta_n[c0]       += mesh->phase_perc_n[p][c0] * etaVE;
              mesh->eta_phys_n[c0]  += mesh->phase_perc_n[p][c0] * eta;
              break;
            case 1 :
              // HARMONIC AVERAGE
              mesh->eta_n[c0]      += mesh->phase_perc_n[p][c0] * 1.0/etaVE;
              mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/eta;
              break;
            case 2 :
              // GEOMETRIC AVERAGE
              mesh->eta_n[c0]      += mesh->phase_perc_n[p][c0] * log(etaVE);
              mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * log(eta);
              break;
          }

      //  if isnan(mesh->eta_n[c0]) {
      //   printf("%2.2e %2.2e %2.2e %2.2e\n",mesh->eta_n[c0],Exx,mesh->sxxd[c0], etaVE);
      //   exit(1);
      // }

          mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
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

          if (mesh->Wdiss[c0]<0.0) {printf("negative dissipation: you crazy! --> Wdiss = %2.2e Wdiss = %2.2e\n", mesh->Wdiss[c0]*scaling->S*scaling->E, Wdiss*scaling->S*scaling->E); }

          mesh->p_corr[c0]      += mesh->phase_perc_n[p][c0] * Pcorr;
          mesh->div_u_el[c0]    += mesh->phase_perc_n[p][c0] * div_el;
          mesh->div_u_pl[c0]    += mesh->phase_perc_n[p][c0] * div_pl;
          mesh->div_u_r[c0]     += mesh->phase_perc_n[p][c0] * div_r;

// printf("%1.2e\n", Xreac);
          mesh->X_n[c0]         += mesh->phase_perc_n[p][c0] * Xreac;
          mesh->OverS_n[c0]     += mesh->phase_perc_n[p][c0] * OverS;

          // Volume changes
          if ( model->density_variations == 1 ) {
            mesh->rho_n[c0]       += mesh->phase_perc_n[p][c0] * rho;
          }
        }
      }

      mesh->d_n[c0]          = 1.0/mesh->d_n[c0];
      if (isinf(mesh->d_n[c0]) || isnan(mesh->d_n[c0]) ) {
        for ( p=0; p<model->Nb_phases; p++) printf("mesh->phase_perc_n[%d][%d] = %lf\n", p, c0, mesh->phase_perc_n[p][c0]);
        printf("Cell went empty!!! Exiting..."); exit(345);
      }
      // if ( model->density_variations == 1 ) mesh->rho_n[c0]        = 1.0/mesh->rho_n[c0];

      // HARMONIC AVERAGE
      if ( average == 1 ) {
        mesh->eta_n[c0]      = 1.0/mesh->eta_n[c0];
        mesh->eta_phys_n[c0] = 1.0/mesh->eta_phys_n[c0];

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
          printf("chemical_diffusion %d\n", model->chemical_diffusion);
          for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
          printf("eta=%2.2e G=%2.2e T=%2.2e P=%2.2e d=%2.2e phi=%2.2e %2.2e %2.2e %2.2e %2.2e\n", eta*scaling->eta, mesh->mu_n[c0]*scaling->S, mesh->T[c0]*scaling->T, mesh->p_in[c0]*scaling->S, mesh->d0_n[c0]*scaling->L, mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
          printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
          printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
          exit(1);
        }
      }
      // GEOMETRIC AVERAGE
      if ( average == 2 ) {
        mesh->eta_n[c0]      = exp(mesh->eta_n[c0]);
        mesh->eta_phys_n[c0] = exp(mesh->eta_phys_n[c0]);
      }

      // Final stress update
      mesh->sxxd[c0]  = 2.0*mesh->eta_n[c0]*Exx;
      mesh->szzd[c0]  = 2.0*mesh->eta_n[c0]*Ezz;
      mesh->sxz_n[c0] = 2.0*mesh->eta_n[c0]*Exz;

    }
  }

// Calculate vertices viscosity
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Wtot, Wel, Wdiss, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e ) firstprivate( el, materials, scaling, average, model, Nx, Nz )
  for ( k1=0; k1<Nx*Nz; k1++ ) {

    k  = mesh->kn[k1];
    l  = mesh->ln[k1];
    c1 = k + l*Nx;

    mesh->sxxd_s[c1]     = 0.0;
    mesh->szzd_s[c1]     = 0.0;
    mesh->sxz[c1]        = 0.0;
    mesh->eta_phys_s[c1] = 0.0;
    mesh->eta_s[c1]      = 0.0;
    mesh->X_s[c1]        = 0.0;

    if ( mesh->BCg.type[c1] != 30 ) {

      if ( model->elastic==1   ) eta_e      = model->dt*mesh->mu_s[c1];
      else           eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
      Exx = mesh->exxd_s[c1] + el*mesh->sxxd0_s[c1]/eta_e/2.0;
      Ezz = mesh->ezzd_s[c1] + el*mesh->szzd0_s[c1]/eta_e/2.0;
      Exz = mesh->exz[c1]    + el*mesh->sxz0[c1]   /eta_e/2.0;
      // printf("-----\n");
      // printf("%2.4e %2.4e  %2.4e \n", mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1]);
      // printf("%2.4e %2.4e  %2.4e \n", Exx, Ezz, Exz);

      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_s[p][c1])>min_fraction ) is_phase_active = true;

        if ( is_phase_active==true ) {

          // printf("P0 before vertices = %.4e\n", mesh->p0_s[c1]);

          eta =  ViscosityConcise( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], Exx, Ezz, Exz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], mesh->p0_s[c1], 0.0, &Xreac, &OverS, &Pcorr, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 0, final_update, c1 );
          mesh->phase_eta_s[p][c1] = etaVE;

          switch ( average ) {
            case 0 :
              mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * etaVE;
              mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * eta;
              break;
            case 1:
              mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * 1.0/etaVE;
              mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * 1.0/eta;
              break;
            case 2:
              mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * log(etaVE);
              mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * log(eta);
              break;
          }
          mesh->X_s[c1]        += mesh->phase_perc_s[p][c1] * Xreac;
        }
      }
      // HARMONIC AVERAGE
      if (average == 1) {
        mesh->eta_s[c1]      = 1.0/mesh->eta_s[c1];
        mesh->eta_phys_s[c1] = 1.0/mesh->eta_phys_s[c1];
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
        mesh->eta_s[c1]       = exp(mesh->eta_s[c1]);
        mesh->eta_phys_s[c1]  = exp(mesh->eta_phys_s[c1]);
      }

      // Final stress update
      mesh->sxxd_s[c1] = 2.0*mesh->eta_s[c1]*Exx;
      mesh->szzd_s[c1] = 2.0*mesh->eta_s[c1]*Ezz;
      mesh->sxz[c1]    = 2.0*mesh->eta_s[c1]*Exz;
    }
  }
  // printf("Txz:\n");
  //  Print2DArrayDouble( mesh->sxz,  mesh->Nx, mesh->Nz, scaling->S );
  //    printf("eta_s:\n");
  //  Print2DArrayDouble( mesh->eta_s,  mesh->Nx, mesh->Nz, scaling->eta );

  //  printf("eta_n:\n");
  //  Print2DArrayDouble( mesh->eta_n,  mesh->Nx-1, mesh->Nz-1, scaling->eta );
  //     printf("Txx:\n");
  //  Print2DArrayDouble( mesh->sxxd,  mesh->Nx-1, mesh->Nz-1, scaling->S );
  //    printf("Tzz:\n");
  //  Print2DArrayDouble( mesh->szzd,  mesh->Nx-1, mesh->Nz-1, scaling->S );
  // printf("Txz:\n");
  //  Print2DArrayDouble( mesh->sxz,  mesh->Nx, mesh->Nz, scaling->S );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Softening(int c0, double** phase_perc, double* dil_arr, double* fric_arr, double* C_arr, double strain_acc, params model, mat_prop materials, int style, int average ) {

  double fric, dil, C, Hc;
  double dfric, fric0, ddil, dil0, dcoh, C0, mu_strain, dstrain;

  // Loop on phases
  for ( int p=0; p<model.Nb_phases; p++) {

    fric = materials.phi[p];
    dil  = materials.psi[p];
    C    = materials.C[p];
    Hc   = materials.H_c[p];
    

    // Apply strain softening
    dstrain   = materials.pls_end[p] - materials.pls_start[p];

    if ( style == 1) {
      // Non linear error function softening

      mu_strain = 0.5*(materials.pls_end[p] + materials.pls_start[p]);

      if (materials.phi_soft[p] == 1) {
        dfric     = materials.phi[p]     - materials.phi_end[p];
        fric      = materials.phi[p]     - dfric/2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.phi_end[p];
        fric0     = materials.phi[p]     - dfric/2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.phi_end[p];
        fric      = fric * dfric / fric0 + materials.phi_end[p];
      }

      if (materials.psi_soft[p] == 1) {
        ddil      = materials.psi[p]     - materials.psi_end[p];
        dil       = materials.psi[p]     - ddil /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.psi_end[p];
        dil0      = materials.psi[p]     - ddil /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.psi_end[p];
        dil       = dil * ddil  / dil0  + materials.psi_end[p];
      }

      if (materials.coh_soft[p] == 1) {
        dcoh      = materials.C[p]       - materials.C_end[p];
        C         = materials.C[p]       - dcoh /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.C_end[p];
        C0        = materials.C[p]       - dcoh /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.C_end[p];
        C         = C * dcoh / C0 + materials.C_end[p];
      }
    }
    else if (style == 2)
    {
      C = fmax(materials.C[p] + Hc * strain_acc * materials.C_end[p], materials.C_end[p]);
    }
    

    else {
      // Piecewise linear function softening

      // If we are below the lower strain limit
      if (strain_acc < materials.pls_start[p]) {
        if (materials.phi_soft[p] == 1) fric = materials.phi[p];
        if (materials.psi_soft[p] == 1) dil  = materials.psi[p];
        if (materials.coh_soft[p] == 1) C    = materials.C[p];
      }
      // If we are above the upper strain limit
      if (strain_acc >= materials.pls_end[p]) {
        if (materials.phi_soft[p] == 1) fric = materials.phi_end[p];
        if (materials.psi_soft[p] == 1) dil  = materials.psi_end[p];
        if (materials.coh_soft[p] == 1) C    = materials.C_end[p];
      }
      // If we are in the softening strain range
      if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
        if (materials.phi_soft[p] == 1) fric = materials.phi[p] + (materials.phi_end[p] - materials.phi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
        if (materials.psi_soft[p] == 1) dil  = materials.psi[p] + (materials.psi_end[p] - materials.psi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
        if (materials.coh_soft[p] == 1) C    = materials.C[p]   + (materials.C_end[p]   - materials.C[p]  ) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
      }

    }

    // Arithmetic
    if (average ==0) {
      fric_arr[c0] += phase_perc[p][c0] * fric;
      dil_arr[c0]  += phase_perc[p][c0] * dil;
      C_arr[c0]    += phase_perc[p][c0] * C;
    }
    // Harmonic
    if (average == 1) {
      fric_arr[c0] += phase_perc[p][c0] *  1.0/fric;
      dil_arr[c0]  += phase_perc[p][c0] *  1.0/dil;
      C_arr[c0]    += phase_perc[p][c0] *  1.0/C;
    }
    // Geometric
    if (average == 2) {
      fric_arr[c0] += phase_perc[p][c0] *  log(fric);
      dil_arr[c0]  += phase_perc[p][c0] *  log(dil);
      C_arr[c0]    += phase_perc[p][c0] *  log(C);
    }
  }
  // Post-process for geometric/harmonic averages
  if ( average==1 ) fric_arr[c0] = 1.0/fric_arr[c0];
  if ( average==2 ) fric_arr[c0] = exp(fric_arr[c0]);
  if ( average==1 ) dil_arr[c0]  = 1.0/dil_arr[c0];
  if ( average==2 ) dil_arr[c0]  = exp(dil_arr[c0]);
  if ( average==1 ) C_arr[c0]    = 1.0/C_arr[c0];
  if ( average==2 ) C_arr[c0]    = exp(C_arr[c0]);

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CohesionFrictionDilationGrid( grid* mesh, markers* particles, mat_prop materials, params model, scale scaling ) {

  int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
  int average = 0;
  double *strain_pl, lap_strain_pl = 0., pl_str_W = 0., pl_str_E = 0., pl_str_S = 0., pl_str_N = 0.,  g= 5e2;
  int style = 0;
  int cent = 1, vert = 0, prop = 1, interp = 0;

  if (model.smooth_softening==1) style = 1;
  if (model.hardening_modulus==1) style = 2;
  

  Nx = mesh->Nx;
  Nz = mesh->Nz;
  Ncx = Nx-1;
  Ncz = Nz-1;

  // Plastic strain
  strain_pl  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
  P2Mastah( &model, *particles,  particles->strain_pl,     mesh, strain_pl,   mesh->BCp.type,  1, 0, interp, cent, model.interp_stencil);

  // Calculate cell centers cohesion and friction
#pragma omp parallel for shared( mesh, strain_pl ) private( p, c0 ) firstprivate( model, materials, average, style, Ncx, Ncz )
  for ( c0=0; c0<Ncx*Ncz; c0++ ) {

    // First - initialize to 0
    mesh->fric_n[c0] = 0.0;
    mesh->dil_n[c0]  = 0.0;
    mesh->C_n[c0]    = 0.0;

    // Compute only if below free surface
    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
      // Call softening routine
      Softening(c0, mesh->phase_perc_n, mesh->dil_n, mesh->fric_n, mesh->C_n, strain_pl[c0], model, materials, style, average );
      // Include random noise
      if (model.marker_noise == 1) mesh->C_n[c0] += mesh->noise_n[c0]*mesh->C_n[c0];

    //   if (mesh->kp[c0] == 0 ) pl_str_W = strain_pl[c0];
    //   else pl_str_W = strain_pl[c0-1];
    //   if (mesh->kp[Ncx-1] == 0 ) pl_str_E = strain_pl[c0];
    //   else pl_str_E = strain_pl[c0+1];
    //   if (mesh->lp[c0] == 0 ) pl_str_S = strain_pl[c0];
    //   else pl_str_S = strain_pl[c0-Ncx];
    //   if (mesh->lp[Ncz-1] == 0 ) pl_str_N = strain_pl[c0];
    //   else pl_str_N = strain_pl[c0+Ncx];
    // //  printf("%2.2e %2.2e\n", mesh->C_n[c0], g*(pl_str_W + pl_str_E + pl_str_S + pl_str_N - 4.0*strain_pl[c0]));
    //   mesh->C_n[c0] += g*( (pl_str_W + pl_str_E - 2.0*strain_pl[c0])/mesh->dx*mesh->dx + (pl_str_S + pl_str_N - 2.0*strain_pl[c0])/mesh->dz*mesh->dz );
    }
  }

  // Freedom
  DoodzFree( strain_pl );

  // Plastic strain
  strain_pl  = DoodzCalloc((model.Nx-0)*(model.Nz-0), sizeof(double));
  // Interp_P2N ( *particles,  particles->strain_pl, mesh, strain_pl, mesh->xg_coord, mesh->zg_coord, 1, 0, &model );
  P2Mastah( &model, *particles, particles->strain_pl, mesh, strain_pl, mesh->BCg.type,  1, 0, interp, vert, (&model)->interp_stencil);

#pragma omp parallel for shared( mesh, strain_pl ) private( p, c1 ) firstprivate( model, materials, average, style, Nx, Nz )
  // Calculate vertices cohesion and friction
  for ( c1=0; c1<Nx*Nz; c1++ ) {

    // First - initialize to 0
    mesh->fric_s[c1] = 0.0;
    mesh->dil_s[c1]  = 0.0;
    mesh->C_s[c1]    = 0.0;

    // Compute only if below free surface
    if ( mesh->BCg.type[c1] != 30 ) {
      // Call softening routine
      Softening(c1, mesh->phase_perc_s, mesh->dil_s, mesh->fric_s, mesh->C_s, strain_pl[c1], model, materials, style, average );
      // Include random noise
      if (model.marker_noise == 1) mesh->C_s[c1] += mesh->noise_s[c1]*mesh->C_s[c1];
      
      // if (mesh->kn[c0] == 0 ) pl_str_W = strain_pl[c1];
      // else pl_str_W = strain_pl[c1-1];
      // if (mesh->kn[Nx-1] == 0 ) pl_str_E = strain_pl[c1];
      // else pl_str_E = strain_pl[c1+1];
      // if (mesh->ln[c0] == 0 ) pl_str_S = strain_pl[c1];
      // else pl_str_S = strain_pl[c1-Nx];
      // if (mesh->ln[Nz-1] == 0 ) pl_str_N = strain_pl[c1];
      // else pl_str_N = strain_pl[c1+Nx];
      // mesh->C_s[c1] += g*( (pl_str_W + pl_str_E - 2.0*strain_pl[c1])/mesh->dx*mesh->dx + (pl_str_S + pl_str_N - 2.0*strain_pl[c1])/mesh->dz*mesh->dz );
    }
  }

  // Freedom
  DoodzFree( strain_pl );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ShearModCompExpGrid( grid* mesh, mat_prop *materials, params *model, scale scaling ) {

  int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
  int average = 1;//%model.eta_average; // SHOULD NOT BE ALLOWED TO BE ELSE THAN 1

  Nx = mesh->Nx;
  Nz = mesh->Nz;
  Ncx = Nx-1;
  Ncz = Nz-1;

  // Calculate cell centers shear modulus
  for ( l=0; l<Ncz; l++ ) {
    for ( k=0; k<Ncx; k++ ) {

      // Cell center index
      c0 = k  + l*(Ncx);

      // First - initialize to 0
      mesh->mu_n[c0]  = 0.0;
      mesh->bet_n[c0] = 0.0;
      mesh->alp[c0]   = 0.0;

      // Compute only if below free surface
      if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {

          // Arithmetic
          if (average == 0) {
            mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] * materials->G[p];
            mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] * materials->bet[p];
          }
          // Harmonic
          if (average == 1) {
            mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] *  1.0/materials->G[p];
            mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/materials->bet[p];
          }
          // Geometric
          if (average == 2) {
            mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] *  log(materials->G[p]);
            mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] *  log(materials->bet[p]);
          }

          // Standard arithmetic interpolation
          mesh->alp  [c0] += mesh->phase_perc_n[p][c0] * materials->alp[p];

        }
        // Post-process for geometric/harmonic averages
        if ( average==1 ) mesh->mu_n[c0]  = 1.0/mesh->mu_n[c0];
        if ( average==2 ) mesh->mu_n[c0]  = exp(mesh->mu_n[c0]);
        if ( average==1 ) mesh->bet_n[c0] = 1.0/mesh->bet_n[c0];
        if ( average==2 ) mesh->bet_n[c0] = exp(mesh->bet_n[c0]);
      }
      if (isinf( mesh->mu_n[c0])) {
        printf("%f x = %2.6e z = %2.6e - BCp %d BCt %d\n", mesh->mu_n[c0], mesh->xc_coord[k]*scaling.L, mesh->zc_coord[l]*scaling.L, mesh->BCp.type[c0], mesh->BCt.type[c0]); exit(32);
      }
    }
  }


  // Calculate vertices shear modulus
  for ( l=0; l<Nz; l++ ) {
    for ( k=0; k<Nx; k++ ) {

      // Vertex index
      c1 = k + l*Nx;

      // First - initialize to 0
      mesh->mu_s[c1]  = 0.0;
      mesh->bet_s[c1] = 0.0;

      // Compute only if below free surface
      if ( mesh->BCg.type[c1] != 30 ) {

        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {

          // Arithmetic
          if (average == 0) {
            mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] * materials->G[p];
            mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] * materials->bet[p];
          }
          // Harmonic
          if (average == 1) {
            mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] *  1.0/materials->G[p];
            mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/materials->bet[p];
          }
          // Geometric
          if (average == 2) {
            mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] *  log(materials->G[p]);
            mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] *  log(materials->bet[p]);
          }

        }

        if ( isinf(1.0/mesh->mu_s[c1]) ) {
          printf("Aaaaargh...!! %2.2e %2.2e ----> ShearModulusCompressibilityExpansivityGrid\n", mesh->phase_perc_s[0][c1], mesh->phase_perc_s[1][c1]);
        }

        // Post-process for geometric/harmonic averages
        if ( average==1 ) mesh->mu_s[c1]  = 1.0/mesh->mu_s[c1];
        if ( average==2 ) mesh->mu_s[c1]  = exp(mesh->mu_s[c1]);
        if ( average==1 ) mesh->bet_s[c1] = 1.0/mesh->bet_s[c1];
        if ( average==2 ) mesh->bet_s[c1] = exp(mesh->bet_s[c1]);
      }
    }
  }

  // Periodic
  double av;
  if (model->periodic_x==1) {
    for( l=0; l<Nz; l++) {
      c1 = l*Nx + Nx-1;
      av = 0.5*(mesh->mu_s[c1] + mesh->mu_s[l*Nx]);
      mesh->mu_s[c1] = av; mesh->mu_s[l*Nx] = av;
      av = 0.5*(mesh->bet_s[c1] + mesh->bet_s[l*Nx]);
      mesh->bet_s[c1] = av; mesh->bet_s[l*Nx] = av;      
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double EvaluateDensity( int p, double T, double P, double X, double phi, params *model, mat_prop *materials ) {

  double rho;
  
  // Constant density
  if ( materials->density_model[p] == 0 ) {
    const double rho_ref = materials->rho[p];
    rho     = rho_ref;
  }

  // T, P, X dependent density based on EOS
  if ( materials->density_model[p] == 1 ) {
    const double rho_ref = materials->rho[p];
    const double drho    = materials->drho[p];
    const double T0      = materials->T0 [p];
    const double alpha   = materials->alp[p];
    const double P0      = materials->P0 [p];
    const double beta    = materials->bet[p];
    rho     = (1.0 -  alpha * (T - T0) ) * (1.0 +  beta * (P - P0) ); // EOS general
    rho     = ((1.0-X)*rho_ref + X*(rho_ref+drho))*rho;                     // Average density based on X
  }

  // T and P dependent density based on phase diagrams
  if ( materials->density_model[p] == 2 ) {
    int PD = materials->phase_diagram[p]; // PD: index of phase diagram
    rho    = Interpolate2Ddata( T, P, model->PDMTmin[PD], model->PDMTmax[PD], model->PDMPmin[PD], model->PDMPmax[PD], model->PDMnT[PD], model->PDMnP[PD], model->PDMrho[PD] );
  }

  // P-T dependent density
  if ( materials->density_model[p] == 3 && model->density_variations == 0 ) {
    const double rho_ref = materials->rho[p];
    const double beta    = materials->bet[p];
    const double alpha   = materials->alp[p];
    rho     = rho_ref*exp(beta*P  - alpha*T);
  }

  // P-T dependent density: models used for Yamato et al. (2022)
  if ( materials->density_model[p] == 3 && model->density_variations == 1) {
    const double rho1    = materials->rho[p];
    const double rho2    = materials->rho[materials->reac_phase[p]];
    const double beta    = materials->bet[p];
    const double alpha   = materials->alp[p];
    const double rho_ref = (1.0-X)*rho1 + X*rho2;
    rho     = rho_ref * exp(beta*P - alpha*T);
  }

  // P-T dependent density: models used for Poh et al. (2021)
  if ( materials->density_model[p] == 5  ) {
    const double rho_sol = materials->rho[p];
    const double rho_liq = materials->rho[p] + materials->drho[p]; // !!!! Achtung: drho should be NEGATIVE to make liquid lighter 
    const double beta    = materials->bet[p];
    const double alpha   = materials->alp[p];
    const double rhos    = rho_sol * exp(beta*P - alpha*T);
    const double rhol    = rho_liq * exp(beta*P - alpha*T);
    // rho     = rhos * (1 - phi + phi*rho_liq/rho_sol );
    rho     = (1-phi)*rhos + phi*rhol;
    // rho = rhos;
    // printf("rho_sol = %2.6e rho_liq = %2.6e drho = %2.6e\n", rho_sol, rho_liq, materials->drho[p]);
  }

  // P dependent density read from the 1D table
  if ( materials->density_model[p] == 4 ) {
    rho     = ItpRho1D( P, model, p );
  }
  return rho;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
void UpdateDensity( grid* mesh, markers* particles, mat_prop *materials, params *model, scale *scaling ) {

  int k, p, c0, Ncx=mesh->Nx-1, Ncz=mesh->Nz-1;
  int    old = 0;
  double rho, rho0, epsi = 1e-13, mean_rho, mean_rho0;
  // printf("Update density fields on mesh\n");

  // mean_rho = SumArray( mesh->rho_n, scaling->rho, Ncx*Ncz, "mean rho" ) /(Ncx*Ncz);

#pragma omp parallel for shared( mesh, materials ) private( rho, rho0, c0, p) firstprivate(old, Ncx, Ncz, model, epsi)
  for ( c0=0; c0<Ncx*Ncz; c0++ ) {

    // Initialise
    mesh->rho_n[c0]  = 0.0;
    if (old==1) mesh->rho0_n[c0] = 0.0;

    // Loop on phases
    for ( p=0; p<model->Nb_phases; p++) {

      if ( fabs(mesh->phase_perc_n[p][c0])>epsi) {
        // Call density evaluation
        rho              = EvaluateDensity( p, mesh->T[c0],    mesh->p_in[c0], mesh->X_n[c0],  mesh->phi_n[c0],  model, materials );
        if (old==1) rho0 = EvaluateDensity( p, mesh->T0_n[c0], mesh->p0_n[c0], mesh->X0_n[c0], mesh->phi0_n[c0], model, materials );

        // Average density base on phase density and phase volume fraction
        if ( mesh->BCp.type[c0] != 30 )             mesh->rho_n[c0]  += mesh->phase_perc_n[p][c0] * rho;
        if (old==1) if ( mesh->BCp.type[c0] != 30 ) mesh->rho0_n[c0] += mesh->phase_perc_n[p][c0] * rho0;
        // if ( mesh->BCp.type[c0] != 30 )             mesh->rho_n[c0]  += mesh->phase_perc_n[p][c0] * 1./rho;
        // if (old==1) if ( mesh->BCp.type[c0] != 30 ) mesh->rho0_n[c0] += mesh->phase_perc_n[p][c0] * 1./rho0;

        // if (mesh->BCp.type[c0] == -1 && mesh->BCp.type[c0+Ncx] == 31 ) {
        //   printf("At the surface: rho = %2.2lf, phase = %d rho = %2.2lf X = %2.2lf density_model = %d\n", mesh->rho_n[c0]*scaling->rho, p, rho, mesh->phase_perc_n[p][c0], materials->density_model[p]);
        // }
      }
    }
    // mesh->rho_n[c0]              = 1.0/mesh->rho_n[c0]; 
    // if (old==1) mesh->rho0_n[c0] = 1.0/mesh->rho0_n[c0];
  }
  
  // Interpolate to vertices
  InterpCentroidsToVerticesDouble( mesh->rho_n, mesh->rho_s, mesh, model );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void HuetAveragingModel( double *B_pwl, double *C_pwl, double *n_pwl1, int phase, double R, double T, int t_pwl, double X, double pre_factor, mat_prop* materials ) {
  // Parameters of end-members
  const double ndis1  = materials->npwl[phase];   const double ndis2  = materials->npwl[materials->reac_phase[phase]];
  const double Adis1  = materials->Apwl[phase];   const double Adis2  = materials->Apwl[materials->reac_phase[phase]];
  const double Qdis1  = materials->Qpwl[phase];   const double Qdis2  = materials->Qpwl[materials->reac_phase[phase]];
  double sum_up, sum_down;
  // Huet et al 2014 ---------------
  const double f1 = 1.0-X;
  const double f2 = X;

  // (1) Calcul des ai
  const double a1 = ndis2 + 1.0;
  const double a2 = ndis1 + 1.0;

  // (2) Calcul de n bulk:
  sum_up   = f1*a1*ndis1 + f2*a2*ndis2;
  sum_down = f1*a1+f2*a2;
  const double n_pwl     = sum_up/sum_down;

  // (3) Calcul de Q bulk:
  sum_up    = f1*a1*Qdis1 + f2*a2*Qdis2;
  sum_down  = f1*a1+f2*a2;
  const double Ea_pwl    = sum_up/sum_down;

  // (4) Calcul de A bulk:
  const double Prod1 = pow(Adis1,(f1*a1/sum_down)) * pow(Adis2,(f2*a2/sum_down));
  const double sum_n = f1*ndis1/(ndis1 + 1.0) + f2*ndis2/(ndis2 + 1.0);
  const double Prod2 = pow(ndis1/(ndis1 + 1.0),f1*a1*ndis1/sum_down) * pow(ndis2/(ndis2+1.0),f2*a2*ndis2/sum_down);
  const double A_pwl = Prod1 * pow(sum_n,-n_pwl) * Prod2;

  // Proper choices of corrections factors
  double F_pwl;
  if ( (int)t_pwl == 0 ) {
    F_pwl = 1.0;
  }
  if ( (int)t_pwl == 1 ) {
    F_pwl = 1.0/6.0*pow(2.0,1.0/n_pwl) * pow(3.0,(n_pwl-1.0)/2.0/n_pwl);
  }
  if ( (int)t_pwl == 2 ) {
    F_pwl = 1.0/4.0*pow(2,1.0/n_pwl);
  }

  // Override power-law flow law parameters
  *B_pwl  = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl)/R/n_pwl/T );
  *C_pwl  = pow(2.0*( (*B_pwl) ), -n_pwl);
  *n_pwl1 = n_pwl;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ItpRho1D( double Pgrid, params* model, int k ) {
  // Declarations
  const    int NP = model->PD1DnP[k];
  const double dP = ( model->PD1Dmax[k] - model->PD1Dmin[k]) / (NP -1);
  if (Pgrid<model->PD1Dmin[k]) Pgrid = model->PD1Dmin[k] + 0.01*dP;
  if (Pgrid>model->PD1Dmax[k]) Pgrid = model->PD1Dmax[k] - 0.01*dP;
  // Find index of minimum/west pressure node
  double dstP = Pgrid - model->PD1Dmin[k];
  int    iP   = ceil( dstP/dP ) - 1;
  // Calculate weighting coefficients for linear interpolant
  double PW   = model->PD1Dmin[k] + iP*dP;
  double wW   = 1.0 - (Pgrid - PW)/dP;
  // Interpolate from 2 neighbours
  double rho   = wW * model->PD1Drho[k][iP] + (1.0-wW) * model->PD1Drho[k][iP+1];
  return rho;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double Interpolate2Ddata( double Tgrid, double Pgrid, double Tmin, double Tmax, double Pmin, double Pmax, int nT, int nP, double* data ) {
  // Declarations
  double data_itp, dstT, dstP, TW, PW;
  double wT, wP;
  int iT, iP, iSW, iSE, iNW, iNE;
  // Determine T and P increments in the current phase diagram
  const double dT = (Tmax - Tmin)/(nT-1.0);
  const double dP = (Pmax - Pmin)/(nP-1.0);
  // Pressure and temperature + correct to remain within the database bounds
  if ( Tgrid < Tmin ) Tgrid = Tmin + 0.01*dT;
  if ( Tgrid > Tmax ) Tgrid = Tmax - 0.01*dT;
  if ( Pgrid < Pmin ) Pgrid = Pmin + 0.01*dP;
  if ( Pgrid > Pmax ) Pgrid = Pmax - 0.01*dP;
  // Find index of minimum/west temperature node
  dstT  = Tgrid - Tmin ;
  iT    = ceil( dstT/dT ) - 1;
  // Find index of minimum/west pressure node
  dstP  = Pgrid - Pmin;
  iP    = ceil( dstP/dP ) - 1;
  // Calculate weights for bilinear interpolant
  TW    = Tmin + iT*dT;
  PW    = Pmin + iP*dP;
  wT    = 1.0 - (Tgrid - TW )/dT;
  wP    = 1.0 - (Pgrid - PW )/dP;
  // Indices of neigbours
  iSW   = iT + iP*nT;
  iSE   = iT + iP*nT+1;
  iNW   = iT + (iP+1)*nT;
  iNE   = iT + (iP+1)*nT+1;
  // Interpolate from 4 neighbours
  data_itp  = 0.0;
  data_itp +=  (1.0 - wT)* (1.0 - wP) * data[iNE];
  data_itp +=  (      wT)* (1.0 - wP) * data[iNW];
  data_itp +=  (1.0 - wT)* (      wP) * data[iSE];
  data_itp +=  (      wT)* (      wP) * data[iSW];
  return data_itp;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Strain rate
void StrainRateComponents( grid* mesh, scale scaling, params* model ) {

  int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
  double dx, dz;//, Eyyt = 0.0*model->bkg_div_rate/3.0;
  double dvxdx, dvydy, dvzdz;
  double out_of_plane = 0.0;
  if (model->out_of_plane==1) out_of_plane = 1.0;

  Nx = mesh->Nx;
  Nz = mesh->Nz;
  Ncx = Nx-1;
  Ncz = Nz-1;
  dx = mesh->dx;
  dz = mesh->dz;

#pragma omp parallel for shared( mesh ) private( k, k1, l, c0, c1, c2, dvxdx, dvydy, dvzdz  ) firstprivate( dx, dz, Nx, Ncx, Ncz, out_of_plane )
  for ( k1=0; k1<Ncx*Ncz; k1++ ) {
    k  = mesh->kp[k1];
    l  = mesh->lp[k1];
    c0 = k  + l*(Nx-1);
    c1 = k  + l*(Nx);
    c2 = k  + l*(Nx+1);

    mesh->div_u[c0] = 0.0;
    mesh->exxd[c0]  = 0.0;
    mesh->ezzd[c0]  = 0.0;

    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

      // Velocity gradients - Total normal strain rates
      dvxdx = (mesh->u_in[c1+1+Nx]   - mesh->u_in[c1+Nx])/dx;
      dvzdz = (mesh->v_in[c2+Nx+1+1] - mesh->v_in[c2+1] )/dz;
      dvydy = out_of_plane*0.5*( dvxdx + dvzdz );

      // Velocity divergence
      mesh->div_u[c0] = dvxdx + dvzdz + dvydy;

      // Normal strain rates
      mesh->exxd[c0]  = dvxdx - model->compressible*1.0/3.0*mesh->div_u[c0]; // SURE? 
      mesh->ezzd[c0]  = dvzdz - model->compressible*1.0/3.0*mesh->div_u[c0]; // SURE?
    }
  }

  // Shear components we only calculate sxz because sxz = szx (stress tensor symmetry)
#pragma omp parallel for shared( mesh ) private( k, k1, l, c1, c2  ) firstprivate( dx, dz, Nx, Nz )
  for ( k1=0; k1<Nx*Nz; k1++ ) {
    k  = mesh->kn[k1];
    l  = mesh->ln[k1];
    c1 = k  + l*(Nx);
    c2 = k  + l*(Nx+1);

    mesh->exz[c1] = 0.0;

    if ( mesh->BCg.type[c1] != 30 ) {
      //            // Original implementation from Duretz et al., 2016
      //            if (mesh->BCu.type[c1] != 30 && mesh->BCu.type[c1+Nx] != 30) mesh->exz[c1] +=  0.5 * (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
      //            if (mesh->BCv.type[c2] != 30 && mesh->BCv.type[c2+1]  != 30) mesh->exz[c1] +=  0.5 * (mesh->v_in[c2+1]  - mesh->v_in[c2])/dx;
      // Simplified implementation
      const double dvxdz = (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
      const double dvzdx = (mesh->v_in[c2+1]  - mesh->v_in[c2])/dx;
      mesh->exz[c1] =  0.5 * (dvxdz + dvzdx);
      mesh->wxz[c1] =  0.5 * (dvzdx - dvxdz);
    }
  }

// TODO: DELETE 
#pragma omp parallel for shared( mesh ) private( k1 ) firstprivate( Nx, Ncx, Ncz )
  for ( k1=0; k1<Ncx*Ncz; k1++ ) {

    int k  = mesh->kp[k1];
    int l  = mesh->lp[k1];
    int c0 = k  + l*(Nx-1);
    int c1 = k  + l*(Nx);

    mesh->exz_n[c0] = 0.0;
    mesh->wxz_n[c0] = 0.0;

    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {
      if (mesh->BCg.type[c1]      != 30 ) { mesh->exz_n[c0] += 0.25*mesh->exz[c1];      }
      if (mesh->BCg.type[c1+1]    != 30 ) { mesh->exz_n[c0] += 0.25*mesh->exz[c1+1];    }
      if (mesh->BCg.type[c1+Nx]   != 30 ) { mesh->exz_n[c0] += 0.25*mesh->exz[c1+Nx];   }
      if (mesh->BCg.type[c1+Nx+1] != 30 ) { mesh->exz_n[c0] += 0.25*mesh->exz[c1+Nx+1]; }
      if (mesh->BCg.type[c1]      != 30 ) { mesh->wxz_n[c0] += 0.25*mesh->wxz[c1];      }
      if (mesh->BCg.type[c1+1]    != 30 ) { mesh->wxz_n[c0] += 0.25*mesh->wxz[c1+1];    }
      if (mesh->BCg.type[c1+Nx]   != 30 ) { mesh->wxz_n[c0] += 0.25*mesh->wxz[c1+Nx];   }
      if (mesh->BCg.type[c1+Nx+1] != 30 ) { mesh->wxz_n[c0] += 0.25*mesh->wxz[c1+Nx+1]; }
    }
  }

  // Interpolate normal strain rate on vertices
  InterpCentroidsToVerticesDouble( mesh->exxd,  mesh->exxd_s,  mesh, model );
  InterpCentroidsToVerticesDouble( mesh->ezzd,  mesh->ezzd_s,  mesh, model );
  InterpCentroidsToVerticesDouble( mesh->div_u, mesh->div_u_s, mesh, model );

  // TODO: DELETE 

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GenerateDeformationMaps( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {

  // This functions generates deformation maps and writem to disk

  // Definition of parameters and allocation of memory
  int nT = model->nT, nE = model->nE, nd = model->nd, ix, iy, iz;
  double stepT, stepE, stepd;
  double Tmin = model->Tmin;
  double Tmax = model->Tmax;
  double Emin = model->Emin;
  double Emax = model->Emax;
  double dmin = model->dmin;
  double dmax = model->dmax;

  double *T, *E, *d, *dlog, *Elog, gs_ref;
  double txx1, tzz1, txz1, etaVE, VEcoeff, eII_el=0.0, eII_pl=0.0, eII_pwl=0.0, eII_exp=0.0, eII_lin=0.0, eII_gbs=0.0, eII_cst=0.0, d1;
  double Wtot, Wel, Wdiss, div_el, div_pl, div_r;
  double Pn = model->Pn , eta;
  int    *dom_mech, ind, mech, loud=0, k;
  double *stress, *visco, t_omp;
  double Xreac;
  double OverS;
  double Pcorr, rho;

  T    = malloc(nT*sizeof(double));
  Elog = malloc(nE*sizeof(double));
  E    = malloc(nE*sizeof(double));
  dlog = malloc(nd*sizeof(double));
  d    = malloc(nd*sizeof(double));

  //Initialization of arrays
  stepT = (Tmax-Tmin)/(nT-1);
  stepE = (Emax-Emin)/(nE-1);
  stepd = (dmax-dmin)/(nd-1);

  // Temperature vector
  for ( ix=0; ix<nT; ix++) {
    if (ix==0) T[ix] = Tmin;
    else       T[ix] = T[ix-1] + stepT;
  }
  // Strain rate vector
  for ( iy=0; iy<nE; iy++) {
    if (iy==0) Elog[iy] = Emin;
    else       Elog[iy] = Elog[iy-1] + stepE;
    E[iy] = pow(10.0,Elog[iy])/scaling->E;
  }
  // Grain size vector
  for (iz=0; iz<nd; iz++) {
    if (iz==0) dlog[iz] = dmin;
    else       dlog[iz] = dlog[iz-1] + stepd;
    d[iz] = pow(10.0,dlog[iz])/scaling->L;
  }

  double Flin, Fgbs, Fpwl, texp;

  // Loop on all phases
  for (k=0;k<model->Nb_phases; k++) {

    printf("Generating deformation map of phase %2d\n", k);

    // Save real values
    Flin   = materials->Flin[k];
    Fgbs   = materials->Fgbs[k];
    Fpwl   = materials->Fpwl[k];
    texp   = materials->texp[k];
    gs_ref = materials->gs_ref[k];

    // Set no correction for deformation maps
    materials->Flin[k]   = 1.0;
    materials->Fgbs[k]   = 1.0;
    materials->Fpwl[k]   = 1.0;
    materials->texp[k]   = 0.0;

    // Allocate maps
    dom_mech = malloc(nT*nE*nd*sizeof(int));
    stress   = malloc(nT*nE*nd*sizeof(double));
    visco    = malloc(nT*nE*nd*sizeof(double));

    t_omp = (double)omp_get_wtime();

    // Boucles spatiales 2D pour créations des carte de déformation

    for ( iz=0; iz<nd; iz++) {

      // Force grain size
      materials->gs_ref[k] = d[iz];

      for ( ix=0; ix<nT; ix++) {
        for ( iy=0; iy<nE; iy++) {

          // Evaluate viscosity and stress
          eta =  ViscosityConcise( k, 0.0, T[ix], Pn, d[iz], 0.0, 0.0, E[iy], E[iy], 0.0, 0.0, 0.0, 0.0, materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst,  &d1, 0.0, materials->psi[k], materials->phi[k], materials->C[k], 0.0, 0.0, &Xreac, &OverS, &Pcorr, &rho, 0.0, 0.0, &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 0, 0, 0, 0 );

          // Select mechanism
          if (eII_pwl>eII_el && eII_pwl>eII_pl  && eII_pwl>eII_exp && eII_pwl>eII_lin && eII_pwl>eII_gbs && eII_pwl>eII_cst) mech = 1; // dislocation
          if (eII_lin>eII_el && eII_lin>eII_pl  && eII_lin>eII_exp && eII_lin>eII_pwl && eII_lin>eII_gbs && eII_lin>eII_cst) mech = 2; // diffusion
          if (eII_gbs>eII_el && eII_gbs>eII_pl  && eII_gbs>eII_lin && eII_gbs>eII_pwl && eII_gbs>eII_exp && eII_gbs>eII_cst) mech = 3; // gbs
          if (eII_exp>eII_el && eII_exp>eII_pl  && eII_exp>eII_lin && eII_exp>eII_pwl && eII_exp>eII_gbs && eII_exp>eII_cst) mech = 4; // peierls
          if (eII_pl >eII_el && eII_pl >eII_pwl && eII_pl >eII_exp && eII_pl >eII_lin && eII_pl >eII_gbs && eII_pl >eII_cst) mech = 5; // plastic
          if (eII_cst>eII_pl && eII_cst>eII_pwl && eII_cst>eII_exp && eII_cst>eII_lin && eII_cst>eII_gbs && eII_cst>eII_el ) mech = 6; // constant viscosity
          if (eII_el >eII_pl && eII_el >eII_pwl && eII_el >eII_exp && eII_el >eII_lin && eII_el >eII_gbs && eII_el >eII_cst) mech = 7; // elastic

          // Global index for 3D matrix flattened in a 1D vector
          ind           = ix + iy*nT + iz*nE*nT;
          dom_mech[ind] = mech;
          stress[ind]   = txx1;
          visco[ind]    = eta;
        }
      }
    }

    printf("** Deformation map %2d ---> %lf sec\n", k, (double)((double)omp_get_wtime() - t_omp));


    // Print deformation maps to screen
    if (loud ==1) {
      for ( iz=0; iz<nd; iz++) {
        printf("GS = %2.2e\n", d[iz]*scaling->L);
        for ( ix=0; ix<nT; ix++) {
          if (ix==0) printf("           ");
          printf("%2.2lf ", T[ix]*scaling->T);
        }
        printf("\n");
        for ( iy=0; iy<nE; iy++) {
          printf("%2.2e   ", E[iy]*scaling->E);
          for ( ix=0; ix<nT; ix++) {
            ind = ix + iy*nT + iz*nd*nT;

            printf("%6d ", dom_mech[ind] );

          }
          printf("\n");
        }
        printf("\n");
      }
    }

    // Print deformation maps to file

    char *filename;
    asprintf( &filename, "DefMap%02d.gzip.h5", k );
    double *CT, *Cd, *CE, *Cstress, *Cvisco;

    CT = DoodzMalloc( sizeof(double)*nT);
    ArrayEqualArray( CT, T, nT );
    //        DoubleToFloat( T, CT, nT );
    ScaleBackD( CT, scaling->T, nT );

    CE = DoodzMalloc( sizeof(double)*nE);
    ArrayEqualArray( CE, E, nE );
    //        DoubleToFloat( E, CE, nE );
    ScaleBackD( CE, scaling->E, nE );

    Cd = DoodzMalloc( sizeof(double)*nd);
    ArrayEqualArray( Cd, d, nd );
    //        DoubleToFloat( d, Cd, nd );
    ScaleBackD( Cd, scaling->L, nd );

    Cstress = DoodzMalloc( sizeof(double)*nT*nE*nd);
    ArrayEqualArray( Cstress, stress, nT*nE*nd );
    //        DoubleToFloat( stress, Cstress, nT*nE*nd );
    ScaleBackD( Cstress, scaling->S, nT*nE*nd );

    Cvisco= DoodzMalloc( sizeof(double)*nT*nE*nd);
    ArrayEqualArray( Cvisco, visco, nT*nE*nd );
    //        DoubleToFloat( visco, Cvisco, nT*nE*nd );
    ScaleBackD( Cvisco, scaling->eta, nT*nE*nd );

    // Fill in DD data structure
    double params[3];
    params[0] = nT;
    params[1] = nE;
    params[2] = nd;

    // Send data to file
    CreateOutputHDF5( filename );
    AddGroupToHDF5( filename, "model" );
    AddGroupToHDF5( filename, "arrays" );

    AddFieldToGroup( filename, "model", "params" , 'd',  3, params,  1 );
    AddFieldToGroup( filename, "arrays", "T"     , 'd', nT, CT,  1 );
    AddFieldToGroup( filename, "arrays", "E"     , 'd', nE, CE,  1 );
    AddFieldToGroup( filename, "arrays", "d"     , 'd', nd, Cd,  1 );
    AddFieldToGroup( filename, "arrays", "stress", 'd', nT*nE*nd, Cstress,  1 );
    AddFieldToGroup( filename, "arrays", "visco" , 'd', nT*nE*nd, Cvisco,  1 );
    AddFieldToGroup( filename, "arrays", "map"   , 'i', nT*nE*nd, dom_mech,  1 );


    free(filename);
    DoodzFree(CT);
    DoodzFree(CE);
    DoodzFree(Cd);
    DoodzFree(Cstress);
    DoodzFree(Cvisco);

    // Free memory 1
    free(dom_mech);
    free(stress);
    free(visco);

    // Set no correction for deformation maps
    materials->Flin[k]   = Flin;
    materials->Fgbs[k]   = Fgbs;
    materials->Fpwl[k]   = Fpwl;
    materials->texp[k]   = texp;
    materials->gs_ref[k] = gs_ref;
  }
  // Free memory 2
  free(T);
  free(E);
  free(d);
  free(dlog);
  free(Elog);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void LogTimeSeries( grid* mesh, params model, scale scaling ) {

  mesh->Uthermal_time[model.step]  = mesh->Uthermal*scaling.S*scaling.L*scaling.L;
  mesh->Uelastic_time[model.step]  = mesh->Uelastic*scaling.S*scaling.L*scaling.L;
  mesh->Work_time[model.step]      = mesh->Work*scaling.S*scaling.L*scaling.L;
  mesh->Time_time[model.step]      = model.time*scaling.t;
  mesh->Short_time[model.step]     = (model.L0-(model.xmax - model.xmin))/model.L0 * 100.0;
  mesh->P_mean_time[model.step]    = mesh->P_mean*scaling.S;
  mesh->T_mean_time[model.step]    = mesh->T_mean*scaling.T;
  mesh->sxxd_mean_time[model.step] = mesh->sxxd_mean*scaling.S;
  mesh->szzd_mean_time[model.step] = mesh->szzd_mean*scaling.S;
  mesh->sxz_mean_time[model.step]  = mesh->sxz_mean*scaling.S;
  mesh->Tii_mean_time[model.step]  = mesh->Tii_mean*scaling.S;
  mesh->exxd_mean_time[model.step] = mesh->exxd_mean*scaling.E;
  mesh->ezzd_mean_time[model.step] = mesh->ezzd_mean*scaling.E;
  mesh->exz_mean_time[model.step]  = mesh->exz_mean*scaling.E;
  mesh->Eii_mean_time[model.step]  = mesh->Eii_mean*scaling.E;

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ComputeMeanQuantitesForTimeSeries( grid *mesh ) {

  int k1, k, l, cell, N=0;
  double T_mean=0.0, P_mean=0.0;
  double exxd_mean=0.0, ezzd_mean=0.0, exz_mean=0.0, Eii_mean=0.0;
  double sxxd_mean=0.0, szzd_mean=0.0, sxz_mean=0.0, Tii_mean=0.0;
  double Txzc, Exzc, Tiic, Eiic; // Interpolated from vertices to cell centers
  int vertSW, vertSE, vertNW, vertNE;
  int Nvx = mesh->Nx;
  int Nvz = mesh->Nz;
  int Ncx = mesh->Nx-1;
  int Ncz = mesh->Nz-1;

#pragma omp parallel for shared( mesh ) firstprivate( Nvx,Ncx,Ncz ) private( k1,k,l,cell,vertSW,vertSE,vertNW,vertNE,Exzc,Txzc,Eiic,Tiic) reduction( +:sxxd_mean,szzd_mean,sxz_mean,Tii_mean,exxd_mean,ezzd_mean,exz_mean,Eii_mean,T_mean,P_mean,N )
  for ( k1=0; k1<Ncx*Ncz; k1++ ) {

    k      = mesh->kp[k1];
    l      = mesh->lp[k1];
    cell   = k  + l*(Ncx);
    vertSW = k  + l*(Nvx);
    vertSE = k  + l*(Nvx) + 1;
    vertNW = k  + (l+1)*(Nvx);
    vertNE = k  + (l+1)*(Nvx) + 1;

    // If cell is active (below free surface)
    if ( mesh->BCp.type[cell] != 30 && mesh->BCp.type[cell] != 31 ) {
      // Interpolate from vertex to centers
      Exzc = 0.25*(mesh->exz[vertSW] + mesh->exz[vertSE] + mesh->exz[vertNW] + mesh->exz[vertNE]);
      Txzc = 0.25*(mesh->sxz[vertSW] + mesh->sxz[vertSE] + mesh->sxz[vertNW] + mesh->sxz[vertNE]);
      // Compute invariants
      Eiic = sqrt( 0.5*(pow(mesh->exxd[cell], 2) + pow(mesh->ezzd[cell], 2) ) + pow( Exzc, 2) );
      Tiic = sqrt( 0.5*(pow(mesh->sxxd[cell], 2) + pow(mesh->szzd[cell], 2) ) + pow( Txzc, 2) );
      // Reduce
      sxxd_mean += mesh->sxxd[cell];
      szzd_mean += mesh->szzd[cell];
      sxz_mean  += Txzc;
      Tii_mean  += Tiic;
      exxd_mean += mesh->exxd[cell];
      ezzd_mean += mesh->ezzd[cell];
      exz_mean  += Exzc;
      Eii_mean  += Eiic;
      T_mean    += mesh->T[cell];
      P_mean    += mesh->p_in[cell];
      N++;
    }
  }

  // Divide by number of active cells
  mesh->T_mean    = T_mean    /N;
  mesh->P_mean    = P_mean    /N;
  mesh->exxd_mean = exxd_mean /N;
  mesh->ezzd_mean = ezzd_mean /N;
  mesh->exz_mean  = exz_mean  /N;
  mesh->Eii_mean  = Eii_mean  /N;
  mesh->sxxd_mean = sxxd_mean /N;
  mesh->szzd_mean = szzd_mean /N;
  mesh->sxz_mean  = sxz_mean  /N;
  mesh->Tii_mean  = Tii_mean  /N;

}

