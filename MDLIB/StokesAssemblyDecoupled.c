// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2019  MDOODZ Developper team
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

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddCoeff3( int* J, double*A, int eqn, int jeq, int *nnzc, double coeff, int NODE_TYPE, double NODE_VAL, double* RHS ) {

    if (NODE_TYPE == 0 || NODE_TYPE == 31 || NODE_TYPE == 11 || NODE_TYPE == 13) {
        RHS[eqn]  -= coeff * NODE_VAL;
    }
    else {
        J[(*nnzc)] = jeq;
        A[(*nnzc)] = coeff;
        (*nnzc)++;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int ix, int iz, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int nz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    const int nxvz = nx+1, ncx = nx-1;
    const double dx = mesh->dx, dz = mesh->dz;
    double OOP = 1.0, uC_corr = 0.0;
    if ( model.out_of_plane==1 ) OOP = 3.0/2.0;
    char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type, *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    // Regular or stress BC nodes?
    double SxxBCW = 0., SxxBCE = 0.;
    if (ix==0    && BCu_t[c1]==2) SxxBCW = 1.0;
    if (ix==nx-1 && BCu_t[c1]==2) SxxBCE = 1.0;
    
    // Stencil
    int iVxC   = c1;
    int iVxW   = iVxC-1;
    int iVxE   = iVxC+1;
    int iVxS   = iVxC-nx;
    int iVxN   = iVxC+nx;
    int iVzSW  = c3-nxvz;
    int iVzSE  = c3-nxvz+1;
    int iVzNW  = c3;
    int iVzNE  = c3+1;
    int iPrW   = c2;
    int iPrE   = c2+1;
    int ixyN   = c1;
    int ixyS   = c1-nx;

    // Linear equation coefficient
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;
    
    // Periodic stencil
    if ( BCu_t[iVxC]==-2 ) {
        iVxW  = c1+nx-2;
        iPrW  = c2+ncx;
        iVzSW = c3-2;
        iVzNW = c3+nxvz-2;
    }

    // Variable PDE coefficient
    const double D11E = mesh->D11_n[iPrE];
    const double D11W = mesh->D11_n[iPrW];
    const double D33N = mesh->D33_s[ixyN];
    const double D33S = mesh->D33_s[ixyS];
    
    // Flags
    double inE=0.0, inW=0.0;
    if (BCp_t[iPrW] == -1) inW = 1.0;
    if (BCp_t[iPrE] == -1) inE = 1.0;
    
    double inS=1.0, inN = 1.0, RegS = 1.0, RegN = 1.0, NeuS = 0.0, NeuN = 0.0, DirS = 0.0, DirN = 0.0;
    if (BCg_t[ixyS] == 30) inS = 0.0;
    if (BCg_t[ixyN] == 30) inN = 0.0;
    if (BCu_t[iVxS] == 13) {NeuS = 1.0; RegS = 0.0;}
    if (BCu_t[iVxN] == 13) {NeuN = 1.0; RegN = 0.0;}
    if (BCu_t[iVxS] == 11) {DirS = 1.0; RegS = 0.0;}
    if (BCu_t[iVxN] == 11) {DirN = 1.0; RegN = 0.0;}

    double DirSW = 0.0, DirNW = 0.0;
    if (BCv_t[c3-nxvz]   == 11) DirSW = 1.0;
    if (BCv_t[c3]        == 11) DirNW = 1.0;
    double DirSE = 0.0, DirNE = 0.0;
    if (BCv_t[c3-nxvz+1] == 11) DirSE = 1.0;
    if (BCv_t[c3+1]      == 11) DirNE = 1.0;

    // Linear equation coefficient (AssembleGeneralStiffness_MDOODZ_7.0_v4_StressBC.ipynb)
    uW = (1.0/3.0)*D11W*inW*(1 - SxxBCW)*(OOP*comp*inW - 3)/pow(dx, 2);
    uC = -1.0/6.0*(SxxBCE*(2*D11W*pow(dz, 2)*inW*(OOP*comp*inW - 3) - 3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS))) + SxxBCW*(2*D11E*pow(dz, 2)*inE*(OOP*comp*inE - 3) - 3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS))) + 2*(3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS)) - pow(dz, 2)*(D11E*inE*(OOP*comp*inE - 3) + D11W*inW*(OOP*comp*inW - 3)))*(SxxBCE + SxxBCW - 1))/(pow(dx, 2)*pow(dz, 2));
    uE = (1.0/3.0)*D11E*inE*(1 - SxxBCE)*(OOP*comp*inE - 3)/pow(dx, 2);
    uS = (1.0/2.0)*D33S*RegS*inS*(SxxBCE + SxxBCW - 2)/pow(dz, 2);
    uN = (1.0/2.0)*D33N*RegN*inN*(SxxBCE + SxxBCW - 2)/pow(dz, 2);
    vSW = (1.0/6.0)*(-3*D33S*SxxBCW*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1) + SxxBCE*(2*D11W*OOP*comp*pow(inW, 2) - 3*D33S*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1)) - 2*(D11W*OOP*comp*pow(inW, 2) - 3*D33S*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    vSE = (1.0/6.0)*(3*D33S*SxxBCE*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1) - SxxBCW*(2*D11E*OOP*comp*pow(inE, 2) - 3*D33S*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1)) + 2*(D11E*OOP*comp*pow(inE, 2) - 3*D33S*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    vNW = (1.0/6.0)*(3*D33N*SxxBCW*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1) - SxxBCE*(2*D11W*OOP*comp*pow(inW, 2) - 3*D33N*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1)) + 2*(D11W*OOP*comp*pow(inW, 2) - 3*D33N*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    vNE = (1.0/6.0)*(-3*D33N*SxxBCE*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1) + SxxBCW*(2*D11E*OOP*comp*pow(inE, 2) - 3*D33N*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1)) - 2*(D11E*OOP*comp*pow(inE, 2) - 3*D33N*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    pW = (SxxBCW - 1)/dx;
    pE = (1 - SxxBCE)/dx;

    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodx  = (mesh->rho_n[c2+1] - mesh->rho_n[c2])*one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        // Importante trique, voire meme gigantesque!
        if (uC+uC_corr<0.0) uC_corr = 0.0;
        uC += uC_corr;
    }
    
    if ( Assemble == 1 ) {
        const bool SetVxS  = BCu_t[iVxS] != 30;
        const bool SetVxW  = BCu_t[iVxW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVxC  = BCu_t[iVxC] != 30;
        const bool SetVxE  = BCu_t[iVxE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        const bool SetVxN  = BCu_t[iVxN] != 30;
        //--------------------
        const bool SetVzSW = BCv_t[iVzSW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVzSE = BCv_t[iVzSE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        const bool SetVzNW = BCv_t[iVzNW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVzNE = BCv_t[iVzNE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        //--------------------
        const bool SetPrW  = BCp_t[iPrW]  != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetPrE  = BCp_t[iPrE]  != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        //--------------------
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (SetVxS)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS],  &(nnzc2A[ith]), uS*celvol, BCu_t[iVxS],     mesh->BCu.val[iVxS], StokesA->bbc);        
        if (SetVxW)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxW],  &(nnzc2A[ith]), uW*celvol, BCu_t[iVxW],  mesh->BCu.val[iVxW],  StokesA->bbc);
        if (SetVxC)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, BCu_t[iVxC],    mesh->BCu.val[iVxC],    StokesA->bbc);
        if (SetVxE)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxE],  &(nnzc2A[ith]), uE*celvol, BCu_t[iVxE],  mesh->BCu.val[iVxE],  StokesA->bbc);
        if (SetVxN)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN],  &(nnzc2A[ith]), uN*celvol, BCu_t[iVxN],   mesh->BCu.val[iVxN], StokesA->bbc);
        //--------------------
        if (SetVzSW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW], &(nnzc2A[ith]), vSW*celvol, BCv_t[iVzSW], mesh->BCv.val[iVzSW], StokesA->bbc);
        if (SetVzSE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE], &(nnzc2A[ith]), vSE*celvol, BCv_t[iVzSE], mesh->BCv.val[iVzSE], StokesA->bbc);
        if (SetVzNW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW], &(nnzc2A[ith]), vNW*celvol, BCv_t[iVzNW], mesh->BCv.val[iVzNW], StokesA->bbc);
        if (SetVzNE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE], &(nnzc2A[ith]), vNE*celvol, BCv_t[iVzNE], mesh->BCv.val[iVzNE], StokesA->bbc);
        //--------------------
        // if ( BCp_t[iPrE] != 30 && BCp_t[iPrW] != 30 ) {
        if (SetPrW)  AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrW] - Stokes->neq_mom, &(nnzc2B[ith]), pW*celvol, BCp_t[iPrW], mesh->BCp.val[iPrW], StokesB->bbc);
        if (SetPrE)  AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrE] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, BCp_t[iPrE], mesh->BCp.val[iPrE], StokesB->bbc);
        // }
    }
    else {
        // Residual function
        double F_Reg = 0., F_W =0., F_E =0., SxxE, SxxW;
        
        // Residual: regular nodes
        F_Reg += (mesh->sxxd[iPrE]   - mesh->sxxd[iPrW]  )/dx;
        F_Reg -= (mesh->p_corr[iPrE] - mesh->p_corr[iPrW])/dx;
        F_Reg += (mesh->sxz[ixyN]    - mesh->sxz[ixyS]   )/dz;

        // Residual: stress BC west
        SxxE = mesh->sxxd[iPrE] - mesh->p_corr[iPrE];
        SxxW = 2.0*mesh->BCu.val[iVxC] - SxxE;
        F_W += (SxxE - SxxW)/dx;
        F_W += (mesh->sxz[ixyN] - mesh->sxz[ixyS])/dz;
        F_W *= 0.5;

        // Residual: stress BC east
        SxxW = mesh->sxxd[iPrW] - mesh->p_corr[iPrW];
        SxxE = 2.0*mesh->BCu.val[iVxC] - SxxW;
        F_E += (SxxE - SxxW)/dx;
        F_E += (mesh->sxz[ixyN] - mesh->sxz[ixyS])/dz;
        F_E *= 0.5;

        // Residual
        StokesA->F[eqn]  = (1.0 - SxxBCE - SxxBCW)*F_Reg + SxxBCW*F_W + SxxBCE*F_E;
        StokesA->F[eqn] += StokesA->b[eqn] - uC_corr*u[iVxC]; // no body force
        StokesA->F[eqn] *= -celvol;
    }
} 

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
void Xjacobian_InnerNodesDecoupled3( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int ix, int iz, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B, int k, int l  ) {
    
    int Newton = 1;
    double out_of_plane = 1.0, uC_corr=0.0;
    if (model.out_of_plane==1) out_of_plane = 3.0/2.0;
    double dx = mesh->dx;
    double dz = mesh->dz;
    int nzvx = model.Nz+1, nz = model.Nz;

    // const int nxvz = nx+1, ncx = nx-1;
    // const double dx = mesh->dx, dz = mesh->dz;
    // double OOP = 1.0, uC_corr = 0.0;
    // if ( model.out_of_plane==1 ) OOP = 3.0/2.0;
    char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type, *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    // Regular or stress BC nodes?
    double SxxBCW = 0., SxxBCE = 0.;
    if (ix==0    && BCu_t[c1]==2) SxxBCW = 1.0;
    if (ix==nx-1 && BCu_t[c1]==2) SxxBCE = 1.0;
    
    int iVxC   = c1;
    int iVxW   = iVxC-1;
    int iVxE   = iVxC+1;
    int iVxS   = iVxC-nx;
    int iVxN   = iVxC+nx;
    int iVxSW  = iVxS-1;
    int iVxSE  = iVxS+1;
    int iVxNW  = iVxN-1;
    int iVxNE  = iVxN+1;
    int iVzSW  = c3-nxvz;
    int iVzSE  = c3-nxvz+1;
    int iVzNW  = c3;
    int iVzNE  = c3+1;
    int iVzSWW = iVzSW-1;
    int iVzSEE = iVzSE+1;
    int iVzNWW = iVzNW-1;
    int iVzNEE = iVzNE+1;
    int iVzNNW = iVzNW+nxvz, iVzSSW = iVzSW-nxvz;
    int iVzNNE = iVzNE+nxvz, iVzSSE = iVzSE-nxvz;
    
    int iPrW  = c2;
    int iPrE  = c2+1;
    int ixyN  = c1;
    int ixyS  = c1-nx;
    int ixySW = ixyS-1;
    int ixySE = ixyS+1;
    int ixyNW = ixyN-1;
    int ixyNE = ixyN+1;
    int iPrSW = iPrW - ncx;
    int iPrSE = iPrE - ncx;
    int iPrNW = iPrW + ncx;
    int iPrNE = iPrE + ncx;
    
    // Periodic ends of stencil for inner points
    if (l>1)    if (BCv_t[iVzSWW] == -12)  iVzSWW = iVzSW + (nx-2);
    if (l<nz-1) if (BCv_t[iVzNWW] == -12)  iVzNWW = iVzNW + (nx-2);
    if (l>1)    if (BCv_t[iVzSEE] == -12)  iVzSEE = iVzSE - (nx-2);
    if (l<nz-1) if (BCv_t[iVzNEE] == -12)  iVzNEE = iVzNE - (nx-2);
    
    if (BCu_t[iVxC] == -2) {
        iVxW   = c1+nx-2;      iPrW = c2+ncx; iVzSW = c3-2; iVzNW = c3+nxvz-2;
        iPrNW  = iPrW+ncx;    iPrSW = iPrW-ncx;
        iVxSW  = iVxW-nx;     iVxNW = iVxW+nx;
        iVzNNW = iVzNW+nxvz; iVzSSW = iVzSW-nxvz;
        iVzSWW = iVzSW-1;    iVzNWW = iVzNW-1;
        // Valgrind / Audresselles 28/12/21
        ixySW  = ixyS + nx;
        ixyNW  = ixyN + nx;
    }

    if (k==0) { // Valgrind / Audresselles 28/12/21
        ixySW  = ixyS;
        ixyNW  = ixyN;
    }

    if (k==nx-1) { // Valgrind / Audresselles 28/12/21
        ixySE  = ixyS; 
        ixyNE  = ixyN; 
    }
            
    // The computation of FD coefficients is only useful for the purpose of the stiffness/Jacobian matrix assembly
    double  uC=0.0;
    double  uS=0.0,  uN=0.0,  uW=0.0,  uE=0.0,  vSW=0.0,  vSE=0.0,  vNW=0.0,  vNE=0.0, pE=0.0, pW=0.0;
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vSWW=0.0, vSEE=0.0, vNWW=0.0, vNEE=0.0, vSSW=0.0, vSSE = 0.0, vNNW=0.0, vNNE=0.0;
    double pSW=0.0, pSE=0.0, pNW=0.0, pNE=0.0;
    
    double D11E = mesh->D11_n[iPrE];
    double D12E = mesh->D12_n[iPrE];
    double D13E = mesh->D13_n[iPrE];
    double D14E = mesh->D14_n[iPrE];
    
    double D11W = mesh->D11_n[iPrW];
    double D12W = mesh->D12_n[iPrW];
    double D13W = mesh->D13_n[iPrW];
    double D14W = mesh->D14_n[iPrW];
    
    double D31N = mesh->D31_s[ixyN];
    double D32N = mesh->D32_s[ixyN];
    double D33N = mesh->D33_s[ixyN];
    double D34N = mesh->D34_s[ixyN];
    
    double D31S = mesh->D31_s[ixyS];
    double D32S = mesh->D32_s[ixyS];
    double D33S = mesh->D33_s[ixyS];
    double D34S = mesh->D34_s[ixyS];
    
    double inE=0.0, inW=0.0, inS=0.0, inN = 0.0, inSv = 0.0, inNv = 0.0;
  
    if (BCp_t[iPrW] == -1) inW = 1.0;
    if (BCp_t[iPrE] == -1) inE = 1.0;
    
    //        if (mesh->BCg.type[ixyS] != 30 && BCu_t[iVxS] != 13) inS = 1.0; // !!!!!!!!!!!!
    //        if (mesh->BCg.type[ixyN] != 30 && BCu_t[iVxN] != 13) inN = 1.0;
    if (BCu_t[iVxS] != 30 && BCu_t[iVxS] != 13 )  inS  = 1.0;
    if (BCu_t[iVxN] != 30 && BCu_t[iVxN] != 13 )  inN  = 1.0;
    if (BCg_t[ixyS] != 30 ) inSv = 1.0;
    if (BCg_t[ixyN] != 30 ) inNv = 1.0;

    // X-tra
    double wE=0.0, wW=0.0, wS=0.0, wN = 0.0;
    double inSWc=0.0,inSEc=0.0,inNWc=0.0,inNEc=0.0;
    double inSWv=0.0,inSEv=0.0,inNWv=0.0,inNEv=0.0;
    
    if ( l>1 ){// || (l==1 && BCu_t[iVxS] == 11 ) ) {  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Apparently incorrect
        if (BCp_t[iPrSW] == -1) inSWc = 1.0;
        if (BCp_t[iPrSE] == -1) inSEc = 1.0;
    }
    
    if ( l<nzvx-2 ){// || (l==nzvx-2 && BCu_t[iVxN] == 11 )) { //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Apparently incorrect
        if (BCp_t[iPrNW] == -1) inNWc = 1.0;
        if (BCp_t[iPrNE] == -1) inNEc = 1.0;
    }

    if ( (k>0)  || (k==0 && (BCv_t[iVzSW] == -1 || BCv_t[iVzSW] == 0)) ) {
        if (BCg_t[ixySW] != 30) inSWv = 1.0;   // modify for periodic
        if (BCg_t[ixyNW] != 30) inNWv = 1.0;   // modify for periodic
    }
    if ( (k<nx-1) || (k==nx-1 && (BCv_t[iVzSE] == -1 || BCv_t[iVzSE] == 0) ) ) {
        if (BCg_t[ixySE] != 30) inSEv = 1.0;   // modify for periodic
        if (BCg_t[ixyNE] != 30) inNEv = 1.0;   // modify for periodic
    }
    
    // New stuff
    double inSWW = 0.0, inNWW = 0.0;
    if (BCv_t[iVzSWW] == -1) inSWW = 1.0;
    if (BCv_t[iVzNWW] == -1) inNWW = 1.0;
    
    double inSEE = 0.0, inNEE = 0.0;
    if (BCv_t[iVzSEE] == -1) inSEE = 1.0;
    if (BCv_t[iVzNEE] == -1) inNEE = 1.0;
    
    wE = inN + inS + inNEv + inSEv;
    wW = inN + inS + inNWv + inSWv;
    wS = inW + inE + inSWc + inSEc;
    wN = inW + inE + inNWc + inNEc;
    
    if (wW>1.0) wW = 1.0/wW;
    if (wE>1.0) wE = 1.0/wE;
    if (wS>1.0) wS = 1.0/wS;
    if (wN>1.0) wN = 1.0/wN;
    
    // Audresselles 29/12/21
    uW = (1.0/3.0)*(-0.75*D13W*dx*(inN*inNWv - inS*inSWv) + dx*inW*(-D31N*wN*(comp*out_of_plane - 3) + D31S*wS*(comp*out_of_plane - 3) - D32N*comp*out_of_plane*wN + D32S*comp*out_of_plane*wS) + dz*inW*(D11W*(comp*out_of_plane - 3) + D12W*comp*out_of_plane))/(pow(dx, 2)*dz);
    uC = (1.0/3.0)*(dx*(3*dx*(D33N*pow(inN, 2)*inNv + D33S*pow(inS, 2)*inSv) + dz*(inE - inW)*(-D31N*wN*(comp*out_of_plane - 3) + D31S*wS*(comp*out_of_plane - 3) - D32N*comp*out_of_plane*wN + D32S*comp*out_of_plane*wS)) - dz*(0.75*dx*(-D13E + D13W)*(pow(inN, 2)*inNv - pow(inS, 2)*inSv) + dz*(D11E*inE*(comp*out_of_plane - 3) + D11W*inW*(comp*out_of_plane - 3) + D12E*comp*inE*out_of_plane + D12W*comp*inW*out_of_plane)))/(pow(dx, 2)*pow(dz, 2));
    uE = (1.0/3.0)*(0.75*D13E*dx*(inN*inNEv - inS*inSEv) + dx*inE*(D31N*wN*(comp*out_of_plane - 3) - D31S*wS*(comp*out_of_plane - 3) + D32N*comp*out_of_plane*wN - D32S*comp*out_of_plane*wS) + dz*inE*(D11E*(comp*out_of_plane - 3) + D12E*comp*out_of_plane))/(pow(dx, 2)*dz);
    uS = (-D33S*dx*pow(inS, 2)*inSv + 0.25*dz*pow(inS, 2)*inSv*(D13E - D13W) + (1.0/3.0)*dz*wS*(inSEc - inSWc)*(D31S*(comp*out_of_plane - 3) + D32S*comp*out_of_plane))/(dx*pow(dz, 2));
    uN = (-D33N*dx*pow(inN, 2)*inNv + 0.25*dz*pow(inN, 2)*inNv*(-D13E + D13W) - 1.0/3.0*dz*wN*(inNEc - inNWc)*(D31N*(comp*out_of_plane - 3) + D32N*comp*out_of_plane))/(dx*pow(dz, 2));
    vSW = (1.0/3.0)*(-dx*(3*D33S*dz*inS*inSv + dx*(D31N*comp*inW*out_of_plane*wN + D31S*comp*out_of_plane*wS*(inSWc - inW) + D32N*inW*wN*(comp*out_of_plane - 3) + D32S*wS*(inSWc - inW)*(comp*out_of_plane - 3))) + dz*(dx*inW*(D11W*comp*out_of_plane + D12W*(comp*out_of_plane - 3)) + 0.75*dz*(D13E*inS*inSv - D13W*(inS*inSv - inSWW*inSWv))))/(pow(dx, 2)*pow(dz, 2));
    vSE = (1.0/3.0)*(dx*(3*D33S*dz*inS*inSv + dx*(-D31N*comp*inE*out_of_plane*wN + D31S*comp*out_of_plane*wS*(inE - inSEc) - D32N*inE*wN*(comp*out_of_plane - 3) + D32S*wS*(inE - inSEc)*(comp*out_of_plane - 3))) - dz*(dx*inE*(D11E*comp*out_of_plane + D12E*(comp*out_of_plane - 3)) + 0.75*dz*(D13E*(inS*inSv - inSEE*inSEv) - D13W*inS*inSv)))/(pow(dx, 2)*pow(dz, 2));
    vNW = (1.0/3.0)*(dx*(3*D33N*dz*inN*inNv - dx*(D31N*comp*out_of_plane*wN*(inNWc - inW) + D31S*comp*inW*out_of_plane*wS + D32N*wN*(inNWc - inW)*(comp*out_of_plane - 3) + D32S*inW*wS*(comp*out_of_plane - 3))) - dz*(dx*inW*(D11W*comp*out_of_plane + D12W*(comp*out_of_plane - 3)) + 0.75*dz*(-D13E*inN*inNv + D13W*(inN*inNv - inNWW*inNWv))))/(pow(dx, 2)*pow(dz, 2));
    vNE = (1.0/3.0)*(-dx*(3*D33N*dz*inN*inNv + dx*(-D31N*comp*out_of_plane*wN*(inE - inNEc) + D31S*comp*inE*out_of_plane*wS - D32N*wN*(inE - inNEc)*(comp*out_of_plane - 3) + D32S*inE*wS*(comp*out_of_plane - 3))) + dz*(dx*inE*(D11E*comp*out_of_plane + D12E*(comp*out_of_plane - 3)) + 0.75*dz*(-D13E*(inN*inNv - inNEE*inNEv) + D13W*inN*inNv)))/(pow(dx, 2)*pow(dz, 2));
    uSW = (1.0/3.0)*(-0.75*D13W*inS*inSWv + D31S*inSWc*wS*(comp*out_of_plane - 3) + D32S*comp*inSWc*out_of_plane*wS)/(dx*dz);
    uSE = (1.0/3.0)*(0.75*D13E*inS*inSEv - D31S*inSEc*wS*(comp*out_of_plane - 3) - D32S*comp*inSEc*out_of_plane*wS)/(dx*dz);
    uNW = (1.0/3.0)*(0.75*D13W*inN*inNWv - D31N*inNWc*wN*(comp*out_of_plane - 3) - D32N*comp*inNWc*out_of_plane*wN)/(dx*dz);
    uNE = (1.0/3.0)*(-0.75*D13E*inN*inNEv + D31N*inNEc*wN*(comp*out_of_plane - 3) + D32N*comp*inNEc*out_of_plane*wN)/(dx*dz);
    vSWW = -0.25*D13W*inSWW*inSWv/pow(dx, 2);
    vSEE = -0.25*D13E*inSEE*inSEv/pow(dx, 2);
    vNWW = -0.25*D13W*inNWW*inNWv/pow(dx, 2);
    vNEE = -0.25*D13E*inNEE*inNEv/pow(dx, 2);
    vSSW = (1.0/3.0)*inSWc*wS*(D31S*comp*out_of_plane + D32S*(comp*out_of_plane - 3))/pow(dz, 2);
    vSSE = (1.0/3.0)*inSEc*wS*(D31S*comp*out_of_plane + D32S*(comp*out_of_plane - 3))/pow(dz, 2);
    vNNW = (1.0/3.0)*inNWc*wN*(D31N*comp*out_of_plane + D32N*(comp*out_of_plane - 3))/pow(dz, 2);
    vNNE = (1.0/3.0)*inNEc*wN*(D31N*comp*out_of_plane + D32N*(comp*out_of_plane - 3))/pow(dz, 2);      
    pW   = -inW*one_dx + inW*(D14W*dz + dx*(-D34N*inN*wN + D34S*inS*wS))/(dx*dz);
    pE   =  inE*one_dx + inE*(-D14E*dz + dx*(-D34N*inN*wN + D34S*inS*wS))/(dx*dz);
    pSW  = D34S*inS*inSWc*wS/dz;
    pSE  = D34S*inS*inSEc*wS/dz;
    pNW  = -D34N*inN*inNWc*wN/dz;
    pNE  = -D34N*inN*inNEc*wN/dz;

    //---------------------------------------------------//

    // Stabilisation with density gradients
    if (stab==1) {
        double drhodx  = (mesh->rho_n[c2+1] - mesh->rho_n[c2])*one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        // Importante trique, voire meme gigantesque!
        if (uC+uC_corr>0.0) uC += uC_corr;
    }
    
    // Add contribution from non-conforming Dirichlets
    if ( BCu_t[iVxS]   == 11 ) uC  -=  uS ;
    if ( BCu_t[iVxN]   == 11 ) uC  -=  uN ;
    if ( BCu_t[iVxSW]  == 11 ) uW  -=  uSW;
    if ( BCu_t[iVxSE]  == 11 ) uE  -=  uSE;
    if ( BCu_t[iVxNW]  == 11 ) uW  -=  uNW;
    if ( BCu_t[iVxNE]  == 11 ) uE  -=  uNE;
    if ( BCv_t[iVzNWW] == 11 ) vNW -= vNWW;
    if ( BCv_t[iVzNEE] == 11 ) vNE -= vNEE;
    if ( BCv_t[iVzSWW] == 11 ) vSW -= vSWW;
    if ( BCv_t[iVzSEE] == 11 ) vSE -= vSEE;

    if ( Assemble == 1 ) {

        const bool SetVxSW = BCu_t[iVxSW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVxS  = BCu_t[iVxS]  != 30;
        const bool SetVxSE = BCu_t[iVxSE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        const bool SetVxW  = BCu_t[iVxW]  != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVxC  = BCu_t[iVxC]  != 30;
        const bool SetVxE  = BCu_t[iVxE]  != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        const bool SetVxNW = BCu_t[iVxNW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVxN  = BCu_t[iVxN]  != 30;
        const bool SetVxNE = BCu_t[iVxNE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        //--------------------
        const bool SetVzSW = BCv_t[iVzSW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVzSE = BCv_t[iVzSE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        const bool SetVzNW = BCv_t[iVzNW] != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetVzNE = BCv_t[iVzNE] != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);
        //--------------------
        const bool SetPrW  = BCp_t[iPrW]  != 30 && (ix>0    || BCu_t[iVxC]==-2);
        const bool SetPrE  = BCp_t[iPrE]  != 30 && (ix<nx-1 || BCu_t[iVxC]==-2);

        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;

        //--------------------
        // dsxx/dx - normal stencil
        if (SetVxSW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, BCu_t[iVxSW],     mesh->BCu.val[iVxSW], StokesA->bbc);        
        if (SetVxS)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS], &(nnzc2A[ith]), uS*celvol, BCu_t[iVxS],     mesh->BCu.val[iVxS], StokesA->bbc);
        if (SetVxSE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, BCu_t[iVxSE],     mesh->BCu.val[iVxSE], StokesA->bbc);
        if (SetVxW)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxW],  &(nnzc2A[ith]), uW*celvol, BCu_t[iVxW],  mesh->BCu.val[iVxW],  StokesA->bbc);
        if (SetVxC)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, BCu_t[iVxC],    mesh->BCu.val[iVxC],    StokesA->bbc);
        if (SetVxE)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxE],  &(nnzc2A[ith]), uE*celvol, BCu_t[iVxE],  mesh->BCu.val[iVxE],  StokesA->bbc);
        if (SetVxNW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, BCu_t[iVxNW],     mesh->BCu.val[iVxNW], StokesA->bbc);
        if (SetVxN)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN], &(nnzc2A[ith]), uN*celvol, BCu_t[iVxN],   mesh->BCu.val[iVxN], StokesA->bbc);
        if (SetVxNE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, BCu_t[iVxNE],     mesh->BCu.val[iVxNE], StokesA->bbc);
        
        //--------------------
        
        // vSSW
        if ( l>1 ) if ( BCv_t[iVzSSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSSW],      &(nnzc2A[ith]), vSSW*celvol, BCv_t[iVzSSW],      mesh->BCv.val[iVzSSW],      StokesA->bbc);
        if ( l>1 ) if ( BCv_t[iVzSSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSSE],      &(nnzc2A[ith]), vSSE*celvol, BCv_t[iVzSSE],      mesh->BCv.val[iVzSSE],      StokesA->bbc);
        
        //--------------------
        
        
        // vSWW (Newton)
        if ( (k > 1 ) ||  ( (k==0 || k==1) &&  BCv_t[iVzSWW] == -1) ) {
            if (BCv_t[iVzSWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSWW],   &(nnzc2A[ith]), vSWW*celvol, BCv_t[iVzSWW],   mesh->BCv.val[iVzSWW],   StokesA->bbc);
        }
        
        // vSW && vSE
        if (SetVzSW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW], &(nnzc2A[ith]), vSW*celvol, BCv_t[iVzSW], mesh->BCv.val[iVzSW],   StokesA->bbc);
        if (SetVzSE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE], &(nnzc2A[ith]), vSE*celvol, BCv_t[iVzSE], mesh->BCv.val[iVzSE], StokesA->bbc);
        
        // vSEE (Newton)
        if ( (k < nx-2 ) ||  ( (k==nx-2 || k==nx-1) &&  BCv_t[iVzSEE] == -1) ) {
            if ( BCv_t[iVzSEE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSEE], &(nnzc2A[ith]), vSEE*celvol, BCv_t[iVzSEE], mesh->BCv.val[iVzSEE], StokesA->bbc);
        }
        
        // vNWW (Newton)
        if ( (k > 1 ) ||  ( (k==0 || k==1) &&  BCv_t[iVzNWW] == -1) ) {
            if ( BCv_t[iVzNWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNWW],        &(nnzc2A[ith]), vNWW*celvol, BCv_t[iVzNWW],        mesh->BCv.val[iVzNWW],        StokesA->bbc);
        }
        
        // vNW && vNE
        if (SetVzNW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],      &(nnzc2A[ith]), vNW*celvol, BCv_t[iVzNW],      mesh->BCv.val[iVzNW],      StokesA->bbc);
        if (SetVzNE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE],      &(nnzc2A[ith]), vNE*celvol, BCv_t[iVzNE],      mesh->BCv.val[iVzNE],      StokesA->bbc);
        
        // vNEE (Newton)
        if ( (k < nx-2 ) ||  ((k==nx-2 || k==nx-1) &&  BCv_t[iVzNEE] == -1) ) {
            if ( BCv_t[iVzNEE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNEE],      &(nnzc2A[ith]), vNEE*celvol, BCv_t[iVzNEE],      mesh->BCv.val[iVzNEE],      StokesA->bbc);
        }
        
        //--------------------

        if ( l<nz-1 ) if ( BCv_t[iVzNNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNNW],      &(nnzc2A[ith]), vNNW*celvol, BCv_t[iVzNNW],      mesh->BCv.val[iVzNNW],      StokesA->bbc);
        if ( l<nz-1 ) if ( BCv_t[iVzNNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNNE],      &(nnzc2A[ith]), vNNE*celvol, BCv_t[iVzNNE],      mesh->BCv.val[iVzNNE],      StokesA->bbc);
        
        //--------------------
        
        if ( Newton==1 && l>1) {
            if ( BCp_t[iPrSW] != 30 ) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSW] - Stokes->neq_mom,   &(nnzc2B[ith]), pSW*celvol, BCp_t[iPrSW],   mesh->BCp.val[iPrSW],   StokesB->bbc);
            if ( BCp_t[iPrSE] != 30 ) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSE] - Stokes->neq_mom,   &(nnzc2B[ith]), pSE*celvol, BCp_t[iPrSE],   mesh->BCp.val[iPrSE],   StokesB->bbc);
        }
        
        // pE && pW -- Valgrind / Audresselles 28/12/21
        // if ( BCp_t[iPrE] != 30 && BCp_t[iPrW] != 30 ) {
            if (SetPrW)  AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrW] - Stokes->neq_mom, &(nnzc2B[ith]), pW*celvol, BCp_t[iPrW], mesh->BCp.val[iPrW], StokesB->bbc);
            if (SetPrE)  AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrE] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, BCp_t[iPrE], mesh->BCp.val[iPrE], StokesB->bbc);
        // }
        
        if ( Newton==1 && l<nz-1 ) {
            if ( BCp_t[iPrNW] != 30 ) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNW] - Stokes->neq_mom,   &(nnzc2B[ith]), pNW*celvol, BCp_t[iPrNW],   mesh->BCp.val[iPrNW],   StokesB->bbc);
            if ( BCp_t[iPrNE] != 30 ) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNE] - Stokes->neq_mom,   &(nnzc2B[ith]), pNE*celvol, mesh->BCp.type[iPrNE],   mesh->BCp.val[iPrNE],   StokesB->bbc);
        }
    }
    else {
        
        // Residual function
        double F_Reg = 0., F_W =0., F_E =0., SxxE, SxxW;
        
        // Residual: regular nodes
        F_Reg += (mesh->sxxd[iPrE]   - mesh->sxxd[iPrW]  )/dx;
        F_Reg -= (mesh->p_corr[iPrE] - mesh->p_corr[iPrW])/dx;
        F_Reg += (mesh->sxz[ixyN]    - mesh->sxz[ixyS]   )/dz;

        // Residual: stress BC west
        SxxE = mesh->sxxd[iPrE] - mesh->p_corr[iPrE];
        SxxW = 2.0*mesh->BCu.val[iVxC] - SxxE;
        F_W += (SxxE - SxxW)/dx;
        F_W += (mesh->sxz[ixyN] - mesh->sxz[ixyS])/dz;
        F_W *= 0.5;

        // Residual: stress BC east
        SxxW = mesh->sxxd[iPrW] - mesh->p_corr[iPrW];
        SxxE = 2.0*mesh->BCu.val[iVxC] - SxxW;
        F_E += (SxxE - SxxW)/dx;
        F_E += (mesh->sxz[ixyN] - mesh->sxz[ixyS])/dz;
        F_E *= 0.5;

        // Residual
        StokesA->F[eqn]  = (1.0 - SxxBCE - SxxBCW)*F_Reg + SxxBCW*F_W + SxxBCE*F_E;
        StokesA->F[eqn] += StokesA->b[eqn] - uC_corr*u[iVxC]; // no body force
        StokesA->F[eqn] *= -celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int ix, int iz, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int nz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {

    const int nxvz = nx+1, ncx = nx-1;
    const double dx = mesh->dx, dz = mesh->dz;
    double OOP = 1.0, vC_corr = 0.0;
    if ( model.out_of_plane==1 ) OOP = 3.0/2.0;
    char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type, *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    // Regular or stress BC nodes?
    double SzzBCS = 0., SzzBCN = 0.;
    if (iz==0    && BCv_t[c3]==2) SzzBCS = 1.0;
    if (iz==nz-1 && BCv_t[c3]==2) SzzBCN = 1.0;
    
    // Stencil
    int iVzC   = c3;
    int iVzW   = iVzC-1;
    int iVzE   = iVzC+1;
    int iVzS   = iVzC-nxvz;
    int iVzN   = iVzC+nxvz;
    // ------------------
    int iVxNW  = c1+nx-1;
    int iVxNE  = c1+nx;
    int iVxSW  = c1-1;
    int iVxSE  = c1;
    // ------------------
    int iPrS   = c2;
    int iPrN   = c2+ncx;
    int ixyE   = c1;
    int ixyW   = c1-1;

    if ( BCv_t[iVzW] == -12 && ix==1      ) iVzW = c3+nxvz-3;
    if ( BCv_t[iVzE] == -12 && ix==nxvz-2 ) iVzE = c3-nxvz+3;

    // Linear equation coefficient
    double vS=0.0, vN=0.0, vW=0.0, vE=0.0, vC=0.0, uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, pN=0.0, pS=0.0;

    // Variable PDE coefficient
    const double D22S  = mesh->D22_n[iPrS];
    const double D22N  = mesh->D22_n[iPrN];
    const double D33E  = mesh->D33_s[ixyE];
    const double D33W  = mesh->D33_s[ixyW];

    // Flags
    double inN=0.0, inS=0.0;
    if (BCp_t[iPrS] == -1) inS = 1.0;
    if (BCp_t[iPrN] == -1) inN = 1.0;
    
    double inW=1.0, inE = 1.0, RegW = 1., RegE = 1., NeuW = 0., NeuE = 0, DirW = 0., DirE = 0.;
    if (BCg_t[ixyW] == 30) inW = 0.0;
    if (BCg_t[ixyE] == 30) inE = 0.0;
    if (BCv_t[iVzW] == 13) {NeuW = 1.0; RegW = 0.0;}
    if (BCv_t[iVzE] == 13) {NeuE = 1.0; RegE = 0.0;}
    if (BCv_t[iVzW] == 11) {DirW = 1.0; RegW = 0.0;}
    if (BCv_t[iVzE] == 11) {DirE = 1.0; RegE = 0.0;}

    double DirSW =0., DirNW = 0. ; 
    if (mesh->BCu.type[c1-1]    == 11) DirSW = 1.0;
    if (mesh->BCu.type[c1+nx-1] == 11) DirNW = 1.0;

    double DirSE = 0., DirNE = 0.;    
    if (mesh->BCu.type[c1]      == 11) DirSE = 1.0;
    if (mesh->BCu.type[c1+nx]   == 11) DirNE = 1.0; 

    // Simpler
    vW = (1.0/2.0)*D33W*RegW*inW*(SzzBCN + SzzBCS - 2)/pow(dx, 2);
    vC = (1.0/6.0)*(-SzzBCN*(2*D22S*pow(dx, 2)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW))) - SzzBCS*(2*D22N*pow(dx, 2)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW))) + 2*(pow(dx, 2)*(D22N + D22S)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW)))*(SzzBCN + SzzBCS - 1))/(pow(dx, 2)*pow(dz, 2));
    vE = (1.0/2.0)*D33E*RegE*inE*(SzzBCN + SzzBCS - 2)/pow(dx, 2);
    vS = (1.0/3.0)*D22S*(1 - SzzBCS)*(OOP*comp - 3)/pow(dz, 2);
    vN = (1.0/3.0)*D22N*(1 - SzzBCN)*(OOP*comp - 3)/pow(dz, 2);
    uSW = (1.0/6.0)*(-3*D33W*SzzBCS*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1) + SzzBCN*(2*D22S*OOP*comp - 3*D33W*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1)) - 2*(D22S*OOP*comp - 3*D33W*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    uSE = (1.0/6.0)*(3*D33E*SzzBCS*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1) - SzzBCN*(2*D22S*OOP*comp - 3*D33E*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1)) + 2*(D22S*OOP*comp - 3*D33E*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    uNW = (1.0/6.0)*(3*D33W*SzzBCN*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1) - SzzBCS*(2*D22N*OOP*comp - 3*D33W*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1)) + 2*(D22N*OOP*comp - 3*D33W*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    uNE = (1.0/6.0)*(-3*D33E*SzzBCN*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1) + SzzBCS*(2*D22N*OOP*comp - 3*D33E*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1)) - 2*(D22N*OOP*comp - 3*D33E*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    pS = (SzzBCS - 1)/dz;
    pN = (1 - SzzBCN)/dz;

    // if (eqn==100) printf("%d %d %d\n", eqn, mesh->BCu.type[iVxSW], mesh->BCu.type[iVxNW]);

    // if (mesh->BCu.type[iVxSW] == 2 && BCv_t[iVzW] == 13) uSW = 0.;
    // if (mesh->BCu.type[iVxNW] == 2 && BCv_t[iVzW] == 13) uNW = 0.;
    //  if (mesh->BCu.type[iVxSW] == 2) uSE = 0.;
    //  if (mesh->BCu.type[iVxNW] == 2) uNE = 0.;

    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodz  = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
        vC_corr = 1.00 * om * model.dt * mesh->gz[c3] * drhodz;
        // Importante trique, voire meme gigantesque!
        if ((vC+vC_corr)<0.0 || (vC+vC_corr)<vW || (vC+vC_corr)<vE) vC_corr = 0.0;
        vC += vC_corr;
    }
    
    if ( Assemble == 1 ) {
        const bool SetVzS  = BCv_t[iVzS] != 30 && iz>0;
        const bool SetVzW  = BCv_t[iVzW] != 30 ;
        const bool SetVzC  = BCv_t[iVzC] != 30;
        const bool SetVzE  = BCv_t[iVzE] != 30 ;
        const bool SetVzN  = BCv_t[iVzN] != 30 && iz<nz-1;
        //--------------------
        const bool SetVxSW = BCu_t[iVxSW] != 30 && iz>0   ;
        const bool SetVxSE = BCu_t[iVxSE] != 30 && iz>0   ;
        const bool SetVxNW = BCu_t[iVxNW] != 30 && iz<nz-1;
        const bool SetVxNE = BCu_t[iVxNE] != 30 && iz<nz-1;
        //--------------------
        const bool SetPrS  = BCp_t[iPrS] != 30 && iz>0;
        const bool SetPrN  = BCp_t[iPrN] != 30 && iz<nz-1;
        //--------------------
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------        
        if (SetVxSW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, BCu_t[iVxSW], mesh->BCu.val[iVxSW], StokesA->bbc );
        if (SetVxSE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, BCu_t[iVxSE], mesh->BCu.val[iVxSE], StokesA->bbc );
        if (SetVxNW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, BCu_t[iVxNW], mesh->BCu.val[iVxNW], StokesA->bbc );
        if (SetVxNE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, BCu_t[iVxNE], mesh->BCu.val[iVxNE], StokesA->bbc );
        //--------------------
        if (SetVzS)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzS],  &(nnzc2A[ith]), vS*celvol, BCv_t[iVzS], mesh->BCv.val[iVzS], StokesA->bbc );
        if (SetVzW)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzW],  &(nnzc2A[ith]), vW*celvol, BCv_t[iVzW],    mesh->BCv.val[iVzW],    StokesA->bbc );
        if (SetVzC)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), vC*celvol, BCv_t[iVzC],      mesh->BCv.val[iVzC],      StokesA->bbc );
        if (SetVzE)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzE],  &(nnzc2A[ith]), vE*celvol, BCv_t[iVzE],    mesh->BCv.val[iVzE],    StokesA->bbc );
        if (SetVzN)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzN],  &(nnzc2A[ith]), vN*celvol, BCv_t[iVzN], mesh->BCv.val[iVzN], StokesA->bbc );
        //--------------------
        // if ( BCp_t[iPrS] != 30 && BCp_t[iPrN] != 30) {
            if (SetPrS) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrS] - Stokes->neq_mom,  &(nnzc2B[ith]), pS*celvol, BCp_t[iPrS], mesh->BCp.val[iPrS], StokesB->bbc );
            if (SetPrN) AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrN] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, BCp_t[iPrN], mesh->BCp.val[iPrN], StokesB->bbc );
        // }
    }
    else {

        // Residual function
        double F_Reg = 0., F_S =0., F_N =0., SzzS, SzzN;

        // Residual: regular nodes
        F_Reg += (mesh->szzd[iPrN]   - mesh->szzd[iPrS])/dz;
        F_Reg -= (mesh->p_corr[iPrN] - mesh->p_corr[iPrS])/dz;
        F_Reg += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;

        // Residual: stress BC south
        SzzN = mesh->szzd[iPrN] - mesh->p_corr[iPrN];
        SzzS = 2.0*mesh->BCv.val[iVzC] - SzzN;
        F_S += (SzzN - SzzS)/dz;
        F_S += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        F_S *= 0.5;

        // Residual: stress BC north
        SzzS = mesh->szzd[iPrS] - mesh->p_corr[iPrS];
        SzzN = 2.0*mesh->BCv.val[iVzC] - SzzS;
        F_N += (SzzN - SzzS)/dz;
        F_N += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        F_N *= 0.5;
        
        // Residual
        StokesA->F[eqn]  = (1.0 - SzzBCS - SzzBCN)*F_Reg + SzzBCS*F_S + SzzBCN*F_N;
        StokesA->F[eqn] += StokesA->b[eqn] - vC_corr*v[iVzC]; 
        StokesA->F[eqn] *= -celvol;

        // StokesA->F[eqn]  = 0.0;
        // StokesA->F[eqn] += (mesh->szzd[iPrN]   - mesh->szzd[iPrS])/dz;
        // StokesA->F[eqn] -= (mesh->p_corr[iPrN] - mesh->p_corr[iPrS])/dz;
        // StokesA->F[eqn] += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        // StokesA->F[eqn] *= -1.0;
        // StokesA->F[eqn] -= StokesA->b[eqn] - vC_corr*v[iVzC]; 
        // StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
void Zjacobian_InnerNodesDecoupled3( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int ix, int iz, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B, int k, int l ) {
    
    double dx = mesh->dx;
    double dz = mesh->dz;
    int Newton = 1;
    double out_of_plane = 1.0;
    if (model.out_of_plane==1) out_of_plane = 3.0/2.0;
    int periodix = model.periodic_x;
    double vC_corr = 0.0;

    int nzvx = model.Nz+1;
    int nz   = model.Nz;
    // const int nxvz = nx+1, ncx = nx-1;
    // const double dx = mesh->dx, dz = mesh->dz;
    // double OOP = 1.0, vC_corr = 0.0;
    // if ( model.out_of_plane==1 ) OOP = 3.0/2.0;
    char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type, *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;


      // Regular or stress BC nodes?
    double SzzBCS = 0., SzzBCN = 0.;
    if (iz==0    && BCv_t[c3]==2) SzzBCS = 1.0;
    if (iz==nz-1 && BCv_t[c3]==2) SzzBCN = 1.0;
    
    int iVzC   = c3;
    int iVzW   = iVzC-1,  iVzE   = iVzC+1;
    int iVzS   = c3-nxvz, iVzN   = c3+nxvz;
    int iVzSW  = iVzS-1,  iVzSE  = iVzS+1;
    int iVzNW  = iVzN-1,  iVzNE  = iVzN+1;
    
    int iVxNW  = c1+nx-1;
    int iVxNE  = c1+nx;
    int iVxSW  = c1-1;
    int iVxSE  = c1;
    int iVxNWW = iVxNW-1;
    int iVxNEE = iVxNE+1;
    int iVxSWW = iVxSW-1;
    int iVxSEE = iVxSE+1;
    int iVxNNW = iVxNW+nx;
    int iVxNNE = iVxNE+nx;
    int iVxSSW = iVxSW-nx;
    int iVxSSE = iVxSE-nx;
    
    int iPrS  = c2;
    int iPrSW = iPrS-1;
    int iPrSE = iPrS+1;
    int iPrN  = c2+ncx;
    int iPrNW = iPrN-1;
    int iPrNE = iPrN+1;
    
    int ixyW  = c1-1;
    int ixySW = ixyW-nx;
    int ixyNW = ixyW+nx;
    int ixyE  = c1;
    int ixySE = ixyE-nx;
    int ixyNE = ixyE+nx;
    
    if (k==1    && BCv_t[iVzW]==-12) {
        iVxNWW = iVxNW   + (nx-2);
        iVxSWW = iVxSW   + (nx-2);
        iPrNW  = iPrN    + (ncx-1);
        iPrSW  = iPrS    + (ncx-1);
        iVzW   = iVzC    + nxvz - 3;// Nx(nx-2); //+ (nxvz-1-1)
        iVzSW  = iVzW    - nxvz;
        iVzNW  = iVzW    + nxvz;
        //        ixyW   = ixyW    + (nx-1); // just added
        //        printf("periodic left\n");
    }
    if (k==nx-1 && BCv_t[iVzE]==-12) {
        iVxNEE = iVxNE - (nx-2);
        iVxSEE = iVxSE - (nx-2);
        iPrNE  = iPrN  - (ncx-1);
        iPrSE  = iPrS  - (ncx-1);
        iVzE   = iVzC  - nxvz + 3;
        iVzNE  = iVzE  + nxvz;
        iVzSE  = iVzE  - nxvz;
        //        ixyE   = ixyE  - (nx-1);
        //        printf("periodic right\n");
    }
   
    // The computation of FD coefficients is only useful for the purpose of the stiffness/Jacobian matrix assembly
    double  vC=0.0;
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0,   vS=0.0,   vW=0.0,   vE=0.0,   vN=0.0, pN=0.0, pS=0.0;
    double vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, uSSW=0.0, uSSE=0.0, uNNW=0.0, uNNE=0.0;
    double uSWW=0.0, uSEE=0.0, uNWW=0.0, uNEE=0.0;
    double pSW=0.0, pSE=0.0, pNW=0.0, pNE=0.0;
    
    double D31W  = mesh->D31_s[ixyW];
    double D32W  = mesh->D32_s[ixyW];
    double D33W  = mesh->D33_s[ixyW];
    double D34W  = mesh->D34_s[ixyW];
    
    double D31E  = mesh->D31_s[ixyE];
    double D32E  = mesh->D32_s[ixyE];
    double D33E  = mesh->D33_s[ixyE];
    double D34E  = mesh->D34_s[ixyE];
    
    double D21S  = mesh->D21_n[iPrS];
    double D22S  = mesh->D22_n[iPrS];
    double D23S  = mesh->D23_n[iPrS];
    double D24S  = mesh->D24_n[iPrS];
    
    double D21N  = mesh->D21_n[iPrN];
    double D22N  = mesh->D22_n[iPrN];
    double D23N  = mesh->D23_n[iPrN];
    double D24N  = mesh->D24_n[iPrN];
    
    double inS=0.0, inN=0.0, inW=0.0, inE = 0.0, inWv = 0.0, inEv = 0.0;
    
    if (BCp_t[iPrS] == -1) inS = 1.0;
    if (BCp_t[iPrN] == -1) inN = 1.0;
    
    //    if (BCg_t[ixyW] != 30 && BCv_t[iVzW] != 13 ) inW = 1.0; // !!!!!!!!!!!!!!!!!!!!
    //    if (BCg_t[ixyE] != 30 && BCv_t[iVzE] != 13 ) inE = 1.0;
    if ( BCv_t[iVzW] != 30 && BCv_t[iVzW] != 13 ) inW  = 1.0;
    if ( BCv_t[iVzE] != 30 && BCv_t[iVzE] != 13 ) inE  = 1.0;
    if ( BCg_t[ixyW] != 30 ) inWv = 1.0;
    if ( BCg_t[ixyE] != 30 ) inEv = 1.0;
    
    // XTRA
    double wS=0.0, wN=0.0, wE=0.0, wW=0.0;
    double inSWc=0.0,inSEc=0.0,inNWc=0.0,inNEc=0.0;
    double inSWv=0.0,inSEv=0.0,inNWv=0.0,inNEv=0.0;
  
    if ( (k>1) || (k==1 && BCv_t[iVzW]==-1) ) {
        if (BCp_t[iPrSW] == -1) inSWc = 1.0;
        if (BCp_t[iPrNW] == -1) inNWc = 1.0;
    }
    
    if ( (k<nxvz-2) || (k==nxvz-2 && BCv_t[iVzE]==-1 ) ) {
        if (BCp_t[iPrSE] == -1) inSEc = 1.0;
        if (BCp_t[iPrNE] == -1) inNEc = 1.0;
    }
    
    if (l>0) {
        if( BCg_t[ixySW] != 30) inSWv = 1.0;
        if( BCg_t[ixySE] != 30) inSEv = 1.0;
    }
    
    if (l<nz-1) {
        if( BCg_t[ixyNW] != 30) inNWv = 1.0;
        if( BCg_t[ixyNE] != 30) inNEv = 1.0;
    }
    
    // New stuff
    double inSSW = 0.0, inSSE = 0.0;
    if (BCu_t[iVxSSW] == -1 || periodix==1) inSSW = 1.0;
    if (BCu_t[iVxSSE] == -1 || periodix==1) inSSE = 1.0;
    
    double inNNW = 0.0, inNNE = 0.0;
    if (BCu_t[iVxNNW] == -1 || periodix==1) inNNW = 1.0;
    if (BCu_t[iVxNNE] == -1 || periodix==1) inNNE = 1.0;
    
    wE = inN + inS + inNEc + inSEc;
    wW = inN + inS + inNWc + inSWc;
    wS = inW + inE + inSWv + inSEv;
    wN = inW + inE + inNWv + inNEv;
    
    if (wW>1.0) wW = 1.0/wW;
    if (wE>1.0) wE = 1.0/wE;
    if (wS>1.0) wS = 1.0/wS;
    if (wN>1.0) wN = 1.0/wN;
            
    // FD Coefficients obtained using AssembleGeneralStiffness_MDOODZ_6.0-simpler-plastic.ipynb
    vW = (-D33W*dz*pow(inW, 2)*inWv + 0.25*dx*pow(inW, 2)*inWv*(D23N*inN - D23S*inS) + (1.0/3.0)*dx*wW*(inNWc - inSWc)*(D31W*comp*out_of_plane + D32W*(comp*out_of_plane - 3)))/(pow(dx, 2)*dz);
    vC = (1.0/3.0)*(-dx*(dx*(D21N*comp*inN*out_of_plane + D21S*comp*inS*out_of_plane + D22N*inN*(comp*out_of_plane - 3) + D22S*inS*(comp*out_of_plane - 3)) + 0.75*dz*(-D23N*inN + D23S*inS)*(pow(inE, 2)*inEv - pow(inW, 2)*inWv)) + dz*(dx*(inN - inS)*(-D31E*comp*out_of_plane*wE + D31W*comp*out_of_plane*wW - D32E*wE*(comp*out_of_plane - 3) + D32W*wW*(comp*out_of_plane - 3)) + 3*dz*(D33E*pow(inE, 2)*inEv + D33W*pow(inW, 2)*inWv)))/(pow(dx, 2)*pow(dz, 2));
    vE = (-D33E*dz*pow(inE, 2)*inEv + 0.25*dx*pow(inE, 2)*inEv*(-D23N*inN + D23S*inS) - 1.0/3.0*dx*wE*(inNEc - inSEc)*(D31E*comp*out_of_plane + D32E*(comp*out_of_plane - 3)))/(pow(dx, 2)*dz);
    vS = (1.0/3.0)*inS*(-0.75*D23S*dz*(inE*inSEv - inSWv*inW) + dx*(D21S*comp*out_of_plane + D22S*(comp*out_of_plane - 3)) + dz*(-D31E*comp*out_of_plane*wE + D31W*comp*out_of_plane*wW - D32E*wE*(comp*out_of_plane - 3) + D32W*wW*(comp*out_of_plane - 3)))/(dx*pow(dz, 2));
    vN = (1.0/3.0)*inN*(0.75*D23N*dz*(inE*inNEv - inNWv*inW) + dx*(D21N*comp*out_of_plane + D22N*(comp*out_of_plane - 3)) + dz*(D31E*comp*out_of_plane*wE - D31W*comp*out_of_plane*wW + D32E*wE*(comp*out_of_plane - 3) - D32W*wW*(comp*out_of_plane - 3)))/(dx*pow(dz, 2));
    uSW = (1.0/3.0)*(dx*(0.75*dx*(D23N*inN*inW*inWv + D23S*inS*(inSSW*inSWv - inW*inWv)) + dz*inS*(D21S*(comp*out_of_plane - 3) + D22S*comp*out_of_plane)) - dz*(3*D33W*dx*inW*inWv + dz*(D31E*inS*wE*(comp*out_of_plane - 3) - D31W*wW*(inS - inSWc)*(comp*out_of_plane - 3) + D32E*comp*inS*out_of_plane*wE - D32W*comp*out_of_plane*wW*(inS - inSWc))))/(pow(dx, 2)*pow(dz, 2));
    uSE = (1.0/3.0)*(-dx*(0.75*dx*(-D23N*inE*inEv*inN + D23S*inS*(inE*inEv - inSEv*inSSE)) + dz*inS*(D21S*(comp*out_of_plane - 3) + D22S*comp*out_of_plane)) + dz*(3*D33E*dx*inE*inEv + dz*(D31E*wE*(inS - inSEc)*(comp*out_of_plane - 3) - D31W*inS*wW*(comp*out_of_plane - 3) + D32E*comp*out_of_plane*wE*(inS - inSEc) - D32W*comp*inS*out_of_plane*wW)))/(pow(dx, 2)*pow(dz, 2));
    uNW = (1.0/3.0)*(dx*(0.75*dx*(D23N*inN*(inNNW*inNWv - inW*inWv) + D23S*inS*inW*inWv) - dz*inN*(D21N*(comp*out_of_plane - 3) + D22N*comp*out_of_plane)) + dz*(3*D33W*dx*inW*inWv + dz*(-D31E*inN*wE*(comp*out_of_plane - 3) + D31W*wW*(inN - inNWc)*(comp*out_of_plane - 3) - D32E*comp*inN*out_of_plane*wE + D32W*comp*out_of_plane*wW*(inN - inNWc))))/(pow(dx, 2)*pow(dz, 2));
    uNE = (1.0/3.0)*(dx*(0.75*dx*(-D23N*inN*(inE*inEv - inNEv*inNNE) + D23S*inE*inEv*inS) + dz*inN*(D21N*(comp*out_of_plane - 3) + D22N*comp*out_of_plane)) - dz*(3*D33E*dx*inE*inEv + dz*(-D31E*wE*(inN - inNEc)*(comp*out_of_plane - 3) + D31W*inN*wW*(comp*out_of_plane - 3) - D32E*comp*out_of_plane*wE*(inN - inNEc) + D32W*comp*inN*out_of_plane*wW)))/(pow(dx, 2)*pow(dz, 2));
    vSW = (1.0/3.0)*(-0.75*D23S*inS*inSWv*inW + D31W*comp*inSWc*out_of_plane*wW + D32W*inSWc*wW*(comp*out_of_plane - 3))/(dx*dz);
    vSE = (1.0/3.0)*(0.75*D23S*inE*inS*inSEv - D31E*comp*inSEc*out_of_plane*wE - D32E*inSEc*wE*(comp*out_of_plane - 3))/(dx*dz);
    vNW = (1.0/3.0)*(0.75*D23N*inN*inNWv*inW - D31W*comp*inNWc*out_of_plane*wW - D32W*inNWc*wW*(comp*out_of_plane - 3))/(dx*dz);
    vNE = (1.0/3.0)*(-0.75*D23N*inE*inN*inNEv + D31E*comp*inNEc*out_of_plane*wE + D32E*inNEc*wE*(comp*out_of_plane - 3))/(dx*dz);
    uSWW = (1.0/3.0)*inSWc*wW*(D31W*(comp*out_of_plane - 3) + D32W*comp*out_of_plane)/pow(dx, 2);
    uSEE = (1.0/3.0)*inSEc*wE*(D31E*(comp*out_of_plane - 3) + D32E*comp*out_of_plane)/pow(dx, 2);
    uNWW = (1.0/3.0)*inNWc*wW*(D31W*(comp*out_of_plane - 3) + D32W*comp*out_of_plane)/pow(dx, 2);
    uNEE = (1.0/3.0)*inNEc*wE*(D31E*(comp*out_of_plane - 3) + D32E*comp*out_of_plane)/pow(dx, 2);
    uSSW = -0.25*D23S*inS*inSSW*inSWv/pow(dz, 2);
    uSSE = -0.25*D23S*inS*inSEv*inSSE/pow(dz, 2);
    uNNW = -0.25*D23N*inN*inNNW*inNWv/pow(dz, 2);
    uNNE = -0.25*D23N*inN*inNEv*inNNE/pow(dz, 2);
    pS  = -inS*one_dz + (inS*(D24S*dx + dz*(-D34E*inE*wE + D34W*inW*wW))/(dx*dz));
    pN  =  inN*one_dz + (inN*(-D24N*dx + dz*(-D34E*inE*wE + D34W*inW*wW))/(dx*dz));
    pSW = (D34W*inSWc*inW*wW/dx);
    pSE = (-D34E*inE*inSEc*wE/dx);
    pNW = (D34W*inNWc*inW*wW/dx);
    pNE = (-D34E*inE*inNEc*wE/dx);
    
    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodz  = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
        vC_corr = 1.00 * om * model.dt * mesh->gz[c3] * drhodz;
        // Importante trique, voire meme gigantesque!
        if (vC+vC_corr<0.0) vC_corr = 0.0;
        vC += vC_corr;
    }
    
    // Add contribution from non-conforming Dirichlets
    if ( BCv_t[iVzW]   == 11 ) vC  -=  vW ;
    if ( BCv_t[iVzE]   == 11 ) vC  -=  vE ;
    if ( BCv_t[iVzSW]  == 11 ) vW  -=  vSW;
    if ( BCv_t[iVzSE]  == 11 ) vE  -=  vSE;
    if ( BCv_t[iVzNW]  == 11 ) vW  -=  vNW;
    if ( BCv_t[iVzNE]  == 11 ) vE  -=  vNE;
    if ( BCu_t[iVxNNW] == 11 ) uNW -= uNNW;
    if ( BCu_t[iVxNNE] == 11 ) uNE -= uNNE;
    if ( BCu_t[iVxSSW] == 11 ) uSW -= uSSW;
    if ( BCu_t[iVxSSE] == 11 ) uSE -= uSSE;
        
    //--------------------
        
    if ( Assemble == 1 ) {
        const bool SetVzSW = BCv_t[iVzSW] != 30 && iz>0;
        const bool SetVzS  = BCv_t[iVzS]  != 30 && iz>0;
        const bool SetVzSE = BCv_t[iVzSE] != 30 && iz>0;
        const bool SetVzW  = BCv_t[iVzW]  != 30 ;
        const bool SetVzC  = BCv_t[iVzC]  != 30;
        const bool SetVzE  = BCv_t[iVzE]  != 30 ;
        const bool SetVzNW = BCv_t[iVzNW] != 30 && iz<nz-1;
        const bool SetVzN  = BCv_t[iVzN]  != 30 && iz<nz-1;
        const bool SetVzNE = BCv_t[iVzNE] != 30 && iz<nz-1;
        //--------------------
        const bool SetVxSW = BCu_t[iVxSW] != 30 && iz>0   ;
        const bool SetVxSE = BCu_t[iVxSE] != 30 && iz>0   ;
        const bool SetVxNW = BCu_t[iVxNW] != 30 && iz<nz-1;
        const bool SetVxNE = BCu_t[iVxNE] != 30 && iz<nz-1;
        //--------------------
        const bool SetPrS  = BCp_t[iPrS] != 30 && iz>0;
        const bool SetPrN  = BCp_t[iPrN] != 30 && iz<nz-1;
        //--------------------
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;

        if ( BCu_t[iVxSSW]  != 30   ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSSW],   &(nnzc2A[ith]), uSSW*celvol, BCu_t[iVxSSW],    mesh->BCu.val[iVxSSW],    StokesA->bbc );
        if ( BCu_t[iVxSSE]  != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSSE],     &(nnzc2A[ith]), uSSE*celvol, BCu_t[iVxSSE],    mesh->BCu.val[iVxSSE],    StokesA->bbc );

        // uSWW
        if ( (k>1 ) || (k==1 && BCu_t[iVxSWW]==-1) ) {
            if( BCu_t[iVxSWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSWW],     &(nnzc2A[ith]), uSWW*celvol, BCu_t[iVxSWW],    mesh->BCu.val[iVxSWW],    StokesA->bbc );
        }
        //
        // uSW && uSE
        if ( BCu_t[iVxSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, BCu_t[iVxSW], mesh->BCu.val[iVxSW], StokesA->bbc );
        if ( BCu_t[iVxSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, BCu_t[iVxSE], mesh->BCu.val[iVxSE], StokesA->bbc );
        //
        // uSEE
        if ( (k<nx-1 ) || (k==nx-1 && BCu_t[iVxSEE]==-1) ) {
            if( BCu_t[c1+1] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSEE],     &(nnzc2A[ith]), uSEE*celvol, BCu_t[iVxSEE],    mesh->BCu.val[iVxSEE],    StokesA->bbc );
        }
        
        // uNWW
        if ( (k>1    ) || (k==1   && BCu_t[iVxNWW]==-1) ) {
            if( BCu_t[c1+nx-2] != 30)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNWW],     &(nnzc2A[ith]), uNWW*celvol, BCu_t[iVxNWW],    mesh->BCu.val[iVxNWW],    StokesA->bbc );
        }
        
        if ( BCu_t[iVxNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, BCu_t[iVxNW], mesh->BCu.val[iVxNW], StokesA->bbc );
        if ( BCu_t[iVxNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, BCu_t[iVxNE], mesh->BCu.val[iVxNE], StokesA->bbc );
        
        // uNEE
        if ( (k<nx-1 ) || (k==nx-1 && BCu_t[iVxNEE]==-1) ) {
            if( BCu_t[iVxNEE] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNEE],     &(nnzc2A[ith]), uNEE*celvol, BCu_t[iVxNEE],    mesh->BCu.val[iVxNEE],    StokesA->bbc );
        }
     
        if ( BCu_t[iVxNNW] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNNW],  &(nnzc2A[ith]), uNNW*celvol, BCu_t[iVxNNW], mesh->BCu.val[iVxNNW], StokesA->bbc );    
        if ( BCu_t[iVxNNE] != 30 )  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNNE],  &(nnzc2A[ith]), uNNE*celvol, BCu_t[iVxNNE], mesh->BCu.val[iVxNNE], StokesA->bbc );
    
        
        //--------------------
        
        if (SetVzSW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW],     &(nnzc2A[ith]), vSW*celvol, BCv_t[iVzSW],      mesh->BCv.val[iVzSW],    StokesA->bbc );
        if (SetVzS)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzS],  &(nnzc2A[ith]), vS*celvol, BCv_t[iVzS], mesh->BCv.val[iVzS], StokesA->bbc );
        if (SetVzSE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE],     &(nnzc2A[ith]), vSE*celvol, BCv_t[iVzSE],      mesh->BCv.val[iVzSE],    StokesA->bbc );
        if (SetVzW)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzW],     &(nnzc2A[ith]), vW*celvol, BCv_t[iVzW],    mesh->BCv.val[iVzW],    StokesA->bbc );        
        if (SetVzC)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, BCv_t[iVzC],      mesh->BCv.val[iVzC],      StokesA->bbc );
        if (SetVzE)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzE],     &(nnzc2A[ith]), vE*celvol, BCv_t[iVzE],    mesh->BCv.val[iVzE],    StokesA->bbc );
        if (SetVzNW) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],     &(nnzc2A[ith]), vNW*celvol, BCv_t[iVzNW],      mesh->BCv.val[iVzNW],    StokesA->bbc );
        if (SetVzN)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzN],  &(nnzc2A[ith]), vN*celvol, BCv_t[iVzN], mesh->BCv.val[iVzN], StokesA->bbc );
        if (SetVzNE) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE],     &(nnzc2A[ith]), vNE*celvol, BCv_t[iVzNE],      mesh->BCv.val[iVzNE],    StokesA->bbc );
        
        //--------------------
        
        // pSW
        if ( (Newton==1 && k>1) || (Newton==1 && k==1 && periodix==1) ) {
            if ( BCp_t[iPrSW] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSW] - Stokes->neq_mom,      &(nnzc2B[ith]), pSW*celvol, BCp_t[iPrSW],     mesh->BCp.val[iPrSW],     StokesB->bbc );
            }
        }
        
        // pS
        if ( BCp_t[iPrS] != 30 && BCp_t[iPrN] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrS] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, BCp_t[iPrS],     mesh->BCp.val[iPrS],     StokesB->bbc );
        }
        
        // pSE
        if ( (Newton==1 && k<nx-1) || (Newton==1 && k==nx-1 && periodix==1) ) {
            if ( BCp_t[iPrSE] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSE] - Stokes->neq_mom,      &(nnzc2B[ith]), pSE*celvol, BCp_t[iPrSE],     mesh->BCp.val[iPrSE],     StokesB->bbc );
            }
        }
        
        // pNW
        if ( (Newton==1 && k>1) || (Newton==1 && k==1 && periodix==1)  ) {
            if ( BCp_t[iPrNW] != 30) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNW] - Stokes->neq_mom,  &(nnzc2B[ith]), pNW*celvol, BCp_t[iPrNW], mesh->BCp.val[iPrNW], StokesB->bbc );
            }
        }
        
        // pN
        if ( BCp_t[iPrS] != 30 && BCp_t[iPrN] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrN] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, BCp_t[iPrN], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
        
        if ( (Newton==1 && k<nx-1) || (Newton==1 && k==nx-1 && periodix==1) ) {
            // pNE
            if ( BCp_t[iPrNE] != 30) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNE] - Stokes->neq_mom,  &(nnzc2B[ith]), pNE*celvol, BCp_t[iPrNE], mesh->BCp.val[iPrNE], StokesB->bbc );
            }
        }
    }
    else {
        
        // Residual function
        double F_Reg = 0., F_S =0., F_N =0., SzzS, SzzN;

        // Residual: regular nodes
        F_Reg += (mesh->szzd[iPrN]   - mesh->szzd[iPrS])/dz;
        F_Reg -= (mesh->p_corr[iPrN] - mesh->p_corr[iPrS])/dz;
        F_Reg += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;

        // Residual: stress BC south
        SzzN = mesh->szzd[iPrN] - mesh->p_corr[iPrN];
        SzzS = 2.0*mesh->BCv.val[iVzC] - SzzN;
        F_S += (SzzN - SzzS)/dz;
        F_S += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        F_S *= 0.5;

        // Residual: stress BC north
        SzzS = mesh->szzd[iPrS] - mesh->p_corr[iPrS];
        SzzN = 2.0*mesh->BCv.val[iVzC] - SzzS;
        F_N += (SzzN - SzzS)/dz;
        F_N += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        F_N *= 0.5;
        
        // Residual
        StokesA->F[eqn]  = (1.0 - SzzBCS - SzzBCN)*F_Reg + SzzBCS*F_S + SzzBCN*F_N;
        StokesA->F[eqn] += StokesA->b[eqn] - vC_corr*v[iVzC]; 
        StokesA->F[eqn] *= -celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Continuity_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesC, SparseMat *StokesD, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempC, double **AtempC, int *nnzc2C, int **JtempD, double **AtempD, int *nnzc2D, int i, int j ) {
    
    double pc=0.0, uW=0.0, uE=0.0, vN=0.0, vS=0.0, oop_fact = 1.0, comp_fact = 0.0;

    // const double rhoW = 0.5*(mesh->rho_n[c1]   + mesh->rho_n[c1-nx]);
    // const double rhoE = 0.5*(mesh->rho_n[c1+1] + mesh->rho_n[c1-nx+1]);
    // const double rhoS = 0.5*(mesh->rho_n[c1-nx]+ mesh->rho_n[c1-nx+1]);
    // const double rhoN = 0.5*(mesh->rho_n[c1]   + mesh->rho_n[c1+1]);
    // const double drhodx = (rhoE-rhoW)*one_dx;
    // const double drhody = (rhoN-rhoS)*one_dz;
    // const double VxC = 0.5*(mesh->u_in[c1]   + mesh->u_in[c1+1]);
    // const double VzC = 0.5*(mesh->v_in[c3]   + mesh->v_in[c3+nxvz]);
    const double adv_rho = 1.0; // + model.dt*VxC*drhodx/mesh->rho_n[c2] + model.dt*VzC*drhodx/mesh->rho_n[c2];

    // if (fabs(adv_rho)>1+1e-10) printf("%2.2e %2.2e\n ", model.dt*VxC*drhodx/mesh->rho_n[c2] + model.dt*VzC*drhodx/mesh->rho_n[c2], drhodx);
    
    // Non-zero out-of-plane strain
    if ( model.out_of_plane == 1 ) oop_fact = 3.0/2.0;

    // Compressibility
    if ( comp == 1 ) comp_fact = 1.0;
    pc = comp_fact*mesh->bet_n[c2]/model.dt;
    
    // div u
    if ( mesh->BCu.type[c1     ] != 13 ) uW = -adv_rho*one_dx * oop_fact;
    if ( mesh->BCu.type[c1+1   ] != 13 ) uE =  adv_rho*one_dx * oop_fact;
    if ( mesh->BCv.type[c3     ] != 13 ) vS = -adv_rho*one_dz * oop_fact;
    if ( mesh->BCv.type[c3+nxvz] != 13 ) vN =  adv_rho*one_dz * oop_fact;
    
    // Stencil assembly / residual
    if ( Assemble == 1 ) {
        StokesC->b[eqn] *= celvol;
        StokesD->b[eqn] *= celvol;
        if ( mesh->BCu.type[c1     ] != 13 )    AddCoeff3( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_u[c1],      &(nnzc2C[ith]), uW*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesC->bbc );
        if ( mesh->BCu.type[c1+1   ] != 13 )    AddCoeff3( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_u[c1+1],    &(nnzc2C[ith]), uE*celvol, mesh->BCu.type[c1+1],    mesh->BCu.val[c1+1],    StokesC->bbc );
        if ( mesh->BCv.type[c3     ] != 13 )    AddCoeff3( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_v[c3],      &(nnzc2C[ith]), vS*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesC->bbc );
        if ( mesh->BCv.type[c3+nxvz] != 13 )    AddCoeff3( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_v[c3+nxvz], &(nnzc2C[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesC->bbc );
    }
    else {
        
        
        if ( model.density_variations == 0 )  {
            // beta * dP / dt
            StokesC->F[eqn] = comp_fact*mesh->bet_n[c2]*(mesh->p_in[c2] -  mesh->p0_n[c2] ) / model.dt + mesh->div_u[c2];
        }
        if ( model.density_variations == 1 )  {
            // d (ln(rho)) / dt
            if ( model.density_variations == 1 )  StokesC->F[eqn] = comp_fact*( log(mesh->rho_n[c2]) -  log(mesh->rho0_n[c2]) ) / model.dt + mesh->div_u[c2];
            // if ( model.density_variations == 1 )  StokesC->F[eqn] = comp_fact*( mesh->rho_n[c2] -  mesh->rho0_n[c2] ) / model.dt +  mesh->rho_n[c2]*mesh->div_u[c2];

        }
        StokesC->F[eqn]-= StokesC->b[eqn];
        StokesC->F[eqn] *= celvol;
        StokesD->F[eqn] = StokesC->F[eqn] ;        
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MergeParallelMatrixDecoupled( SparseMat *Stokes, double **Atemp, int **Jtemp, int **Itemp, grid* mesh, int* estart, int* eend, int *nnzc, int *nnzc2, int *last_eqn, int n_th, char* BC_type, int *eqn_number ) {
    
    int eqn, ith, begin[n_th], end[n_th], c, last_eqn_serial=0, k, l;
    
    for (ith=0; ith<n_th; ith++) {
        
        // Get last equation number
        if (last_eqn[ith]>last_eqn_serial) {
            last_eqn_serial = last_eqn[ith];
        }
        
        if (ith==0) {
            begin[ith] = *(nnzc);
            end[ith]   = *(nnzc) + nnzc2[ith]-1;
        }
        else {
            begin[ith] = end[ith-1]+1;
            end[ith]   = end[ith-1]+1+nnzc2[ith]-1;
        }
    }
    
    // Recombine A and J
#pragma omp parallel private(ith, eqn, c, k, l ) shared( Stokes, Atemp, Jtemp, Itemp, begin, BC_type, eqn_number )
    {
        ith = omp_get_thread_num();
        
        for( c=estart[ith]; c<eend[ith]+1; c++) {
            
            // Get equation numbering
            //            if ( BC_type[c] != 0 && BC_type[c] != 30 && BC_type[c] != 31 && BC_type[c] != 13 && BC_type[c] != 11 ) {
            if ( BC_type[c] != 0 && BC_type[c] != 30 && BC_type[c] != 31 && BC_type[c] != 13 && BC_type[c] != 11 && BC_type[c] != -12  ) {
                eqn = eqn_number[c] -  Stokes->neq_mom;
                Stokes->Ic[eqn] = Itemp[ith][eqn] + begin[ith];
            }
        }
        
        for (k=0; k<nnzc2[ith]; k++) {
            l = begin[ith] + k;
            Stokes->A[l] = Atemp[ith][k];
            Stokes->J[l] = Jtemp[ith][k];
        }
    }
    
    // Update total non-zero number
    for (ith=0; ith<n_th; ith++) {
        *(nnzc) += nnzc2[ith];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeDecompositionDecoupled( int *estart, int *eend, int *DD, int *last_eqn  ) {
    DoodzFree( DD );
    DoodzFree( estart );
    DoodzFree( eend );
    DoodzFree( last_eqn );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateDecompositionDecoupled( int n_th, int **estart, int **eend, int **DD, int **last_eqn ) {
    
    /* Array containing the size of the block for each thread */
    *DD       = DoodzMalloc( sizeof(int) * n_th );
    *estart   = DoodzMalloc( sizeof(int) * n_th );
    *eend     = DoodzMalloc( sizeof(int) * n_th );
    *last_eqn = DoodzCalloc( n_th, sizeof(int) );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DomainDecompositionDecoupled( int N, int n_th, int *estart, int *eend, int *DD ) {
    
    int ith, eps, block_size, x0 = 0, x1 = 0;
    eps = N%n_th;
    block_size  = (N - eps) / n_th;
    
    // Divide domain into blocks
    for (ith=0; ith<n_th; ith++) {
        if (ith<n_th-1)
            DD[ith] = block_size;
        else {
            DD[ith] = block_size + eps;
        }
        
        x1 = x0 + DD[ith];
        estart[ith] = x0;
        eend[ith]   = x1-1;
        x0          = x1;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateTempMatArraysDecoupled(double ***Atemp, int ***Itemp, int ***Jtemp, int n_th, int nnz, int neq, int* DD, int **nnzc2 ) {
    
    int ith;
    
    // Temporary parallel Matrix arrays
    *Atemp = DoodzMalloc( sizeof(double*) * n_th );
    *Itemp = DoodzMalloc( sizeof(int*) * n_th );
    *Jtemp = DoodzMalloc( sizeof(int*) * n_th );
    *nnzc2 = DoodzCalloc( n_th, sizeof(int) );
    
    for (ith=0; ith<n_th; ith++) {
        (*Atemp)[ith] = DoodzMalloc( sizeof(double) * (int)(nnz/n_th) );
        (*Itemp)[ith] = DoodzCalloc( (neq+1), sizeof(int) );
        (*Jtemp)[ith] = DoodzMalloc( sizeof(int) * (int)(nnz/n_th) );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeTempMatArraysDecoupled(double **Atemp, int **Itemp, int **Jtemp, int n_th, int *nnzc2 ) {
    
    int ith;
    
    // Temporary parallel Matrix arrays
    for (ith=0; ith<n_th; ith++) {
        DoodzFree( Atemp[ith] );
        DoodzFree( Itemp[ith] );
        DoodzFree( Jtemp[ith] );
    }
    DoodzFree( Atemp );
    DoodzFree( Itemp );
    DoodzFree( Jtemp );
    DoodzFree( nnzc2 );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildStokesOperatorDecoupled( grid *mesh, params model, int lev, double *p_corr, double *p, double *u, double *v, SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, int Assemble ) {
    
    // FLAG
    // Assemble=1: Build the linear system of discrete equations
    // Assemble=0: Evaluate Stokes residuals
    
    int    cc, k, l, c1, c2, c3, nx=mesh->Nx, nz=mesh->Nz, nxvz=nx+1, nzvx=nz+1, ncx=nx-1, ncz=nz-1;
    
    // Pre-calculate FD coefs
    double celvol    = mesh->dx*mesh->dz;
    double one_dx    = 1.0/mesh->dx;
    double one_dz    = 1.0/mesh->dz;
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    double one_dx_dz = 1.0/mesh->dx/mesh->dz;
    
    // Switches
    int eqn,  sign = 1, comp = 0, stab = 0;
    if (model.free_surface_stab>0) stab = 1;
    if (model.compressible==1  ) comp = 1;
    double theta = model.free_surface_stab;
        
    // Decompose domain
    int n_th, N, ith;
    int *DD, *estart, *eend, *last_eqn;
    
    // Temporary parallel matrix
    double **AtempA;
    int **JtempA, **ItempA;
    double **AtempB;
    int **JtempB, **ItempB;
    double **AtempC;
    int **JtempC, **ItempC;
    double **AtempD;
    int **JtempD, **ItempD;
    
    int nnzA, nnzB, nnzC, nnzD, nnzcA=0, nnzcB=0, nnzcC=0, nnzcD=0;
    int *nnzc2A, *nnzc2B, *nnzc2C, *nnzc2D;
    
#pragma omp parallel shared(n_th)
    {
        n_th = omp_get_num_threads();
    }
#pragma omp barrier

    if (Assemble==1) printf("Assemble Picard operator on %d threads\n", n_th);
    
    // Matrix initialisation
    if ( Assemble == 1 ) {
        nnzA  = 9*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzB  = 5*((mesh->Nx-1) * (mesh->Nz-1));
        nnzC  = nnzB;
        nnzD  = 1*((mesh->Nx-1) * (mesh->Nz-1));
        StokesA->neq = Stokes->neq_mom;
        StokesB->neq = Stokes->neq_mom;
        StokesC->neq = Stokes->neq_cont; StokesC->neq_mom = Stokes->neq_mom;
        StokesD->neq = Stokes->neq_cont; StokesD->neq_mom = Stokes->neq_mom;
        //printf("Assembling  decoupled Stokes matrix...\n");
        AllocMat( StokesA, nnzA );
        AllocMat( StokesB, nnzB );
        AllocMat( StokesC, nnzC );
        AllocMat( StokesD, nnzD );
    }
    
    // Build velocity block RHS
    int inc=0;
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
                StokesA->b[inc] = mesh->roger_x[cc];
                inc++;
            }
        }
    }
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13  && mesh->BCv.type[cc] != -12  ) {
                StokesA->b[inc] = mesh->roger_z[cc];
                inc++;
            }
        }
    }
    
    // Build pressure block RHS
    inc=0;
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            cc = k + l*ncx;
            if ( mesh->BCp.type[cc] != 31 && mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30 ) {
                StokesC->b[inc] = mesh->rhs_p[cc];
                inc++;
            }
        }
    }
    
    //------------------------------------------------------------------//
    //------------------------- U-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nx*nzvx;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }

    // double SUM=0.;
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p_corr, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, comp )
    {
        ith = omp_get_thread_num();
        
        for( c1=estart[ith]; c1<eend[ith]+1; c1++) {
            
            k      = mesh->kvx[c1];
            l      = mesh->lvx[c1];
            c2     = k-1 + (l-1)*ncx;
            c3     = k   + l*nxvz;
            
            // Get equation numbering
            if ( mesh->BCu.type[c1]==-1 || mesh->BCu.type[c1]==2 || mesh->BCu.type[c1]==-2 ) {
                eqn = Stokes->eqn_u[c1];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if (  k>=0 && k<=nx-1 && l>0 && l<nzvx-1 ) {
                    // Maybe one should avoid sides (k>0 && k<nxvz-1) if  mesh->BCu.type[c1] != 2 ?? 
                    Xmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, k, l, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, nz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
            }
        }
    }

    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- V-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nxvz*nz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p_corr, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c3=estart[ith]; c3<eend[ith]+1; c3++) {
            
            k  = mesh->kvz[c3];
            l  = mesh->lvz[c3];
            c1 = k   + l*nx;
            c2 = k-1 + (l-1)*ncx;
            
            // Get equation numbering
            if ( mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3]==2 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>=0 && l<=nz-1 ) {
                    // Maybe one should avoid sides (l>0 && l<nz-1) if  mesh->BCv.type[c3] != 2 ?? 
                    Zmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, k, l, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, nz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- Continuity -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = ncx*ncz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempC, &ItempC, &JtempC, n_th, nnzC, Stokes->neq_cont, DD, &nnzc2C  );
        AllocateTempMatArraysDecoupled( &AtempD, &ItempD, &JtempD, n_th, nnzD, Stokes->neq_cont, DD, &nnzc2D  );
    }
       
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesC, StokesD, u, v, p, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn, comp ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, nx, ncx, nxvz, nzvx )
    {   
        ith = omp_get_thread_num();
        
        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
            
            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c1  = k   + (l+1)*nx;
            c3  = k   + l*nxvz + 1;
            
            //--------------------- INNER NODES ---------------------//
            if ( mesh->BCp.type[c2] == -1) {
                
                comp = mesh->comp_cells[c2];
                
                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
                last_eqn[ith]   = eqn ;
                
                if ( Assemble == 1 ) {
                    ItempC[ith][eqn] = nnzc2C[ith];
                    ItempD[ith][eqn] = nnzc2D[ith]; //printf("%d ",  nnzc2D[ith]);
                }
                //
                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrixDecoupled( StokesC, AtempC, JtempC, ItempC, mesh, estart, eend, &nnzcC, nnzc2C, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
        MergeParallelMatrixDecoupled( StokesD, AtempD, JtempD, ItempD, mesh, estart, eend, &nnzcD, nnzc2D, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
    }
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempC, ItempC, JtempC, n_th, nnzc2C );
        FreeTempMatArraysDecoupled( AtempD, ItempD, JtempD, n_th, nnzc2D );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //----------------------------- End --------------------------------//
    //------------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        
        // Add contribution from the BC's
        for (k=0; k<StokesA->neq; k++) {
            StokesA->b[k] += StokesA->bbc[k];
            StokesB->b[k] = StokesA->b[k];
        }
        
        //            MinMaxArray(StokesA->b, 1, StokesA->neq, "rhs_mom_here" );
        //            MinMaxArray(StokesC->b, 1, StokesC->neq, "rhs_cont_here" );
        
        
        // Add contribution from the BC's
        for (k=0; k<StokesC->neq; k++) {
            StokesC->b[k] += StokesC->bbc[k];
            StokesD->b[k] = StokesC->b[k];
        }
        
        // Final index
        StokesA->Ic[StokesA->neq] = nnzcA;
        StokesB->Ic[StokesB->neq] = nnzcB;
        StokesC->Ic[StokesC->neq] = nnzcC;
        StokesD->Ic[StokesD->neq] = nnzcD;
        
        //        for (ieq=0; ieq<Stokes->neq_cont+1; ieq++) printf( "%d ", StokesC->I[ieq] );
        
        StokesA->nnz = nnzcA;
        StokesB->nnz = nnzcB;
        StokesC->nnz = nnzcC;
        StokesD->nnz = nnzcD;
        
        // Resize arrays to proper number of non-zeros
        double *bufd;
        int *bufi;
        bufi      = DoodzRealloc(StokesA->J, nnzcA*sizeof(int));
        bufd      = DoodzRealloc(StokesA->A, nnzcA*sizeof(double));
        StokesA->J = bufi;
        StokesA->A = bufd;
        
        bufi      = DoodzRealloc(StokesB->J, nnzcB*sizeof(int));
        bufd      = DoodzRealloc(StokesB->A, nnzcB*sizeof(double));
        StokesB->J = bufi;
        StokesB->A = bufd;
        
        bufi      = DoodzRealloc(StokesC->J, nnzcC*sizeof(int));
        bufd      = DoodzRealloc(StokesC->A, nnzcC*sizeof(double));
        StokesC->J = bufi;
        StokesC->A = bufd;
        
        // ACHTUNG: not sure why this had to be commented - this maybe affect other things
        // bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        // bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        // StokesD->J = bufi;
        // StokesD->A = bufd;
        
        //        printf("System size: ndof = %d, nzA = %d nzB = %d nzC = %d nzD = %d\n", Stokes->neq, nnzcA, nnzcB, nnzcC, nnzcD);
        
        printf("Number of momentum equations: %d\n", StokesA->neq);
        
        //        // Extract Diagonal of A - Viscous block
        //        int i, j, locNNZ;
        //        int I1, J1;
        //
        //        if (model.diag_scaling == 1) {
        //
        //#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
        //        for (i=0;i<StokesA->neq; i++) {
        //            I1     = StokesA->Ic[i];
        //            locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
        //            for (J1=0;J1<locNNZ; J1++) {
        //                 j = StokesA->J[I1 + J1];
        //                if (i==j) StokesA->d[i] = StokesA->A[I1 + J1];
        //             }
        //        }
        //        // Extract Diagonal of D - Pressure block
        //#pragma omp parallel for shared(StokesD) private(i)
        //        for (i=0;i<StokesD->neq; i++) {
        //            StokesC->d[i] = 1.0;
        //        }
        //
        //        // Scale A
        //#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
        //        for (i=0;i<StokesA->neq; i++) {
        //
        //            StokesA->b[i] *= StokesA->d[i];  // scale RHS
        //
        //            I1     = StokesA->Ic[i];
        //            locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
        //            for (J1=0;J1<locNNZ; J1++) {
        //                j = StokesA->J[I1 + J1];
        //                StokesA->A[I1 + J1] *= StokesA->d[i]*StokesA->d[j];
        //            }
        //        }
        //
        //        // Scale B
        //#pragma omp parallel for shared(StokesB,StokesA) private(I1,J1, i, j, locNNZ )
        //        for (i=0;i<StokesB->neq; i++) {
        //            I1     = StokesB->Ic[i];
        //            locNNZ = StokesB->Ic[i+1] - StokesB->Ic[i];
        //            for (J1=0;J1<locNNZ; J1++) {
        //                j = StokesB->J[I1 + J1];
        //                StokesB->A[I1 + J1] *= StokesA->d[i]*StokesC->d[j];
        //            }
        //        }
        //
        //        // Scale C
        //#pragma omp parallel for shared(StokesC,StokesA) private(I1,J1, i, j, locNNZ )
        //        for (i=0;i<StokesC->neq; i++) {
        //
        //            StokesC->b[i] *= StokesC->d[i]; // scale RHS
        //
        //            I1     = StokesC->Ic[i];
        //            locNNZ = StokesC->Ic[i+1] - StokesC->Ic[i];
        //            for (J1=0;J1<locNNZ; J1++) {
        //                j = StokesC->J[I1 + J1];
        //                StokesC->A[I1 + J1] *= StokesC->d[i]*StokesA->d[j];
        //            }
        //        }
        //
        //    }
        //
        //        MinMaxArray(StokesA->d, 1, StokesA->neq, "diag. A" );
        
        //        printf("\n");
        //         MinMaxArray(StokesA->A, 1, StokesA->nnz, "A" );
        //        printf("\n");
        
        
        //--------------------------------------//
        
        if ( model.writer_debug == 1 ) {
            
            char *filename;
            asprintf( &filename, "Stokes_%02dcpu_step%02d_iter%02d.gzip.h5", n_th, model.step, model.nit );
            printf("Writing Stokes matrix file: %s to disk...\n", filename);
            
            // Fill in DD data structure
            OutputSparseMatrix OutputDDA, OutputDDB, OutputDDC, OutputDDD;
             
            OutputDDA.V         = StokesA->A;
            OutputDDA.Ic        = StokesA->Ic;
            OutputDDA.J         = StokesA->J;
            OutputDDA.b         = StokesA->b;
            OutputDDA.F         = StokesA->F;
            OutputDDB.V         = StokesB->A;
            OutputDDB.Ic        = StokesB->Ic;
            OutputDDB.J         = StokesB->J;
            OutputDDB.b         = StokesB->b;
            OutputDDC.V         = StokesC->A;
            OutputDDC.Ic        = StokesC->Ic;
            OutputDDC.J         = StokesC->J;
            OutputDDC.b         = StokesC->b;
            OutputDDC.F         = StokesC->F;
            OutputDDD.V         = StokesD->A;
            OutputDDD.Ic        = StokesD->Ic;
            OutputDDD.J         = StokesD->J;
            OutputDDD.b         = StokesD->b;
            OutputDDA.eta_cell  = mesh->eta_n;
            OutputDDA.params[0] = nx;
            OutputDDA.params[1] = nz;
            OutputDDA.params[2] = mesh->dx;
            OutputDDA.params[3] = mesh->dz;
            OutputDDA.eqn_u     = DoodzMalloc(nx*nzvx*sizeof(int));
            OutputDDA.eqn_v     = DoodzMalloc(nxvz*nz*sizeof(int));
            OutputDDA.eqn_p     = DoodzMalloc(ncx*ncz*sizeof(int));
            
            for( l=0; l<nzvx; l++) {
                for( k=0; k<nx; k++) {
                    cc = k + l*nx;
                    c1 = k + l*nx;
                    OutputDDA.eqn_u[cc]=Stokes->eqn_u[cc];
                }
            }
            for( l=0; l<nz; l++) {
                for( k=0; k<nxvz; k++) {
                    cc = k + l*nxvz;
                    c3 = nx*nzvx + k + l*nxvz;
                    OutputDDA.eqn_v[cc]=Stokes->eqn_v[cc];
                }
            }
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    OutputDDA.eqn_p[cc]=Stokes->eqn_p[cc];
                }
            }
            
            // Send data to file
            CreateOutputHDF5( filename );
            AddGroupToHDF5( filename, "model" );
            AddGroupToHDF5( filename, "matrix" );
            AddGroupToHDF5( filename, "numbering" );
            AddGroupToHDF5( filename, "fields" );
            AddFieldToGroup( filename, "model", "params" , 'd', 4, OutputDDA.params,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_u" , 'i', nx*nzvx, OutputDDA.eqn_u,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_v" , 'i', nxvz*nz, OutputDDA.eqn_v,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDDA.eqn_p,  1 );
            AddFieldToGroup( filename, "fields", "eta_n" , 'd', ncx*ncz, mesh->eta_n,  1 );
            AddFieldToGroup( filename, "fields", "eta_s" , 'd', nx*nz, mesh->eta_s,  1 );
            AddFieldToGroup( filename, "fields", "tag_s" , 'c', nx*nz, mesh->BCg.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_u" , 'c', nx*nzvx, mesh->BCu.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_v" , 'c', nxvz*nz, mesh->BCv.type,  1 );
            
            AddFieldToGroup( filename, "matrix", "IA" , 'i', StokesA->neq+1, OutputDDA.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JA" , 'i', nnzcA,  OutputDDA.J,  1 );
            AddFieldToGroup( filename, "matrix", "VA" , 'd', nnzcA,  OutputDDA.V,  1 );
            AddFieldToGroup( filename, "matrix", "IB" , 'i', StokesB->neq+1, OutputDDB.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JB" , 'i', nnzcB,  OutputDDB.J,  1 );
            AddFieldToGroup( filename, "matrix", "VB" , 'd', nnzcB,  OutputDDB.V,  1 );
            AddFieldToGroup( filename, "matrix", "IC" , 'i', StokesC->neq+1, OutputDDC.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JC" , 'i', nnzcC,  OutputDDC.J,  1 );
            AddFieldToGroup( filename, "matrix", "VC" , 'd', nnzcC,  OutputDDC.V,  1 );
            AddFieldToGroup( filename, "matrix", "ID" , 'i', StokesD->neq+1, OutputDDD.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JD" , 'i', nnzcD,  OutputDDD.J,  1 );
            AddFieldToGroup( filename, "matrix", "VD" , 'd', nnzcD,  OutputDDD.V,  1 );
            AddFieldToGroup( filename, "matrix", "eta_cell", 'd', ncx*ncz, OutputDDA.eta_cell,  1 );
            AddFieldToGroup( filename, "matrix", "rhs_mom" , 'd', StokesA->neq, OutputDDA.b,  1 );
            AddFieldToGroup( filename, "matrix", "rhs_cont", 'd', StokesC->neq, OutputDDC.b,  1 );
            AddFieldToGroup( filename, "matrix", "F_mom" , 'd', StokesA->neq, OutputDDA.F,  1 );
            AddFieldToGroup( filename, "matrix", "F_cont", 'd', StokesC->neq, OutputDDC.F,  1 );
            //            AddFieldToGroup( filename, "matrix", "rhs", 'd', Stokes->neq, OutputDD.b,  1 );
            DoodzFree(OutputDDA.eqn_u);
            DoodzFree(OutputDDA.eqn_v);
            DoodzFree(OutputDDA.eqn_p);
            free(filename);
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildJacobianOperatorDecoupled( grid *mesh, params model, int lev, double *p_corr, double *p, double *u, double *v, SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, int Assemble ) {
    
    // FLAG
    // Assemble=1: Build the linear system of discrete equations
    // Assemble=0: Evaluate Stokes residuals
    
    int    cc, k, l, c1, c2, c3, nx=mesh->Nx, nz=mesh->Nz, nxvz=nx+1, nzvx=nz+1, ncx=nx-1, ncz=nz-1;
    
    // Pre-calculate FD coefs
    double celvol    = mesh->dx*mesh->dz;
    double one_dx    = 1.0/mesh->dx;
    double one_dz    = 1.0/mesh->dz;
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    double one_dx_dz = 1.0/mesh->dx/mesh->dz;
    
    // Switches
    int eqn,  sign = 1, comp = 1, stab = 0;
    if (model.free_surface_stab>0) stab = 1;
    if (model.compressible==1  ) comp = 1;
    double theta = model.free_surface_stab;
    
    // Decompose domain
    int n_th, N, ith;
    int *DD, *estart, *eend, *last_eqn;
    
    // Temporary parallel matrix
    double **AtempA;
    int **JtempA, **ItempA;
    double **AtempB;
    int **JtempB, **ItempB;
    double **AtempC;
    int **JtempC, **ItempC;
    double **AtempD;
    int **JtempD, **ItempD;
    
    int nnzA, nnzB, nnzC, nnzD, nnzcA=0, nnzcB=0, nnzcC=0, nnzcD=0;
    int *nnzc2A, *nnzc2B, *nnzc2C, *nnzc2D;
    
#pragma omp parallel shared(n_th)
    {
        n_th = omp_get_num_threads();
    }
#pragma omp barrier

    if (Assemble==1) printf("Assemble Jacobian on %d threads\n", n_th);
    
    StokesA->neq = Stokes->neq_mom;
    StokesB->neq = Stokes->neq_mom;
    StokesC->neq = Stokes->neq_cont; StokesC->neq_mom = Stokes->neq_mom;
    StokesD->neq = Stokes->neq_cont; StokesD->neq_mom = Stokes->neq_mom;
    
    // Matrix initialisation
    if ( Assemble == 1 ) {
        
        nnzA  = 21*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzB  = 6*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzC  = 2*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzD  = 1*((mesh->Nx-1) * (mesh->Nz-1));
        
        //printf("Assembling  decoupled Jacobian matrix...\n");
        AllocMat( StokesA, nnzA );
        AllocMat( StokesB, nnzB );
        AllocMat( StokesC, nnzC );
        AllocMat( StokesD, nnzD );
    }
    
    // Build velocity block RHS
    int inc=0;
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13  && mesh->BCu.type[cc] != -12 ) {
                StokesA->b[inc] = mesh->roger_x[cc];
                inc++;
            }
        }
    }
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12  ) {
                StokesA->b[inc] = mesh->roger_z[cc];
                inc++;
            }
        }
    }
    
    // Build pressure block RHS
    inc=0;
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            cc = k + l*ncx;
            if ( mesh->BCp.type[cc] != 31 && mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30 ) {
                StokesC->b[inc] = mesh->rhs_p[cc] ;
                inc++;
            }
        }
    }
    
    //------------------------------------------------------------------//
    //------------------------- U-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nx*nzvx;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p_corr, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c1=estart[ith]; c1<eend[ith]+1; c1++) {
            
            k      = mesh->kvx[c1];
            l      = mesh->lvx[c1];
            c2     = k-1 + (l-1)*ncx;
            c3     = k   + l*nxvz;
            
            // Get equation numbering
            if ( mesh->BCu.type[c1] == -1 || mesh->BCu.type[c1] == 2  || mesh->BCu.type[c1] == -2 ) {
                eqn = Stokes->eqn_u[c1];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( l>0 && l<nzvx-1 && k>0 && k<nx-1 ) {
                    Xjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, k, l, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );                    
                }
                
                // Periodic
                if ( k==0  && mesh->BCu.type[c1]==-2  ) {
                    Xjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, k, l, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );
                }
                
                // // NOT CODED YET
                // if ( k==0    && mesh->BCu.type[c1]==2 ) {
                //     Xmomentum_WestNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                // }
                
                // if ( k==nx-1 && mesh->BCu.type[c1]==2 ) {
                //     Xmomentum_EastNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                // }
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- V-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nxvz*nz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p_corr, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c3=estart[ith]; c3<eend[ith]+1; c3++) {
            
            k  = mesh->kvz[c3];
            l  = mesh->lvz[c3];
            c1 = k   + l*nx;
            c2 = k-1 + (l-1)*ncx;
            
            // Get equation numbering
            if ( mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3] == 2 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                    Zjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, k, l, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );
                }
                
                //                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                //
                //                                  Zmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                //                              }
                
                // //--------------------- INNER NODES ---------------------//
                
                // if ( l==0 && mesh->BCv.type[c3] == 2 ) {
                //     Zmomentum_SouthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                // }
                
                // //--------------------- INNER NODES ---------------------//
                
                // if ( l==nz-1 && mesh->BCv.type[c3] == 2 ) {
                //     Zmomentum_NorthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                // }
                
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- Continuity -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = ncx*ncz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempC, &ItempC, &JtempC, n_th, nnzC, Stokes->neq_cont, DD, &nnzc2C  );
        AllocateTempMatArraysDecoupled( &AtempD, &ItempD, &JtempD, n_th, nnzD, Stokes->neq_cont, DD, &nnzc2D  );
    }
    
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesC, StokesD, u, v, p, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn, comp ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, nx, ncx, nxvz, nzvx )
    {
        
        ith = omp_get_thread_num();
        
        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
            
            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c1  = k   + (l+1)*nx;
            c3  = k   + l*nxvz + 1;
            
            //--------------------- INNER NODES ---------------------//
            if ( mesh->BCp.type[c2] == -1) {
                
                comp = mesh->comp_cells[c1];
                
                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
                last_eqn[ith]   = eqn ;
                
                if ( Assemble == 1 ) {
                    ItempC[ith][eqn] = nnzc2C[ith];
                    ItempD[ith][eqn] = nnzc2D[ith];
                }
                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrixDecoupled( StokesC, AtempC, JtempC, ItempC, mesh, estart, eend, &nnzcC, nnzc2C, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
        MergeParallelMatrixDecoupled( StokesD, AtempD, JtempD, ItempD, mesh, estart, eend, &nnzcD, nnzc2D, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
    }
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempC, ItempC, JtempC, n_th, nnzc2C );
        FreeTempMatArraysDecoupled( AtempD, ItempD, JtempD, n_th, nnzc2D );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //----------------------------- End --------------------------------//
    //------------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        
        // Add contribution from the BC's
        for (k=0; k<StokesA->neq; k++) {
            StokesA->b[k] += StokesA->bbc[k];
            StokesB->b[k] = StokesA->b[k];
        }
        
        // Add contribution from the BC's
        for (k=0; k<StokesC->neq; k++) {
            StokesC->b[k] += StokesC->bbc[k];
            StokesD->b[k] = StokesC->b[k];
        }
        
        // Final index
        StokesA->Ic[StokesA->neq] = nnzcA;
        StokesB->Ic[StokesB->neq] = nnzcB;
        StokesC->Ic[StokesC->neq] = nnzcC;
        StokesD->Ic[StokesD->neq] = nnzcD;
        
        StokesA->nnz = nnzcA;
        StokesB->nnz = nnzcB;
        StokesC->nnz = nnzcC;
        StokesD->nnz = nnzcD;
        
        // Resize arrays to proper number of non-zeros
        double *bufd;
        int *bufi;
        bufi      = DoodzRealloc(StokesA->J, nnzcA*sizeof(int));
        bufd      = DoodzRealloc(StokesA->A, nnzcA*sizeof(double));
        StokesA->J = bufi;
        StokesA->A = bufd;
        
        bufi      = DoodzRealloc(StokesB->J, nnzcB*sizeof(int));
        bufd      = DoodzRealloc(StokesB->A, nnzcB*sizeof(double));
        StokesB->J = bufi;
        StokesB->A = bufd;
        
        bufi      = DoodzRealloc(StokesC->J, nnzcC*sizeof(int));
        bufd      = DoodzRealloc(StokesC->A, nnzcC*sizeof(double));
        StokesC->J = bufi;
        StokesC->A = bufd;
        
        // bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        // bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        // StokesD->J = bufi;
        // StokesD->A = bufd;
        
        //        printf("System size: ndof = %d, nzA = %d nzB = %d nzC = %d nzD = %d\n", Stokes->neq, nnzcA, nnzcB, nnzcC, nnzcD);
        
        //        // Extract Diagonal of A - Viscous block
        //        int i, j, locNNZ;
        //        int I1, J1;
        //
        //        if (model.diag_scaling == 1) {
        //
        //#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
        //            for (i=0;i<StokesA->neq; i++) {
        //                I1     = StokesA->Ic[i];
        //                locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
        //                for (J1=0;J1<locNNZ; J1++) {
        //                    j = StokesA->J[I1 + J1];
        //                    if (i==j) StokesA->d[i] = StokesA->A[I1 + J1];
        //                }
        //            }
        //            // Extract Diagonal of D - Pressure block
        //#pragma omp parallel for shared(StokesD) private(i)
        //            for (i=0;i<StokesD->neq; i++) {
        //                StokesC->d[i] = 1.0;
        //            }
        //
        //            // Scale A
        //#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
        //            for (i=0;i<StokesA->neq; i++) {
        //
        //                StokesA->b[i] *= StokesA->d[i];  // scale RHS
        //
        //                I1     = StokesA->Ic[i];
        //                locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
        //                for (J1=0;J1<locNNZ; J1++) {
        //                    j = StokesA->J[I1 + J1];
        //                    StokesA->A[I1 + J1] *= StokesA->d[i]*StokesA->d[j];
        //                }
        //            }
        //
        //            // Scale B
        //#pragma omp parallel for shared(StokesB,StokesA) private(I1,J1, i, j, locNNZ )
        //            for (i=0;i<StokesB->neq; i++) {
        //                I1     = StokesB->Ic[i];
        //                locNNZ = StokesB->Ic[i+1] - StokesB->Ic[i];
        //                for (J1=0;J1<locNNZ; J1++) {
        //                    j = StokesB->J[I1 + J1];
        //                    StokesB->A[I1 + J1] *= StokesA->d[i]*StokesC->d[j];
        //                }
        //            }
        //
        //            // Scale C
        //#pragma omp parallel for shared(StokesC,StokesA) private(I1,J1, i, j, locNNZ )
        //            for (i=0;i<StokesC->neq; i++) {
        //
        //                StokesC->b[i] *= StokesC->d[i]; // scale RHS
        //
        //                I1     = StokesC->Ic[i];
        //                locNNZ = StokesC->Ic[i+1] - StokesC->Ic[i];
        //                for (J1=0;J1<locNNZ; J1++) {
        //                    j = StokesC->J[I1 + J1];
        //                    StokesC->A[I1 + J1] *= StokesC->d[i]*StokesA->d[j];
        //                }
        //            }
        //
        //        }
        //
        //        MinMaxArray(StokesA->d, 1, StokesA->neq, "diag. A" );
        //
        //        printf("\n");
        //
        //            MinMaxArray(StokesA->A, 1, StokesA->nnz, "A" );
        //        printf("\n");
        
        
        //--------------------------------------//
#ifndef _VG_
        if ( model.writer_debug == 1 ) {
            
            char *filename;
            asprintf( &filename, "Jacobian_%02dcpu_step%02d_iter%02d.gzip.h5", n_th, model.step, model.nit );
            printf("Writing Jacobian matrix file: %s to disk...\n", filename);
            MinMaxArray(StokesA->F, 1, StokesA->neq, "Fu" );
            MinMaxArray(StokesC->F, 1, StokesC->neq, "Fp" );
            
            // Fill in DD data structure
            OutputSparseMatrix OutputDDA, OutputDDB, OutputDDC, OutputDDD;
            
            OutputDDA.V = StokesA->A;
            OutputDDA.Ic = StokesA->Ic;
            OutputDDA.J = StokesA->J;
            OutputDDA.b = StokesA->b;
            OutputDDA.F = StokesA->F;
            OutputDDB.V = StokesB->A;
            OutputDDB.Ic = StokesB->Ic;
            OutputDDB.J = StokesB->J;
            OutputDDB.b = StokesB->b;
            OutputDDC.V = StokesC->A;
            OutputDDC.Ic = StokesC->Ic;
            OutputDDC.J = StokesC->J;
            OutputDDC.b = StokesC->b;
            OutputDDC.F = StokesC->F;
            OutputDDD.V = StokesD->A;
            OutputDDD.Ic = StokesD->Ic;
            OutputDDD.J = StokesD->J;
            OutputDDD.b = StokesD->b;
            OutputDDA.eta_cell = mesh->eta_n;
            OutputDDA.params[0] = nx;
            OutputDDA.params[1] = nz;
            OutputDDA.params[2] = mesh->dx;
            OutputDDA.params[3] = mesh->dz;
            OutputDDA.eqn_u = DoodzMalloc(nx*nzvx*sizeof(int));
            OutputDDA.eqn_v = DoodzMalloc(nxvz*nz*sizeof(int));
            OutputDDA.eqn_p = DoodzMalloc(ncx*ncz*sizeof(int));
            
            for( l=0; l<nzvx; l++) {
                for( k=0; k<nx; k++) {
                    cc = k + l*nx;
                    c1 = k + l*nx;
                    OutputDDA.eqn_u[cc]=Stokes->eqn_u[cc];
                }
            }
            for( l=0; l<nz; l++) {
                for( k=0; k<nxvz; k++) {
                    cc = k + l*nxvz;
                    c3 = nx*nzvx + k + l*nxvz;
                    OutputDDA.eqn_v[cc]=Stokes->eqn_v[cc];
                }
            }
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    OutputDDA.eqn_p[cc]=Stokes->eqn_p[cc];
                }
            }
            
            // Send data to file
            CreateOutputHDF5( filename );
            AddGroupToHDF5( filename, "model" );
            AddGroupToHDF5( filename, "matrix" );
            AddGroupToHDF5( filename, "numbering" );
            AddGroupToHDF5( filename, "fields" );
            AddFieldToGroup( filename, "model", "params" , 'd', 4, OutputDDA.params,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_u" , 'i', nx*nzvx, OutputDDA.eqn_u,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_v" , 'i', nxvz*nz, OutputDDA.eqn_v,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDDA.eqn_p,  1 );
            AddFieldToGroup( filename, "fields", "eta_n" , 'd', ncx*ncz, mesh->eta_n,  1 );
            AddFieldToGroup( filename, "fields", "eta_s" , 'd', nx*nz, mesh->eta_s,  1 );
            AddFieldToGroup( filename, "fields", "tag_s" , 'c', nx*nz, mesh->BCg.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_u" , 'c', nx*nzvx, mesh->BCu.type,  1 );
            AddFieldToGroup( filename, "fields", "tag_v" , 'c', nxvz*nz, mesh->BCv.type,  1 );
            
            AddFieldToGroup( filename, "matrix", "IA" , 'i', StokesA->neq+1, OutputDDA.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JA" , 'i', nnzcA,  OutputDDA.J,  1 );
            AddFieldToGroup( filename, "matrix", "VA" , 'd', nnzcA,  OutputDDA.V,  1 );
            AddFieldToGroup( filename, "matrix", "IB" , 'i', StokesB->neq+1, OutputDDB.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JB" , 'i', nnzcB,  OutputDDB.J,  1 );
            AddFieldToGroup( filename, "matrix", "VB" , 'd', nnzcB,  OutputDDB.V,  1 );
            AddFieldToGroup( filename, "matrix", "IC" , 'i', StokesC->neq+1, OutputDDC.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JC" , 'i', nnzcC,  OutputDDC.J,  1 );
            AddFieldToGroup( filename, "matrix", "VC" , 'd', nnzcC,  OutputDDC.V,  1 );
            AddFieldToGroup( filename, "matrix", "ID" , 'i', StokesD->neq+1, OutputDDD.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "JD" , 'i', nnzcD,  OutputDDD.J,  1 );
            AddFieldToGroup( filename, "matrix", "VD" , 'd', nnzcD,  OutputDDD.V,  1 );
            AddFieldToGroup( filename, "matrix", "eta_cell", 'd', ncx*ncz, OutputDDA.eta_cell,  1 );
            AddFieldToGroup( filename, "matrix", "rhs_mom" , 'd', StokesA->neq, OutputDDA.b,  1 );
            AddFieldToGroup( filename, "matrix", "rhs_cont", 'd', StokesC->neq, OutputDDC.b,  1 );
            AddFieldToGroup( filename, "matrix", "F_mom" , 'd', StokesA->neq, OutputDDA.F,  1 );
            AddFieldToGroup( filename, "matrix", "F_cont", 'd', StokesC->neq, OutputDDC.F,  1 );
            
            DoodzFree(OutputDDA.eqn_u);
            DoodzFree(OutputDDA.eqn_v);
            DoodzFree(OutputDDA.eqn_p);
            free(filename);
        }
#endif
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MergeParallelMatrix( SparseMat *Stokes, double **Atemp, int **Jtemp, int **Itemp, grid* mesh, int* estart, int* eend, int *nnzc, int *nnzc2, int *last_eqn, int n_th, char* BC_type, int *eqn_number ) {

    int eqn, ith, begin[n_th], end[n_th], c, last_eqn_serial=0, k, l;

    for (ith=0; ith<n_th; ith++) {

        // Get last equation number
        if (last_eqn[ith]>last_eqn_serial) {
            last_eqn_serial = last_eqn[ith];
        }

        if (ith==0) {
            begin[ith] = *(nnzc);
            end[ith]   = *(nnzc) + nnzc2[ith]-1;
        }
        else {
            begin[ith] = end[ith-1]+1;
            end[ith]   = end[ith-1]+1+nnzc2[ith]-1;
        }
    }

    // Recombine A and J
#pragma omp parallel private(ith, eqn, c, k, l ) shared( Stokes, Atemp, Jtemp, Itemp, begin, BC_type, eqn_number )
    {
        ith = omp_get_thread_num();

        for( c=estart[ith]; c<eend[ith]+1; c++) {

            // Get equation numbering
            if ( BC_type[c] != 0 && BC_type[c] != 30 && BC_type[c] != 31 && BC_type[c] != 13 && BC_type[c] != 11 && BC_type[c] != -12 ) {
                eqn = eqn_number[c];
                Stokes->Ic[eqn] = Itemp[ith][eqn] + begin[ith];
            }
        }

        for (k=0; k<nnzc2[ith]; k++) {
            l = begin[ith] + k;
            Stokes->A[l] = Atemp[ith][k];
            Stokes->J[l] = Jtemp[ith][k];
        }
    }

    // Update total non-zero number
    for (ith=0; ith<n_th; ith++) {
        *(nnzc) += nnzc2[ith];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
