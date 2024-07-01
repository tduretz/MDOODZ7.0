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

// MD6
void Xjacobian_InnerNodesDecoupled3( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B, int k, int l  ) {
    
    int Newton = 1;
    double out_of_plane = 1.0, uC_corr=0.0;
    if (model.out_of_plane==1) out_of_plane = 3.0/2.0;
    double dx = mesh->dx;
    double dz = mesh->dz;
    int nzvx = model.Nz+1, nz = model.Nz;
    
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
    if (l>1)    if (mesh->BCv.type[iVzSWW] == -12)  iVzSWW = iVzSW + (nx-2);
    if (l<nz-1) if (mesh->BCv.type[iVzNWW] == -12)  iVzNWW = iVzNW + (nx-2);
    if (l>1)    if (mesh->BCv.type[iVzSEE] == -12)  iVzSEE = iVzSE - (nx-2);
    if (l<nz-1) if (mesh->BCv.type[iVzNEE] == -12)  iVzNEE = iVzNE - (nx-2);
    
    if (mesh->BCu.type[iVxC] == -2) {
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
  
    if (mesh->BCp.type[iPrW] == -1) inW = 1.0;
    if (mesh->BCp.type[iPrE] == -1) inE = 1.0;
    
    //        if (mesh->BCg.type[ixyS] != 30 && mesh->BCu.type[iVxS] != 13) inS = 1.0; // !!!!!!!!!!!!
    //        if (mesh->BCg.type[ixyN] != 30 && mesh->BCu.type[iVxN] != 13) inN = 1.0;
    if (mesh->BCu.type[iVxS] != 30 && mesh->BCu.type[iVxS] != 13 )  inS  = 1.0;
    if (mesh->BCu.type[iVxN] != 30 && mesh->BCu.type[iVxN] != 13 )  inN  = 1.0;
    if (mesh->BCg.type[ixyS] != 30 ) inSv = 1.0;
    if (mesh->BCg.type[ixyN] != 30 ) inNv = 1.0;

    // X-tra
    double wE=0.0, wW=0.0, wS=0.0, wN = 0.0;
    double inSWc=0.0,inSEc=0.0,inNWc=0.0,inNEc=0.0;
    double inSWv=0.0,inSEv=0.0,inNWv=0.0,inNEv=0.0;
    
    if ( l>1 ){// || (l==1 && mesh->BCu.type[iVxS] == 11 ) ) {  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Apparently incorrect
        if (mesh->BCp.type[iPrSW] == -1) inSWc = 1.0;
        if (mesh->BCp.type[iPrSE] == -1) inSEc = 1.0;
    }
    
    if ( l<nzvx-2 ){// || (l==nzvx-2 && mesh->BCu.type[iVxN] == 11 )) { //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Apparently incorrect
        if (mesh->BCp.type[iPrNW] == -1) inNWc = 1.0;
        if (mesh->BCp.type[iPrNE] == -1) inNEc = 1.0;
    }

    if ( (k>0)  || (k==0 && (mesh->BCv.type[iVzSW] == -1 || mesh->BCv.type[iVzSW] == 0)) ) {
        if (mesh->BCg.type[ixySW] != 30) inSWv = 1.0;   // modify for periodic
        if (mesh->BCg.type[ixyNW] != 30) inNWv = 1.0;   // modify for periodic
    }
    if ( (k<nx-1) || (k==nx-1 && (mesh->BCv.type[iVzSE] == -1 || mesh->BCv.type[iVzSE] == 0) ) ) {
        if (mesh->BCg.type[ixySE] != 30) inSEv = 1.0;   // modify for periodic
        if (mesh->BCg.type[ixyNE] != 30) inNEv = 1.0;   // modify for periodic
    }
    
    // New stuff
    double inSWW = 0.0, inNWW = 0.0;
    if (mesh->BCv.type[iVzSWW] == -1) inSWW = 1.0;
    if (mesh->BCv.type[iVzNWW] == -1) inNWW = 1.0;
    
    double inSEE = 0.0, inNEE = 0.0;
    if (mesh->BCv.type[iVzSEE] == -1) inSEE = 1.0;
    if (mesh->BCv.type[iVzNEE] == -1) inNEE = 1.0;
    
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

    // //---------------------------------------------------//

    // // double 
    // inW = 0., inE = 0., inS = 0., inN = 0.;
    // if (mesh->BCp.type[iPrW] != 30 && mesh->BCp.type[iPrW] != 31) inW = 1.0;
    // if (mesh->BCp.type[iPrE] != 30 && mesh->BCp.type[iPrE] != 31) inE = 1.0;
    // if (mesh->BCg.type[ixyS] != 30) inS = 1.0;
    // if (mesh->BCg.type[ixyN] != 30) inN = 1.0;

    // // double 
    // inNWc = 0., inNEc = 0., inSWc = 0., inSEc = 0.; 
    // if (mesh->BCp.type[iPrNW] != 30 && mesh->BCp.type[iPrNW] != 31) inNWc = 1.0;
    // if (mesh->BCp.type[iPrNE] != 30 && mesh->BCp.type[iPrNE] != 31) inNEc = 1.0;
    // if (mesh->BCp.type[iPrSW] != 30 && mesh->BCp.type[iPrSW] != 31) inSWc = 1.0;
    // if (mesh->BCp.type[iPrSE] != 30 && mesh->BCp.type[iPrSE] != 31) inSEc = 1.0;

    // // double 
    // inNWv = 0., inNEv = 0., inSWv = 0., inSEv = 0.; 
    // if (mesh->BCg.type[ixyNW] != 30 && mesh->BCg.type[ixyNW] != 31) inNWv = 1.0;
    // if (mesh->BCg.type[ixyNE] != 30 && mesh->BCg.type[ixyNE] != 31) inNEv = 1.0;
    // if (mesh->BCg.type[ixySW] != 30 && mesh->BCg.type[ixySW] != 31) inSWv = 1.0;
    // if (mesh->BCg.type[ixySE] != 30 && mesh->BCg.type[ixySE] != 31) inSEv = 1.0;

    // uW = (1.0/3.0)*(-0.75*D13W*dx*(pow(inNWv, 2) - pow(inSWv, 2)) + 0.25*dx*pow(inW, 2)*(-D31N*(comp*inW*out_of_plane - 3) + D31S*(comp*inW*out_of_plane - 3) - D32N*comp*inW*out_of_plane + D32S*comp*inW*out_of_plane) + dz*inW*(D11W*(comp*inW*out_of_plane - 3) + D12W*comp*inW*out_of_plane))/(pow(dx, 2)*dz);
    // uC = (dx*(dx*(D33N*inN + D33S*inS) + dz*(-D31N*(pow(inE, 2)*(0.083333333333333329*comp*inE*out_of_plane - 0.25) + pow(inW, 2)*(-0.083333333333333329*comp*inW*out_of_plane + 0.25)) + D31S*(pow(inE, 2)*(0.083333333333333329*comp*inE*out_of_plane - 0.25) + pow(inW, 2)*(-0.083333333333333329*comp*inW*out_of_plane + 0.25)) - 0.083333333333333329*D32N*comp*out_of_plane*(pow(inE, 3) - pow(inW, 3)) + 0.083333333333333329*D32S*comp*out_of_plane*(pow(inE, 3) - pow(inW, 3)))) - 1.0/3.0*dz*(0.75*dx*(-D13E + D13W)*(pow(inN, 2) - pow(inS, 2)) + dz*(D11E*inE*(comp*inE*out_of_plane - 3) + D11W*inW*(comp*inW*out_of_plane - 3) + D12E*comp*pow(inE, 2)*out_of_plane + D12W*comp*pow(inW, 2)*out_of_plane)))/(pow(dx, 2)*pow(dz, 2));
    // uE = (1.0/3.0)*(0.75*D13E*dx*(pow(inNEv, 2) - pow(inSEv, 2)) + 0.25*dx*pow(inE, 2)*(D31N*(comp*inE*out_of_plane - 3) - D31S*(comp*inE*out_of_plane - 3) + D32N*comp*inE*out_of_plane - D32S*comp*inE*out_of_plane) + dz*inE*(D11E*(comp*inE*out_of_plane - 3) + D12E*comp*inE*out_of_plane))/(pow(dx, 2)*dz);
    // uS = (-D33S*dx*inS + 0.25*dz*pow(inS, 2)*(D13E - D13W) + dz*(D31S*(pow(inSEc, 2)*(0.083333333333333329*comp*inSEc*out_of_plane - 0.25) + pow(inSWc, 2)*(-0.083333333333333329*comp*inSWc*out_of_plane + 0.25)) + 0.083333333333333329*D32S*comp*out_of_plane*(pow(inSEc, 3) - pow(inSWc, 3))))/(dx*pow(dz, 2));
    // uN = (-D33N*dx*inN + 0.25*dz*pow(inN, 2)*(-D13E + D13W) - dz*(D31N*(pow(inNEc, 2)*(0.083333333333333329*comp*inNEc*out_of_plane - 0.25) + pow(inNWc, 2)*(-0.083333333333333329*comp*inNWc*out_of_plane + 0.25)) + 0.083333333333333329*D32N*comp*out_of_plane*(pow(inNEc, 3) - pow(inNWc, 3))))/(dx*pow(dz, 2));
    // vSW = (-dx*(D33S*dz*inS + dx*(0.083333333333333329*D31N*comp*pow(inW, 3)*out_of_plane + 0.083333333333333329*D31S*comp*out_of_plane*(pow(inSWc, 3) - pow(inW, 3)) + 0.083333333333333329*D32N*pow(inW, 2)*(comp*inW*out_of_plane - 3) + D32S*(pow(inSWc, 2)*(0.083333333333333329*comp*inSWc*out_of_plane - 0.25) + pow(inW, 2)*(-0.083333333333333329*comp*inW*out_of_plane + 0.25)))) + (1.0/3.0)*dz*(dx*inW*(D11W*comp*inW*out_of_plane + D12W*(comp*inW*out_of_plane - 3)) + 0.75*dz*(D13E*pow(inS, 2) - D13W*(pow(inS, 2) - pow(inSWv, 2)))))/(pow(dx, 2)*pow(dz, 2));
    // vSE = (dx*(D33S*dz*inS + dx*(-0.083333333333333329*D31N*comp*pow(inE, 3)*out_of_plane + 0.083333333333333329*D31S*comp*out_of_plane*(pow(inE, 3) - pow(inSEc, 3)) - 0.083333333333333329*D32N*pow(inE, 2)*(comp*inE*out_of_plane - 3) + D32S*(pow(inE, 2)*(0.083333333333333329*comp*inE*out_of_plane - 0.25) + pow(inSEc, 2)*(-0.083333333333333329*comp*inSEc*out_of_plane + 0.25)))) - 1.0/3.0*dz*(dx*inE*(D11E*comp*inE*out_of_plane + D12E*(comp*inE*out_of_plane - 3)) + 0.75*dz*(D13E*(pow(inS, 2) - pow(inSEv, 2)) - D13W*pow(inS, 2))))/(pow(dx, 2)*pow(dz, 2));
    // vNW = (dx*(D33N*dz*inN - dx*(0.083333333333333329*D31N*comp*out_of_plane*(pow(inNWc, 3) - pow(inW, 3)) + 0.083333333333333329*D31S*comp*pow(inW, 3)*out_of_plane + D32N*(pow(inNWc, 2)*(0.083333333333333329*comp*inNWc*out_of_plane - 0.25) + pow(inW, 2)*(-0.083333333333333329*comp*inW*out_of_plane + 0.25)) + 0.083333333333333329*D32S*pow(inW, 2)*(comp*inW*out_of_plane - 3))) - 1.0/3.0*dz*(dx*inW*(D11W*comp*inW*out_of_plane + D12W*(comp*inW*out_of_plane - 3)) + 0.75*dz*(-D13E*pow(inN, 2) + D13W*(pow(inN, 2) - pow(inNWv, 2)))))/(pow(dx, 2)*pow(dz, 2));
    // vNE = (-dx*(D33N*dz*inN + dx*(-0.083333333333333329*D31N*comp*out_of_plane*(pow(inE, 3) - pow(inNEc, 3)) + 0.083333333333333329*D31S*comp*pow(inE, 3)*out_of_plane - D32N*(pow(inE, 2)*(0.083333333333333329*comp*inE*out_of_plane - 0.25) + pow(inNEc, 2)*(-0.083333333333333329*comp*inNEc*out_of_plane + 0.25)) + 0.083333333333333329*D32S*pow(inE, 2)*(comp*inE*out_of_plane - 3))) + (1.0/3.0)*dz*(dx*inE*(D11E*comp*inE*out_of_plane + D12E*(comp*inE*out_of_plane - 3)) + 0.75*dz*(-D13E*(pow(inN, 2) - pow(inNEv, 2)) + D13W*pow(inN, 2))))/(pow(dx, 2)*pow(dz, 2));
    // uSW = (-0.25*D13W*pow(inSWv, 2) + 0.083333333333333329*D31S*pow(inSWc, 2)*(comp*inSWc*out_of_plane - 3) + 0.083333333333333329*D32S*comp*pow(inSWc, 3)*out_of_plane)/(dx*dz);
    // uSE = (0.25*D13E*pow(inSEv, 2) - 0.083333333333333329*D31S*pow(inSEc, 2)*(comp*inSEc*out_of_plane - 3) - 0.083333333333333329*D32S*comp*pow(inSEc, 3)*out_of_plane)/(dx*dz);
    // uNW = (0.25*D13W*pow(inNWv, 2) - 0.083333333333333329*D31N*pow(inNWc, 2)*(comp*inNWc*out_of_plane - 3) - 0.083333333333333329*D32N*comp*pow(inNWc, 3)*out_of_plane)/(dx*dz);
    // uNE = (-0.25*D13E*pow(inNEv, 2) + 0.083333333333333329*D31N*pow(inNEc, 2)*(comp*inNEc*out_of_plane - 3) + 0.083333333333333329*D32N*comp*pow(inNEc, 3)*out_of_plane)/(dx*dz);
    // vSWW = -0.25*D13W*pow(inSWv, 2)/pow(dx, 2);
    // vSEE = -0.25*D13E*pow(inSEv, 2)/pow(dx, 2);
    // vNWW = -0.25*D13W*pow(inNWv, 2)/pow(dx, 2);
    // vNEE = -0.25*D13E*pow(inNEv, 2)/pow(dx, 2);
    // vSSW = 0.083333333333333329*pow(inSWc, 2)*(D31S*comp*inSWc*out_of_plane + D32S*(comp*inSWc*out_of_plane - 3))/pow(dz, 2);
    // vSSE = 0.083333333333333329*pow(inSEc, 2)*(D31S*comp*inSEc*out_of_plane + D32S*(comp*inSEc*out_of_plane - 3))/pow(dz, 2);
    // vNNW = 0.083333333333333329*pow(inNWc, 2)*(D31N*comp*inNWc*out_of_plane + D32N*(comp*inNWc*out_of_plane - 3))/pow(dz, 2);
    // vNNE = 0.083333333333333329*pow(inNEc, 2)*(D31N*comp*inNEc*out_of_plane + D32N*(comp*inNEc*out_of_plane - 3))/pow(dz, 2);
    // pW = (0.25*dx*inW*(-D34N + D34S) + dz*(D14W - 1))/(dx*dz);
    // pE = (0.25*dx*inE*(-D34N + D34S) + dz*(1 - D14E))/(dx*dz);
    // pSW = 0.25*D34S*inSWc/dz;
    // pSE = 0.25*D34S*inSEc/dz;
    // pNW = -0.25*D34N*inSWc/dz;
    // pNE = -0.25*D34N*inSEc/dz;

    // Stabilisation with density gradients
    if (stab==1) {
        double drhodx  = (mesh->rho_n[c2+1] - mesh->rho_n[c2])*one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        // Importante trique, voire meme gigantesque!
        if (uC+uC_corr>0.0) uC += uC_corr;
    }
    
    // Add contribution from non-conforming Dirichlets
    if ( mesh->BCu.type[iVxS]   == 11 ) uC  -=  uS ;
    if ( mesh->BCu.type[iVxN]   == 11 ) uC  -=  uN ;
    if ( mesh->BCu.type[iVxSW]  == 11 ) uW  -=  uSW;
    if ( mesh->BCu.type[iVxSE]  == 11 ) uE  -=  uSE;
    if ( mesh->BCu.type[iVxNW]  == 11 ) uW  -=  uNW;
    if ( mesh->BCu.type[iVxNE]  == 11 ) uE  -=  uNE;
    if ( mesh->BCv.type[iVzNWW] == 11 ) vNW -= vNWW;
    if ( mesh->BCv.type[iVzNEE] == 11 ) vNE -= vNEE;
    if ( mesh->BCv.type[iVzSWW] == 11 ) vSW -= vSWW;
    if ( mesh->BCv.type[iVzSEE] == 11 ) vSE -= vSEE;

    if ( Assemble == 1 ) {

        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;

        //--------------------
        // dsxx/dx - normal stencil
        
        // uSW (Newton)
        if (mesh->BCu.type[iVxSW] != 30) {
            if (mesh->BCu.type[iVxSW] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[iVxSW],     mesh->BCu.val[iVxSW], StokesA->bbc);
            else                             AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[iVxSW], 2.0*mesh->BCu.val[iVxSW], StokesA->bbc);
        }
        
        // uS
        if (mesh->BCu.type[iVxS] != 30) {
            if (mesh->BCu.type[iVxS] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[iVxS],     mesh->BCu.val[iVxS], StokesA->bbc);
            else                            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[iVxS], 2.0*mesh->BCu.val[iVxS], StokesA->bbc);
        }
        
        // uSE (Newton)
        if (mesh->BCu.type[iVxSE] != 30) {
            if (mesh->BCu.type[iVxSE] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[iVxSE],     mesh->BCu.val[iVxSE], StokesA->bbc);
            else                             AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[iVxSE], 2.0*mesh->BCu.val[iVxSE], StokesA->bbc);
        }
        
        // uW
        if (mesh->BCu.type[iVxW]  != 30)     AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxW],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[iVxW],  mesh->BCu.val[iVxW],  StokesA->bbc);
        
        // uC
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[iVxC],    mesh->BCu.val[iVxC],    StokesA->bbc);
        
        // uE
        if (mesh->BCu.type[iVxE]  != 30)     AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxE],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[iVxE],  mesh->BCu.val[iVxE],  StokesA->bbc);
        
        // uNW (Newton)
        if (mesh->BCu.type[iVxNW] != 30) {
            if (mesh->BCu.type[iVxNW] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[iVxNW],     mesh->BCu.val[iVxNW], StokesA->bbc);
            else                             AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[iVxNW], 2.0*mesh->BCu.val[iVxNW], StokesA->bbc);
        }
        
        // uN
        if (mesh->BCu.type[iVxN] != 30) {
            if (mesh->BCu.type[iVxN] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[iVxN],   mesh->BCu.val[iVxN], StokesA->bbc);
            else                            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[iVxN], 2*mesh->BCu.val[iVxN], StokesA->bbc);
        }
        
        // uNE (Newton)
        if (mesh->BCu.type[iVxNE] != 30) {
            if (mesh->BCu.type[iVxNE] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[iVxNE],     mesh->BCu.val[iVxNE], StokesA->bbc);
            else                             AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[iVxNE], 2.0*mesh->BCu.val[iVxNE], StokesA->bbc);
        }
        
        //--------------------
        
        // vSSW
        if ( l>1 ) if ( mesh->BCv.type[iVzSSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSSW],      &(nnzc2A[ith]), vSSW*celvol, mesh->BCv.type[iVzSSW],      mesh->BCv.val[iVzSSW],      StokesA->bbc);
        
        // vSSE
        if ( l>1 ) if ( mesh->BCv.type[iVzSSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSSE],      &(nnzc2A[ith]), vSSE*celvol, mesh->BCv.type[iVzSSE],      mesh->BCv.val[iVzSSE],      StokesA->bbc);
        
        //--------------------
        
        
        // vSWW (Newton)
        if ( (k > 1 ) ||  ( (k==0 || k==1) &&  mesh->BCv.type[iVzSWW] == -1) ) {
            if (mesh->BCv.type[iVzSWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSWW],   &(nnzc2A[ith]), vSWW*celvol, mesh->BCv.type[iVzSWW],   mesh->BCv.val[iVzSWW],   StokesA->bbc);
        }
        
        // vSW && vSE
        if ( mesh->BCv.type[iVzSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW], &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW], mesh->BCv.val[iVzSW],   StokesA->bbc);
        if ( mesh->BCv.type[iVzSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[iVzSE], mesh->BCv.val[iVzSE], StokesA->bbc);
        
        // vSEE (Newton)
        if ( (k < nx-2 ) ||  ( (k==nx-2 || k==nx-1) &&  mesh->BCv.type[iVzSEE] == -1) ) {
            if ( mesh->BCv.type[iVzSEE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSEE], &(nnzc2A[ith]), vSEE*celvol, mesh->BCv.type[iVzSEE], mesh->BCv.val[iVzSEE], StokesA->bbc);
        }
        
        // vNWW (Newton)
        if ( (k > 1 ) ||  ( (k==0 || k==1) &&  mesh->BCv.type[iVzNWW] == -1) ) {
            if ( mesh->BCv.type[iVzNWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNWW],        &(nnzc2A[ith]), vNWW*celvol, mesh->BCv.type[iVzNWW],        mesh->BCv.val[iVzNWW],        StokesA->bbc);
        }
        
        // vNW && vNE
        if ( mesh->BCv.type[iVzNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],      &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW],      mesh->BCv.val[iVzNW],      StokesA->bbc);
        if ( mesh->BCv.type[iVzNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[iVzNE],      mesh->BCv.val[iVzNE],      StokesA->bbc);
        
        // vNEE (Newton)
        if ( (k < nx-2 ) ||  ((k==nx-2 || k==nx-1) &&  mesh->BCv.type[iVzNEE] == -1) ) {
            if ( mesh->BCv.type[iVzNEE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNEE],      &(nnzc2A[ith]), vNEE*celvol, mesh->BCv.type[iVzNEE],      mesh->BCv.val[iVzNEE],      StokesA->bbc);
        }
        
        //--------------------
        
        // vNNW
        if ( l<nz-1 ) if ( mesh->BCv.type[iVzNNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNNW],      &(nnzc2A[ith]), vNNW*celvol, mesh->BCv.type[iVzNNW],      mesh->BCv.val[iVzNNW],      StokesA->bbc);
        
        // vNNE
        if ( l<nz-1 ) if ( mesh->BCv.type[iVzNNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNNE],      &(nnzc2A[ith]), vNNE*celvol, mesh->BCv.type[iVzNNE],      mesh->BCv.val[iVzNNE],      StokesA->bbc);
        
        //--------------------
        
        if ( Newton==1 && l>1) {
            // pSW
            if ( mesh->BCp.type[iPrSW] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSW] - Stokes->neq_mom,   &(nnzc2B[ith]), pSW*celvol, mesh->BCp.type[iPrSW],   mesh->BCp.val[iPrSW],   StokesB->bbc);
            }
            
            // pSE
            if ( mesh->BCp.type[iPrSE] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSE] - Stokes->neq_mom,   &(nnzc2B[ith]), pSE*celvol, mesh->BCp.type[iPrSE],   mesh->BCp.val[iPrSE],   StokesB->bbc);
            }
        }
        
        // pE && pW -- Valgrind / Audresselles 28/12/21
        if ( mesh->BCp.type[iPrE] != 30 && mesh->BCp.type[iPrW] != 30 ) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrW] - Stokes->neq_mom, &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[iPrW], mesh->BCp.val[iPrW], StokesB->bbc);
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrE] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[iPrE], mesh->BCp.val[iPrE], StokesB->bbc);
        }
        
        if ( Newton==1 && l<nz-1 ) {
            // pNW
            if ( mesh->BCp.type[iPrNW] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNW] - Stokes->neq_mom,   &(nnzc2B[ith]), pNW*celvol, mesh->BCp.type[iPrNW],   mesh->BCp.val[iPrNW],   StokesB->bbc);
            }
            
            // pNW
            if ( mesh->BCp.type[iPrNE] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNE] - Stokes->neq_mom,   &(nnzc2B[ith]), pNE*celvol, mesh->BCp.type[iPrNE],   mesh->BCp.val[iPrNE],   StokesB->bbc);
            }
        }
        
    }
    else {
        
        // Residual function
        StokesA->F[eqn]  = 0.0;
        StokesA->F[eqn] += (mesh->sxxd[iPrE] - mesh->sxxd[iPrW])/dx;
        StokesA->F[eqn] -= (mesh->p_corr[iPrE] - mesh->p_corr[iPrW])/dx;
        StokesA->F[eqn] += (mesh->sxz[ixyN]  - mesh->sxz[ixyS]) /dz;
        StokesA->F[eqn] *= -1.0;
        StokesA->F[eqn] -= (StokesA->b[eqn]);
        StokesA->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
void Zjacobian_InnerNodesDecoupled3( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B, int k, int l ) {
    
    double dx = mesh->dx;
    double dz = mesh->dz;
    int Newton = 1;
    double out_of_plane = 1.0;
    if (model.out_of_plane==1) out_of_plane = 3.0/2.0;
    int periodix = model.periodic_x;
    double vC_corr = 0.0;
    
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
    
    if (k==1    && mesh->BCv.type[iVzW]==-12) {
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
    if (k==nx-1 && mesh->BCv.type[iVzE]==-12) {
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
    
    if (mesh->BCp.type[iPrS] == -1) inS = 1.0;
    if (mesh->BCp.type[iPrN] == -1) inN = 1.0;
    
    //    if (mesh->BCg.type[ixyW] != 30 && mesh->BCv.type[iVzW] != 13 ) inW = 1.0; // !!!!!!!!!!!!!!!!!!!!
    //    if (mesh->BCg.type[ixyE] != 30 && mesh->BCv.type[iVzE] != 13 ) inE = 1.0;
    if ( mesh->BCv.type[iVzW] != 30 && mesh->BCv.type[iVzW] != 13 ) inW  = 1.0;
    if ( mesh->BCv.type[iVzE] != 30 && mesh->BCv.type[iVzE] != 13 ) inE  = 1.0;
    if ( mesh->BCg.type[ixyW] != 30 ) inWv = 1.0;
    if ( mesh->BCg.type[ixyE] != 30 ) inEv = 1.0;
    
    // XTRA
    double wS=0.0, wN=0.0, wE=0.0, wW=0.0;
    double inSWc=0.0,inSEc=0.0,inNWc=0.0,inNEc=0.0;
    double inSWv=0.0,inSEv=0.0,inNWv=0.0,inNEv=0.0;
    int nzvx = model.Nz+1;
    int nz   = model.Nz;
    
    if ( (k>1) || (k==1 && mesh->BCv.type[iVzW]==-1) ) {
        if (mesh->BCp.type[iPrSW] == -1) inSWc = 1.0;
        if (mesh->BCp.type[iPrNW] == -1) inNWc = 1.0;
    }
    
    if ( (k<nxvz-2) || (k==nxvz-2 && mesh->BCv.type[iVzE]==-1 ) ) {
        if (mesh->BCp.type[iPrSE] == -1) inSEc = 1.0;
        if (mesh->BCp.type[iPrNE] == -1) inNEc = 1.0;
    }
    
    if (l>0) {
        if( mesh->BCg.type[ixySW] != 30) inSWv = 1.0;
        if( mesh->BCg.type[ixySE] != 30) inSEv = 1.0;
    }
    
    if (l<nz-1) {
        if( mesh->BCg.type[ixyNW] != 30) inNWv = 1.0;
        if( mesh->BCg.type[ixyNE] != 30) inNEv = 1.0;
    }
    
    // New stuff
    double inSSW = 0.0, inSSE = 0.0;
    if (mesh->BCu.type[iVxSSW] == -1 || periodix==1) inSSW = 1.0;
    if (mesh->BCu.type[iVxSSE] == -1 || periodix==1) inSSE = 1.0;
    
    double inNNW = 0.0, inNNE = 0.0;
    if (mesh->BCu.type[iVxNNW] == -1 || periodix==1) inNNW = 1.0;
    if (mesh->BCu.type[iVxNNE] == -1 || periodix==1) inNNE = 1.0;
    
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
    if ( mesh->BCv.type[iVzW]   == 11 ) vC  -=  vW ;
    if ( mesh->BCv.type[iVzE]   == 11 ) vC  -=  vE ;
    if ( mesh->BCv.type[iVzSW]  == 11 ) vW  -=  vSW;
    if ( mesh->BCv.type[iVzSE]  == 11 ) vE  -=  vSE;
    if ( mesh->BCv.type[iVzNW]  == 11 ) vW  -=  vNW;
    if ( mesh->BCv.type[iVzNE]  == 11 ) vE  -=  vNE;
    if ( mesh->BCu.type[iVxNNW] == 11 ) uNW -= uNNW;
    if ( mesh->BCu.type[iVxNNE] == 11 ) uNE -= uNNE;
    if ( mesh->BCu.type[iVxSSW] == 11 ) uSW -= uSSW;
    if ( mesh->BCu.type[iVxSSE] == 11 ) uSE -= uSSE;
        

    //--------------------
        
    if ( Assemble == 1 ) {

        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;

        // uSSW (Newton)
        if ( mesh->BCu.type[iVxSSW]  != 30   ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSSW],   &(nnzc2A[ith]), uSSW*celvol, mesh->BCu.type[iVxSSW],    mesh->BCu.val[iVxSSW],    StokesA->bbc );
        }
        
        // uSSE (Newton)
        if ( mesh->BCu.type[iVxSSE]   != 30 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSSE],     &(nnzc2A[ith]), uSSE*celvol, mesh->BCu.type[iVxSSE],    mesh->BCu.val[iVxSSE],    StokesA->bbc );
        }
        //
        // uSWW
        if ( (k>1 ) || (k==1 && mesh->BCu.type[iVxSWW]==-1) ) {
            if( mesh->BCu.type[iVxSWW] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSWW],     &(nnzc2A[ith]), uSWW*celvol, mesh->BCu.type[iVxSWW],    mesh->BCu.val[iVxSWW],    StokesA->bbc );
        }
        //
        // uSW && uSE
        if ( mesh->BCu.type[iVxSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSW], &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[iVxSW], mesh->BCu.val[iVxSW], StokesA->bbc );
        if ( mesh->BCu.type[iVxSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSE], &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[iVxSE], mesh->BCu.val[iVxSE], StokesA->bbc );
        //
        // uSEE
        if ( (k<nx-1 ) || (k==nx-1 && mesh->BCu.type[iVxSEE]==-1) ) {
            if( mesh->BCu.type[c1+1] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxSEE],     &(nnzc2A[ith]), uSEE*celvol, mesh->BCu.type[iVxSEE],    mesh->BCu.val[iVxSEE],    StokesA->bbc );
        }
        
        // uNWW
        if ( (k>1    ) || (k==1   && mesh->BCu.type[iVxNWW]==-1) ) {
            if( mesh->BCu.type[c1+nx-2] != 30)  AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNWW],     &(nnzc2A[ith]), uNWW*celvol, mesh->BCu.type[iVxNWW],    mesh->BCu.val[iVxNWW],    StokesA->bbc );
        }
        
        // uNW && uNE
        if ( mesh->BCu.type[iVxNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNW], &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[iVxNW], mesh->BCu.val[iVxNW], StokesA->bbc );
        if ( mesh->BCu.type[iVxNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNE], &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[iVxNE], mesh->BCu.val[iVxNE], StokesA->bbc );
        
        // uNEE
        if ( (k<nx-1 ) || (k==nx-1 && mesh->BCu.type[iVxNEE]==-1) ) {
            if( mesh->BCu.type[iVxNEE] != 30) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNEE],     &(nnzc2A[ith]), uNEE*celvol, mesh->BCu.type[iVxNEE],    mesh->BCu.val[iVxNEE],    StokesA->bbc );
        }
        
        // uNNW (Newton)
        if ( mesh->BCu.type[iVxNNW]   != 30  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNNW],  &(nnzc2A[ith]), uNNW*celvol, mesh->BCu.type[iVxNNW], mesh->BCu.val[iVxNNW], StokesA->bbc );
        }
        
        // uNNE (Newton)
        if ( mesh->BCu.type[iVxNNE] != 30 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxNNE],  &(nnzc2A[ith]), uNNE*celvol, mesh->BCu.type[iVxNNE], mesh->BCu.val[iVxNNE], StokesA->bbc );
        }
        
        //--------------------
        
        // vSW (Newton)
        if ( mesh->BCv.type[iVzSW] == -1  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW],     &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW],      mesh->BCv.val[iVzSW],    StokesA->bbc );
                    //    } else if (mesh->BCv.type[iVzSW] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW],     &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW],    2*mesh->BCv.val[iVzSW],    StokesA->bbc );
        }
        
        // vS
        if ( mesh->BCv.type[iVzS] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzS],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[iVzS], mesh->BCv.val[iVzS], StokesA->bbc );
        
        // vSE (Newton)
        if ( mesh->BCv.type[iVzSE] == -1  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE],     &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[iVzSE],      mesh->BCv.val[iVzSE],    StokesA->bbc );
                    //    } else if (mesh->BCv.type[iVzSE] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE],     &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[iVzSE],    2*mesh->BCv.val[iVzSE],    StokesA->bbc );
        }
        
        // vW
        if ( mesh->BCv.type[iVzW] == -1  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzW],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[iVzW],    mesh->BCv.val[iVzW],    StokesA->bbc );
                    //    } else if (mesh->BCv.type[iVzW] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzW],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[iVzW],    2*mesh->BCv.val[iVzW],    StokesA->bbc );
        }
        
        // vC
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[iVzC],      mesh->BCv.val[iVzC],      StokesA->bbc );
        
        // vE
        if ( mesh->BCv.type[iVzE] == -1 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzE],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[iVzE],    mesh->BCv.val[iVzE],    StokesA->bbc );
                        // } else if (mesh->BCv.type[iVzE] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzE],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[iVzE],    2*mesh->BCv.val[iVzE],    StokesA->bbc );
        }
        
        // vNW (Newton)
        if ( mesh->BCv.type[iVzNW] == -1  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],     &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW],      mesh->BCv.val[iVzNW],    StokesA->bbc );
                    //    } else if (mesh->BCv.type[iVzNW] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],     &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW],    2*mesh->BCv.val[iVzNW],    StokesA->bbc );
        }
        
        // vN
        if ( mesh->BCv.type[iVzN] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzN],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[iVzN], mesh->BCv.val[iVzN], StokesA->bbc );
        
        // vNE (Newton)
        if ( mesh->BCv.type[iVzNE] == -1  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE],     &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[iVzNE],      mesh->BCv.val[iVzNE],    StokesA->bbc );
                    //    } else if (mesh->BCv.type[iVzNE] !=-12) { AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE],     &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[iVzNE],    2*mesh->BCv.val[iVzNE],    StokesA->bbc );
        }
        
        //--------------------
        // Newton  = 1;
        
        // pSW
        if ( (Newton==1 && k>1) || (Newton==1 && k==1 && periodix==1) ) {
            if ( mesh->BCp.type[iPrSW] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSW] - Stokes->neq_mom,      &(nnzc2B[ith]), pSW*celvol, mesh->BCp.type[iPrSW],     mesh->BCp.val[iPrSW],     StokesB->bbc );
            }
        }
        
        // pS
        if ( mesh->BCp.type[iPrS] != 30 && mesh->BCp.type[iPrN] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrS] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[iPrS],     mesh->BCp.val[iPrS],     StokesB->bbc );
        }
        
        // pSE
        if ( (Newton==1 && k<nx-1) || (Newton==1 && k==nx-1 && periodix==1) ) {
            if ( mesh->BCp.type[iPrSE] != 30 ) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrSE] - Stokes->neq_mom,      &(nnzc2B[ith]), pSE*celvol, mesh->BCp.type[iPrSE],     mesh->BCp.val[iPrSE],     StokesB->bbc );
            }
        }
        
        // pNW
        if ( (Newton==1 && k>1) || (Newton==1 && k==1 && periodix==1)  ) {
            if ( mesh->BCp.type[iPrNW] != 30) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNW] - Stokes->neq_mom,  &(nnzc2B[ith]), pNW*celvol, mesh->BCp.type[iPrNW], mesh->BCp.val[iPrNW], StokesB->bbc );
            }
        }
        
        // pN
        if ( mesh->BCp.type[iPrS] != 30 && mesh->BCp.type[iPrN] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrN] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[iPrN], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
        
        if ( (Newton==1 && k<nx-1) || (Newton==1 && k==nx-1 && periodix==1) ) {
            // pNE
            if ( mesh->BCp.type[iPrNE] != 30) {
                AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrNE] - Stokes->neq_mom,  &(nnzc2B[ith]), pNE*celvol, mesh->BCp.type[iPrNE], mesh->BCp.val[iPrNE], StokesB->bbc );
            }
        }
    }
    else {
        
        // Residual
        StokesA->F[eqn]  = 0.0;
        StokesA->F[eqn] += (mesh->szzd[iPrN] - mesh->szzd[iPrS])/dz;
        StokesA->F[eqn] -= (mesh->p_corr[iPrN] - mesh->p_corr[iPrS])/dz;
        StokesA->F[eqn] += (mesh->sxz[ixyE]  - mesh->sxz[ixyW]) /dx;
        StokesA->F[eqn] *= -1.0;
        StokesA->F[eqn] -= (StokesA->b[eqn]) - vC_corr*v[iVzC];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Continuity_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesC, SparseMat *StokesD, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempC, double **AtempC, int *nnzc2C, int **JtempD, double **AtempD, int *nnzc2D, int i, int j ) {
    
    double pc=0.0, uW=0.0, uE=0.0, vN=0.0, vS=0.0, oop_fact = 1.0, comp_fact = 0.0;

    const double rhoW = 0.5*(mesh->rho_n[c1]   + mesh->rho_n[c1-nx]);
    const double rhoE = 0.5*(mesh->rho_n[c1+1] + mesh->rho_n[c1-nx+1]);
    const double rhoS = 0.5*(mesh->rho_n[c1-nx]+ mesh->rho_n[c1-nx+1]);
    const double rhoN = 0.5*(mesh->rho_n[c1]   + mesh->rho_n[c1+1]);
    const double drhodx = (rhoE-rhoW)*one_dx;
    const double drhody = (rhoN-rhoS)*one_dz;
    const double VxC = 0.5*(mesh->u_in[c1]   + mesh->u_in[c1+1]);
    const double VzC = 0.5*(mesh->v_in[c3]   + mesh->v_in[c3+nxvz]);
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

void Xmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double dx = mesh->dx;
    double dz = mesh->dz;
    double out_of_plane = 1.0;
    if ( model.out_of_plane==1 ) out_of_plane = 3.0/2.0;
    double uC_corr = 0.0;

    comp  = 0;
    
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
    
    // Periodic stencil
    if ( mesh->BCu.type[iVxC]==-2 ) {
        iVxW  = c1+nx-2;
        iPrW  = c2+ncx;
        iVzSW = c3-2;
        iVzNW = c3+nxvz-2;
    }
    
    double D11E = mesh->D11_n[iPrE];
    double D11W = mesh->D11_n[iPrW];
    double D33N = mesh->D33_s[ixyN];
    double D33S = mesh->D33_s[ixyS];
    comp = 1;
    
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;

    // Flags
    double inE=0.0, inW=0.0;
    if (mesh->BCp.type[iPrW] == -1) inW = 1.0;
    if (mesh->BCp.type[iPrE] == -1) inE = 1.0;
    
    double inS=1.0, inN = 1.0, NeuS = 0.0, NeuN = 0.0, DirS = 0.0, DirN = 0.0;
    if (mesh->BCg.type[ixyS] == 30 || mesh->BCu.type[iVxS] == 13) {NeuS = 1.0; inS = 0.0;}
    if (mesh->BCg.type[ixyN] == 30 || mesh->BCu.type[iVxN] == 13) {NeuN = 1.0; inN = 0.0;}
    if (                              mesh->BCu.type[iVxS] == 11) {DirS = 1.0; inS = 0.0;}
    if (                              mesh->BCu.type[iVxN] == 11) {DirN = 1.0; inN = 0.0;}

    // Simpler
    uW = D11W*inW*(0.33333333333333331*inW*out_of_plane - 1)/pow(dx, 2);
    uC = (pow(dx, 2)*(D33N*(2*DirS + inN) + D33S*(2*DirN + inS)) - pow(dz, 2)*(D11E*inE*(0.33333333333333331*inE*out_of_plane - 1) + D11W*inW*(0.33333333333333331*inW*out_of_plane - 1)))/(pow(dx, 2)*pow(dz, 2));
    uE = D11E*inE*(0.33333333333333331*inE*out_of_plane - 1)/pow(dx, 2);
    uS = -D33S*inS/pow(dz, 2);
    uN = -D33N*inN/pow(dz, 2);
    vSW = (0.33333333333333331*D11W*pow(inW, 2)*out_of_plane - D33S)/(dx*dz);
    vSE = (-0.33333333333333331*D11E*pow(inE, 2)*out_of_plane + D33S)/(dx*dz);
    vNW = (-0.33333333333333331*D11W*pow(inW, 2)*out_of_plane + D33N)/(dx*dz);
    vNE = (0.33333333333333331*D11E*pow(inE, 2)*out_of_plane - D33N)/(dx*dz);
    pW  =  -one_dx;
    pE  =   one_dx;

    // if (eqn== 1) printf("uN = %1.4f\n", uN);
    // if (eqn==10) printf("uS = %1.4f\n", uS);

    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodx  = (mesh->rho_n[c2+1] - mesh->rho_n[c2])*one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        // Importante trique, voire meme gigantesque!
        if (uC+uC_corr<0.0) uC_corr = 0.0;
        uC += uC_corr;
    }
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[iVxS] != 30)    AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[iVxS],     mesh->BCu.val[iVxS], StokesA->bbc);        
        if (mesh->BCu.type[iVxW] != 30)    AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxW],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[iVxW],  mesh->BCu.val[iVxW],  StokesA->bbc);
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[iVxC],    mesh->BCu.val[iVxC],    StokesA->bbc);
        if (mesh->BCu.type[iVxE] != 30)   AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxE],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[iVxE],  mesh->BCu.val[iVxE],  StokesA->bbc);
        if (mesh->BCu.type[iVxN] != 30)   AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[iVxN],   mesh->BCu.val[iVxN], StokesA->bbc);
        //--------------------
        if ( mesh->BCv.type[iVzSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW], &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW], mesh->BCv.val[iVzSW], StokesA->bbc);
        if ( mesh->BCv.type[iVzSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[iVzSE], mesh->BCv.val[iVzSE], StokesA->bbc);
        if ( mesh->BCv.type[iVzNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW], &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW], mesh->BCv.val[iVzNW], StokesA->bbc);
        if ( mesh->BCv.type[iVzNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE], &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[iVzNE], mesh->BCv.val[iVzNE], StokesA->bbc);
        //--------------------
        if ( mesh->BCp.type[iPrE] != 30 && mesh->BCp.type[iPrW] != 30 ) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrW] - Stokes->neq_mom, &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[iPrW], mesh->BCp.val[iPrW], StokesB->bbc);
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrE] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[iPrE], mesh->BCp.val[iPrE], StokesB->bbc);
        }
    }
    else {
        // Residual function
        StokesA->F[eqn]  = 0.0;
        // StokesA->F[eqn] += (inE*mesh->sxxd[iPrE]   - inW*mesh->sxxd[iPrW]  )/dx;
        // StokesA->F[eqn] -= (inE*mesh->p_corr[iPrE] - inW*mesh->p_corr[iPrW])/dx;
        // StokesA->F[eqn] += (inN*mesh->sxz[ixyN]    - inS*mesh->sxz[ixyS]   )/dz;
        StokesA->F[eqn] += (mesh->sxxd[iPrE]   - mesh->sxxd[iPrW]  )/dx;
        StokesA->F[eqn] -= (mesh->p_corr[iPrE] - mesh->p_corr[iPrW])/dx;
        StokesA->F[eqn] += (mesh->sxz[ixyN]    - mesh->sxz[ixyS]   )/dz;
        StokesA->F[eqn] *= -1.0;
        StokesA->F[eqn] -= StokesA->b[eqn] - uC_corr*u[iVxC]; // no body force
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_WestNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double dx = mesh->dx;
    double dz = mesh->dz;
    double out_of_plane = 1.0;
    if ( model.out_of_plane==1 ) out_of_plane = 3.0/2.0;
    double uC_corr = 0.0;
    
    int iVxC   = c1;
    int iVxE   = iVxC+1;
    int iVxS   = iVxC-nx;
    int iVxN   = iVxC+nx;
    int iVzSW  = c3-nxvz;
    int iVzSE  = c3-nxvz+1;
    int iVzNW  = c3;
    int iVzNE  = c3+1;
    int iPrE   = c2+1;
    int ixyN   = c1;
    int ixyS   = c1-nx;
    
    double D11E = 2*mesh->eta_n[iPrE];
    double D33N =   mesh->eta_s[ixyN];
    double D33S =   mesh->eta_s[ixyS];
    comp = 1.0;
    
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;

    // Flags
    double inE=0.0, inW=0.0;
    if (mesh->BCp.type[iPrE] == -1) inE = 1.0;
    
    double inS=1.0, inN = 1.0, NeuS = 0.0, NeuN = 0.0, DirS = 0.0, DirN = 0.0;
    if (mesh->BCg.type[ixyS] == 30 || mesh->BCu.type[iVxS] == 13) {NeuS = 1.0; inS = 0.0;}
    if (mesh->BCg.type[ixyN] == 30 || mesh->BCu.type[iVxN] == 13) {NeuN = 1.0; inN = 0.0;}
    if (                              mesh->BCu.type[iVxS] == 11) {DirS = 1.0; inS = 0.0;}
    if (                              mesh->BCu.type[iVxN] == 11) {DirN = 1.0; inN = 0.0;}


    // Simpler
    uC = (-D11E*pow(dz, 2)*inE*(0.33333333333333331*inE*out_of_plane - 1) + pow(dx, 2)*(D33N*(2*DirS + inN) + D33S*(2*DirN + inS)))/(pow(dx, 2)*pow(dz, 2));
    uE = D11E*inE*(0.33333333333333331*inE*out_of_plane - 1)/pow(dx, 2);
    uS = -D33S*inS/pow(dz, 2);
    uN = -D33N*inN/pow(dz, 2);
    vSW = -D33S/(dx*dz);
    vSE = (-0.33333333333333331*D11E*pow(inE, 2)*out_of_plane + D33S)/(dx*dz);
    vNW = D33N/(dx*dz);
    vNE = (0.33333333333333331*D11E*pow(inE, 2)*out_of_plane - D33N)/(dx*dz);
    pE  =   inE*one_dx;

    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodx  = (mesh->rho_n[c2+1] - mesh->rho_n[c2])*one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        // Importante trique, voire meme gigantesque!
        if (uC+uC_corr<0.0) uC_corr = 0.0;
        uC += uC_corr;
    }
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[iVxS] != 30)    AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxS], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[iVxS],     mesh->BCu.val[iVxS], StokesA->bbc);        
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[iVxC],    mesh->BCu.val[iVxC],    StokesA->bbc);
        if (mesh->BCu.type[iVxE]  != 30)   AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxE],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[iVxE],  mesh->BCu.val[iVxE],  StokesA->bbc);
        if (mesh->BCu.type[iVxN] != 30)   AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxN], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[iVxN],   mesh->BCu.val[iVxN], StokesA->bbc);
        //--------------------
        if ( mesh->BCv.type[iVzSW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW], &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW], mesh->BCv.val[iVzSW], StokesA->bbc);
        if ( mesh->BCv.type[iVzSE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSE], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[iVzSE], mesh->BCv.val[iVzSE], StokesA->bbc);
        if ( mesh->BCv.type[iVzNW] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW], &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW], mesh->BCv.val[iVzNW], StokesA->bbc);
        if ( mesh->BCv.type[iVzNE] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNE], &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[iVzNE], mesh->BCv.val[iVzNE], StokesA->bbc);
        // printf("%d\n", Stokes->eqn_v[iVzSW]);
        // printf("%d %d\n", eqn, Stokes->eqn_v[iVzSE]);
        // printf("%d\n", Stokes->eqn_v[iVzNW]);
        // printf("%d %d\n", eqn, Stokes->eqn_v[iVzNE]);



        //--------------------
        if ( mesh->BCp.type[iPrE] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrE] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[iPrE], mesh->BCp.val[iPrE], StokesB->bbc);
        }   
    }
    else {
        // Residual function
        const double SxxE = mesh->sxxd[iPrE] - mesh->p_corr[iPrE];
        const double SxxW = 2*mesh->BCu.val[iVxC] - SxxE;
        StokesA->F[eqn]  = 0.0;
        StokesA->F[eqn] += (SxxE - SxxW )/dx;
        StokesA->F[eqn] += (mesh->sxz[ixyN] - mesh->sxz[ixyS])/dz*2; // dodgy factor 2 make it symmetric
        StokesA->F[eqn] *= -1/2;
        StokesA->F[eqn] -= StokesA->b[eqn] - uC_corr*u[iVxC]; // no body force
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_EastNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double etaW, etaN, etaS, etaE;
    etaW = mesh->eta_n[c2];
    etaE = etaW;
    etaN = mesh->eta_s[c1];
    etaS = mesh->eta_s[c1-nx];
    
    //    double rhoVx = 0.5*(mesh->rho_s[c1] + mesh->rho_s[c1-nx]);
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;
    
    StokesA->b[eqn] += mesh->BCu.val[c1]*one_dx;
    
    // dsxx/dx
    if ( mesh->BCu.type[c1-1] != 30  ) {
        uW  =  2.0*one_dx_dx * etaW         - comp*2.0/3.0*one_dx_dx * etaW;
        uC  = -2.0*one_dx_dx * (etaW )      + comp*2.0/3.0*one_dx_dx * etaW; //!!!!!!!!
        //        uE  =  2*one_dx_dx * etaE         - comp*2.0/3.0*one_dx_dx * etaE;
    }
    
    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {uS   =  one_dz_dz * etaS; uC  +=  -one_dz_dz * etaS;}
        if ( mesh->BCu.type[c1+nx] != 13 ) {uN   =  one_dz_dz * etaN; uC  +=  -one_dz_dz * etaN;}
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) uC  +=  -one_dz_dz * etaN;
        if ( mesh->BCu.type[c1+nx] != 13 ) uN   =   one_dz_dz * etaN;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) uC  +=  -one_dz_dz * etaS;
        if ( mesh->BCu.type[c1-nx] != 13 ) uS   =   one_dz_dz * etaS;
    }
    
    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * etaS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * etaN;
    
    // eta*dvz/dx
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE = -one_dx_dz * etaS + comp*2.0/3.0*one_dx_dz * etaE;
        vSW =  one_dx_dz * etaS - comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        vNE =  one_dx_dz * etaN - comp*2.0/3.0*one_dx_dz * etaE;
        vNW = -one_dx_dz * etaN + comp*2.0/3.0*one_dx_dz * etaW;
    }
    // Pressure gradient
    if ( mesh->BCp.type[c2] != 30) {
        pW  =   one_dx;
    }
    
    //    // Inertia
    //    if ( model.isinertial == 1 || model.isinertial == 2 ) {
    //        uC -= sign * rhoVx/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
    //        //        uW += rhoVx/2*one_dx*mesh->u_in[c1];
    //        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    //    }
    
    // Stabilisation with density gradients
    if ( stab == 1 ) {
        //        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_n[c2+1] - mesh->rho_n[c2]) * one_dx;
        //        uC += correction;
        //        correction = - om*0.5*model.dt * model.gx * (mesh->rho_s[c1] - mesh->rho_s[c1-nx]) * one_dz;
        //        vSW += 0.25*correction;
        //        vNE += 0.25*correction;
        //        vNW += 0.25*correction;
        //        vSE += 0.25*correction;
    }
    
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    
    if ( Assemble == 1 ) {
        
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        if (mesh->BCu.type[c1-1]  != 30)     AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  StokesA->bbc);
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2.0*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        //--------------------
        
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   StokesA->bbc);
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        StokesA->bbc);
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        
        //--------------------
        if ( mesh->BCp.type[c2] != 30  ) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   StokesB->bbc);
        }
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += uC*u[c1];
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) StokesA->F[eqn] += uS*u[c1-nx];
        }
        if (mesh->BCu.type[c1-1]  != 30)     StokesA->F[eqn] += uW*u[c1-1];
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) StokesA->F[eqn] += uN*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            if (mesh->BCv.type[c3-nxvz+1] != 11 && mesh->BCv.type[c3-nxvz+1] != 0) StokesA->F[eqn] += vSE*v[c3-nxvz+1];
            if (mesh->BCv.type[c3-nxvz  ] != 11 && mesh->BCv.type[c3-nxvz  ] != 0) StokesA->F[eqn] += vSW*v[c3-nxvz];
        }
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            if (mesh->BCv.type[c3+1] != 11 && mesh->BCv.type[c3+1] != 0) StokesA->F[eqn] += vNE*v[c3+1];
            if (mesh->BCv.type[c3  ] != 11 && mesh->BCv.type[c3  ] != 0) StokesA->F[eqn] += vNW*v[c3];
        }
        //--------------------
        if ( mesh->BCp.type[c2] != 30  ) StokesA->F[eqn]  +=  pW*p[c2];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {

    double dx = mesh->dx;
    double dz = mesh->dz;
    double out_of_plane = 1.0;
    if ( model.out_of_plane==1 ) out_of_plane = 3.0/2.0;
    double vC_corr = 0.0;
    
    int iVzC   = c3;
    int iVzW   = iVzC-1;
    int iVzE   = iVzC+1;
    int iVzS   = iVzC-nxvz;
    int iVzN   = iVzC+nxvz;
    int iVxSW  = c1-nx;
    int iVxSE  = c1-nx+1;
    int iVxNW  = c1;
    int iVxNE  = c1+1;
    int iPrS   = c2;
    int iPrN   = c2+ncx;
    int ixyE   = c1;
    int ixyW   = c1-1;
    
    double D22S  = mesh->D22_n[iPrS];
    double D22N  = mesh->D22_n[iPrN];
    double D33E  = mesh->D33_s[ixyE];
    double D33W  = mesh->D33_s[ixyW];
    comp = 1;
    
    double vS=0.0, vN=0.0, vW=0.0, vE=0.0, vC=0.0, uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, pN=0.0, pS=0.0;

    // Flags
    double inN=0.0, inS=0.0;
    if (mesh->BCp.type[iPrS] == -1) inS = 1.0;
    if (mesh->BCp.type[iPrN] == -1) inN = 1.0;
    
    double inW=1.0, inE = 1.0, NeuW = 0., NeuE = 0, DirW = 0., DirE = 0.;
    if (mesh->BCg.type[ixyW] == 30 || mesh->BCv.type[iVzW] == 13) {NeuW = 1.0; inW = 0.0;}
    if (mesh->BCg.type[ixyE] == 30 || mesh->BCv.type[iVzE] == 13) {NeuE = 1.0; inE = 0.0;}
    if (                              mesh->BCv.type[iVzW] == 11) {DirW = 1.0; inW = 0.0;}
    if (                              mesh->BCv.type[iVzE] == 11) {DirE = 1.0; inE = 0.0;}      

    // Simpler
    vW = -D33W*inW/pow(dx, 2);
    vC = (-pow(dx, 2)*(D22N + D22S)*(0.33333333333333331*out_of_plane - 1) + pow(dz, 2)*(D33E*(2*DirE + inE) + D33W*(2*DirW + inW)))/(pow(dx, 2)*pow(dz, 2));
    vE = -D33E*inE/pow(dx, 2);
    vS = D22S*(0.33333333333333331*out_of_plane - 1)/pow(dz, 2);
    vN = D22N*(0.33333333333333331*out_of_plane - 1)/pow(dz, 2);
    uSW = (0.33333333333333331*D22S*out_of_plane - D33W)/(dx*dz);
    uSE = (-0.33333333333333331*D22S*out_of_plane + D33E)/(dx*dz);
    uNW = (-0.33333333333333331*D22N*out_of_plane + D33W)/(dx*dz);
    uNE = (0.33333333333333331*D22N*out_of_plane - D33E)/(dx*dz);
    pS  =  -one_dz;
    pN  =   one_dz;

    // if (eqn==100) printf("uSW = %2.2lf\n", uSW);

    // Stabilisation with density gradients
    if ( stab==1 ) {
        double drhodz  = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
        vC_corr = 1.00 * om * model.dt * mesh->gz[c3] * drhodz;
        // Importante trique, voire meme gigantesque!
        if ((vC+vC_corr)<0.0 || (vC+vC_corr)<vW || (vC+vC_corr)<vE) vC_corr = 0.0;
        vC += vC_corr;
    }
    
    if ( Assemble == 1 ) {

        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------        
        if ( mesh->BCu.type[c1+nx]   != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        if ( mesh->BCu.type[c1]      != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        if ( mesh->BCu.type[c1+nx-1] != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
        if ( mesh->BCu.type[c1-1]    != 30 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], StokesA->bbc );
        if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3-1] != -12  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        // Periodic
        if ( mesh->BCv.type[c3-1] == -12 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz-3],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3+nxvz-3],    mesh->BCv.val[c3+nxvz-3],    StokesA->bbc );
        }
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3+1] != -12  ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        // Periodic
        if ( mesh->BCv.type[c3+1] == -12 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+3],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3-nxvz+3],    mesh->BCv.val[c3-nxvz+3],    StokesA->bbc );
        }
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesA->bbc );
        //--------------------
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     StokesB->bbc );
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+ncx] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
    }
    else {
        // Residual function
        StokesA->F[eqn]  = 0.0;
        // StokesA->F[eqn] += (inN*mesh->szzd[iPrN] - inS*mesh->szzd[iPrS])/dz;
        // StokesA->F[eqn] -= (inN*mesh->p_corr[iPrN] - inS*mesh->p_corr[iPrS])/dz;
        // StokesA->F[eqn] += (inE*mesh->sxz[ixyE]  - inW*mesh->sxz[ixyW]) /dx;
        StokesA->F[eqn] += (mesh->szzd[iPrN]   - mesh->szzd[iPrS])/dz;
        StokesA->F[eqn] -= (mesh->p_corr[iPrN] - mesh->p_corr[iPrS])/dz;
        StokesA->F[eqn] += (mesh->sxz[ixyE] - mesh->sxz[ixyW]) /dx;
        StokesA->F[eqn] *= -1.0;
        StokesA->F[eqn] -= StokesA->b[eqn] - vC_corr*v[iVzC]; 
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_SouthNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vS=0.0, vW=0.0, vC=0.0, vE=0.0, vN=0.0, pN=0.0, pS=0.0;
    
    // Coefficients
    double etaS, etaN, etaW, etaE;
    etaN  = mesh->eta_n[c2+ncx];
    etaS  = etaN;
    etaE  = mesh->eta_s[c1];
    etaW  = mesh->eta_s[c1-1];
    
    StokesA->b[eqn] -= mesh->BCv.val[c3]*one_dz;
    
    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
    }
    
    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
    }
    
    // dsyy/dz
    if (  mesh->BCv.type[c3+nxvz] !=30 ) { //mesh->BCv.type[c3-nxvz] != 30 &&
        vC  = -2.0*one_dz_dz * ( etaN ) + comp*2.0/3.0*etaS*one_dz_dz;
        vN  =  2.0*one_dz_dz *   etaN   - comp*2.0/3.0*etaN*one_dz_dz;
    }
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * etaW; vC += -one_dx_dx * etaW;}
        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * etaE; vC += -one_dx_dx * etaE;}
        
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * etaE;
        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * etaE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * etaW;
        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * etaW;
    }
    
    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * etaW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * etaE;
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
        //        double drhodx = (mesh->rho_s[c1]   - mesh->rho_s[c1-1])*one_dx;
        vC  +=  1.00 * om * model.dt * model.gz * drhodz;
        //        // Non-symmetric contibution to the system of equation
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }
    
    // Pressure gradient
    if ( mesh->BCp.type[c2+ncx] != 30 ) {
        pN  =  -one_dz;
    }
    
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    
    //    printf("%1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e\n", vC, vW,  vE, vN, uSW, uSE, uNW, uNE, pN);
    
    
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        }
        //--------------------
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesA->bbc );
        //--------------------
        if ( mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+ncx] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
        //        printf("Vz South residual = %2.2e b= %2.2e sum.c=%2.2e\n", StokesA->F[eqn], StokesA->b[eqn], uSW+uSE+uNW+uNE+vS+vW+vC+vE+vN+pN+pS);
        
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += vC*v[c3];
        if (mesh->BCv.type[c3-1] != 30) {
            if (mesh->BCv.type[c3-1] != 11) StokesA->F[eqn] += vW*v[c3-1];
        }
        if (mesh->BCv.type[c3+nxvz]  != 30) StokesA->F[eqn] += vN*v[c3+nxvz];
        if (mesh->BCv.type[c3+1] != 30) {
            if (mesh->BCv.type[c3+1] != 11) StokesA->F[eqn] += vE*v[c3+1];
        }
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            if (mesh->BCu.type[c1-1] != 11 && mesh->BCu.type[c1-1] != 0) StokesA->F[eqn] += uSW*u[c1-1];
            if (mesh->BCu.type[c1  ] != 11 && mesh->BCu.type[c1  ] != 0) StokesA->F[eqn] += uSE*u[c1  ];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            if (mesh->BCu.type[c1+nx-1] != 11 && mesh->BCu.type[c1+nx-1] != 0) StokesA->F[eqn] += uNW*u[c1+nx-1];
            if (mesh->BCu.type[c1+nx  ] != 11 && mesh->BCu.type[c1+nx  ] != 0) StokesA->F[eqn] += uNE*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCp.type[c2+ncx] != 30  ) StokesA->F[eqn]  +=  pN*p[c2+ncx];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_NorthNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vS=0.0, vW=0.0, vC=0.0, vE=0.0, vN=0.0, pN=0.0, pS=0.0;
    
    // Coefficients
    double etaS, etaW, etaE, etaN;
    etaS  = mesh->eta_n[c2];
    etaE  = mesh->eta_s[c1];
    etaW  = mesh->eta_s[c1-1];
    etaN  = etaS;
    
    StokesA->b[eqn] += mesh->BCv.val[c3]*one_dz;
    
    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
    }
    
    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
    }
    
    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 ) {
        vS  =  2.0*one_dz_dz * etaS   - comp*2.0/3.0*etaS*one_dz_dz;
        vC  = -2.0*one_dz_dz * (etaS) + comp*2.0/3.0*etaS*one_dz_dz;
    }
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * etaW; vC += -one_dx_dx * etaW;}
        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * etaE; vC += -one_dx_dx * etaE;}
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * etaE;
        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * etaE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * etaW;
        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * etaW;
    }
    
    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * etaW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * etaE;
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
        //        double drhodx = (mesh->rho_s[c1]   - mesh->rho_s[c1-1])*one_dx;
        vC  += - 1.00 * om * model.dt * model.gz * drhodz;
        //        // Non-symmetric contibution to the system of equation
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }
    
    // Pressure gradient
    if (  mesh->BCp.type[c2] != 30 ) {
        pS  =   one_dz;
    }
    
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
            AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], StokesA->bbc );
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2.0*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        AddCoeff3( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff3( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2.0*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        //--------------------
        if (  mesh->BCp.type[c2] != 30) {
            AddCoeff3( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     StokesB->bbc );
        }
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += vC*v[c3];
        if (mesh->BCv.type[c3-1] != 30) {
            if (mesh->BCv.type[c3-1] != 11) StokesA->F[eqn] += vW*v[c3-1];
        }
        if (mesh->BCv.type[c3-nxvz]  != 30) StokesA->F[eqn] += vS*v[c3-nxvz];
        if (mesh->BCv.type[c3+1] != 30) {
            if (mesh->BCv.type[c3+1] != 11) StokesA->F[eqn] += vE*v[c3+1];
        }
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            if (mesh->BCu.type[c1-1] != 11 && mesh->BCu.type[c1-1] != 0) StokesA->F[eqn] += uSW*u[c1-1];
            if (mesh->BCu.type[c1  ] != 11 && mesh->BCu.type[c1  ] != 0) StokesA->F[eqn] += uSE*u[c1  ];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            if (mesh->BCu.type[c1+nx-1] != 11 && mesh->BCu.type[c1+nx-1] != 0) StokesA->F[eqn] += uNW*u[c1+nx-1];
            if (mesh->BCu.type[c1+nx  ] != 11 && mesh->BCu.type[c1+nx  ] != 0) StokesA->F[eqn] += uNE*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCp.type[c2-ncx] != 30  ) StokesA->F[eqn]  +=  pS*p[c2-ncx];
        StokesA->F[eqn] *= celvol;
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
    
    //    if ( stab==1 ) printf("It tastes like Bo'\n");
    
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
                if ( l>0 && l<nzvx-1 && k>0 && k<nx-1 ) {
                    Xmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                // Periodic W
                if ( k==0    && mesh->BCu.type[c1]==-2 ) {
                    Xmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                // Neumann W
                if ( k==0    && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_WestNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                // Neumann E
                if ( k==nx-1 && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_EastNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
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
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                    Zmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==0 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_SouthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==nz-1 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_NorthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
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
    
    //    printf("NC = %d %d %d %d %d\n", Stokes->neq_cont, estart[0], eend[0], DD[0], last_eqn[0]);
    
    //#pragma omp parallel shared( eend, estart, mesh, Stokes, u, v, p, nx, ncx, nzvx, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    //    {
    //
    //        ith = omp_get_thread_num();
    //
    //        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
    //
    //            k   = mesh->kp[c2];
    //            l   = mesh->lp[c2];
    //            c1  = k   + (l+1)*nx;
    //            c3  = k   + l*nxvz + 1;
    //
    //            //--------------------- INNER NODES ---------------------//
    //            if ( mesh->BCp.type[c2] == -1) {
    //
    //
    //                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
    //                last_eqn[ith]   = eqn ;
    //
    //                if ( Assemble == 1 ) {
    //                    ItempC[ith][eqn] = nnzc2C[ith];
    //                    ItempD[ith][eqn] = nnzc2D[ith]; //printf("%d ",  nnzc2D[ith]);
    //                }
    //                //
    //                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
    //            }
    //        }
    //    }
    
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
        
        bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        StokesD->J = bufi;
        StokesD->A = bufd;
        
        
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
                    Xjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );                    
                }
                
                // Periodic
                if ( k==0  && mesh->BCu.type[c1]==-2  ) {
                    Xjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );
                }
                
                // NOT CODED YET
                if ( k==0    && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_WestNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                if ( k==nx-1 && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_EastNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
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
            if ( mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3] == 2 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                    Zjacobian_InnerNodesDecoupled3( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B, k, l );
                }
                
                //                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                //
                //                                  Zmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                //                              }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==0 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_SouthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==nz-1 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_NorthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p_corr, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
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
        
        bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        StokesD->J = bufi;
        StokesD->A = bufd;
        
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
