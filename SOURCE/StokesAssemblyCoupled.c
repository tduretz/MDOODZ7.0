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
#include "header_MDOODZ.h"
#include "time.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddCoeff2( int* J, double*A, int eqn, int jeq, int *nnzc, double coeff, int NODE_TYPE, double NODE_VAL, double* RHS ) {

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

void Continuity_InnerNodes( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int i, int j ) {

    double pc, uW, uE, vN, vS;

    // Compressibility
    if ( comp == 0 ) {
        pc = 0;
    }
    else {
        pc =  1;//-mesh->bet[c2]/model.dt;
    }

    // -div u
    if (comp == 0) {
        uW = one_dx;
        uE = -uW;
        vS = one_dz;
        vN = -vS;
    }
    // dt/Beta * div u
    else {
        uW = -model.dt/mesh->bet_n[c2]*one_dx;
        uE = -uW;
        vS = -model.dt/mesh->bet_n[c2]*one_dz;
        vN = -vS;
    }

    // Stencil assembly / residual
    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //        if ( (mesh->BCu.type[c1]!=30 && Stokes->eqn_u[c1]==-1) || (mesh->BCu.type[c1+1]!=30 && Stokes->eqn_u[c1+1]==-1) || (mesh->BCv.type[c3]!=30 && Stokes->eqn_v[c3]==-1) || (mesh->BCv.type[c3+nxvz]!=30 && Stokes->eqn_v[c3+nxvz]==-1) ) {
        //            printf("i=%d j =%d uW=%d uE=%d vS=%d vN=%d\n", i, j, Stokes->eqn_u[c1],Stokes->eqn_u[c1+1],Stokes->eqn_v[c3],Stokes->eqn_v[c3+nxvz]);
        //        }
        if ( mesh->BCu.type[c1] != 13 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],      &(nnzc2[ith]), uW*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc );
        if ( mesh->BCu.type[c1+1] != 13 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+1],    &(nnzc2[ith]), uE*celvol, mesh->BCu.type[c1+1],    mesh->BCu.val[c1+1],    Stokes->bbc );
        if ( mesh->BCv.type[c3] != 13 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],      &(nnzc2[ith]), vS*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      Stokes->bbc );
        if ( mesh->BCv.type[c3+nxvz] != 13 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+nxvz], &(nnzc2[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn                   , &(nnzc2[ith]), pc*celvol, mesh->BCp.type[c2],      mesh->BCp.val[c2],      Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = pc*p[c2] + uW*u[c1] + uE*u[c1+1] + vS*v[c3] + vN*v[c3+nxvz];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_WestNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    Stokes->b[eqn] += 2*mesh->BCu.val[c1] *one_dx;

    double AE = mesh->eta_n[c2+1];
    double AN = mesh->eta_s[c1];
    double AS = mesh->eta_s[c1-nx];

    // Poisson
    double uS  =  one_dz_dz * AS;
    double uC  = -2*one_dx_dx * (2*AE) - one_dz_dz * (AN + AS);//  + comp*2.0/3.0*one_dx_dx * 2*AE;
    double uE  =  2*one_dx_dx * 2*AE;//                                 - comp*2.0/3.0*one_dx_dx * 2*AE;
    double uN  =  one_dz_dz * AN;

    // Shear terms
    double vSW =  one_dx_dz * AS;// - comp*2.0/3.0*one_dx_dz * AW;
    double vSE = -one_dx_dz * AS;// + comp*2.0/3.0*one_dx_dz * AE;
    double vNW = -one_dx_dz * AN;// + comp*2.0/3.0*one_dx_dz * AW;
    double vNE =  one_dx_dz * AN;// - comp*2.0/3.0*one_dx_dz * AE;

    // Pressure gradient
    double pW  =   one_dx;
    double pE  =  -pW;

    double rhoVx = 0.5*(mesh->rho_s[c1] + mesh->rho_s[c1-nx]);

    // Inertia
    if (model.isinertial == 1 || model.isinertial == 2 ) {
        //        uC -= sign * mesh->rhoVx[c1]/model.dt;
        uC -= sign * rhoVx/model.dt;
    }
    if (model.isinertial == 2) {
        //        uN -= mesh->rhoVx[c1]/2*one_dz*mesh->VzVx[c1];
        //        uS += mesh->rhoVx[c1]/2*one_dz*mesh->VzVx[c1];
        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
        uS += rhoVx/2*one_dz*mesh->VzVx[c1];

    }

    // Stencil assembly / residual
    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+1], &(nnzc2[ith]), 2*pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], Stokes->bbc );
    }
    else {
        Stokes->F[eqn]  = 2*pE*p[c2+1];
        Stokes->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1];
        Stokes->F[eqn] += uS*u[c1-nx] + uC*u[c1] + uE*u[c1+1] + uN*u[c1+nx];
        Stokes->F[eqn] -= (Stokes->b[eqn]);
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_EastNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    Stokes->b[eqn] -= 2*mesh->BCu.val[c1] *one_dx;

    double AW = mesh->eta_n[c2];
    double AN = mesh->eta_s[c1];
    double AS = mesh->eta_s[c1-nx];

    // Poisson
    double uS  =  one_dz_dz * AS;
    double uW  =  2*one_dx_dx * 2*AW;// - comp*2.0/3.0*one_dx_dx * 2*AW;
    double uC  = -2*one_dx_dx * (2*AW) - one_dz_dz * (AN + AS);// + comp*2.0/3.0*one_dx_dx * 2*AW;
    double uN  =  one_dz_dz * AN;

    // Shear terms
    double vSW =  one_dx_dz * AS;// - comp*2.0/3.0*one_dx_dz * AW;
    double vSE = -one_dx_dz * AS;// + comp*2.0/3.0*one_dx_dz * AE;
    double vNW = -one_dx_dz * AN;// + comp*2.0/3.0*one_dx_dz * AW;
    double vNE =  one_dx_dz * AN;// - comp*2.0/3.0*one_dx_dz * AE;

    // Pressure gradient
    double pW  =   one_dx;

    double rhoVx = 0.5*(mesh->rho_s[c1] + mesh->rho_s[c1-nx]);

    // Inertia
    if (model.isinertial == 1 || model.isinertial == 2 ) {
        uC -= sign * rhoVx/model.dt;
    }
    if (model.isinertial == 2) {
        uN -= rhoVx/2.0*one_dz*mesh->VzVx[c1];
        uS += rhoVx/2.0*one_dz*mesh->VzVx[c1];
    }

    // Stencil assembly / residual
    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //-------------------- No contribution from east guy
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]),  uS*celvol, mesh->BCu.type[c1-nx],   mesh->BCu.val[c1-nx], Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2[ith]),  uW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],  Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],    &(nnzc2[ith]),  uC*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],    Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]),  uN*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx], Stokes->bbc  );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],    &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1],  &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],         &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],       &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],   &(nnzc2[ith]), 2*pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   Stokes->bbc );
    }
    else {
        Stokes->F[eqn]  = 2*pW*p[c2];
        Stokes->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1];
        Stokes->F[eqn] += uS*u[c1-nx] + uW*u[c1-1] + uC*u[c1] + uN*u[c1+nx];
        Stokes->F[eqn] -= (Stokes->b[eqn]);
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_NorthNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    double uS = -1 * mesh->p_scale * one_dx_dz;
    double uC = -uS;
    Stokes->b[eqn] = 0.0;

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = uS*u[c1-nx] + uC*u[c1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_NorthDirichlet( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    double uS = 1 * mesh->p_scale * one_dx_dz;
    double uC = uS;
    Stokes->b[eqn] = 2*mesh->BCu.val[c1] * mesh->p_scale * one_dx_dz;

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]),  uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]),  uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = uS*u[c1-nx] + uC*u[c1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_SouthNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    double uN = 1 * mesh->p_scale * one_dx_dz;
    double uC = -uN;
    Stokes->b[eqn] = 0;

    if ( Assemble == 1 ) {
        //        printf("%lf %lf\n", mesh->p_scale,  one_dx_dz);

        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = uN*u[c1+nx] + uC*u[c1];
        //        printf("%lf %lf\n", mesh->p_scale,  one_dx_dz);
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_SouthDirichlet( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    double uN = 1 * mesh->p_scale * one_dx_dz;
    double uC = uN;
    Stokes->b[eqn] = 2*mesh->BCu.val[c1] * mesh->p_scale * one_dx_dz;

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]),  uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]),  uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = uN*u[c1+nx] + uC*u[c1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void Xmomentum_WestPeriodic( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    // Full stencil
    double AE = mesh->eta_n[c2 + 1];
    double AW = mesh->eta_n[c2 + ncx];
    double AN = mesh->eta_s[l*nx];
    double AS = mesh->eta_s[l*nx - nx];

    // Poisson
    double uS  =  one_dz_dz * AS;
    double uW  =  2*one_dx_dx * AW - comp*2.0/3.0*one_dx_dx * AW;
    double uC  = -2*one_dx_dx * (AE + AW) - one_dz_dz * (AN + AS)  + comp*2.0/3.0*one_dx_dx * AE + comp*2.0/3.0*one_dx_dx * AW;
    double uE  =  2*one_dx_dx * AE - comp*2.0/3.0*one_dx_dx * AE;
    double uN  =  one_dz_dz * AN;

    // Shear terms
    double vSW =  one_dx_dz * AS - comp*2.0/3.0*one_dx_dz * AW;
    double vSE = -one_dx_dz * AS + comp*2.0/3.0*one_dx_dz * AE;
    double vNW = -one_dx_dz * AN + comp*2.0/3.0*one_dx_dz * AW;
    double vNE =  one_dx_dz * AN - comp*2.0/3.0*one_dx_dz * AE;

    // Pressure gradient
    double pW  =   one_dx;
    double pE  =  -pW;

    double rhoVx = 0.5*(mesh->rho_s[c1] + mesh->rho_s[c1-nx]);


    // Inertia
    if (model.isinertial == 1 || model.isinertial == 2 ) {
        uC -= sign * rhoVx/model.dt;
    }
    if (model.isinertial == 2) {
        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
        uW += rhoVx/2*one_dx*mesh->u_in[c1];
        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    }

    // Stabilisation with density gradients
    if (stab == 1) {
        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_n[c2+1] - mesh->rho_n[c2]) * one_dx;
        uC += correction;
        correction = - om*0.5*model.dt * model.gx * (mesh->rho_s[l*nx] - mesh->rho_s[l*nx-nx]) * one_dz;
        vSW += 0.25*correction;
        vNE += 0.25*correction;
        vNW += 0.25*correction;
        vSE += 0.25*correction;
    }

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx],  &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx],   mesh->BCu.val[c1-nx],   Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                   &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+1],   &(nnzc2[ith]), uE*celvol, mesh->BCu.type[c1+1],    mesh->BCu.val[c1+1],    Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx-2],&(nnzc2[ith]), uW*celvol, mesh->BCu.type[c1+nx-2], mesh->BCu.val[c1+nx-2], Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx],  &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   Stokes->bbc  );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1],    &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1],    mesh->BCv.val[c3-nxvz+1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+nx-1], &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz+nx-1], mesh->BCv.val[c3-nxvz+nx-1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],         &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],         mesh->BCv.val[c3+1],         Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+nx-1],      &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3+nx-1],      mesh->BCv.val[c3+nx-1],      Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+1],    &(nnzc2[ith]), pE*celvol, mesh->BCp.type[c2+1],   mesh->BCp.val[c2+1],   Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+ncx],  &(nnzc2[ith]), pW*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], Stokes->bbc );
    }
    else {
        Stokes->F[eqn]  = pE*p[c2+1] + pW*p[c2+ncx];
        Stokes->F[eqn] += vSW*v[c3-nxvz+nx-1] + vSE*v[c3-nxvz+1] + vNW*v[c3+nx-1] + vNE*v[c3+1];
        Stokes->F[eqn] += uS*u[c1-nx] + uW*u[c1+nx-2] + uC*u[c1] + uE*u[c1+1] + uN*u[c1+nx];
        Stokes->F[eqn] -= (Stokes->b[eqn]);
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_EastPeriodic( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {


    double uW = 1 * mesh->p_scale * one_dx_dz;
    double uC = -uW;

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-ncx], &(nnzc2[ith]), uW*celvol, mesh->BCu.type[c1-ncx], mesh->BCu.val[c1-ncx], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                   &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],     mesh->BCu.val[c1],     Stokes->bbc  );
    }
    else {
        Stokes->F[eqn] = uW*u[c1-ncx] + uC*u[c1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void Xmomentum_InnerNodes( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, int om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {
//
//    double AE = mesh->eta_n[c2+1];
//    double AW = mesh->eta_n[c2];
//    double AN = mesh->eta_s[c1];
//    double AS = mesh->eta_s[c1-nx];
//
//    // Poisson
//    double uS  =  one_dz_dz * AS;
//    double uW  =  2*one_dx_dx * AW                                 - comp*2.0/3.0*one_dx_dx * AW;
//    double uC  = -2*one_dx_dx * (AE + AW) - one_dz_dz * (AN + AS)  + comp*2.0/3.0*one_dx_dx * AE + comp*2.0/3.0*one_dx_dx * AW;
//    double uE  =  2*one_dx_dx * AE                                 - comp*2.0/3.0*one_dx_dx * AE;
//    double uN  =  one_dz_dz * AN;
//
//    // Shear terms
//    double vSW =  one_dx_dz * AS - comp*2.0/3.0*one_dx_dz * AW;
//    double vSE = -one_dx_dz * AS + comp*2.0/3.0*one_dx_dz * AE;
//    double vNW = -one_dx_dz * AN + comp*2.0/3.0*one_dx_dz * AW;
//    double vNE =  one_dx_dz * AN - comp*2.0/3.0*one_dx_dz * AE;
//
//    // Pressure gradient
//    double pW  =   one_dx;
//    double pE  =  -pW;
//
//    // Inertia
//    if ( model.isinertial == 1 || model.isinertial == 2 ) {
//        uC -= sign * mesh->rhoVx[c1]/model.dt;
//    }
//    if ( model.isinertial == 2 ) {
//        uN -= mesh->rhoVx[c1]/2*one_dz*mesh->VzVx[c1];
//        uS += mesh->rhoVx[c1]/2*one_dz*mesh->VzVx[c1];
//        uW += mesh->rhoVx[c1]/2*one_dx*mesh->u_in[c1];
//        uE -= mesh->rhoVx[c1]/2*one_dx*mesh->u_in[c1];
//    }
//
//    // Stabilisation with density gradients
//    if ( stab == 1 ) {
//        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_n[c2+1] - mesh->rho_n[c2]) * one_dx;
//        uC += correction;
//        correction = - om*0.5*model.dt * model.gx * (mesh->rho_s[c1] - mesh->rho_s[c1-nx]) * one_dz;
//        vSW += 0.25*correction;
//        vNE += 0.25*correction;
//        vNW += 0.25*correction;
//        vSE += 0.25*correction;
//    }
//
//    if ( Assemble == 1 ) {
//        Stokes->b[eqn] *= celvol;
//        //--------------------
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], Stokes->bbc );
//        //--------------------
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      Stokes->bbc );
//        //--------------------
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],   &(nnzc2[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   Stokes->bbc );
//        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+1], &(nnzc2[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], Stokes->bbc );
//    }
//    else {
//        Stokes->F[eqn]  = pW*p[c2] + pE*p[c2+1];
//        Stokes->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1];
//        Stokes->F[eqn] += uS*u[c1-nx] + uW*u[c1-1] + uC*u[c1] + uE*u[c1+1] + uN*u[c1+nx];
//        Stokes->F[eqn] -= (Stokes->b[eqn]);
//    }
//
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_InnerNodes( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double AE = mesh->eta_n[c2+1];
    double AW = mesh->eta_n[c2];
    double AN = mesh->eta_s[c1];
    double AS = mesh->eta_s[c1-nx];

    double rhoVx = 0.5*(mesh->rho_s[c1] + mesh->rho_s[c1-nx]);


    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;

    // dsxx/dx
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] != 30 ) {
        uW  =  2*one_dx_dx * AW         - comp*2.0/3.0*one_dx_dx * AW;
        uC  = -2*one_dx_dx * (AE + AW)  + comp*2.0/3.0*one_dx_dx * AE + comp*2.0/3.0*one_dx_dx * AW;
        uE  =  2*one_dx_dx * AE         - comp*2.0/3.0*one_dx_dx * AE;
    }
    if ( mesh->BCu.type[c1-1] == 30 && mesh->BCu.type[c1+1] != 30 ) {
        uC  = -2*one_dx_dx * (AE) + comp*2.0/3.0*one_dx_dx * AE;
        uE  =  2*one_dx_dx * AE - comp*2.0/3.0*one_dx_dx * AE;
    }
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] == 30 ) {
        uW  =  2*one_dx_dx * AW - comp*2.0/3.0*one_dx_dx * AW;
        uC  = -2*one_dx_dx * (AW) + comp*2.0/3.0*one_dx_dx * AW;
    }

    //    // eta*dvx/dz
    //    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
    //        {uS   =  one_dz_dz * AS; uC  +=  -one_dz_dz * AS;}
    //         {uN   =  one_dz_dz * AN; uC  +=  -one_dz_dz * AN;}
    //    }
    //    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
    //         uC  +=  -one_dz_dz * AN;
    //         uN   =   one_dz_dz * AN;
    //    }
    //    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
    //         uC  +=  -one_dz_dz * AS;
    //         uS   =   one_dz_dz * AS;
    //    }

    //    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * AS;
    //    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * AN;


    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {uS   =  one_dz_dz * AS; uC  +=  -one_dz_dz * AS;}
        if ( mesh->BCu.type[c1+nx] != 13 ) {uN   =  one_dz_dz * AN; uC  +=  -one_dz_dz * AN;}
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) uC  +=  -one_dz_dz * AN;
        if ( mesh->BCu.type[c1+nx] != 13 ) uN   =   one_dz_dz * AN;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) uC  +=  -one_dz_dz * AS;
        if ( mesh->BCu.type[c1-nx] != 13 ) uS   =   one_dz_dz * AS;
    }

    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * AS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * AN;

    // eta*dvz/dx
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE = -one_dx_dz * AS + comp*2.0/3.0*one_dx_dz * AE;
        vSW =  one_dx_dz * AS - comp*2.0/3.0*one_dx_dz * AW;
    }

    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        vNE =  one_dx_dz * AN - comp*2.0/3.0*one_dx_dz * AE;
        vNW = -one_dx_dz * AN + comp*2.0/3.0*one_dx_dz * AW;
    }

    //    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3+nxvz+1] != 30  ) {
    //        vSE +=  comp*2.0/3.0*one_dx_dz * AE;
    //        vNE += -comp*2.0/3.0*one_dx_dz * AE;
    //    }
    //
    //    if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+nxvz] != 30  ) {
    //        vSW +=  -comp*2.0/3.0*one_dx_dz * AW;
    //        vNW +=   comp*2.0/3.0*one_dx_dz * AW;
    //    }



    //    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
    //        vNE =  one_dx_dz * AN - comp*2.0/3.0*one_dx_dz * AE;
    //        vSE = -one_dx_dz * AS + comp*2.0/3.0*one_dx_dz * AE;
    //    }
    //
    //    if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3-nxvz] != 30 )  {
    //        vNW = -one_dx_dz * AN + comp*2.0/3.0*one_dx_dz * AW;
    //        vSW =  one_dx_dz * AS - comp*2.0/3.0*one_dx_dz * AW;
    //    }


    // Pressure gradient
    if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
        pW  =   one_dx;
        pE  =  -pW;
    }


    // Inertia
    if ( model.isinertial == 1 || model.isinertial == 2 ) {
        uC -= sign * rhoVx/model.dt;
    }
    if ( model.isinertial == 2 ) {
        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
        uW += rhoVx/2*one_dx*mesh->u_in[c1];
        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    }

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

        //        if (isnan(uC)>0 || isnan(uE)>0 || isnan(uW)>0 || isnan(uN)>0 || isnan(uS)>0 || isnan(vSE)>0 || isnan(vNE)>0 || isnan(vSW)>0 || isnan(vNW)>0 || isnan(pE)>0 ||  isnan(pW)>0)
        //        {
        //        printf("%2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf %2.2lf \n", uC, uE, uW, uN, uS, vSE, vNE, vSW, vNW, pE, pW);
        //            printf("%2.2lf %2.2lf %2.2lf %2.2lf\n", AE, AW, AN, AS);
        //            printf("%d %d\n", c2+1, ncx*(mesh->Nz-1));
        //            IsNanArray2DFP(mesh->eta_n, (mesh->Nx-1)*(mesh->Nz-1));
        //        }

        Stokes->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], Stokes->bbc );
        }
        if (mesh->BCu.type[c1-1]  != 30) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    Stokes->bbc );
        if (mesh->BCu.type[c1+1]  != 30) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  Stokes->bbc );
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2*mesh->BCu.val[c1+nx], Stokes->bbc );
        }
        //--------------------

        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], Stokes->bbc );
        }

        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      Stokes->bbc );
        }

        //        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3] != 30  ) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   Stokes->bbc );
        //        }
        //
        //        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], Stokes->bbc );
        //        }
        //        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3] != 30  ) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc );
        //        }
        //
        //        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      Stokes->bbc );
        //        }

        //--------------------
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],   &(nnzc2[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+1], &(nnzc2[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], Stokes->bbc );
        }

        //        if ( mesh->BCp.type[c2+1] == -1 && mesh->BCp.type[c2] == -1 ) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],   &(nnzc2[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   Stokes->bbc );
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+1], &(nnzc2[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], Stokes->bbc );
        //        }
        //
        //        // Pressure gradient
        //        if ( mesh->BCp.type[c2+1] == 31 && mesh->BCp.type[c2] == -1 ) {
        //             AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],   &(nnzc2[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   Stokes->bbc );
        //        }
        //
        //        // Pressure gradient
        //        if ( mesh->BCp.type[c2-1] == 31 && mesh->BCp.type[c2] == -1 ) {
        //           AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2-1], &(nnzc2[ith]), pE*celvol, mesh->BCp.type[c2-1], mesh->BCp.val[c2-1], Stokes->bbc );
        //        }

        Stokes->F[eqn] = u[c1];

    }
    else {

        Stokes->F[eqn] = uC*u[c1];
        if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1-nx] != 11 ) Stokes->F[eqn] += uS*u[c1-nx];
        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1+nx] != 11 ) Stokes->F[eqn] += uN*u[c1+nx];
        if ( mesh->BCu.type[c1-1]  != 30 ) Stokes->F[eqn] += uW*u[c1-1];
        if ( mesh->BCu.type[c1+1]  != 30 ) Stokes->F[eqn] += uE*u[c1+1];
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            Stokes->F[eqn]  += pW*p[c2] + pE*p[c2+1];
        }
        if ( mesh->BCv.type[c3-nxvz+1] != 30 && mesh->BCv.type[c3-nxvz] != 30  ) {
            Stokes->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1];
        }
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30 ) {
            Stokes->F[eqn] += vNW*v[c3] + vNE*v[c3+1];
        }
//        if ( mesh->BCu.type[c1-nx] == 11 )   Stokes->F[eqn] += -2*AS*one_dz_dz*mesh->BCu.val[c1-nx];
//        if ( mesh->BCu.type[c1+nx] == 11 )   Stokes->F[eqn] += -2*AN*one_dz_dz*mesh->BCu.val[c1+nx];
        Stokes->F[eqn] -= (Stokes->b[eqn])+ 0.0*Stokes->bbc[eqn];
        Stokes->F[eqn] *= celvol;


        if (eqn==277) {printf("\nFx = %2.2e %d\n",Stokes->F[eqn], mesh->BCu.type[c1] );

            ////        double F1 = pW*p[c2] + pE*p[c2+1] + vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1] + uS*u[c1-nx] + uN*u[c1+nx] + uW*u[c1-1] + uE*u[c1+1] + uC*u[c1] + (StokesA->b[eqn] - StokesA->bbc[eqn]);
            ////        if (fabs(F1)>1e-6) {
            ////        printf("F0 = %2.2e %2.2e\n", StokesA->F[eqn], F1*celvol);
            printf("x momentum --  AW = %2.2e AE = %2.2e AS = %2.2e AN = %2.2e\n", AW, AE, AS, AN);
            printf("uC = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1+0],  uC, u[c1] );
            printf("uW = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1-1],  uW, u[c1-1] );
            printf("uE = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1+1],  uE, u[c1+1] );
            printf("uS = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1-nx], uS, u[c1-nx] );
            printf("uN = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1+nx], uN, u[c1+nx] );
            printf("vSE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz+1], vSE, v[c3-nxvz+1] );
            printf("vNE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3+1], vNE, v[c3+1] );
            printf("vSW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz], vSW, v[c3-nxvz] );
            printf("vNW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3], vNW, v[c3] );
            printf("pW = %02d %2.2e %2.2e \n",     mesh->BCp.type[c2], pW, p[c2] );
            printf("pE = %02d %2.2e %2.2e %2.2e %2.2e bbc = %2.2e b = %2.2e\n",mesh->BCp.type[c2+1], pE, p[c2+1], Stokes->b[eqn] + Stokes->bbc[eqn], mesh->BCu.val[c1], Stokes->bbc[eqn] , Stokes->b[eqn]);
            printf("%2.2e EQN=%d\n", (uC*u[c1] + uW*u[c1-1])*4e-5* celvol, eqn);
        }



    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_WestNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double vC = -1 * mesh->p_scale * one_dx_dz;
    double vE = -vC;
    Stokes->b[eqn] = 0 * one_dx_dz;

    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],   mesh->BCv.val[c3],   Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],  &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1], mesh->BCv.val[c3+1], Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vE*v[c3+1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_WestDirichlet( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double vC = 1 * mesh->p_scale * one_dx_dz;
    double vE = vC;
    Stokes->b[eqn] = 2*mesh->BCv.val[c3] * mesh->p_scale * one_dx_dz;

    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],   mesh->BCv.val[c3],   Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],  &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1], mesh->BCv.val[c3+1], Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vE*v[c3+1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_WestPeriodic( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    //    jeq = nzvx*nx + l*nxvz + nxvz-2;

    double vC = 1 * mesh->p_scale * one_dx_dz;
    double vE = -vC;
    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                             &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],              mesh->BCv.val[c3],              Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[l*nxvz + nxvz-2],  &(nnzc2[ith]), vE*celvol, mesh->BCv.type[l*nxvz + nxvz-2], mesh->BCv.val[l*nxvz + nxvz-2], Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vE*v[l*nxvz + nxvz-2];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_EastNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double vC = 1 * mesh->p_scale * one_dx_dz;
    double vW = -vC;
    Stokes->b[eqn] = 0 * one_dx_dz;

    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],  &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1], mesh->BCv.val[c3-1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],   mesh->BCv.val[c3],   Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vW*v[c3-1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_EastDirichlet( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double vC = 1 * mesh->p_scale * one_dx_dz;
    double vW = vC;
    Stokes->b[eqn] = 2*mesh->BCv.val[c3] * mesh->p_scale * one_dx_dz;

    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],  &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1], mesh->BCv.val[c3-1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                  &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],   mesh->BCv.val[c3],   Stokes->bbc );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vW*v[c3-1];
        Stokes->F[eqn] -= Stokes->b[eqn];
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_EastPeriodic( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2, int l ) {

    //    jeq = nzvx*nx + l*nxvz +1;
    double vC = 1 * mesh->p_scale * one_dx_dz;
    double vW = -vC;
    Stokes->b[eqn] = 0 ;

    if ( Assemble == 1) {
        Stokes->b[eqn] *= celvol;
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[l*nxvz +1],  &(nnzc2[ith]), vW*celvol, mesh->BCv.type[l*nxvz +1], mesh->BCv.val[l*nxvz +1], Stokes->bbc  );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                       &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        Stokes->bbc  );
    }
    else {
        Stokes->F[eqn] = vC*v[c3] + vW*v[l*nxvz +1];
        Stokes->F[eqn] *= celvol;

    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_NorthNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    Stokes->b[eqn] += -2*mesh->BCv.val[c3] *one_dz;

    // Coefficients
    double AS  = mesh->eta_n[c2];
    double AE  = mesh->eta_s[c1];
    double AW  = mesh->eta_s[c1-1];

    // Shear terms
    double uSW =  one_dx_dz * AW;// - comp*2.0/3.0*AS*one_dx_dz;
    double uSE = -one_dx_dz * AE;// + comp*2.0/3.0*AS*one_dx_dz;
    double uNW = -one_dx_dz * AW;// + comp*2.0/3.0*AN*one_dx_dz;
    double uNE =  one_dx_dz * AE;// - comp*2.0/3.0*AN*one_dx_dz;

    // Poisson
    double vS  =  2*one_dz_dz * 2*AS;// - comp*2.0/3.0*2*AS*one_dz_dz;
    double vW  =  one_dx_dx * AW;
    double vC  = -2*one_dz_dz * (2*AS) - one_dx_dx * (AE + AW);// + comp*2.0/3.0*2*AS*one_dz_dz;
    double vE  =  one_dx_dx * AE;

    // Pressure gradient
    double pS  =   one_dz;

    double rhoVz = 0.5 * (mesh->rho_s[c1] + mesh->rho_s[c1+1]);

    // Inertial terms
    if (model.isinertial == 1 || model.isinertial == 2 ) {
        vC -= sign * rhoVz/model.dt;
    }

    if ( model.isinertial == 2 ) {
        //        vE -= mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
        //        vW += mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
        vE -= rhoVz/2*one_dx*mesh->VxVz[c3];
        vW += rhoVz/2*one_dx*mesh->VxVz[c3];

    }

    // Stabilisation with density gradients
    if (stab==1) {
        double correction = - om*0.5*model.dt * model.gz * (mesh->rho_n[c2+ncx] - mesh->rho_n[c2]) * one_dz;
        vC += correction;
        correction = - om*0.5*model.dt * model.gz * (mesh->rho_s[c1] - mesh->rho_s[c1-1]) * one_dx;
        uNW += 0.25*correction;
        uNE += 0.25*correction;
        uSW += 0.25*correction;
        uSE += 0.25*correction;
    }

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                     &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      Stokes->bbc );
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    Stokes->bbc );
        //--------------------
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],      &(nnzc2[ith]), 2*pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     Stokes->bbc );
    }
    else {
        Stokes->F[eqn]  = 2*pS*p[c2];
        Stokes->F[eqn] += vC*v[c3] + vS*v[c3-nxvz] + vW*v[c3-1] + vE*v[c3+1];
        Stokes->F[eqn] += uSW*u[c1-1] + uSE*u[c1] + uNW*u[c1+nx-1] + uNE*u[c1+nx];
        //        Stokes->F[eqn] -= (Stokes->b[eqn] + Stokes->bbc[eqn]);
        Stokes->F[eqn] -= (Stokes->b[eqn]);
        Stokes->F[eqn] *= celvol;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_SouthNeumann( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    Stokes->b[eqn] += mesh->BCv.val[c3] *one_dz;

    double uSW=0, uSE=0, uNW=0, uNE=0, vS=0, vW=0, vC=0, vE=0, vN=0, pN=0, pS=0;

    // Coefficients
    double AS  = 0;
    double AN  = mesh->eta_n[c2+ncx];
    double AE  = mesh->eta_s[c1];
    double AW  = mesh->eta_s[c1-1];

    //    // (eta du/dx) S
    //    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
    //        uSW =  one_dx_dz * AW - comp*2.0/3.0*AS*one_dx_dz;
    //        uSE = -one_dx_dz * AE + comp*2.0/3.0*AS*one_dx_dz;
    //    }

    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * AW + comp*2.0/3.0*AN*one_dx_dz;
        uNE =  one_dx_dz * AE - comp*2.0/3.0*AN*one_dx_dz;
    }

    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3+nxvz] !=30 ) {
        vS  =  2*one_dz_dz * AS        - comp*2.0/3.0*AS*one_dz_dz;
        vC  = -2*one_dz_dz * (AN + AS) + comp*2.0/3.0*AN*one_dz_dz + comp*2.0/3.0*AS*one_dz_dz;
        vN  =  2*one_dz_dz * AN        - comp*2.0/3.0*AN*one_dz_dz;
    }
    if ( mesh->BCv.type[c3-nxvz] == 30 && mesh->BCv.type[c3+nxvz] != 30 ) {
        vC  = -2*one_dz_dz * (AN) + comp*2.0/3.0*AN*one_dz_dz;
        vN  =  2*one_dz_dz * AN   - comp*2.0/3.0*AN*one_dz_dz;
    }
    if ( mesh->BCv.type[c3-nxvz] != 30  && mesh->BCv.type[c3+nxvz] == 30 ) {
        vS  =  2*one_dz_dz * AS   - comp*2.0/3.0*AS*one_dz_dz;
        vC  = -2*one_dz_dz * (AS) + comp*2.0/3.0*AS*one_dz_dz;
    }

    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        {vW  =  one_dx_dx * AW; vC += -one_dx_dx * AW;}
        {vE  =  one_dx_dx * AE; vC += -one_dx_dx * AE;}
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        vE  =  one_dx_dx * AE;
        vC += -one_dx_dx * AE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        vW  =  one_dx_dx * AW;
        vC += -one_dx_dx * AW;
    }

    //    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * AW;
    //    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * AE;


    //    // d/dx eta*dvzdx
    //    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
    //        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * AW; vC += -one_dx_dx * AW;}
    //        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * AE; vC += -one_dx_dx * AE;}
    //    }
    //    // d/dx eta*dvzdx
    //    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
    //        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * AE;
    //        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * AE;
    //    }
    //    // d/dx eta*dvzdx
    //    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
    //        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * AW;
    //        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * AW;
    //    }
    //
    //    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * AW;
    //    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * AE;

    //    // Inertial terms
    //    if (model.isinertial == 1 || model.isinertial == 2 ) {
    //        vC -= sign * mesh->rhoVz[c3]/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        vN -= mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        vS += mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        vE -= mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
    //        vW += mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
    //    }


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
    if ( mesh->BCp.type[c2+ncx] != 30 ) {
        pN  =  -one_dz;
    }


    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;

    if ( Assemble == 1 ) {
        Stokes->b[eqn] *= celvol;
        //--------------------
        //        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    Stokes->bbc );
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc );
        //        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   Stokes->bbc );
        }
        //--------------------
        //        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], Stokes->bbc );
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    2*mesh->BCv.val[c3-1],    Stokes->bbc );
        }
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                     &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      Stokes->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    2*mesh->BCv.val[c3+1],    Stokes->bbc );
        }
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], Stokes->bbc );
        //--------------------
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30) {
            //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],      &(nnzc2[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+ncx],  &(nnzc2[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], Stokes->bbc );
        }

    }
    else {
        if ( mesh->BCp.type[c2+ncx] != 30 ) {
            //            printf("v=%2.2e p=%2.2e\n", pN*celvol, p[c2+ncx]);

            Stokes->F[eqn]  = pN*p[c2+ncx];

        }
        Stokes->F[eqn] += vC*v[c3];
        //        if ( mesh->BCv.type[c3-nxvz] != 30 ) Stokes->F[eqn] += vS*v[c3-nxvz];
        if ( mesh->BCv.type[c3-1]    != 30 ) Stokes->F[eqn] += vW*v[c3-1];
        if ( mesh->BCv.type[c3+1]    != 30 ) Stokes->F[eqn] += vE*v[c3+1];
        if ( mesh->BCv.type[c3+nxvz] != 30 ) Stokes->F[eqn] += vN*v[c3+nxvz];
        //        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        //            Stokes->F[eqn] += uSW*u[c1-1] + uSE*u[c1];
        ////            printf("F -> uSW = %2.2e uSE = %2.2e\n", uSW, uSE );
        //        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            Stokes->F[eqn] += uNW*u[c1+nx-1] + uNE*u[c1+nx];
        }
        //        if ( mesh->BCv.type[c3-1] == 11 )   Stokes->F[eqn] += -2*AW*one_dx_dx*mesh->BCv.val[c3-1];
        //        if ( mesh->BCv.type[c3+1] == 11 )   Stokes->F[eqn] += -2*AE*one_dx_dx*mesh->BCv.val[c3+1];
        Stokes->F[eqn] -= (Stokes->b[eqn]); //
        Stokes->F[eqn] *= celvol;
        //        printf("%2.2e\n", Stokes->F[eqn]*1e-8);


    }

    //    // Coefficients
    //    double AN  = mesh->eta_n[c2+ncx];
    //    double AE  = mesh->eta_s[c1];
    //    double AW  = mesh->eta_s[c1-1];
    //
    //    // Shear terms
    //    double uSW =  one_dx_dz * AW;// - comp*2.0/3.0*AS*one_dx_dz;
    //    double uSE = -one_dx_dz * AE;// + comp*2.0/3.0*AS*one_dx_dz;
    //    double uNW = -one_dx_dz * AW;// + comp*2.0/3.0*AN*one_dx_dz;
    //    double uNE =  one_dx_dz * AE;// - comp*2.0/3.0*AN*one_dx_dz;
    //
    //    // Poisson
    //    double vW  =  one_dx_dx * AW;
    //    double vC  = -2*one_dz_dz * (2*AN) - one_dx_dx * (AE + AW);// + comp*2.0/3.0*2*AN*one_dz_dz;
    //    double vE  =  one_dx_dx * AE;
    //    double vN  =  2*one_dz_dz * 2*AN;// - comp*2.0/3.0*2*AN*one_dz_dz;
    //
    //    double rhoVz = 0.5 * (mesh->rho_s[c1] + mesh->rho_s[c1+1]);
    //
    //    vE -= rhoVz/2*one_dx*mesh->VxVz[c3];
    //    vW += rhoVz/2*one_dx*mesh->VxVz[c3];
    //
    //    // Inertial terms
    //    if (model.isinertial == 1 || model.isinertial == 2 ) {
    //        vC -= sign * rhoVz/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        //                        vN -= mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        //                        vS += mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        vE -= rhoVz/2*one_dx*mesh->VxVz[c3];
    //        vW += rhoVz/2*one_dx*mesh->VxVz[c3];
    //    }
    //
    //    // Stabilisation with density gradients
    //    if (stab==1) {
    //        double correction = - om*0.5*model.dt * model.gz * (mesh->rho_n[c2+ncx] - mesh->rho_n[c2]) * one_dz;
    //        vC += correction;
    //        correction = - om*0.5*model.dt * model.gz * (mesh->rho_s[c1] - mesh->rho_s[c1-1]) * one_dx;
    //        uNW += 0.25*correction;
    //        uNE += 0.25*correction;
    //        uSW += 0.25*correction;
    //        uSE += 0.25*correction;
    //    }
    //
    //    // Pressure gradient
    //    double pS  =   one_dz;
    //    double pN  =  -pS;
    //
    //    if ( Assemble == 1 ) {
    //        Stokes->b[eqn] *= celvol;
    //        //--------------------
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   Stokes->bbc );
    //        //--------------------
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                     &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    Stokes->bbc );
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], Stokes->bbc );
    //        //--------------------
    //        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+ncx],  &(nnzc2[ith]), 1*pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], Stokes->bbc );
    //    }
    //    else {
    //        Stokes->F[eqn]  = 1* pN*p[c2+ncx];
    //        Stokes->F[eqn] += vC*v[c3] + vW*v[c3-1] + vE*v[c3+1] + vN*v[c3+nxvz];
    //        Stokes->F[eqn] += uSW*u[c1-1] + uSE*u[c1] + uNW*u[c1+nx-1] + uNE*u[c1+nx];
    //        Stokes->F[eqn] -= (Stokes->b[eqn]);
    ////        Stokes->F[eqn] -= (Stokes->b[eqn] + Stokes->bbc[eqn]);
    //        Stokes->F[eqn] *= celvol;
    //
    //
    //    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void Zmomentum_InnerNodes( SparseMat *Stokes, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **Jtemp, double **Atemp, int *nnzc2 ) {

    double uSW=0, uSE=0, uNW=0, uNE=0, vS=0, vW=0, vC=0, vE=0, vN=0, pN=0, pS=0;

    // Coefficients
    double AS  = mesh->eta_n[c2];
    double AN  = mesh->eta_n[c2+ncx];
    double AE  = mesh->eta_s[c1];
    double AW  = mesh->eta_s[c1-1];

    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * AW - comp*2.0/3.0*AS*one_dx_dz;
        uSE = -one_dx_dz * AE + comp*2.0/3.0*AS*one_dx_dz;
    }

    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * AW + comp*2.0/3.0*AN*one_dx_dz;
        uNE =  one_dx_dz * AE - comp*2.0/3.0*AN*one_dx_dz;
    }

    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3+nxvz] !=30 ) {
        vS  =  2*one_dz_dz * AS        - comp*2.0/3.0*AS*one_dz_dz;
        vC  = -2*one_dz_dz * (AN + AS) + comp*2.0/3.0*AN*one_dz_dz + comp*2.0/3.0*AS*one_dz_dz;
        vN  =  2*one_dz_dz * AN        - comp*2.0/3.0*AN*one_dz_dz;
    }
    if ( mesh->BCv.type[c3-nxvz] == 30 && mesh->BCv.type[c3+nxvz] != 30 ) {
        vC  = -2*one_dz_dz * (AN) + comp*2.0/3.0*AN*one_dz_dz;
        vN  =  2*one_dz_dz * AN   - comp*2.0/3.0*AN*one_dz_dz;
    }
    if ( mesh->BCv.type[c3-nxvz] != 30  && mesh->BCv.type[c3+nxvz] == 30 ) {
        vS  =  2*one_dz_dz * AS   - comp*2.0/3.0*AS*one_dz_dz;
        vC  = -2*one_dz_dz * (AS) + comp*2.0/3.0*AS*one_dz_dz;
    }

    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * AW; vC += -one_dx_dx * AW;}
        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * AE; vC += -one_dx_dx * AE;}
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * AE;
        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * AE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * AW;
        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * AW;
    }

    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * AW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * AE;

    //    // Inertial terms
    //    if (model.isinertial == 1 || model.isinertial == 2 ) {
    //        vC -= sign * mesh->rhoVz[c3]/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        vN -= mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        vS += mesh->rhoVz[c3]/2*one_dz*mesh->v_in[c3];
    //        vE -= mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
    //        vW += mesh->rhoVz[c3]/2*one_dx*mesh->VxVz[c3];
    //    }

    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2]);
        double drhodx = (mesh->rho_s[c1+1]   - mesh->rho_s[c1]);
        double correction = - 1.0/om * 0.5 * model.dt * model.gz * drhodz * one_dz;
        vC += correction;
//        double drhodz  = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
//        double vC_corr = - 1.00 * om * model.dt * model.gz * drhodz;

//        om = 0.15;
//        double vC_corr = -  1.0/om * 0.5 * model.dt * model.gz * drhodz;

        // Importante trique, voire meme gigantesque!
//        if (vC+vC_corr<0.0)  vC += vC_corr;

//        printf("1.00 * om * model.dt * model.gz = %2.2e\n", 1.00 * om * model.dt * model.gz);
//        printf(" om  = %2.2e  model.dt  = %2.2e  model.gz  = %2.2e\n",  om ,  model.dt, model.gz );


//        correction = - om*0.25*model.dt * model.gz * drhodx * one_dx;
//        uNW += correction;
//        uNE += correction;
//        uSW += correction;
//        uSE += correction;
    }


//    // Stabilisation with density gradients
//    if (stab==1) {
//        double drhodz = (mesh->rho_n[c2+ncx] - mesh->rho_n[c2])*one_dz;
////        double drhodx = (mesh->rho_s[c1]   - mesh->rho_s[c1-1])*one_dx;
//        vC  += - 1.00 * om * model.dt * model.gz * drhodz;
//        //        // Non-symmetric contibution to the system of equation
//        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
//        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
//        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
//        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
//    }

    // Pressure gradient
    if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
        pS  =   one_dz;
        pN  =  -one_dz;
    }

    //    // Pressure gradient
    //    if ( mesh->BCp.type[c2] == -1 && mesh->BCp.type[c2+ncx] == -1 ) {
    //        pS  =   one_dz;
    //        pN  =  -one_dz;
    //    }
    //
    //    // Pressure gradient
    //    if ( mesh->BCp.type[c2+ncx] == 31 && mesh->BCp.type[c2] == -1 ) {
    //        pS  =   0.1*one_dz;
    //    }

    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;

    if ( Assemble == 1 ) {

        Stokes->b[eqn] *= celvol;
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      Stokes->bbc );
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   Stokes->bbc );
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], Stokes->bbc );
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2[ith]), vW*celvol, mesh->BCv.type[c3-1],    2*mesh->BCv.val[c3-1],    Stokes->bbc );
        }
        AddCoeff2( Jtemp[ith], Atemp[ith], eqn, eqn,                     &(nnzc2[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      Stokes->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    Stokes->bbc );
            else AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2[ith]), vE*celvol, mesh->BCv.type[c3+1],    2*mesh->BCv.val[c3+1],    Stokes->bbc );
        }
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], Stokes->bbc );
        //--------------------
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],      &(nnzc2[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     Stokes->bbc );
            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+ncx],  &(nnzc2[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], Stokes->bbc );
        }

        //        if ( mesh->BCp.type[c2] == -1 && mesh->BCp.type[c2+ncx] == -1) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],      &(nnzc2[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     Stokes->bbc );
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2+ncx],  &(nnzc2[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], Stokes->bbc );
        //        }
        //
        //        if ( mesh->BCp.type[c2+ncx] == 31 && mesh->BCp.type[c2] == -1 ) {
        //            AddCoeff2( Jtemp[ith], Atemp[ith], eqn, Stokes->eqn_p[c2],      &(nnzc2[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     Stokes->bbc );
        //        }

        Stokes->F[eqn] = v[c3];
    }
    else {
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
            Stokes->F[eqn]  = pS*p[c2] + pN*p[c2+ncx];
        }
        //        if ( mesh->BCp.type[c2+ncx] == 31 && mesh->BCp.type[c2] == -1 ) {
        //            Stokes->F[eqn]  = pS*p[c2];
        //        }
        Stokes->F[eqn] += vC*v[c3];
        if ( mesh->BCv.type[c3-nxvz] != 30 ) Stokes->F[eqn] += vS*v[c3-nxvz];
        if ( mesh->BCv.type[c3-1]    != 30 && mesh->BCv.type[c3-1]    != 11  ) Stokes->F[eqn] += vW*v[c3-1];
        if ( mesh->BCv.type[c3+1]    != 30 && mesh->BCv.type[c3+1]    != 11  ) Stokes->F[eqn] += vE*v[c3+1];
        if ( mesh->BCv.type[c3+nxvz] != 30 ) Stokes->F[eqn] += vN*v[c3+nxvz];
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            Stokes->F[eqn] += uSW*u[c1-1] + uSE*u[c1];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            Stokes->F[eqn] += uNW*u[c1+nx-1] + uNE*u[c1+nx];
        }
//        if ( mesh->BCv.type[c3-1] == 11 )   Stokes->F[eqn] += -2*AW*one_dx_dx*mesh->BCv.val[c3-1];
//        if ( mesh->BCv.type[c3+1] == 11 )   Stokes->F[eqn] += -2*AE*one_dx_dx*mesh->BCv.val[c3+1];
        Stokes->F[eqn] -= (Stokes->b[eqn]); //
        Stokes->F[eqn] *= celvol;
        //        printf("%2.2e\n", Stokes->F[eqn]*1e-8);


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

void FreeDecomposition( int *estart, int *eend, int *DD, int *nnzc2, int *last_eqn  ) {
    DoodzFree( DD );
    DoodzFree( estart );
    DoodzFree( eend );
    DoodzFree( nnzc2 );
    DoodzFree( last_eqn );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateDecomposition( int n_th, int **estart, int **eend, int **DD, int **nnzc2, int **last_eqn ) {

    /* Array containing the size of the block for each thread */
    *DD       = DoodzMalloc( sizeof(int) * n_th );
    *estart   = DoodzMalloc( sizeof(int) * n_th );
    *eend     = DoodzMalloc( sizeof(int) * n_th );
    *nnzc2    = DoodzCalloc( n_th, sizeof(int) );
    *last_eqn = DoodzCalloc( n_th, sizeof(int) );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DomainDecomposition( int N, int n_th, int *estart, int *eend, int *DD, int *nnzc2 ) {

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

void AllocateTempMatArrays(double ***Atemp, int ***Itemp, int ***Jtemp, int n_th, int nnz, int neq, int* DD ) {

    int ith;

    // Temporary parallel Matrix arrays
    *Atemp = DoodzMalloc( sizeof(double*) * n_th );
    *Itemp = DoodzMalloc( sizeof(int*) * n_th );
    *Jtemp = DoodzMalloc( sizeof(int*) * n_th );

    for (ith=0; ith<n_th; ith++) {
        (*Atemp)[ith] = DoodzMalloc( sizeof(double) * (int)(nnz/n_th) );
        (*Itemp)[ith] = DoodzCalloc( (neq+1), sizeof(int) );
        (*Jtemp)[ith] = DoodzMalloc( sizeof(int) * (int)(nnz/n_th) );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeTempMatArrays(double **Atemp, int **Itemp, int **Jtemp, int n_th ) {

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
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void BuildStokesOperator( grid *mesh, params model, int lev, double *p, double *u, double *v, SparseMat *Stokes, int Assemble ) {

    // FLAG
    // Assemble=1: Build the linear system of discrete equations
    // Assemble=0: Evaluate Stokes residuals

    int    cc, k, l, c1, c2, c3, nx=mesh->Nx, nz=mesh->Nz, nxvz=nx+1, nzvx=nz+1, ncx=nx-1, ncz=nz-1;
    int    nnz = 0;

    // Pre-calculate FD coefs
    double celvol    = mesh->dx*mesh->dz;
    double one_dx    = 1.0/mesh->dx;
    double one_dz    = 1.0/mesh->dz;
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    double one_dx_dz = 1.0/mesh->dx/mesh->dz;

    // Switches
    int eqn, nnzc = 0, sign = 1, comp = 0, stab = model.free_surf_stab;
//    double theta = model.free_surf_stab;

    if (model.free_surf_stab>0.0) stab = 1;
    double theta = model.free_surf_stab;

    // Decompose domain
    int n_th, N, ith;
    int *DD, *estart, *eend, *nnzc2, *last_eqn;

    // Temporary parallel matrix
    double **Atemp;
    int **Jtemp, **Itemp;

#pragma omp parallel shared(n_th)
    {
        n_th = omp_get_num_threads();
    }
#pragma omp barrier

    // Matrix initialisation
    if ( Assemble == 1 ) {

        nnz  = 5*((mesh->Nx-1) * (mesh->Nz-1)) + 2*11*(mesh->Nx * mesh->Nz);

        printf("Assembling Stokes matrix, assumming nnz_max = %d...\n", nnz);
        AllocMat( Stokes, nnz );

        for (k=0; k<Stokes->neq+1; k++) {
            Stokes->Ic[k] = 0;
        }

        mesh->p_scale = 1e-10;//1/one_dx_dz;//mesh->p_scale * mesh->p_scale* mesh->p_scale;
    }


    // Build right-hand side
    int inc=0;
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 ) {
                Stokes->F[inc] = 0.0;
                Stokes->b[inc] = 0.0;
                Stokes->b[inc] = mesh->roger_x[cc];
                inc++;
            }
        }
    }
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 ) {
                Stokes->F[inc] = 0.0;
                Stokes->b[inc] = 0.0;
                Stokes->b[inc] = mesh->roger_z[cc];
                inc++;
            }
        }
    }
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            cc = k + l*ncx;
            if ( mesh->BCp.type[cc] != 31 && mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30 ) {
                Stokes->F[inc] = 0.0;
                Stokes->b[inc] = 0.0;
                Stokes->b[inc] = mesh->rhs_p[cc];
                inc++;
            }
        }
    }

    //------------------------------------------------------------------//
    //------------------------- U-momentum -----------------------------//
    //------------------------------------------------------------------//

    // Parallel decomposition
    N = nx*nzvx;
    AllocateDecomposition(  n_th, &estart, &eend, &DD, &nnzc2, &last_eqn );
    DomainDecomposition( N, n_th, estart, eend, DD, nnzc2 );

    // Allocate parallel matrix

    AllocateTempMatArrays( &Atemp, &Itemp, &Jtemp, n_th, nnz, Stokes->neq, DD );


#pragma omp parallel shared( eend, estart, mesh, Stokes, u, v, p, nx, ncx, nzvx, nnzc2, Atemp, Jtemp, Itemp, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();

        for( c1=estart[ith]; c1<eend[ith]+1; c1++) {

            k      = mesh->kvx[c1];
            l      = mesh->lvx[c1];
            c2     = k-1 + (l-1)*ncx;
            c3     = k   + l*nxvz;

            // Get equation numbering
            if ( mesh->BCu.type[c1] != 0 && mesh->BCu.type[c1] != 30 && mesh->BCu.type[c1] != 11 && mesh->BCu.type[c1] != 13 ) {
                eqn = Stokes->eqn_u[c1];
                last_eqn[ith]   = eqn;

                if ( Assemble == 1 ) {
                    Itemp[ith][eqn] = nnzc2[ith];
                }

                //--------------------- WEST PERIODIC ---------------------//
                if ( (k==0 && (l>0 && l<nzvx-1) ) && mesh->BCu.type[c1] == -2 ) {

                    Xmomentum_WestPeriodic( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- EAST PERIODIC ---------------------//
                if ( ( k==nx-1 && (l>0 && l< nzvx-1) ) && mesh->BCu.type[c1] == -2) {

                    Xmomentum_EastPeriodic( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- SOUTH NEUMANN ---------------------//
                if ( l==0 && mesh->BCu.type[c1] == 3 ) {

                    Xmomentum_SouthNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- SOUTH DIRICHLET ---------------------//
                if ( l==0 && mesh->BCu.type[c1] == 1 ) {

                    Xmomentum_SouthDirichlet( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- NORTH NEUMANN ---------------------//
                if ( l==nzvx-1 && mesh->BCu.type[c1] == 3 ) {

                    Xmomentum_NorthNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- NORTH DIRICHLET ---------------------//
                if ( l==nzvx-1 && mesh->BCu.type[c1] == 1 ) {

                    Xmomentum_NorthDirichlet( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- WEST NEUMANN ---------------------//
                if ( mesh->BCu.type[c1] == 2 && k==0 && (l>0 && l< nzvx-1) ) {

                    Xmomentum_WestNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- EAST NEUMANN ---------------------//
                if ( mesh->BCu.type[c1] == 2 && k==nx-1 && (l>0 && l< nzvx-1) ) {

                    Xmomentum_EastNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- INNER NODES ---------------------//
                if ( l>0 && l<nzvx-1 && k>0 && k<nx-1 ) {

                    Xmomentum_InnerNodes( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );

//                    printf("%.5e ", fabs(Stokes->F[eqn]));
                }
            }
//            if ((c1%nx)>(nx-2)) printf("\n");
        }
    }

    //-----------------------------------------------------------------//

    if ( Assemble == 1 ) {
        MergeParallelMatrix( Stokes, Atemp, Jtemp, Itemp, mesh, estart, eend, &nnzc, nnzc2, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
    }

    FreeTempMatArrays( Atemp, Itemp, Jtemp, n_th );

    FreeDecomposition( estart, eend, DD, nnzc2, last_eqn );

    //------------------------------------------------------------------//
    //------------------------- V-momentum -----------------------------//
    //------------------------------------------------------------------//

    // Parallel decomposition
    N = nxvz*nz;
    AllocateDecomposition(  n_th, &estart, &eend, &DD, &nnzc2, &last_eqn );
    DomainDecomposition( N, n_th, estart, eend, DD, nnzc2 );

    // Allocate parallel matrix
    AllocateTempMatArrays( &Atemp, &Itemp, &Jtemp, n_th, nnz, Stokes->neq, DD );


#pragma omp parallel shared( eend, estart, mesh, Stokes, u, v, p, nx, ncx, nzvx, nnzc2, Atemp, Jtemp, Itemp, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();

        for( c3=estart[ith]; c3<eend[ith]+1; c3++) {

            k  = mesh->kvz[c3];
            l  = mesh->lvz[c3];
            c1 = k   + l*nx;
            c2 = k-1 + (l-1)*ncx;

            // Get equation numbering
            if ( mesh->BCv.type[c3] != 0 && mesh->BCv.type[c3] != 30  && mesh->BCv.type[c3] != 11 && mesh->BCv.type[c3] != 13 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;

                if ( Assemble == 1 ) {
                    Itemp[ith][eqn] = nnzc2[ith];
                }

                //                //--------------------- WEST NEUMANN ---------------------//
                if ( k==0 && mesh->BCv.type[c3] == 3 ) {

                    Zmomentum_WestNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- WEST DIRICHLET ---------------------//
                if ( k==0 && mesh->BCv.type[c3] == 1) {

                    Zmomentum_WestDirichlet( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- WEST PERIODIC ---------------------//
                if ( k==0 && mesh->BCv.type[c3] == -2 ) {

                    Zmomentum_WestPeriodic( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- EAST NEUMANN ---------------------//

                if ( k==nxvz-1 && mesh->BCv.type[c3] == 3 ) {

                    Zmomentum_EastNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- EAST DIRICHLET ---------------------//

                if ( k==nxvz-1 && mesh->BCv.type[c3] == 1 ) {

                    Zmomentum_EastDirichlet( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- EAST PERIODIC ---------------------//
                if ( k==nxvz-1 && mesh->BCv.type[c3] == -2 ) {

                    Zmomentum_EastPeriodic( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, l );
                }

                //--------------------- NORTH NEUMANN ---------------------//
                if (  mesh->BCv.type[c3] == 2 && l==nz-1 && (k>0 && k<nxvz-1)) {

                    Zmomentum_NorthNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- SOUTH NEUMANN ---------------------//
                if (  mesh->BCv.type[c3] == 2 ) { //&& l==0 && (k>0 && k<nxvz-1)

                    //                    printf("Res. South Neumann\n");
                    Zmomentum_SouthNeumann( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );
                }

                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {

                    Zmomentum_InnerNodes( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2 );

//                     printf("%.5e ", fabs(Stokes->F[eqn]));
                }
//                if ((c3%nxvz)>9) printf("\n");
            }
        }
    }

    //-----------------------------------------------------------------//

    if ( Assemble == 1 ) {
        MergeParallelMatrix( Stokes, Atemp, Jtemp, Itemp, mesh, estart, eend, &nnzc, nnzc2, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
    }

    FreeTempMatArrays( Atemp, Itemp, Jtemp, n_th );
    FreeDecomposition( estart, eend, DD, nnzc2, last_eqn );

    //------------------------------------------------------------------//
    //------------------------- Continuity -----------------------------//
    //------------------------------------------------------------------//

    // Parallel decomposition
    N = ncx*ncz;
    AllocateDecomposition(  n_th, &estart, &eend, &DD, &nnzc2, &last_eqn );
    DomainDecomposition( N, n_th, estart, eend, DD, nnzc2 );

    // Allocate parallel matrix
    AllocateTempMatArrays( &Atemp, &Itemp, &Jtemp, n_th, nnz, Stokes->neq, DD );


#pragma omp parallel shared( eend, estart, mesh, Stokes, u, v, p, nx, ncx, nzvx, nnzc2, Atemp, Jtemp, Itemp, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {

        ith = omp_get_thread_num();

        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {

            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c1  = k   + (l+1)*nx;
            c3  = k   + l*nxvz + 1;

            //--------------------- INNER NODES ---------------------//
            if ( mesh->BCp.type[c2] != 31 && mesh->BCp.type[c2] != 0 && mesh->BCp.type[c2] != 30 ) {

                eqn = Stokes->eqn_p[c2];
                last_eqn[ith]   = eqn;

                if ( Assemble == 1 ) {
                    Itemp[ith][eqn] = nnzc2[ith];
                }
                Continuity_InnerNodes( Stokes, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, Jtemp, Atemp, nnzc2, k, l );
            }
        }
    }

    //-----------------------------------------------------------------//

    if ( Assemble == 1 ) {
        MergeParallelMatrix( Stokes, Atemp, Jtemp, Itemp, mesh, estart, eend, &nnzc, nnzc2, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
    }

    FreeTempMatArrays( Atemp, Itemp, Jtemp, n_th );

    FreeDecomposition( estart, eend, DD, nnzc2, last_eqn );


    //------------------------------------------------------------------//
    //----------------------------- End --------------------------------//
    //------------------------------------------------------------------//



    if ( Assemble == 1 ) {

        // Add contribution from the BC's
        for (k=0; k<Stokes->neq; k++) {
            Stokes->b[k] += Stokes->bbc[k];
        }

        // Final index
        Stokes->Ic[Stokes->neq] = nnzc;

        MinMaxArrayI(Stokes->Ic, 1, Stokes->neq+1, "I" );
        MinMaxArrayI(Stokes->J, 1, nnzc, "J" );
        MinMaxArray(Stokes->A, 1, nnzc, "V" );


        // Resize arrays to proper number of non-zeros
        double *bufd;
        int *bufi;
        bufi      = DoodzRealloc(Stokes->J, nnzc*sizeof(int));
        bufd      = DoodzRealloc(Stokes->A, nnzc*sizeof(double));
        Stokes->J = bufi;
        Stokes->A = bufd;

        printf("System size: ndof = %d, nnz = %d\n", Stokes->neq, nnzc);

        //--------------------------------------//

//#ifndef _VG_
//        if ( model.write_debug == 1 ) {

            char *filename;
            asprintf( &filename, "MatrixCoupled.gzip_%dcpu.h5", n_th );

            // Fill in DD data structure
            OutputSparseMatrix OutputDD;
            OutputDD.V = Stokes->A;
            OutputDD.Ic = Stokes->Ic;
            OutputDD.J = Stokes->J;
            OutputDD.b = Stokes->b;
            OutputDD.eta_cell = mesh->eta_n;
            OutputDD.params[0] = nx;
            OutputDD.params[1] = nz;
            OutputDD.params[2] = mesh->dx;
            OutputDD.params[3] = mesh->dz;
            OutputDD.eqn_u = DoodzMalloc(nx*nzvx*sizeof(int));
            OutputDD.eqn_v = DoodzMalloc(nxvz*nz*sizeof(int));
            OutputDD.eqn_p = DoodzMalloc(ncx*ncz*sizeof(int));

            for( l=0; l<nzvx; l++) {
                for( k=0; k<nx; k++) {
                    cc = k + l*nx;
                    c1 = k + l*nx;
                    OutputDD.eqn_u[cc]=Stokes->eqn_u[cc];
                }
            }
            for( l=0; l<nz; l++) {
                for( k=0; k<nxvz; k++) {
                    cc = k + l*nxvz;
                    c3 = nx*nzvx + k + l*nxvz;
                    OutputDD.eqn_v[cc]=Stokes->eqn_v[cc];
                }
            }
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    OutputDD.eqn_p[cc]=Stokes->eqn_p[cc];
                }
            }

            // Send data to file
            create_output_hdf5( filename );
            AddGroup_to_hdf5( filename, "model" );
            AddGroup_to_hdf5( filename, "matrix" );
            AddGroup_to_hdf5( filename, "numbering" );
            AddGroup_to_hdf5( filename, "fields" );
            AddFieldToGroup_generic( _TRUE_, filename, "model", "params" , 'd', 4, OutputDD.params,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_u" , 'i', nx*nzvx, OutputDD.eqn_u,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_v" , 'i', nxvz*nz, OutputDD.eqn_v,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDD.eqn_p,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_n" , 'd', ncx*ncz, mesh->eta_n,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_s" , 'd', nx*nz, mesh->eta_s,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_s" , 'c', nx*nz, mesh->BCg.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_u" , 'c', nx*nzvx, mesh->BCu.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_v" , 'c', nxvz*nz, mesh->BCv.type,  1 );


            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "I" , 'i', Stokes->neq+1, OutputDD.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "J" , 'i', nnzc,  OutputDD.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "V" , 'd', nnzc,  OutputDD.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "eta_cell" , 'd', ncx*ncz, OutputDD.eta_cell,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs" , 'd', Stokes->neq, OutputDD.b,  1 );

            DoodzFree(OutputDD.eqn_u);
            DoodzFree(OutputDD.eqn_v);
            DoodzFree(OutputDD.eqn_p);
            free(filename);
//        }
//#endif


    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
