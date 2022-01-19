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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"
#include "header_MDOODZ.h"

#include "cs.h" 

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double DotProduct( DoodzFP *vec1, DoodzFP *vec2, int size ) {

    int k;
    double result = 0.0;
    double *dummy;
    dummy = DoodzCalloc( size, sizeof(DoodzFP) );
#pragma omp parallel for shared( dummy, vec1, vec2 ) private( k ) //firstprivate( celvol, nx, nzvx ) reduction(+:resx)
    for ( k=0; k<size; k++ ) {
        dummy[k] = vec1[k]*vec2[k];
    }
#pragma omp parallel for shared( dummy ) private( k ) reduction(+:result)
    for ( k=0; k<size; k++ ) {
        result += dummy[k];
    }
    DoodzFree(dummy);
    return result;
}

//double DotProduct( DoodzFP *a, DoodzFP *b, int n ) {
//
//    int i, chunk = 10;
//    double result = 0.0;
//
//
//#pragma omp parallel for      \
//default(shared) private(i)  \
//    schedule(static,chunk)      \
//    reduction(+:result)
//
//    for (i=0; i < n; i++)
//        result += (a[i] * b[i]);
//    return result;
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PermuteVector( DoodzFP *vec1, DoodzFP *vec2, int *P, int size ) {
    
    int k;
    for ( k=0; k<size; k++ ) {
        vec1[k] = vec2[P[k]];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PermuteInverseVector( DoodzFP *vec1, DoodzFP *vec2, int *P, int size ) {
    
    int k;
    for ( k=0; k<size; k++ ) {
        vec1[P[k]] = vec2[k];
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NormResidualCholmod( double *ru, double *rp, cholmod_dense *Fu, cholmod_dense *Fp, int neq_mom, int neq_cont, params model, scale scaling, int quiet ) {
    
    double celvol, Area=0.0, resu, resp;
    int cc;
    
    celvol = model.dx*model.dz;
    resu = 0.0;
    resp = 0.0;
    int mode = 1;
    
    // Integrate residuals
#pragma omp parallel for shared( Fu, mode ) private( cc ) firstprivate( celvol ) reduction(+:resu)
    for( cc=0; cc<neq_mom; cc++) {
        if ( mode == 1 ) resu += ((double*)Fu->x)[cc]*((double*)Fu->x)[cc];//*celvol;
        if ( mode == 0 ) resu += ((double*)Fu->x)[cc];
    }
    
#pragma omp parallel for shared( Fp, mode ) private( cc ) firstprivate( celvol ) reduction(+:resp, Area)
    for( cc=0; cc<neq_cont; cc++) {
        Area += celvol;
        if ( mode == 1 ) resp += ((double*)Fp->x)[cc]*((double*)Fp->x)[cc];//*celvol;
        if ( mode == 0 ) resp += ((double*)Fp->x)[cc];
    }
    
    // Sqrt
    resu =  sqrt(resu/neq_mom);
    resp =  sqrt(resp/neq_cont);
    
    *ru = resu;
    *rp = resp;
    
    if ( quiet == 0 ) {
        printf("Fu = %2.6e\n", fabs(*ru) ); // Units of momentum
        printf("Fp = %2.6e\n", fabs(*rp) ); // Units of velocity gradient
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildInitialSolutions( double* u0, double* p0, grid* mesh ) {
    
    int inc=0, k, l, kk;
    
    // Get Vx
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            if ( mesh->BCu.type[kk] != 0 && mesh->BCu.type[kk] != 30 && mesh->BCu.type[kk] != 13 && mesh->BCu.type[kk] != 11 && mesh->BCu.type[kk] != -12 ) {
                u0[inc] = mesh->u_in[kk];
                inc++;
            }
        }
    }
    
    // Get Vz
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            if ( mesh->BCv.type[kk] != 0  && mesh->BCv.type[kk] != 30 && mesh->BCv.type[kk] != 13  && mesh->BCv.type[kk] != 11 && mesh->BCv.type[kk] != -12 ) {
                u0[inc] = mesh->u_in[kk];
                inc++;
            }
        }
    }
    
    inc = 0;
    // Get P
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            if ( mesh->BCp.type[kk] != 0  && mesh->BCp.type[kk] != 30 && mesh->BCp.type[kk] != 31 ) {
                p0[inc] = mesh->p_in[kk];
                inc++;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BackToSolutionVector( cholmod_dense* u, cholmod_dense* p, double* x, grid* mesh, SparseMat* Stokes ) {
    
    int inc = 0, k, l, kk, inc_tot=0;
    
    // Get Vx
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            if ( mesh->BCu.type[kk] != 0 && mesh->BCu.type[kk] != 30 && mesh->BCu.type[kk] != 13 && mesh->BCu.type[kk] != 11 && mesh->BCu.type[kk] != -12) {
                x[inc_tot] = ((double*)u->x)[inc];
                inc++;
                inc_tot++;
            }
        }
    }
    
    // Get Vz
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            if ( mesh->BCv.type[kk] != 0  && mesh->BCv.type[kk] != 30 && mesh->BCv.type[kk] != 13  && mesh->BCv.type[kk] != 11 && mesh->BCv.type[kk] != -12  ) {
                x[inc_tot] = ((double*)u->x)[inc];
                inc++;
                inc_tot++;
            }
        }
    }
    
    inc = 0;
    // Get P
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            if ( mesh->BCp.type[kk] == -1 ) {
                x[inc_tot] = ((double*)p->x)[inc];
                
                inc++;
                inc_tot++;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void copy_cs_to_cholmod_matrix( cholmod_sparse* Acm, cs* Ac ) {
    
    int k;
#pragma omp parallel for shared( Acm, Ac ) private(k)
    for (k=0; k<Ac->n+1; k++) {
        ((int*)Acm->p)[k] = Ac->p[k];
    }
#pragma omp parallel for shared( Acm, Ac ) private(k)
    for (k=0; k<Ac->nzmax; k++) {
        ((int*)Acm->i)[k]    = Ac->i[k];
        ((double*)Acm->x)[k] = Ac->x[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_vec_to_cholmod_dense( cholmod_dense* xcm, DoodzFP* x ) {
    
    int k;
#pragma omp parallel for shared( xcm, x ) private(k)
    for (k=0; k<xcm->nrow; k++) {
        ((double*)xcm->x)[k] = x[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cholmod_dense_to_cholmod_dense( cholmod_dense* xcm1, cholmod_dense* xcm2 ) {
    
    int k;
#pragma omp parallel for shared( xcm1, xcm2 ) private(k)
    for (k=0; k<xcm1->nrow; k++) {
        ((double*)xcm1->x)[k] = ((double*)xcm2->x)[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void cholmod_dense_plus_cholmod_dense( cholmod_dense* xcm1, cholmod_dense* xcm2 ) {
    
    int k;
#pragma omp parallel for shared( xcm1, xcm2 ) private(k)
    for (k=0; k<xcm1->nrow; k++) {
        ((double*)xcm1->x)[k] += ((double*)xcm2->x)[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cholmod_to_cs_matrix( cholmod_sparse* Acm, cs* Ac ) {
    
    int k;
#pragma omp parallel for shared( Acm, Ac ) private(k)
    for (k=0; k<Ac->n+1; k++) {
        Ac->p[k] = ((int*)Acm->p)[k];
    }
#pragma omp parallel for shared( Acm, Ac ) private(k)
    for (k=0; k<Ac->nzmax; k++) {
        Ac->i[k] = ((int*)Acm->i)[k];
        Ac->x[k] = ((double*)Acm->x)[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DecompressCSRtoTriplets( int size, int* Icsr, int* Itriplets ) {
    
    int m, p, k, l, off ;
    
    // Decompress from CSR to triplets
    Itriplets[0] = 0;
    m = 0;
    p = 0;
    
    for ( k=1; k<size + 1; k++) {
        off = Icsr[k] - Icsr[k-1];
        for (l=0; l<off; l++) {
            Itriplets[m+l] = p;
        }
        m += off;
        p += 1;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
