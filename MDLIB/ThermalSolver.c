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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "header_MDOODZ.h"
#include "umfpack.h"

#include "time.h"
#include "amd.h"
#include "cs.h"
#include "cholmod.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#define error printf

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void DirectEnergy( double *a, int *ia, int *ja, double *x, double *b, int n, int nnz  ) {
//
//    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
//    int status;
//    void *Symbolic = NULL;
//	void *Numeric = NULL;
//    int msglvl = 0;
//
//    if (msglvl==1) printf ("\nUMFPACK V%d.%d (%s) demo: _di_ version\n",  UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;
//
//    /* get the default control parameters */
//    umfpack_di_defaults (Control) ;
//
//    /* change the default print level for this demo */
//    /* (otherwise, nothing will print) */
//    if (msglvl==1) Control [UMFPACK_PRL] = 6 ;
//
//    /* print the license agreement */
//    if (msglvl==1) umfpack_di_report_status (Control, UMFPACK_OK) ;
//    if (msglvl==1) Control [UMFPACK_PRL] = 5 ;
//
//    Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
//    //    Control [UMFPACK_PIVOT_TOLERANCE] = 1;
//    Control [UMFPACK_ORDERING] = 1; // not 2
//    Control [UMFPACK_AGGRESSIVE] = 1;
//    Control [UMFPACK_IRSTEP]   = 0;
//    Control [UMFPACK_FIXQ] = 1;
//
//    /* print the control parameters */
//    if (msglvl==1) umfpack_di_report_control (Control) ;
//
//    //    int n = mat->neq;
//
//    /* ---------------------------------------------------------------------- */
//    /* symbolic factorization */
//    /* ---------------------------------------------------------------------- */
//
//    status = umfpack_di_symbolic (n, n, ia, ja, a, &Symbolic,
//                                  Control, Info) ;
//    if (status < 0)
//    {
//        umfpack_di_report_info (Control, Info) ;
//        umfpack_di_report_status (Control, status) ;
//        error ("umfpack_di_symbolic failed") ;
//    }
//
//    /* ---------------------------------------------------------------------- */
//    /* numeric factorization */
//    /* ---------------------------------------------------------------------- */
//
//    status = umfpack_di_numeric (ia, ja, a, Symbolic, &Numeric,
//                                 Control, Info) ;
//    if (status < 0)
//    {
//        umfpack_di_report_info (Control, Info) ;
//        umfpack_di_report_status (Control, status) ;
//        error ("umfpack_di_numeric failed") ;
//    }
//    umfpack_di_free_symbolic(&Symbolic);
//
//    /* ---------------------------------------------------------------------- */
//    /* solve Ax=b */
//    /* ---------------------------------------------------------------------- */
//
//    status = umfpack_di_solve (UMFPACK_Aat, ia, ja, a, x, b,
//                               Numeric, Control, Info) ;
//    if (msglvl==1)umfpack_di_report_info (Control, Info) ;
//    umfpack_di_report_status (Control, status) ;
//    if (status < 0)
//    {
//        error ("umfpack_di_solve failed") ;
//    }
//	umfpack_di_free_numeric(&Numeric);
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cs_to_cholmod_matrix1( cholmod_sparse* Acm, cs* Ac ) {

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

void copy_vec_to_cholmod_dense1( cholmod_dense* xcm, DoodzFP* x ) {

    int k;
#pragma omp parallel for shared( xcm, x ) private(k)
    for (k=0; k<xcm->nrow; k++) {
        ((double*)xcm->x)[k] = x[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cholmod_dense_to_cholmod_dense1( cholmod_dense* xcm1, cholmod_dense* xcm2 ) {

    int k;
#pragma omp parallel for shared( xcm1, xcm2 ) private(k)
    for (k=0; k<xcm1->nrow; k++) {
        ((double*)xcm1->x)[k] = ((double*)xcm2->x)[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cholmod_to_cs_matrix1( cholmod_sparse* Acm, cs* Ac ) {

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

cs_di* TransposeA( cholmod_common *c, double *a, int *ia, int *ja, int n, int nnz ) {

    cs_di  A, *Ac;
    clock_t t_omp;
    cs_di *At;

    // Prepare A
    A.nzmax = nnz;
    A.nz    = nnz;
    A.m     = n;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, ia, A.i );
    ArrayEqualArrayI( A.p, ja,  A.nzmax );
    ArrayEqualArray(  A.x, a,  A.nzmax );
    Ac  = cs_di_compress( &A );

    cs_di *AAt;
    At  = cs_di_transpose( Ac, 1 );

    DoodzFree( A.p );
    DoodzFree( A.x );
    DoodzFree( A.i );
    cs_spfree(Ac);
    // cs_spfree(At);

return At;
}


cholmod_factor* FactorEnergyCHOLMOD( cholmod_common *c, cs_di *At, double *a, int *ia, int *ja, int n, int nnz, int polar_mode ) {

    cs_di  A, *Ac, *AAt;
    clock_t t_omp;

    // Prepare A
    A.nzmax = nnz;
    A.nz    = nnz;
    A.m     = n;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, ia, A.i );
    ArrayEqualArrayI( A.p, ja,  A.nzmax );
    ArrayEqualArray(  A.x, a,  A.nzmax );
    Ac  = cs_di_compress( &A );

    // CHOLMOD
    cholmod_sparse *Acm, *Acml;
    cholmod_factor *Afact, *Lfact;

    if ( polar_mode == 1 ) {
        // AAt = A*At
        AAt = cs_di_multiply( Ac, At);
        Acm = cholmod_allocate_sparse (AAt->m, AAt->n, AAt->nzmax, 0, 1, 0, 1, c) ;
        copy_cs_to_cholmod_matrix1( Acm, AAt );
    }
    else {
        // Copy matrix to CHOLMOD context
        Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, c) ;
        copy_cs_to_cholmod_matrix1( Acm, Ac );
    }

    // Extract lower triangular part
    Acml = cholmod_copy ( Acm, -1, 1, c );

    // Factorization
    c->nmethods  = 1;
    c->method[0].ordering=CHOLMOD_AMD;
    c->postorder = 1;
//    cholmod_check_sparse ( Acml, c ) ;
    Afact = cholmod_analyze( Acml, c ) ;
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Acml, Afact, c) ;
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

    if ( Afact->is_super == 0 ) cholmod_change_factor( 1, 1, Afact->is_super, 0, 0, Afact, c );

    // free matrices
    cholmod_free_sparse( &Acm, c );
    cholmod_free_sparse( &Acml, c );

    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    if ( polar_mode == 1 ) cs_spfree(AAt);

    Lfact = cholmod_copy_factor( Afact, c );
    cholmod_free_factor ( &Afact, c) ;

    return Lfact;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveEnergyCHOLMOD( cholmod_common *c, cs_di *At, cholmod_factor *Afact, double *x, double *b, int n, int nnz, int polar_mode ) {

    int k;
    cholmod_dense *bcm, *xcm;
    double *y = DoodzCalloc(sizeof(double),n );

    // Copy matrix RHS to CHOLMOD context
    bcm  = cholmod_allocate_dense( n, 1, n, CHOLMOD_REAL, c );
    copy_vec_to_cholmod_dense1( bcm, b );

    // Back-substitutions
    xcm = cholmod_solve (CHOLMOD_A, Afact, bcm, c) ;

    // Retrieve  solution from cholmod array
    for (k=0; k<n; k++) x[k] = ((double*)xcm->x)[k];

    if ( polar_mode == 1 ) {
        cs_di_gaxpy( At, x, y );
        for (k=0; k<n; k++) x[k] = y[k];
    }

    // Free
    cholmod_free_dense( &xcm, c );
    cholmod_free_dense( &bcm, c );
    DoodzFree( y );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
