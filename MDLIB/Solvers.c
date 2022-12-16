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
#include "time.h"
#include "mdoodz-private.h"
#include "umfpack.h"
#include "amd.h"
#include "cs.h"
#include "cholmod.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#define error printf

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void copy_cholmod_matrix_to_cs( SparseMat *matA, cs_di* A, int m, int n, int nnz, int nnz1, bool decomp ) {
    A->nzmax = nnz;
    A->nz    = nnz;
    A->m     = m;
    A->n     = n;
    A->p     = DoodzCalloc(     nnz1, sizeof(int) );
    A->i     = DoodzCalloc( A->nzmax, sizeof(int) );
    A->x     = DoodzCalloc( A->nzmax, sizeof(double) );
    if ( decomp ) {
        DecompressCSRtoTriplets( A->m, matA->Ic, A->i );
        ArrayEqualArrayI( A->p, matA->J,  A->nzmax );
        ArrayEqualArray(  A->x, matA->A,  A->nzmax );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CheckArrays( cholmod_dense* fu, cholmod_dense* fp, cholmod_dense* bu, cholmod_dense* bp, grid* mesh, SparseMat* Stokes ) {
    
    int inc = 0, k, l, kk, inc_tot=0;
    
    // Get Vx
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            if ( mesh->BCu.type[kk] != 0 && mesh->BCu.type[kk] != 30 && mesh->BCu.type[kk] != 13 && mesh->BCu.type[kk] != 11 && mesh->BCu.type[kk] != -12) {
                if ( fabs( ((double*)fu->x)[inc]) > 1e5 ) {
                    printf("CHK X --- %2.2e --- %2.2e\n", ((double*)fu->x)[inc], ((double*)bu->x)[inc]);
                }

//                x[inc_tot] = ((double*)u->x)[inc];
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
                if ( fabs( ((double*)fu->x)[inc]) > 1e5 ) {
                    printf("CHK Y --- %2.2e --- %2.2e\n\n", ((double*)fu->x)[inc], ((double*)bu->x)[inc]);
                }
//                x[inc_tot] = ((double*)u->x)[inc];
                inc++;
                inc_tot++;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ArrayTimesArray( double* arr1, int* scalar, double* arr2, int size ) {
    int k;
#pragma omp parallel for shared( arr1, arr2, scalar) private(k) schedule( static )
    for(k=0;k<size;k++) {
        arr1[k] *= (double)scalar[k]*arr2[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void kspgcr( cholmod_sparse *M, cholmod_dense *b, cholmod_dense *x, cholmod_factor *Lfact, int N, cholmod_common *c, double eps, int noisy, int *its_tot) {
    
    // Initialise KSP
    int restart = 25, cc;//6;
    int max_it  = 1000;
    int ncycles = 0;
    int its     = 0, i1, i2, success=0;
    double  norm_r, rnorm0, fact, r_dot_v, nrm;
    double **VV, **SS;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    cholmod_dense *f, *s, *v, *val, *Y, *E;
 
    f      = cholmod_zeros (N, 1, CHOLMOD_REAL, c );
    v      = cholmod_zeros (       N,     1, CHOLMOD_REAL, c );
    val    = cholmod_zeros ( restart,     1, CHOLMOD_REAL, c );
    s      = cholmod_zeros (N, 1, CHOLMOD_REAL, c );
    Y      = cholmod_zeros (N, 1, CHOLMOD_REAL, c );
    E      = cholmod_zeros (N, 1, CHOLMOD_REAL, c );
       
    // Allocate temp vecs
    VV  = (double**) DoodzCalloc( restart, sizeof(double*));
    SS  = (double**) DoodzCalloc( restart, sizeof(double*));
    for (cc=0; cc<restart; cc++) {
        VV[cc]  = (double*) DoodzCalloc( N, sizeof(double));
        SS[cc]  = (double*) DoodzCalloc( N, sizeof(double));
    }

    copy_cholmod_dense_to_cholmod_dense( f, b );
    cholmod_sdmult ( M, 0, mone, one, x, f, c) ;
    norm_r = cholmod_norm_dense ( f, 2, c );
    
    // Evaluate reference norm
    rnorm0 = norm_r;
    if (noisy > 1) printf("       %1.4d KSP GCR Residual %1.12e %1.12e\n", 0, norm_r, norm_r/rnorm0);
    
    while ( success == 0 && its<max_it ) {
        
        for ( i1=0; i1<restart; i1++ ) {
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            cholmod_solve2(CHOLMOD_A, Lfact, f, NULL, &s, NULL, &Y, &E, c);
            
            // Simplest preconditionning: PC = I
            //            ArrayEqualArray( s->x, f->x, N );
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            
            // Action of Jacobian on s
            cholmod_sdmult ( M, 0, one, zero, s, v, c) ;
            
            // Approximation of the Jv product
            for (i2=0; i2<i1; i2++) {
                ((double*)val->x)[i2] = DotProduct( ((double*)v->x), VV[i2], N );
            }
            for(i2=0; i2<i1+1; i2++) {
                fact = -((double*)val->x)[i2];
                ArrayPlusScalarArray( ((double*)v->x), fact, VV[i2], N );
                ArrayPlusScalarArray( ((double*)s->x), fact, SS[i2], N );
            }
            
            r_dot_v = DotProduct( ((double*)f->x), ((double*)v->x), N );
            nrm     = sqrt(  DotProduct( ((double*)v->x), ((double*)v->x), N )  );
            r_dot_v = r_dot_v / nrm;
            
            fact = 1.0/nrm;
            ArrayTimesScalar( ((double*)v->x), fact, N );
            ArrayTimesScalar( ((double*)s->x), fact, N );
            
            fact = r_dot_v;
            ArrayPlusScalarArray( ((double*)x->x), fact, ((double*)s->x), N );
            fact =-r_dot_v;
            ArrayPlusScalarArray( ((double*)f->x), fact, ((double*)v->x), N );
            
            // Check convergence
            norm_r = cholmod_norm_dense ( f, 2, c );
            if (norm_r < eps * rnorm0 ) { // || norm_r < eps 
                success = 1;
                break;
            }
            if (noisy>1) printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
            
            // Store arrays
            ArrayEqualArray( VV[i1], ((double*)v->x), N );
            ArrayEqualArray( SS[i1], ((double*)s->x), N );

            its++;
        }
        its++;
        ncycles++;
    }
    if (noisy>1) printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
    
    cholmod_free_dense( &Y,   c );
    cholmod_free_dense( &E,   c );
    cholmod_free_dense( &s,   c );
    cholmod_free_dense( &f,   c );
    cholmod_free_dense( &v,   c );
    cholmod_free_dense( &val, c );
    
    // Free temp vecs
    for (cc=0; cc<restart; cc++) {
        DoodzFree(VV[cc]);
        DoodzFree(SS[cc]);
    }
    DoodzFree(VV);
    DoodzFree(SS);
    *its_tot = its;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokes( SparseMat *mat, DirectSolver *pardi, double *rhs, double *sol ) {
    
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
    int status;
    void *Symbolic = NULL;
    void *Numeric = NULL;
    int msglvl = 0;
    
    clock_t t_omp;
    
    t_omp = (double)omp_get_wtime();
    
    if (msglvl==1) printf ("\nUMFPACK V%d.%d (%s) demo: _di_ version\n",  UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;
    
    /* get the default control parameters */
    umfpack_di_defaults (Control) ;
    
    /* change the default print level for this demo */
    /* (otherwise, nothing will print) */
    if (msglvl==1) Control [UMFPACK_PRL] = 6 ;
    
    /* print the license agreement */
    if (msglvl==1) umfpack_di_report_status (Control, UMFPACK_OK) ;
    if (msglvl==1) Control [UMFPACK_PRL] = 5 ;
    
    Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    //    Control [UMFPACK_PIVOT_TOLERANCE] = 1;
    Control [UMFPACK_ORDERING] = 1; // not 2
    Control [UMFPACK_AGGRESSIVE] = 1;
    Control [UMFPACK_IRSTEP]   = 0;
    Control [UMFPACK_FIXQ] = 1;
    
    /* print the control parameters */
    if (msglvl==1) umfpack_di_report_control (Control) ;
    
    int n = mat->neq;
    
    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_symbolic (n, n, mat->Ic, mat->J, mat->A, &Symbolic,
                                  Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
        error ("umfpack_di_symbolic failed") ;
    }
    
    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_numeric (mat->Ic, mat->J, mat->A, Symbolic, &Numeric,
                                 Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
        error ("umfpack_di_numeric failed") ;
    }
    umfpack_di_free_symbolic(&Symbolic);
    
    /* ---------------------------------------------------------------------- */
    /* solve Ax=b */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_solve (UMFPACK_Aat, mat->Ic, mat->J, mat->A, sol, rhs,
                               Numeric, Control, Info) ;
    if (msglvl==1)umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    if (status < 0)
    {
        error ("umfpack_di_solve failed") ;
    }
    umfpack_di_free_numeric(&Numeric);
    
    printf("** Time for UMFPACK = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokesDecoupled( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes ) {
    
    double  celvol=model.dx*model.dz ;
    cs_di  A, B, D, C, *B1, *L, *Ac, *Bc,  *Cc,  *Dc, *L1;
    DoodzFP  *u0, *p0;
    int  noisy=1;
    int nitmax=20, k;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    DoodzFP ru, rp;
    double maxdiv0, mindiv, maxdiv, maxdivit=0, rel_tol_div=model.rel_tol_div;
    double minru0, maxru0, minru, maxru;
    
    cholmod_common c ;
    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm, *Dcm0;
    cholmod_factor *Lfact ;
    
    //    cholmod_start( &c ) ;
    //    pardi->Analyze = 1;
    
    Lfact = pardi->Lfact;
    c     = pardi->c;
    
    // --------------- Pre-process--------------- //
    
    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma  = model.penalty, penalty, min_eta, max_eta, min_G, max_G;
    rsize = Stokes->neq_cont;
    Dcm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_G  , &max_G   );
    penalty = min_G*model.dt *model.dt * celvol;
    //    printf("min eta = %2.2e, max eta = %2.2e, penalty try=%2.2e, penalty =%2.2e, celvol=%2.2e...\n", min_eta, max_eta, penalty, -gamma * celvol, celvol);
    
    if (model.auto_penalty == 1) {
        printf("Trying AUTO_PENALTY...\n");
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
        
        //        penalty = -min_eta*model.dt * celvol;
        //        penalty = -1/(min_G*model.dt *model.dt);
        //        penalty = -min_eta/model.dt / celvol;
        //         printf("min eta = %2.2e, max eta = %2.2e, penalty try=%2.2e, penalty try=%2.2e...\n", min_eta, max_eta, penalty, -gamma * celvol);
    }
    else {
        penalty = gamma / celvol;
        printf("Penalty factor = %2.2e\n", penalty);
    }
    
    // Here Dcm0 is the inverse of the pressure block
    for (k=0;k<Dcm0->nzmax;k++) ((double*)Dcm0->x)[k] *= penalty; // Should be /celvol
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();
    
    // Build initial solution vector
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    BuildInitialSolutions( u0, p0, mesh );
    
    // Prepare A
    A.nzmax = matA->nnz;
    A.nz    = matA->nnz;
    A.m     = matA->neq;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, matA->Ic, A.i );
    ArrayEqualArrayI( A.p, matA->J,  A.nzmax );
    ArrayEqualArray(  A.x, matA->A,  A.nzmax );
    Ac  = cs_di_compress( &A );
    //    cs_droptol( Ac, 1.0e-18 );
    
    // Prepare B
    B.nzmax = matB->nnz;
    B.nz    = matB->nnz;
    B.m     = matB->neq;
    B.n     = matC->neq;
    B.p     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.i     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.x     = DoodzCalloc( B.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( B.m, matB->Ic, B.i );
    ArrayEqualArrayI( B.p, matB->J,  B.nzmax );
    ArrayEqualArray(  B.x, matB->A,  B.nzmax );
    Bc  = cs_di_compress( &B );
    
    // Prepare C
    C.nzmax = matC->nnz;
    C.nz    = matC->nnz;
    C.m     = matC->neq;
    C.n     = matB->neq;
    C.p     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.i     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.x     = DoodzCalloc( C.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( C.m, matC->Ic, C.i );
    ArrayEqualArrayI( C.p, matC->J,  C.nzmax );
    ArrayEqualArray(  C.x, matC->A,  C.nzmax );
    Cc  = cs_di_compress( &C );
    
    // Prepare D
    D.nzmax = Stokes->neq_cont;
    D.nz    = Stokes->neq_cont;
    D.m     = Stokes->neq_cont;
    D.n     = Stokes->neq_cont;
    D.p     = DoodzCalloc( Stokes->neq_cont+1, sizeof(int) );
    D.i     = DoodzCalloc( Stokes->neq_cont, sizeof(int) );
    D.x     = DoodzCalloc( Stokes->neq_cont, sizeof(double) );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    //    cholmod_free_sparse( &Dcm0, &c );
    Dc  = cs_di_compress( &D );
    
    // --------------- Schur complement --------------- //
    
    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );
    
    //----- test - in case C' != B (assume B is deficient)
    B1 = cs_di_transpose( Cc, 1);
    
    // minus sign: B = -C'
    for (k=0;k<B1->nzmax;k++) {
        B1->x[k] *= -1.0;
    }
    //-----
    
    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply( B1, L);
    cs_spfree(L);
    
    // Matrix addition: L = A - B*(D*C)
    L = cs_di_add( Ac, L1, 1, -1);
    cs_spfree(L1);
    
    // --------------- Cholesky --------------- //
    
    // LL' Cholesky
    Lcm = cholmod_allocate_sparse (L->m, L->n, L->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, L );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    Bcm = cholmod_allocate_sparse (B1->m, B1->n, B1->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, B1 );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );

    // Keep lower part
    Lcml = cholmod_copy ( Lcm, -1, 1, &c );
    
    if (Lcml == NULL || Lcml->stype == 0)
    {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }
    
    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        
    }
    
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        
    cholmod_dense *u, *p, *bu, *bp, *pdum, *udum, *fu, *fp, *E, *Y;
    pdum = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    udum = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    u    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    E    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    Y    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    p    = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    fu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    fp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    
    copy_vec_to_cholmod_dense( u, u0 );
    copy_vec_to_cholmod_dense( p, p0 );
    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    printf("Initial residual:\n");
    copy_cholmod_dense_to_cholmod_dense( fu, bu );   // fu = bu
    cholmod_sdmult ( Acm, 0, mone, one, u, fu, &c) ; // fu -= A*u
    cholmod_sdmult ( Bcm, 0, mone, one, p, fu, &c) ; // fu -= B*p
    
    copy_cholmod_dense_to_cholmod_dense( fp, bp );   // fp = bp
    cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ; // fp -= C*u
    MinMaxArrayVal( fu->x, matA->neq, &minru0, &maxru0 );
    NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
    MinMaxArray( bu->x, 1, matA->neq, "u ini");
    MinMaxArray( bp->x, 1, matC->neq, "p ini");
    
    for ( k=0; k<nitmax; k++) {
        
        cholmod_sdmult ( Dcm, 0, one, zero, bp, pdum, &c) ;   // pdum <-- Dcm * bp
        copy_cholmod_dense_to_cholmod_dense( udum, bu );      // udum <-- bu
        cholmod_sdmult ( Bcm, 0, mone, one, p, udum, &c) ;    // udum <-- bu - B*p
        cholmod_sdmult ( Bcm, 0, mone, one, pdum, udum, &c) ; // udum <-- bu - B*p - B*(D*bp)
        
        cholmod_free_dense( &u, &c );
        //        t_omp = (double)omp_get_wtime;
        cholmod_solve2(CHOLMOD_A, Lfact, udum, NULL, &u, NULL, &Y, &E, &c);
        //        printf( "substitution: %lf s\n", (double)((double)omp_get_wtime() - t_omp));
        
        copy_cholmod_dense_to_cholmod_dense( pdum, bp );      // pdum <-- bp
        cholmod_sdmult ( Ccm, 0, mone, one, u, pdum, &c);     // pdum <-- bp - C*u
        cholmod_sdmult ( Dcm, 0, one, one, pdum, p, &c) ;     // p <-- p + D*(bp - C*u)
        
        copy_cholmod_dense_to_cholmod_dense( fp, bp );
        cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ;
        //MinMaxArray(fp->x, scaling.E, matC->neq, "divU" );
        
        if (k>0) maxdivit = maxdiv;
        MinMaxArrayVal( fp->x, matC->neq, &mindiv, &maxdiv );
        MinMaxArrayVal( fu->x, matA->neq, &minru , &maxru  );
        
        if (k==0) maxdiv0 = maxdiv;
        
        if ( noisy == 1 ) {
            //            printf("PH iteration %01d:\n", k+1 );
            copy_cholmod_dense_to_cholmod_dense( fu, bu  );
            cholmod_sdmult ( Acm, 0, mone, one, u, fu, &c);
            cholmod_sdmult ( Bcm, 0, mone, one, p, fu, &c);
            copy_cholmod_dense_to_cholmod_dense( fp, bp );
            cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c);
            printf("PH it. %01d. rel. max. div. = %2.2e -- max. abs. div = %2.2e-- max. rel. div = %2.2e  / max. mom. = %2.2e rel. max. mom. = %2.2e\n", k, fabs(maxdiv/maxdiv0), maxdiv, maxdivit, maxru, maxru/maxru0);
            
            //cholmod_sdmult ( Dcm, 0, mone, one, p, fp, &c) ;
            //            NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
        }
        //        if (fabs(maxdiv/maxdiv0)<rel_tol_div) break;
        //        if (k>0 && fabs(maxdiv)/fabs(maxdivit)>0.75) break;
        if (k>1 && (fabs(maxdiv)<model.abs_tol_div || maxdiv/maxdiv0<rel_tol_div ) ) break;
        
        //        if ( ru<1e-11/(scaling.F/pow(scaling.L,3)) && rp<1e-11/scaling.E) break;
    }
    
    //    if (fabs(maxdiv/maxdiv0)>rel_tol_div || fabs(maxdiv) > rel_tol_div ){
    if (fabs(maxdiv)>model.abs_tol_div && maxdiv/maxdiv0>rel_tol_div) {
        printf("The code has exited since the incompressibility constrain was not satisfied to abs. tol. = %2.2e and rel. tol. = %2.2e\n Try modifying the PENALTY factor or check MIN/MAX viscosities\n Good luck!\n", model.abs_tol_div, rel_tol_div);
        exit(1);
    }
    printf("** PH - iterations = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    SumArray(u->x, 1.0, matC->neq, "u");
    SumArray(p->x, 1.0, matC->neq, "p");
    
    // --------------- Solution vector --------------- //
    BackToSolutionVector( u, p, sol, mesh, Stokes );
    
    // separate residuals
    double *F, Area =0.0, resx=0.0,resz=0.0, resp=0.0;
    int ndofx=0, ndofz=0, ndofp=0, cc;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    F = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    
    BackToSolutionVector( fu, fp, F, mesh, Stokes );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]];
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]];
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]];
        }
    }
    
    // Sqrt
    resx =  sqrt(resx/ndofx);
    resz =  sqrt(resz/ndofz);
    resp =  sqrt(resp/ndofp);
    
    if ( noisy == 1 ) {
        printf("Fu = %2.6e\n", resx ); // Units of momentum
        printf("Fv = %2.6e\n", resz ); // Units of momentum
        printf("Fp = %2.6e\n", resp ); // Units of velocity gradient
    }
    
    // --------------- Free --------------- //
    
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &u, &c );
    cholmod_free_dense( &E, &c );
    cholmod_free_dense( &Y, &c );
    cholmod_free_dense( &p, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );
    cholmod_free_sparse( &Lcm, &c );
    
    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm, &c );
 
    pardi->Lfact = Lfact;
    pardi->c     = c;
    
    cs_spfree(L);
    
    DoodzFree(u0);
    DoodzFree(p0);
    
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    cs_spfree(B1);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);
    
    DoodzFree(F);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokesDecoupledComp( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes ) {
    
    double  celvol = model.dx*model.dz ;
    //    double  celvol = 1.0 ;
    cs_di  A, B, D, C, *B1, *L, *Ac, *Bc,  *Cc,  *Dc, *L1;
    DoodzFP  *u0, *p0;
    int  noisy=1;
    int nitmax=20, k, i;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    DoodzFP ru, rp;
    double maxdiv0, mindiv, maxdiv, maxdivit=0, rel_tol_div=model.rel_tol_div;
    double minru0, maxru0, minru, maxru;
    
    cholmod_common c ;
    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm, *Dcm0, *D1cm, *D1cm0;
    cholmod_factor *Lfact ;
    
    Lfact = pardi->Lfact;
    c     = pardi->c;
    
    // --------------- Pre-process--------------- //
    
    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma  = model.penalty, penalty, min_eta, max_eta, min_G, max_G;
    rsize = Stokes->neq_cont;
    Dcm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    D1cm0 = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_G  , &max_G   );
    penalty = min_G*model.dt *model.dt * celvol;
    
    if (model.auto_penalty == 1) {
        printf("Trying AUTO_PENALTY...\n");
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
    }
    else {
        penalty = gamma / celvol;
        printf("Penalty factor = %2.2e\n", penalty);
    }
    
#pragma omp parallel for shared(D1cm0, Dcm0, mesh, Stokes, matA, matD ) private( i ) firstprivate( model, celvol )
    for( k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
        if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31 ) {
            i = Stokes->eqn_p[k] - matA->neq;
            // Here Dcm0 is the pressure block
            if (mesh->comp_cells[k]==0) ((double*)D1cm0->x)[i] *= 0.0;
            if (mesh->comp_cells[k]==1) ((double*)D1cm0->x)[i]  = mesh->bet_n[k] / model.dt * celvol * matD->d[k]*matD->d[k];
            // Here Dcm0 is the inverse of the pressure block
            if (mesh->comp_cells[k]==0) ((double*)Dcm0->x)[i] *= penalty; // Should be /celvol
            if (mesh->comp_cells[k]==1) ((double*)Dcm0->x)[i]  = 1.0 /  ((double*)D1cm0->x)[k]; // Should be /celvol
        }
    }
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();
    
    // Build initial solution vector
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    BuildInitialSolutions( u0, p0, mesh );
    
    // Prepare A
    A.nzmax = matA->nnz;
    A.nz    = matA->nnz;
    A.m     = matA->neq;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, matA->Ic, A.i );
    ArrayEqualArrayI( A.p, matA->J,  A.nzmax );
    ArrayEqualArray(  A.x, matA->A,  A.nzmax );
    Ac  = cs_di_compress( &A );
    //    cs_droptol( Ac, 1.0e-18 );
    
    // Prepare B
    B.nzmax = matB->nnz;
    B.nz    = matB->nnz;
    B.m     = matB->neq;
    B.n     = matC->neq;
    B.p     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.i     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.x     = DoodzCalloc( B.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( B.m, matB->Ic, B.i );
    ArrayEqualArrayI( B.p, matB->J,  B.nzmax );
    ArrayEqualArray(  B.x, matB->A,  B.nzmax );
    Bc  = cs_di_compress( &B );
    
    // Prepare C
    C.nzmax = matC->nnz;
    C.nz    = matC->nnz;
    C.m     = matC->neq;
    C.n     = matB->neq;
    C.p     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.i     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.x     = DoodzCalloc( C.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( C.m, matC->Ic, C.i );
    ArrayEqualArrayI( C.p, matC->J,  C.nzmax );
    ArrayEqualArray(  C.x, matC->A,  C.nzmax );
    Cc  = cs_di_compress( &C );
    
    // Prepare D
    D.nzmax = Stokes->neq_cont;
    D.nz    = Stokes->neq_cont;
    D.m     = Stokes->neq_cont;
    D.n     = Stokes->neq_cont;
    D.p     = DoodzCalloc( Stokes->neq_cont+1, sizeof(int) );
    D.i     = DoodzCalloc( Stokes->neq_cont, sizeof(int) );
    D.x     = DoodzCalloc( Stokes->neq_cont, sizeof(double) );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    Dc  = cs_di_compress( &D );
    
    // --------------- Schur complement --------------- //
    
    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );
    
    //----- test - in case C' != B (assume B is deficient)
    B1 = cs_di_transpose( Cc, 1);
    
    // minus sign: B = -C'
    for (k=0;k<B1->nzmax;k++) {
        B1->x[k] *= -1.0;
    }
    //-----
    
    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply( B1, L);
    cs_spfree(L);
    
    // Matrix addition: L = A - B*(D*C)
    L = cs_di_add( Ac, L1, 1, -1);
    cs_spfree(L1);
    
    // --------------- Cholesky --------------- //
    
    // LL' Cholesky
    Lcm = cholmod_allocate_sparse (L->m, L->n, L->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, L );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    Bcm = cholmod_allocate_sparse (B1->m, B1->n, B1->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, B1 );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );
    
    // Keep lower part
    Lcml = cholmod_copy ( Lcm, -1, 1, &c );
    
    if (Lcml == NULL || Lcml->stype == 0)
    {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }
    
    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        
    }
    
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    cholmod_dense *u, *p, *du, *dp, *bu, *bp, *pdum, *udum, *fu, *fp, *E, *Y;
    pdum = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    udum = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    u    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    E    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    Y    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    p    = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    du   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    dp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    fu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    fp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    
    copy_vec_to_cholmod_dense( u, u0 );
    copy_vec_to_cholmod_dense( p, p0 );
    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    // Residuals
    copy_cholmod_dense_to_cholmod_dense( fu, bu );    // fu = bu
    cholmod_sdmult( Acm, 0, mone, one, u, fu, &c) ;   // fu -= A*u
    cholmod_sdmult( Bcm, 0, mone, one, p, fu, &c) ;   // fu -= B*p
    copy_cholmod_dense_to_cholmod_dense( fp, bp );    // fp = bp
    cholmod_sdmult( Ccm, 0, mone, one, u, fp, &c) ;   // fp -= C*u
    cholmod_sdmult( D1cm0, 0, mone, one, p, fp, &c) ; // fp -= D*p
    
    // Check: Initial residuals
    printf("Initial residual:\n");
    MinMaxArrayVal( fu->x, matA->neq, &minru0, &maxru0 );
    NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
    
    for ( k=0; k<nitmax; k++) {
        
        // Residuals
        copy_cholmod_dense_to_cholmod_dense( fu, bu );       // fu = bu
        cholmod_sdmult( Acm, 0, mone, one, u, fu, &c) ;      // fu -= A*u
        cholmod_sdmult( Bcm, 0, mone, one, p, fu, &c) ;      // fu -= B*p
        copy_cholmod_dense_to_cholmod_dense( fp, bp );       // fp = bp
        cholmod_sdmult( Ccm, 0, mone, one, u, fp, &c) ;      // fp -= C*u
        cholmod_sdmult( D1cm0, 0, mone, one, p, fp, &c) ;    // fp -= D*p
                
        // Check: Stopping criteria
        MinMaxArrayVal( fp->x, matC->neq, &mindiv, &maxdiv );
        MinMaxArrayVal( fu->x, matA->neq, &minru , &maxru  );
        if (k==0) maxdiv0 = maxdiv;
        if ( noisy == 1 ) printf("PH comp it. %01d. rel. max. div. = %2.2e -- max. abs. div = %2.2e-- max. rel. div = %2.2e  / max. mom. = %2.2e rel. max. mom. = %2.2e\n", k, fabs(maxdiv/maxdiv0), maxdiv, maxdivit, maxru, maxru/maxru0);
        if (k>1 && (fabs(maxdiv)<model.abs_tol_div || maxdiv/maxdiv0<rel_tol_div ) ) break;

        // Schur complement
        cholmod_sdmult( Dcm, 0, one, zero, fp, pdum, &c) ;   // pdum <-- D * fp
        copy_cholmod_dense_to_cholmod_dense( udum, fu );     // udum <-- fu
        cholmod_sdmult( Bcm, 0, mone, one, pdum, udum, &c) ; // udum <-- bu - B*(D*fp)
        
        // Solve for velocity correction
        cholmod_solve2(CHOLMOD_A, Lfact, udum, NULL, &du, NULL, &Y, &E, &c);
        
        // Solve for pressure correction
        copy_cholmod_dense_to_cholmod_dense( pdum, fp );       // pdum <-- bp
        cholmod_sdmult ( Ccm, 0, mone,  one,   du, pdum, &c);  // pdum <-- bp - C*u
        cholmod_sdmult ( Dcm, 0,  one, zero, pdum,   dp, &c);  // dp <-- D*(bp - C*u)
        
        // Velocity and pressure updates
        cholmod_dense_plus_cholmod_dense( u, du );
        cholmod_dense_plus_cholmod_dense( p, dp );
        
        // // ----- Original MATLAB version ----- //
        //        ru       = -(        K*u +        grad*p  - Fu ); resu = norm(ru)/length(ru);
        //        rp       = -(      div*u +   compD.*PP*p  - Fp ); resp = norm(rp)/length(rp);
        //        rusc     = ru - grad*(PPi*rp);
        //        du(sA)   = cs_ltsolve(AsccP,cs_lsolve(AsccP,rusc(sA)));
        //        dp       = PPi*(rp-div*du);
        //        u        = u + du; % Updates
        //        p        = p + dp;
        // // ----- Original MATLAB version ----- //
    }
    if (fabs(maxdiv)>model.abs_tol_div && maxdiv/maxdiv0>rel_tol_div) {
        printf("The code has exited since the incompressibility constrain was not satisfied to abs. tol. = %2.2e and rel. tol. = %2.2e\n Try modifying the PENALTY factor or check MIN/MAX viscosities\n Good luck!\n", model.abs_tol_div, rel_tol_div);
        exit(1);
    }
    printf("** PH - iterations = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    SumArray(u->x, 1.0, matC->neq, "u");
    SumArray(p->x, 1.0, matC->neq, "p");
    
    // --------------- Solution vector --------------- //
    BackToSolutionVector( u, p, sol, mesh, Stokes );
    
    // separate residuals
    double *F, Area =0.0, resx=0.0,resz=0.0, resp=0.0;
    int ndofx=0, ndofz=0, ndofp=0, cc;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    F = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    
    BackToSolutionVector( fu, fp, F, mesh, Stokes );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]];//*celvol;
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]];//*celvol;
            //            if ( mesh->BCv.type[cc] == 2) printf("F=%2.2e\n", F[Stokes->eqn_v[cc]]);
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]];//*celvol;
        }
    }
    
    // Sqrt
    resx =  sqrt(resx/ndofx);
    resz =  sqrt(resz/ndofz);
    resp =  sqrt(resp/ndofp);
    
    if ( noisy == 1 ) {
        printf("Fu = %2.6e\n", resx ); // Units of momentum
        printf("Fv = %2.6e\n", resz ); // Units of momentum
        printf("Fp = %2.6e\n", resp ); // Units of velocity gradient
    }
    
    // --------------- Free --------------- //
    
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &u, &c );
    cholmod_free_dense( &E, &c );
    cholmod_free_dense( &Y, &c );
    cholmod_free_dense( &p, &c );
    cholmod_free_dense( &du, &c );
    cholmod_free_dense( &dp, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );
    cholmod_free_sparse( &Lcm, &c );
    
    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm,  &c );
    cholmod_free_sparse( &Dcm0, &c );
    cholmod_free_sparse( &D1cm0, &c );

    pardi->Lfact = Lfact;
    pardi->c     = c;
    
    cs_spfree(L);
    DoodzFree(u0);
    DoodzFree(p0);
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    cs_spfree(B1);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);
    DoodzFree(F);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void KillerSolver( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes, SparseMat *Jacobian, SparseMat *JmatA,  SparseMat *JmatB,  SparseMat *JmatC ) {

    cs_di  A, B, D, C, *L, *Ac, *Bc,  *Cc,  *Dc, *L1, *L2;
    cs_di  AJ, BJ, CJ, *AJc, *BJc, *CJc;
    cs_di  *PC, *Jt, *Jts, *Js;
    DoodzFP  *u0, *p0, *F;
    int  noisy = 1, pc_type = model.pc_type, nitmax = 20, k, cc, i;
    double celvol = model.dx*model.dz;
    double maxdiv0, mindiv, maxdiv, maxdivit=0;

    cholmod_common c;
    cholmod_sparse *Lcm, *Kcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm; 
    cholmod_sparse *AcmJ, *BcmJ, *CcmJ;
    cholmod_sparse *Dcm0, *D1cm0, *whos_incompressible; 
    cholmod_factor *Lfact ;
    cholmod_sparse *M, *AB, *CD, *Iu, *Ip, *D_zero; 
    cholmod_dense *b, *x, *f, *val, *s, *v;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};

    Lfact = pardi->Lfact;
    c     = pardi->c;

    // --------------- Pre-process--------------- //
    printf("Killer solver...\n");
    printf("Preparing Matrices...\n");

    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma        = model.penalty, penalty;
    int vol_change = model.VolChangeReac;

    rsize               = Stokes->neq_cont;
    Dcm0                = cholmod_speye( rsize, rsize, CHOLMOD_REAL, &c );
    D1cm0               = cholmod_speye( rsize, rsize, CHOLMOD_REAL, &c );
    whos_incompressible = cholmod_speye( rsize, rsize, CHOLMOD_REAL, &c );

    penalty = gamma / celvol;
    printf("Penalty factor = %2.2e\n", penalty);

#pragma omp parallel for shared( whos_incompressible, D1cm0, Dcm0, mesh, Stokes, matA, matD ) private( i ) firstprivate( model, celvol, vol_change )
    for( k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
        if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31 ) {
            i = Stokes->eqn_p[k] - matA->neq;
            // Here Dcm0 is the pressure block - This relates to physics (0 is incompressible, Beta/dt is compressible)
            if (mesh->comp_cells[k]==0) ((double*)D1cm0->x)[i] *= 0.0;
            if (mesh->comp_cells[k]==1 && vol_change == 0 ) ((double*)D1cm0->x)[i]  = mesh->bet_n[k] / model.dt * celvol * matD->d[i]*matD->d[i];
            if (mesh->comp_cells[k]==1 && vol_change == 1 ) ((double*)D1cm0->x)[i]  = mesh->drhodp_n[k] / (mesh->rho_n[k]*model.dt) * celvol * matD->d[i]*matD->d[i];
            // printf("%2.2e %2.2e\n", mesh->bet_n[k], mesh->drhodp_n[k]);
            // Here Dcm0 is the inverse of the pressure block - This relates to numerics in this incompressible case (penalty) or physics in the compressible case (dt/Beta)
            if (mesh->comp_cells[k]==0) ((double*)Dcm0->x)[i]  *= penalty;
            if (mesh->comp_cells[k]==1) ((double*)Dcm0->x)[i]   = 1.0 /  ((double*)D1cm0->x)[i];
            // Detect cell which are compressible
            if (mesh->comp_cells[k]==1) ((double*)whos_incompressible->x)[i] *= 0.0;
        }
    }
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();

    // Build initial solution vector
    F  = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    // BuildInitialSolutions( u0, p0, mesh );

    //------------------------------------------------------------------------------------------------//

    // Prepare A, B, C
    copy_cholmod_matrix_to_cs( matA, &A, matA->neq, matA->neq, matA->nnz, matA->nnz, true );
    Ac = cs_di_compress( &A ); // cs_droptol( Ac, 1.0e-15 );
    copy_cholmod_matrix_to_cs( matB, &B, matB->neq, matC->neq, matB->nnz, matB->nnz, true );
    Bc = cs_di_compress( &B );
    copy_cholmod_matrix_to_cs( matC, &C, matC->neq, matB->neq, matC->nnz, matC->nnz, true );
    Cc = cs_di_compress( &C );
    // Prepare D
    copy_cholmod_matrix_to_cs( matD, &D, matD->neq, matD->neq, Stokes->neq_cont, Stokes->neq_cont+1, false );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    Dc = cs_di_compress( &D );

    //------------------------------------------------------------------------------------------------//

    // Prepare AJ, BJ,  CJ
    copy_cholmod_matrix_to_cs( JmatA, &AJ, JmatA->neq, JmatA->neq, JmatA->nnz, JmatA->nnz, true );
    AJc = cs_di_compress( &AJ );
    copy_cholmod_matrix_to_cs( JmatB, &BJ, JmatB->neq, JmatC->neq, JmatB->nnz, JmatB->nnz, true );
    BJc  = cs_di_compress( &BJ );
    copy_cholmod_matrix_to_cs( JmatC, &CJ, JmatC->neq, JmatB->neq, JmatC->nnz, JmatC->nnz, true );
    CJc  = cs_di_compress( &CJ );

    //------------------------------------------------------------------------------------------------//

    // Construct preconditionner:
    if ( pc_type == 0 ) {
        printf("Construct preconditionner:   PC  = K\n");
        PC = cs_di_add( Ac, Ac, 1.0, 0.0 );
    }

    if ( pc_type == 1 ) {
        printf("Construct preconditionner:   PC  = 1/2 * (J'+ J)\n");
        Jt = cs_di_transpose( AJc, 1);
        PC = cs_di_add( AJc, Jt, 0.5, 0.5);
        cs_spfree(Jt);
    }

    //------------------------------------------------------------------------------------------------//

    // --------------- Schur complements --------------- //
    printf("Compute Schur complement 1:  Jt  = J  - grad*(PPI*div)\n");
    printf("Compute Schur complement 2:  Jts = PC - grad*(PPI*div)\n");

    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );

    //    //----- test - in case C' != B (assume B is deficient)
    //    B1 = cs_di_transpose( Cc, 1);
    //
    //    // minus sign: B = -C'
    //    for (k=0; k<B1->nzmax; k++)  B1->x[k] *= -1.0; // could be implicitly included in the next lines

    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply(  Bc, L);
    L2 = cs_di_multiply( BJc, L);
    cs_spfree(L);

    // Matrix addition: Js = AJ - B*(D*C)
    Jts = cs_di_add(  PC, L1, 1, -1);
    Js  = cs_di_add( AJc, L2, 1, -1);
    cs_spfree(L1);
    cs_spfree(L2);

    //------------------------------------------------------------------------------------------------//
    // Factor Jts
    printf("Cholesky factors of Jts...\n");

    int N = AJc->m;

    // Prepare Jacobian preconditioner blocks
    Lcm = cholmod_allocate_sparse (Jts->m, Jts->n, Jts->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, Jts );
    Kcm = cholmod_allocate_sparse (Js->m, Js->n, Js->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Kcm, Js );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    Bcm = cholmod_allocate_sparse (Bc->m, Bc->n, Bc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, Bc );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );

    // Prepare Jacobian blocks
    AcmJ = cholmod_allocate_sparse (AJc->m, AJc->n, AJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( AcmJ, AJc );
    BcmJ = cholmod_allocate_sparse (BJc->m, BJc->n, BJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( BcmJ, BJc );
    CcmJ = cholmod_allocate_sparse (CJc->m, CJc->n, CJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( CcmJ, CJc );

    // Keep lower part
    Lcml  = cholmod_copy ( Lcm,  -1, 1, &c );

    if (Lcml == NULL || Lcml->stype == 0) {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }

    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    }

    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

    if (c.status==1) { // Here the code should not continue running and start producings Nans!
        printf("CHOLDMOD failed because the matrix is not positive definite...\n");
        printf("MDoodz7.0 has to stop...\n");
        printf("Viscosity or diagonal elements of the tangent matrix are probably negative!\n");
        exit(1);
    }

    //------------------------------------------------------------------------------------------------//

    // Powell-Hestenes iterations
    t_omp = (double)omp_get_wtime();
    printf("Powell-Hestenes iterations, noisy = %d...\n", noisy);

    cholmod_dense  *du, *dp, *bu, *bp, *pdum, *udum, *fu, *fp;
    double minru0, maxru0, minru, maxru, ru, rp;
    int its_KSP, its_KSP_tot=0;
    udum  = cholmod_zeros( matA->neq, 1, CHOLMOD_REAL, &c );
    pdum  = cholmod_zeros( matC->neq, 1, CHOLMOD_REAL, &c );
    bu    = cholmod_zeros( matA->neq, 1, CHOLMOD_REAL, &c );
    bp    = cholmod_zeros( matC->neq, 1, CHOLMOD_REAL, &c );
    du    = cholmod_zeros( matA->neq, 1, CHOLMOD_REAL, &c );
    dp    = cholmod_zeros( matC->neq, 1, CHOLMOD_REAL, &c );
    fu    = cholmod_zeros( matA->neq, 1, CHOLMOD_REAL, &c );
    fp    = cholmod_zeros( matC->neq, 1, CHOLMOD_REAL, &c );

    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    // Residual U 
    copy_cholmod_dense_to_cholmod_dense( fu, bu );          // fu = bu
    cholmod_sdmult(  AcmJ, 0, mone, one, du, fu, &c );      // fu -= A*u
    cholmod_sdmult(  BcmJ, 0, mone, one, dp, fu, &c );      // fu -= B*p
    // Residual P 
    copy_cholmod_dense_to_cholmod_dense( fp, bp );          // fp = bp
    cholmod_sdmult(  CcmJ, 0, mone, one, du, fp, &c );      // fp -= C*u
    cholmod_sdmult( D1cm0, 0, mone, one, dp, fp, &c );      // fp -= D*p
    
    // Check: Initial residuals
    printf("Initial residual:\n");
    MinMaxArrayVal( fu->x, matA->neq, &minru0, &maxru0 );
    NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
    
    for ( k=0; k<nitmax; k++ ) {

        // Schur complement
        cholmod_sdmult ( Dcm, 0, one, zero, bp, pdum, &c);                  // pdum <-- D * fp
        copy_cholmod_dense_to_cholmod_dense( udum, bu );                    // udum <-- bu
        cholmod_sdmult ( BcmJ, 0, mone, one, pdum, udum, &c);               // udum <-- bu - B*(D*fp)
        cholmod_sdmult ( whos_incompressible, 0, one, zero, dp, pdum, &c);  // pdum <-- whos_compressible * dp
        cholmod_sdmult ( BcmJ, 0, mone, one,   pdum, udum, &c);             // udum <-- bu - B*(D*fp) - B*dp          !!!!!!! needed?

        // Solve for velocity correction
        kspgcr( Kcm, udum, du, Lfact, matA->neq, &c, model.rel_tol_KSP, noisy, &its_KSP);
        its_KSP_tot += its_KSP;

        // Solve for pressure correction
        copy_cholmod_dense_to_cholmod_dense( pdum, bp );        // pdum <-- bp
        cholmod_sdmult( CcmJ,  0, mone, one,   du, pdum, &c );  // pdum <-- bp - C*u
        cholmod_sdmult( D1cm0, 0, mone, one,   dp, pdum, &c );  // fp -= D*p
        cholmod_sdmult( Dcm ,  0,  one, one, pdum,   dp, &c );  // dp <-- dp + D*(bp - C*u)

        // Residual U 
        copy_cholmod_dense_to_cholmod_dense( fu, bu );          // fu = bu
        cholmod_sdmult(  AcmJ, 0, mone, one, du, fu, &c );      // fu -= A*u
        cholmod_sdmult(  BcmJ, 0, mone, one, dp, fu, &c );      // fu -= B*p
        // Residual P 
        copy_cholmod_dense_to_cholmod_dense( fp, bp );          // fp = bp
        cholmod_sdmult(  CcmJ, 0, mone, one, du, fp, &c );      // fp -= C*u
        cholmod_sdmult( D1cm0, 0, mone, one, dp, fp, &c );      // fp -= D*p

        // Check: Stopping criteria
        MinMaxArrayVal( fp->x, matC->neq, &mindiv, &maxdiv );
        MinMaxArrayVal( fu->x, matA->neq, &minru , &maxru  );
        if (k==0) { maxdiv0 = maxdiv; maxru0 = maxru; }

        if ( noisy > 0 ) printf("PH comp it. %01d. its_KSP = %02d: max. cont. = %2.2e - rel. max. div. = %2.2e / max. mom. = %2.2e - rel. max. mom. = %2.2e\n", k, its_KSP, maxdiv, fabs(maxdiv/maxdiv0), maxru, maxru/maxru0);
        if ( k>1 && (fabs(maxdiv)<model.abs_tol_div || maxdiv/maxdiv0<model.rel_tol_div )  && (fabs(maxru)<model.abs_tol_mom || maxru/maxru0<model.rel_tol_mom ) ) break;
    }

    // A posteriori checks
    if (fabs(maxdiv)>model.abs_tol_div && maxdiv/maxdiv0>model.rel_tol_div) {
        printf("The code has exited since the incompressibility constrain was not satisfied to abs. tol. = %2.2e and rel. tol. = %2.2e\n Try modifying the PENALTY factor or check MIN/MAX viscosities\n Good luck!\n", model.abs_tol_div, model.rel_tol_div);
        exit(1);
    }
    printf("** PH - iterations = %lf sec - its_KSP_tot = %02d\n", (double)((double)omp_get_wtime() - t_omp), its_KSP_tot);

    // --------------- Solution vector --------------- //
    BackToSolutionVector( du, dp, sol, mesh, Stokes );

    // separate residuals
    double Area =0.0, resx=0.0,resz=0.0, resp=0.0;
    int ndofx=0, ndofz=0, ndofp=0;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;

    BackToSolutionVector( fu, fp, F, mesh, Stokes );

    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]];//*celvol;
        }
    }

#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]];//*celvol;
            //            if ( mesh->BCv.type[cc] == 2) printf("F=%2.2e\n", F[Stokes->eqn_v[cc]]);
        }
    }

#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]];//*celvol;
        }
    }

    // Sqrt
    resx =  sqrt(resx/ndofx);
    resz =  sqrt(resz/ndofz);
    resp =  sqrt(resp/ndofp);

    if ( noisy > 0 ) {
        printf("Fu = %2.6e\n", resx ); // Units of momentum
        printf("Fv = %2.6e\n", resz ); // Units of momentum
        printf("Fp = %2.6e\n", resp ); // Units of velocity gradient
    }

    // Freedom
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &du, &c );
    cholmod_free_dense( &dp, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );

    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Kcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm, &c );
    cholmod_free_sparse( &Dcm0, &c );
    cholmod_free_sparse( &D1cm0, &c );
    cholmod_free_sparse( &whos_incompressible, &c );

    cholmod_free_sparse( &AcmJ, &c );
    cholmod_free_sparse( &BcmJ, &c );
    cholmod_free_sparse( &CcmJ, &c );

    pardi->Lfact = Lfact;
    pardi->c     = c;

    DoodzFree(u0);
    DoodzFree(p0);
    DoodzFree(F);
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);

    DoodzFree(AJ.p);
    DoodzFree(AJ.i);
    DoodzFree(AJ.x);
    cs_spfree(AJc);
    DoodzFree(BJ.p);
    DoodzFree(BJ.i);
    DoodzFree(BJ.x);
    cs_spfree(BJc);
    DoodzFree(CJ.p);
    DoodzFree(CJ.i);
    DoodzFree(CJ.x);
    cs_spfree(CJc);

    cs_spfree(PC);
    cs_spfree(Js);
    cs_spfree(Jts);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void KSPStokesDecoupled( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes, SparseMat *Jacobian, SparseMat *JmatA,  SparseMat *JmatB,  SparseMat *JmatC ) {

    //double Control [AMD_CONTROL], Info [AMD_INFO] ;
    cs_di  A, B, D, C, *B1, *L, *Ac, *Bc,  *Cc,  *Dc, *L1;
    cs_di  AJ, BJ, CJ, *AJc, *BJc, *CJc;
    //int    *P, msglvl = 0;
    DoodzFP  *u0, *p0, *F;
    int  noisy=1;
    int k, cc, i; //nitmax=5
    double celvol = model.dx*model.dz;
    
    cholmod_common c ;
    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm; //, *A1
    cholmod_sparse *AcmJ, *BcmJ, *CcmJ;
    cholmod_sparse *Dcm0, *D1cm0; //, *BDC, *DC, *Lcm0, *Lcml0, *Acm0, *Bcm0, *Ccm0
    cholmod_factor *Lfact ;
    cholmod_sparse *M, *AB, *CD, *Iu, *Ip, *D_zero; // *bs, *bus, *bps, *M0
    cholmod_dense *b, *x, *f, *val, *s, *v;
    //    cholmod_start( &c ) ;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    
    Lfact = pardi->Lfact;
    c     = pardi->c;
    
    // --------------- Pre-process--------------- //
    
    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma  = model.penalty, penalty;//1e12;//1e10*model.Nx*model.Nz;
    int vol_change = model.VolChangeReac;

    rsize = Stokes->neq_cont;
    Dcm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    D1cm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    
    penalty = gamma / celvol;
    printf("Penalty factor = %2.2e\n", penalty);
    //    for (k=0;k<Dcm0->nzmax;k++) ((double*)Dcm0->x)[k] *= gamma*celvol;
    //    printf("-gamma*celvol = %2.2e %2.2e %2.2e %d %d\n", -gamma*celvol, model.dx*scaling.L, model.dz*scaling.L, model.Nx, model.Nz);
    
    
#pragma omp parallel for shared(D1cm0, Dcm0, mesh, Stokes, matA, matD ) private( i ) firstprivate( model, celvol, vol_change )
    for( k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
        if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31 ) {
            i = Stokes->eqn_p[k] - matA->neq;
            
            //            cPC    = drhodPc ./ (rhoc.*dt); %% new diagonal coefficient (= 1/(K*dt) for standard EOS)
            
            
            // Here Dcm0 is the pressure block
            if (mesh->comp_cells[k]==0) ((double*)D1cm0->x)[i] *= 0.0;
            if (mesh->comp_cells[k]==1 && vol_change == 0 ) ((double*)D1cm0->x)[i]  = mesh->bet_n[k] / model.dt * celvol * matD->d[i]*matD->d[i];
            if (mesh->comp_cells[k]==1 && vol_change == 1 ) ((double*)D1cm0->x)[i]  = mesh->drhodp_n[k] / (mesh->rho_n[k]*model.dt) * celvol * matD->d[i]*matD->d[i];

            // Here Dcm0 is the inverse of the pressure block
            if (mesh->comp_cells[k]==0) ((double*)Dcm0->x)[i] *= penalty; // Should be /celvol
            if (mesh->comp_cells[k]==1) ((double*)Dcm0->x)[i]  = 1.0 /  ((double*)D1cm0->x)[k]; // Should be /celvol
            //          printf("%2.2e %2.2e %2.2e %2.2e\n", mesh->bet[k]/model.dt, penalty, mesh->bet[k]*(1/scaling.S), model.dt*scaling.t);
        }
    }
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();
    
    // Build initial solution vector
    F  = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    BuildInitialSolutions( u0, p0, mesh );
    
    //------------------------------------------------------------------------------------------------//
    
    // Prepare A
    A.nzmax = matA->nnz;
    A.nz    = matA->nnz;
    A.m     = matA->neq;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, matA->Ic, A.i );
    ArrayEqualArrayI( A.p, matA->J,  A.nzmax );
    ArrayEqualArray(  A.x, matA->A,  A.nzmax );
    Ac  = cs_di_compress( &A );
    //    cs_droptol( Ac, 1.0e-15 );
    
    // Prepare B
    B.nzmax = matB->nnz;
    B.nz    = matB->nnz;
    B.m     = matB->neq;
    B.n     = matC->neq;
    B.p     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.i     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.x     = DoodzCalloc( B.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( B.m, matB->Ic, B.i );
    ArrayEqualArrayI( B.p, matB->J,  B.nzmax );
    ArrayEqualArray(  B.x, matB->A,  B.nzmax );
    Bc  = cs_di_compress( &B );
    
    // Prepare C
    C.nzmax = matC->nnz;
    C.nz    = matC->nnz;
    C.m     = matC->neq;
    C.n     = matB->neq;
    C.p     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.i     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.x     = DoodzCalloc( C.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( C.m, matC->Ic, C.i );
    ArrayEqualArrayI( C.p, matC->J,  C.nzmax );
    ArrayEqualArray(  C.x, matC->A,  C.nzmax );
    Cc  = cs_di_compress( &C );
    
    // Prepare D
    D.nzmax = Stokes->neq_cont;
    D.nz    = Stokes->neq_cont;
    D.m     = matD->neq;
    D.n     = matD->neq;
    D.p     = DoodzCalloc( D.nzmax+1, sizeof(int) );
    D.i     = DoodzCalloc( D.nzmax, sizeof(int) );
    D.x     = DoodzCalloc( D.nzmax, sizeof(double) );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    Dc  = cs_di_compress( &D );
    
    
    //------------------------------------------------------------------------------------------------//
    
    // Prepare AJ
    AJ.nzmax = JmatA->nnz;
    AJ.nz    = JmatA->nnz;
    AJ.m     = JmatA->neq;
    AJ.n     = AJ.m;
    //    printf( "A.nzmax = %d\n", matA->nnz );
    AJ.p     = DoodzCalloc( AJ.nzmax, sizeof(int) );
    AJ.i     = DoodzCalloc( AJ.nzmax, sizeof(int) );
    AJ.x     = DoodzCalloc( AJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( AJ.m, JmatA->Ic, AJ.i );
    ArrayEqualArrayI( AJ.p, JmatA->J,  AJ.nzmax );
    ArrayEqualArray(  AJ.x, JmatA->A,  AJ.nzmax );
    AJc  = cs_di_compress( &AJ );
    //    cs_droptol( AJc, 1.0e-15 );
    
    // Prepare BJ
    BJ.nzmax = JmatB->nnz;
    BJ.nz    = JmatB->nnz;
    BJ.m     = JmatB->neq;
    BJ.n     = JmatC->neq;
    BJ.p     = DoodzCalloc( BJ.nzmax, sizeof(int) );
    BJ.i     = DoodzCalloc( BJ.nzmax, sizeof(int) );
    BJ.x     = DoodzCalloc( BJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( BJ.m, JmatB->Ic, BJ.i );
    ArrayEqualArrayI( BJ.p, JmatB->J,  BJ.nzmax );
    ArrayEqualArray(  BJ.x, JmatB->A,  BJ.nzmax );
    BJc  = cs_di_compress( &BJ );
    
    // Prepare CJ
    CJ.nzmax = JmatC->nnz;
    CJ.nz    = JmatC->nnz;
    CJ.m     = JmatC->neq;
    CJ.n     = JmatB->neq;
    CJ.p     = DoodzCalloc( CJ.nzmax, sizeof(int) );
    CJ.i     = DoodzCalloc( CJ.nzmax, sizeof(int) );
    CJ.x     = DoodzCalloc( CJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( CJ.m, JmatC->Ic, CJ.i );
    ArrayEqualArrayI( CJ.p, JmatC->J,  CJ.nzmax );
    ArrayEqualArray(  CJ.x, JmatC->A,  CJ.nzmax );
    CJc  = cs_di_compress( &CJ );
    
    //------------------------------------------------------------------------------------------------//
    
    // --------------- Schur complement --------------- //
    
    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );
    
    //----- test - in case C' != B (assume B is deficient)
    B1 = cs_di_transpose( Cc, 1);
    
    // minus sign: B = -C'
    for (k=0;k<B1->nzmax;k++) {
        B1->x[k] *= -1.0;
    }
    //-----
    
    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply( B1, L);
    cs_spfree(L);
    
    // Matrix addition: L = A - B*(D*C)
    L = cs_di_add( Ac, L1, 1, -1);
    cs_spfree(L1);
    
    // --------------- Cholesky --------------- //
    
    int N = Ac->m + Cc->m, cc2;
    
    // Prepare Jacobian preconditioner blocks
    Lcm = cholmod_allocate_sparse (L->m, L->n, L->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, L );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    Bcm = cholmod_allocate_sparse (B1->m, B1->n, B1->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, B1 );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );
    
    // Prepare Jacobian blocks
    AcmJ = cholmod_allocate_sparse (AJc->m, AJc->n, AJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( AcmJ, AJc );
    BcmJ = cholmod_allocate_sparse (BJc->m, BJc->n, BJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( BcmJ, BJc );
    CcmJ = cholmod_allocate_sparse (CJc->m, CJc->n, CJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( CcmJ, CJc );
    
    // Keep lower part
    Lcml  = cholmod_copy ( Lcm,  -1, 1, &c );
    
    if (Lcml == NULL || Lcml->stype == 0)
    {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }
    
    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    }
    
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    cholmod_dense *u, *p, *bu, *bp, *pdum, *udum,  *fu, *fp, *Y, *E;
    pdum = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    udum = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    u    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    E    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    Y    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    p    = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    
    // separate residuals
    double Area =0., resx=0.,resz=0., resp=0.;
    int ndofx=0, ndofz=0, ndofp=0;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    
    //----------------- DEVELOPMENT -----------------//
    
    // Initialise KSP
    int restart = 20;//6;
    int max_it  = 1000;
    int ncycles = 0;
    int its     = 0, i1, i2, success=0;
    double eps  = model.rel_tol_KSP, norm_r, rnorm0, fact, r_dot_v, nrm;
    double **VV, **SS;
    b      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    x      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    f      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    D_zero = cholmod_spzeros ( Dcm->nrow, Dcm->ncol, 0, CHOLMOD_REAL, &c ) ;
    fu     = cholmod_zeros (   Ac->m,     1, CHOLMOD_REAL, &c );
    fp     = cholmod_zeros (   Dc->m,     1, CHOLMOD_REAL, &c );
    val    = cholmod_zeros ( restart,     1, CHOLMOD_REAL, &c );
    s      = cholmod_zeros (       N,     1, CHOLMOD_REAL, &c );
    v      = cholmod_zeros (       N,     1, CHOLMOD_REAL, &c );
    Iu     = cholmod_speye (   Ac->m, Ac->n, CHOLMOD_REAL, &c );
    Ip     = cholmod_speye (   Dc->m, Dc->n, CHOLMOD_REAL, &c );
    
    // Allocate temp vecs
    VV  = (double**) DoodzCalloc( restart, sizeof(double*));
    SS  = (double**) DoodzCalloc( restart, sizeof(double*));
    for (cc=0; cc<restart; cc++) {
        VV[cc]  = (double*) DoodzCalloc( N, sizeof(double));
        SS[cc]  = (double*) DoodzCalloc( N, sizeof(double));
    }
    
    // Initial fields
    copy_vec_to_cholmod_dense( u, u0 );
    copy_vec_to_cholmod_dense( p, p0 );
    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    // Concatenate b = [bu; bp]
    for (cc=0; cc<bu->nrow; cc++) {
        ((double*)b->x)[cc ] = ((double*)bu->x)[cc];
    }
    for (cc=0; cc<bp->nrow; cc++) {
        cc2=bu->nrow+cc;
        ((double*)b->x)[cc2] = ((double*)bp->x)[cc];
    }
    
    // 1. Concatenate M = [A,B;C,D]
    AB  = cholmod_horzcat ( AcmJ, BcmJ,   1, &c );    // AB = [A B]
    CD  = cholmod_horzcat ( CcmJ, D1cm0, 1, &c );     // CD = [C D]
    
    // 2. Concatenate M = [A,B;C,D]
    M   = cholmod_vertcat ( AB,  CD,     1, &c );     // M  = [AB; CD]
    cholmod_free_sparse( &AB,  &c );
    cholmod_free_sparse( &CD,  &c );
    
    // Concatenate s = [u; p]
    for (cc=0; cc<u->nrow; cc++) {
        ((double*)x->x)[cc] =  ((double*)u->x)[cc];
    }
    for (cc=0; cc<p->nrow; cc++) {
        cc2=u->nrow+cc;
        ((double*)x->x)[cc2] = ((double*)p->x)[cc];
    }
    
    // Compute residual f = b - Mx
    copy_cholmod_dense_to_cholmod_dense( f, b );
    cholmod_sdmult ( M, 0, mone, one, x, f, &c) ;
    MinMaxArray(M->x, 1, M->nzmax, "M");
    MinMaxArray(x->x, 1, x->nzmax, "x");
    MinMaxArray(b->x, 1, b->nzmax, "b");
    norm_r = cholmod_norm_dense ( f, 2, &c );
    
    // Evaluate reference norm
    rnorm0 = norm_r;
    printf("       %1.4d KSP GCR Residual %1.12e %1.12e\n", 0, norm_r, norm_r/rnorm0);
    
    while ( success == 0 && its<max_it ) {
        
        for ( i1=0; i1<restart; i1++ ) {
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            
            // extract ru from r = [ru; rp]
#pragma omp parallel for shared( fu, f ) private( cc )
            for (cc=0; cc<fu->nrow; cc++) {
                ((double*)fu->x)[cc] = ((double*)f->x)[cc];
            }
#pragma omp parallel for shared( fp, f ) private( cc, cc2 ) firstprivate( fu )
            for (cc=0; cc<fp->nrow; cc++) {
                cc2=fu->nrow+cc;
                ((double*)fp->x)[cc] = ((double*)f->x)[cc2];
            }
            
            cholmod_sdmult ( Dcm, 0,  one, zero, fp, pdum, &c );
            cholmod_sdmult ( Bcm, 0, mone,  one, pdum, fu, &c );
            cholmod_solve2(CHOLMOD_A, Lfact, fu, NULL, &u, NULL, &Y, &E, &c);
            cholmod_sdmult ( Dcm, 0,  one, zero, fp,    p, &c );
            cholmod_sdmult ( Ccm, 0,  one, zero,     u, pdum, &c );
            cholmod_sdmult ( Dcm, 0, mone,  one,  pdum,    p, &c );
            
            // Concatenate s = [u; p]
#pragma omp parallel for shared( s, u ) private( cc )
            for (cc=0; cc<u->nrow; cc++) {
                ((double*)s->x)[cc] = ((double*)u->x)[cc];
            }
#pragma omp parallel for shared( s, p ) private( cc, cc2 ) firstprivate( u )
            for (cc=0; cc<p->nrow; cc++) {
                cc2=u->nrow+cc;
                ((double*)s->x)[cc2] =((double*)p->x)[cc];
            }
            
            // Simplest preconditionning: PC = I
            //            ArrayEqualArray( s->x, f->x, N );
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            
            // Action of Jacobian on s
            cholmod_sdmult ( M, 0, one, zero, s, v, &c) ;
            
            // Approximation of the Jv product
            for (i2=0; i2<i1; i2++) {
                ((double*)val->x)[i2] = DotProduct( ((double*)v->x), VV[i2], N );
            }
            
            for(i2=0; i2<i1+1; i2++) {
                fact = -((double*)val->x)[i2];
                ArrayPlusScalarArray( ((double*)v->x), fact, VV[i2], N );
                ArrayPlusScalarArray( ((double*)s->x), fact, SS[i2], N );
            }
            
            r_dot_v = DotProduct( ((double*)f->x), ((double*)v->x), N );
            nrm     = sqrt(  DotProduct( ((double*)v->x), ((double*)v->x), N )  );
            r_dot_v = r_dot_v / nrm;
            
            fact = 1.0/nrm;
            ArrayTimesScalar( ((double*)v->x), fact, N );
            ArrayTimesScalar( ((double*)s->x), fact, N );
            
            fact = r_dot_v;
            ArrayPlusScalarArray( ((double*)x->x), fact, ((double*)s->x), N );
            fact =-r_dot_v;
            ArrayPlusScalarArray( ((double*)f->x), fact, ((double*)v->x), N );
            
            // Check convergence
            norm_r = cholmod_norm_dense ( f, 2, &c );
            if (norm_r < eps * rnorm0 ) { // || norm_r < epsa
                success = 1;
                break;
            }
            //            printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
            
            
            // Store arrays
            ArrayEqualArray( VV[i1], ((double*)v->x), N );
            ArrayEqualArray( SS[i1], ((double*)s->x), N );
            its++;
        }

        its++;
        ncycles++;
    }
    printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
    MinMaxArray(fp->x, scaling.E, matC->neq, "divU" );
    
    // Extract ru from r = [ru; rp]
#pragma omp parallel for shared( f, x, u, fu ) private( cc )
    for (cc=0; cc<u->nrow; cc++) {
        ((double*)fu->x)[cc] = ((double*)f->x)[cc];
        ((double*) u->x)[cc] = ((double*)x->x)[cc];
    }
#pragma omp parallel for shared( f, x, p, fp ) private( cc, cc2 ) firstprivate( u )
    for (cc=0; cc<p->nrow; cc++) {
        cc2=u->nrow+cc;
        ((double*) p->x)[cc] = ((double*)x->x)[cc2];
        ((double*)fp->x)[cc] = ((double*)f->x)[cc2];
    }
    
    // --------------- Solution vector --------------- //
    BackToSolutionVector( u, p, sol, mesh, Stokes );
    
    BackToSolutionVector( fu, fp, F, mesh, Stokes );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]];//*celvol;
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]];//*celvol;
            //            if ( mesh->BCv.type[cc] == 2) printf("F=%2.2e\n", F[Stokes->eqn_v[cc]]);
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]];//*celvol;
        }
    }
    
    // Sqrt
    resx =  sqrt(resx)/Area/(ndofx);
    resz =  sqrt(resz)/Area/(ndofz);
    resp =  sqrt(resp)/Area/(ndofp);
    
    if ( noisy == 1 ) {
        printf("Fu = %2.4e\n", resx * (scaling.F/pow(scaling.L,3))); // Units of momentum
        printf("Fv = %2.4e\n", resz * (scaling.F/pow(scaling.L,3))); // Units of momentum
        printf("Fp = %2.4e\n", resp * (scaling.E)); // Units of velocity gradient
    }
    
    // Free temp vecs
    for (cc=0; cc<restart; cc++) {
        DoodzFree(VV[cc]);
        DoodzFree(SS[cc]);
    }
    DoodzFree(VV);
    DoodzFree(SS);
    
    // Free
    cholmod_free_sparse( &M,  &c );
    cholmod_free_dense ( &f,  &c );
    cholmod_free_dense ( &b,  &c );
    cholmod_free_dense ( &x,  &c );
    cholmod_free_dense( &val,  &c );
    cholmod_free_dense( &s,  &c );
    cholmod_free_dense( &v,  &c );
    cholmod_free_sparse( &Iu,  &c );
    cholmod_free_sparse( &Ip,  &c );
    cholmod_free_sparse( &D_zero,  &c );
    
    //----------------- DEVELOPMENT -----------------//
    
    // --------------- Free --------------- //
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &u, &c );
    cholmod_free_dense( &E, &c );
    cholmod_free_dense( &Y, &c );
    cholmod_free_dense( &p, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );
    cholmod_free_sparse( &Lcm, &c );

    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm, &c );
    cholmod_free_sparse( &Dcm0, &c );
    cholmod_free_sparse( &D1cm0, &c );
    
    cholmod_free_sparse( &AcmJ, &c );
    cholmod_free_sparse( &BcmJ, &c );
    cholmod_free_sparse( &CcmJ, &c );

    pardi->Lfact = Lfact;
    pardi->c     = c;
    
    cs_spfree(L);
    
    DoodzFree(u0);
    DoodzFree(p0);
    
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    cs_spfree(B1);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);
    
    DoodzFree(AJ.p);
    DoodzFree(AJ.i);
    DoodzFree(AJ.x);
    cs_spfree(AJc);
    DoodzFree(BJ.p);
    DoodzFree(BJ.i);
    DoodzFree(BJ.x);
    cs_spfree(BJc);
    DoodzFree(CJ.p);
    DoodzFree(CJ.i);
    DoodzFree(CJ.x);
    cs_spfree(CJc);    
    DoodzFree(F);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/