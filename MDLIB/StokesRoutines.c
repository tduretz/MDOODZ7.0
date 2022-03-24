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
#include "header_MDOODZ.h"
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

void ApplyBC( grid* mesh, params* model ) {

    int nx=model->Nx, nz=model->Nz, nzvx=nz+1, nxvz=nx+1;
    int i, j;

    // Vx Neumann
    for( i=0; i<nx; i++) {
        // South
        if ( mesh->BCu.type[i] == 13 ) {
            mesh->u_in[i] = mesh->u_in[i+nx];
        }
        // North
        if ( mesh->BCu.type[i + (nzvx-1)*nx] == 13 ) {
            mesh->u_in[i + (nzvx-1)*nx] = mesh->u_in[i + (nzvx-2)*nx];
        }
    }

    // Vx Dirichlet
    for( i=0; i<nx; i++) {
        // South
        if ( mesh->BCu.type[i] == 11 ) {
            mesh->u_in[i] = 2.0*mesh->BCu.val[i] - mesh->u_in[i+nx];
        }
        // North
        if ( mesh->BCu.type[i + (nzvx-1)*nx] == 11 ) {
            mesh->u_in[i + (nzvx-1)*nx] = 2.0*mesh->BCu.val[i + (nzvx-1)*nx] - mesh->u_in[i + (nzvx-2)*nx];
        }
    }

    // Vz Neumann
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == 13 ) {
            mesh->v_in[j*nxvz] = mesh->v_in[j*nxvz + 1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == 13 ) {
            mesh->v_in[j*nxvz+(nxvz-1)] = mesh->v_in[j*nxvz+(nxvz-1) -1];
        }
    }

    // Vz Dirichlet
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == 11 ) {
            mesh->v_in[j*nxvz] = 2.0*mesh->BCv.val[j*nxvz] - mesh->v_in[j*nxvz + 1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == 11 ) {
            mesh->v_in[j*nxvz + (nxvz-1)] = 2.0*mesh->BCv.val[j*nxvz+(nxvz-1)] - mesh->v_in[j*nxvz + (nxvz-1) -1];
        }
    }

    // Vz Periodic
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == -12 ) {
            mesh->v_in[j*nxvz] = mesh->v_in[j*nxvz + (nxvz-1) -1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == -12 ) {
            mesh->v_in[j*nxvz+(nxvz-1)] = mesh->v_in[j*nxvz+1];
        }
    }

    // Vx Periodic
    for( j=0; j<nzvx; j++) {
        if ( mesh->BCu.type[j*nx + nx-1] == -12 ) {
            mesh->u_in[j*nx + nx-1] = mesh->u_in[j*nx];
        }
    }

    //    printf("Vx\n");
    //    Print2DArrayDouble( mesh->u_in, nx, nzvx, 1.0 );
    //    printf("Vz\n");
    //    Print2DArrayDouble( mesh->v_in, nxvz, nz, 1.0 );
    //    printf("p\n");
    //    Print2DArrayDouble( mesh->p_in, nx-1, nz-1, 1.0 );
}



/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateNonLinearity( grid* mesh, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, params *model, Nparams *Nmodel, scale scaling, int mode, double h_contin ) {

    // Strain rate component evaluation
    StrainRateComponents( mesh, scaling, model );
    //    MinMaxArrayTag( mesh->exxd,      scaling.E, (mesh->Nx-1)*(mesh->Nz-1),     "exx     ", mesh->BCp.type );

    //-----------------------------------------------//

    NonNewtonianViscosityGrid ( mesh, &materials, model, *Nmodel, &scaling );

    //-----------------------------------------------//

    // Evaluate right hand side
    EvaluateRHS( mesh, *model, scaling, materials.rho[0] );

    //-----------------------------------------------//

    // Fill up the rheological matrices arrays
    RheologicalOperators( mesh, model, &scaling, 0 );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DetectCompressibleCells ( grid* mesh, params *model ) {

    int cc, nx=model->Nx, nz=model->Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1, kk=0;

     printf("---> Detecting compressibles cells\n");
     for( cc=0; cc<ncx*ncz; cc++) {

         if ( mesh->BCp.type[cc] != 30 ) {

             if ( mesh->bet_n[cc] > 1e-13 ) {
                 mesh->comp_cells[cc] = 1;
                 kk++;
             }

         }
     }
    printf("---> %04d compressibles cells detected\n", kk);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void ExtractSolutions2( SparseMat *Stokes, grid* mesh, params* model, double* dx, double alpha ) {

    int cc, nx=model->Nx, nz=model->Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1, kk;

    double eps = 1e-13;

// Test relaxed u-v-p solutions
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nzvx, nx )
for( cc=0; cc<nzvx*nx; cc++) {
    if ( mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
        mesh->u_in[cc] = mesh->u_in[cc] + alpha*dx[Stokes->eqn_u[cc]];
    }
}

    // Periodic
    for( cc=0; cc<nzvx; cc++) {
        kk = cc*nx + nx-1;
        if ( mesh->BCu.type[kk] ==-12) {
            mesh->u_in[kk] = mesh->u_in[kk-mesh->Nx+1];
        }
    }

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nz, nxvz )
for( cc=0; cc<nz*nxvz; cc++) {
    if ( mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
        mesh->v_in[cc] = mesh->v_in[cc] + alpha*dx[Stokes->eqn_v[cc]];
    }
}

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, ncx, ncz )
for( cc=0; cc<ncz*ncx; cc++) {

//    mesh->p_in[cc] = 0.0;
//    mesh->dp[cc]   = 0.0;
    if ( mesh->BCp.type[cc] != 30 && mesh->BCp.type[cc] != 0  && mesh->BCp.type[cc] != 31 ) {
        mesh->p_in[cc] =  mesh->p_in[cc] + alpha*dx[Stokes->eqn_p[cc]];
//        mesh->p_in[cc]  = mesh->p_trial[cc] + alpha*dx[Stokes->eqn_p[cc]];
//        mesh->dp[cc]    = alpha*dx[Stokes->eqn_p[cc]];
    }
}

// Apply Bc to Vx and Vz
ApplyBC( mesh, model );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double LineSearch( SparseMat *Stokes, double *dx, grid* mesh, params *model, Nparams* Nmodel, markers* particles, markers *topo_chain, surface *topo, mat_prop materials, scale scaling ) {

    printf("---- Line search for Stokes equations ----\n");

    int ix, kk, k, cc;
    int nx= model->Nx, nz=model->Nz, ncx=nx-1, ncz=nz-1, nzvx=nz+1, nxvz=nx+1;
    double *rx, *rz, *rp, alpha, *alphav, dalpha, minx;
    double *u, *v, *p;
    Nparams residuals;
    clock_t  t_omp;
    int success = 0, niter=0, nitermax=1;
    double minalpha[nitermax], maxalpha[nitermax], frac = 1.0;
    int ntry[nitermax];

    Nmodel->stagnated = 0;

    ntry[0]     = 6;
    minalpha[0] = -2.5;
    maxalpha[0] = -0.0;

    //    ntry[0]     = 11;
    //    minalpha[0] = -2.25;
    //    maxalpha[0] = -0.25;

    //    ntry[0]     = 1;
    //    minalpha[0] = -1.0;
    //    maxalpha[0] = -1.0;

    t_omp = (double)omp_get_wtime();

    // Allocate
    u          = DoodzMalloc( nx*nzvx*sizeof(double) );
    v          = DoodzMalloc( nxvz*nz*sizeof(double) );
    p          = DoodzMalloc( ncx*ncz*sizeof(double) );

    // Save initial fields (before updates)
    ArrayEqualArray( u, mesh->u_in, nx*nzvx );
    ArrayEqualArray( v, mesh->v_in, nxvz*nz );
    ArrayEqualArray( p, mesh->p_in, ncx*ncz );

    // Do line search iterations
    while ( success == 0 && niter < nitermax ) {

        // allocate array
        alphav = DoodzMalloc( ntry[niter]*sizeof(double) );
        rx     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rz     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rp     = DoodzMalloc( ntry[niter]*sizeof(double) );

        alpha    = maxalpha[niter];
        dalpha   = fabs(maxalpha[niter]-minalpha[niter])/(ntry[niter]-1);

        // Search for optimal relaxation parameters
        for( kk=0; kk<ntry[niter]; kk++ ) {

            // Update alpha
            if (kk>0) alpha -= dalpha;
            alphav[kk] = alpha;

            // Start from initial solutions
            ArrayEqualArray( mesh->u_in, u, nx*nzvx );
            ArrayEqualArray( mesh->v_in, v, nxvz*nz );
            ArrayEqualArray( mesh->p_in, p, ncx*ncz );

            // Test relaxed u-v-p solutions
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nzvx, nx )
            for( cc=0; cc<nzvx*nx; cc++) {
                if ( mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12  ) {
                    mesh->u_in[cc] += alpha*dx[Stokes->eqn_u[cc]];
                }
            }

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nz, nxvz )
            for( cc=0; cc<nz*nxvz; cc++) {
                if ( mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
                    mesh->v_in[cc] += alpha*dx[Stokes->eqn_v[cc]];
                }
            }

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, ncx, ncz )
            for( cc=0; cc<ncz*ncx; cc++) {

                if ( mesh->BCp.type[cc] != 30 && mesh->BCp.type[cc] != 0  && mesh->BCp.type[cc] != 31 ) {
                    mesh->p_in[cc] += alpha*dx[Stokes->eqn_p[cc]];
                }
            }

            // Apply Bc to Vx and Vz
            ApplyBC( mesh, model );

            //            Print2DArrayDouble( u,  nx, nzvx, scaling.V );
            //            Print2DArrayDouble( v,  nxvz, nz, scaling.V );
            //            Print2DArrayDouble( p,  ncx, ncz, scaling.S );

            //------------------------------------------------------------------------------------------------------
            // Update non-linearity
            UpdateNonLinearity( mesh, particles, topo_chain, topo, materials, model, Nmodel, scaling, 0, 1.0 );

            //------------------------------------------------------------------------------------------------------

            // Calculate residual
//            model->free_surf_stab = 0;
            EvaluateStokesResidual( Stokes, &residuals, mesh, *model, scaling, 1 );
//            model->free_surf_stab = dummy;

            rx[kk] = residuals.resx;
            rz[kk] = residuals.resz;
            rp[kk] = residuals.resp;
            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.4e rz = %2.4e rp = %2.2e\n", alpha, rx[kk]* (scaling.F/pow(scaling.L,3)), rz[kk]* (scaling.F/pow(scaling.L,3)), rp[kk]*scaling.E );

        }

        // Look for the minimum predicted residuals
        minx  = rz[0];
        ix = 0;
        for( k=1; k<ntry[niter]; k++ ) {
            if(rz[k]<minx) {
                minx = rz[k];
                ix = k;
            }
        }
        alpha = alphav[ix];

        // if the minmimun residuals are lower than starting ones, then success
        if ( rx[ix] < frac*Nmodel->resx_f || rz[ix]<frac*Nmodel->resz_f  ) { //|| rp[ix]<frac*Nmodel->resp
            success = 1;
            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ix], rx[ix]* (scaling.F/pow(scaling.L,3)), rz[ix]* (scaling.F/pow(scaling.L,3)), rp[ix]*scaling.E);
        }

        DoodzFree(rx);
        DoodzFree(rz);
        DoodzFree(rp);
        DoodzFree(alphav);
        niter++;
    }

    if ( fabs(alpha)<1e-13 || success == 0 ) {
        printf( "Found minimum of the function -- iteration cycle stagnates\n...\n" );
        Nmodel->stagnated = 1;
        alpha = 1.0;
    }

    // Reset solutions
    ArrayEqualArray( mesh->u_in, u, nx*nzvx );
    ArrayEqualArray( mesh->v_in, v, nxvz*nz );
    ArrayEqualArray( mesh->p_in, p, ncx*ncz );

    // Free
    DoodzFree( u );
    DoodzFree( v );
    DoodzFree( p );
    printf("** Line search took = %f sec\n", (double)((double)omp_get_wtime() - t_omp) );

    return alpha;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double LineSearchDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, double *dx, grid* mesh, params *model, Nparams* Nmodel, markers* particles, markers *topo_chain, surface *topo, mat_prop materials, scale scaling ) {

    printf("---- Line search for decoupled Stokes equations ----\n");

    int kk, k, cc;
    int nx= model->Nx, nz=model->Nz, ncx=nx-1, ncz=nz-1, nzvx=nz+1, nxvz=nx+1;
    double *rx, *rz, *rp, alpha, *alphav, dalpha, minx;
    double *u, *v, *p;
    Nparams residuals;
    clock_t  t_omp;
    int success = 0, niter=0, nitermax=1;
    double minalpha, maxalpha, frac = 1.0;
    int ntry;
    double min_step = model->line_search_min;

    Nmodel->stagnated = 0;
    residuals.nit = Nmodel->nit;

    ntry     = 6;
    minalpha = -2.5;
    maxalpha = -min_step;

    if (model->Newton==1) minalpha = -1.0;
    if (model->Newton==0) maxalpha = -0.0;

    //    ntry[0]     = 11;
    //    minalpha[0] = -2.25;
    //    maxalpha[0] = -0.25;

    //    ntry[0]     = 1;
    //    minalpha[0] = -1.0;
    //    maxalpha[0] = -1.0;

    //    printf("dalpha %lf %lf %d \n", maxalpha[0], minalpha[0], ntry[0]-1);

    t_omp = (double)omp_get_wtime();

    // Allocate
    u          = DoodzMalloc( nx*nzvx*sizeof(double) );
    v          = DoodzMalloc( nxvz*nz*sizeof(double) );
    p          = DoodzMalloc( ncx*ncz*sizeof(double) );

    // Save initial fields (before updates)
    ArrayEqualArray( u, mesh->u_in, nx*nzvx );
    ArrayEqualArray( v, mesh->v_in, nxvz*nz );
    ArrayEqualArray( p, mesh->p_in, ncx*ncz );

        // allocate array
        alphav = DoodzMalloc( ntry*sizeof(double) );
        rx     = DoodzMalloc( ntry*sizeof(double) );
        rz     = DoodzMalloc( ntry*sizeof(double) );
        rp     = DoodzMalloc( ntry*sizeof(double) );

        alpha    = maxalpha;
        dalpha   = fabs(maxalpha-minalpha)/(ntry-1);

        // Search for optimal relaxation parameters
        for( kk=0; kk<ntry; kk++ ) {

            // Update alpha
            if (kk>0) alpha -= dalpha;
            alphav[kk] = alpha;

            // Start from initial solutions
            ArrayEqualArray( mesh->u_in, u, nx*nzvx );
            ArrayEqualArray( mesh->v_in, v, nxvz*nz );
            ArrayEqualArray( mesh->p_in, p, ncx*ncz );

            // Test relaxed u-v-p solutions
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nzvx, nx )
            for( cc=0; cc<nzvx*nx; cc++) {
                if ( mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
                    mesh->u_in[cc] += alpha*dx[Stokes->eqn_u[cc]];
                }
            }

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nz, nxvz )
            for( cc=0; cc<nz*nxvz; cc++) {
                if ( mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
                    mesh->v_in[cc] += alpha*dx[Stokes->eqn_v[cc]];
                }
            }

#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, ncx, ncz )
            for( cc=0; cc<ncz*ncx; cc++) {

                if ( mesh->BCp.type[cc] != 30 && mesh->BCp.type[cc] != 0  && mesh->BCp.type[cc] != 31 ) {
                    mesh->p_in[cc] += alpha*dx[Stokes->eqn_p[cc]];
                }
            }

            // Apply Bc to Vx and Vz
            ApplyBC( mesh, model );

//            ArrayEqualArray( mesh->p_trial, mesh->p_in, ncx*ncz );

//            MinMaxArray( mesh->u_in, 1, nx*nzvx, "u in");
//            SumArray( mesh->u_in, 1, nx*nzvx, "u in");

            //------------------------------------------------------------------------------------------------------

            // Some stuff to be put on vertices                       < ---------------------- get P from centroids to vertices
            InterpCentroidsToVerticesDouble( mesh->T,    mesh->T_s,   mesh, model );
//            InterpCentroidsToVerticesDouble( mesh->p_in, mesh->P_s,   mesh, model );
            InterpCentroidsToVerticesDouble( mesh->d0_n,   mesh->d0_s,  mesh, model );
            InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s, mesh, model );

            // Update non-linearity
            UpdateNonLinearity( mesh, particles, topo_chain, topo, materials, model, Nmodel, scaling, 0, 1.0 );

            //------------------------------------------------------------------------------------------------------

            // Calculate residual
            //            model->free_surf_stab = 0;
            EvaluateStokesResidualDecoupled( Stokes, StokesA, StokesB, StokesC, StokesD, &residuals, mesh, *model, scaling, 1 );
            //            model->free_surf_stab = dummy;

            rx[kk] = residuals.resx;
            rz[kk] = residuals.resz;
            rp[kk] = residuals.resp;
            //            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.6e rz = %2.6e rp = %2.6e\n", alpha, rx[kk]* (scaling.F/pow(scaling.L,3)), rz[kk]* (scaling.F/pow(scaling.L,3)), rp[kk]*scaling.E );
            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alpha, rx[kk], rz[kk], rp[kk]);


        }

        
        // SEARCH --- 1

        // Look for the minimum predicted residuals
        double r, minxzp, minxz, minz, fxz, fx, fz, fxzp;
        int ixzp, ixz, ix, iz;
        double fxzp0;
        minxzp = sqrt( pow( rx[0],2 ) + pow( rz[0],2 ) + pow( rp[0],2 ) );
        minxz  = sqrt( pow( rx[0],2 ) + pow( rz[0],2 ) );
        minx   = rx[0];
        minz   = rz[0];
        ixz    = 0;
        ix     = 0;
        iz     = 0;
        ixzp   = 0;
        for( k=1; k<ntry; k++ ) {
            fxzp = sqrt( pow( rx[k],2 ) + pow( rz[k],2 ) + pow( rp[k],2 ) );
            fxz  = sqrt( pow( rx[k],2 ) + pow( rz[k],2 ) );
            fx   = rx[k];
            fz   = rz[k];
            if( fxzp < minxzp ) {
                minxzp = fxzp;
                ixzp   = k;
            }
            if( fxz < minxz ) {
                minxz = fxz;
                ixz  = k;
            }
            if( fx < minx ) {
                minx= fx;
                ix  = k;
            }
            if( fz < minz ) {
                minz= fz;
                iz  = k;
            }
        }
        
        // if the minmimun residuals are lower than starting ones, then success
        fxzp  = sqrt( pow( rx[ixzp],2 ) + pow( rz[ixzp],2 ) + pow( rp[ixzp],2 ) );
        fxzp0 = sqrt( pow(Nmodel->resx_f, 2) + pow( Nmodel->resz_f, 2) +  pow(Nmodel->resp_f, 2) );
        
        if ( fxzp < frac*fxzp0   ) { //|| rp[ix]<frac*Nmodel->resp
            alpha = alphav[ixzp];
            success = 1;
            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ixzp], rx[ixzp], rz[ixzp], rp[ixzp]);
        }

        // if the minmimun residuals are lower than starting ones, then success
        if (success==0 && rx[ixz] < frac*Nmodel->resx_f && rz[ixz]<frac*Nmodel->resz_f  ) { //|| rp[ix]<frac*Nmodel->resp
            alpha = alphav[ixz];
            success = 1;
            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ixz], rx[ixz], rz[ixz], rp[ixz]);
        }
        if (success==0 && rx[ix] < frac*Nmodel->resx_f) {
            alpha = alphav[ix];
            success = 1;
        }
        if (success==0 && rz[iz] < frac*Nmodel->resz_f) {
            alpha = alphav[iz];
            success = 1;
        }
            
        DoodzFree(rx);
        DoodzFree(rz);
        DoodzFree(rp);
        DoodzFree(alphav);

    if ( fabs(alpha)<1e-13 || success == 0 ) {
        printf( "Found minimum of the function -- cannot iterate further down\n" );
        if ( Nmodel->let_res_grow == 1 ) {
            printf("Letting residual grow\n!");
            Nmodel->stagnated = 0;
            alpha = min_step;
        }
        else {
            printf("Stagnating...\n!");
            Nmodel->stagnated = 1;
            alpha = 0.0;
        }
    }

    // Reset solutions
    ArrayEqualArray( mesh->u_in, u, nx*nzvx );
    ArrayEqualArray( mesh->v_in, v, nxvz*nz );
    ArrayEqualArray( mesh->p_in, p, ncx*ncz );

    // Free
    DoodzFree( u );
    DoodzFree( v );
    DoodzFree( p );
    printf("** Line search took = %f sec\n", (double)((double)omp_get_wtime() - t_omp) );

    return alpha;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokes( SparseMat *Stokes, DirectSolver* PardisoStokes ) {

    printf("---- Solve Stokes ----\n");

    clock_t t_omp;

    // Call direct solver
    t_omp = (double)omp_get_wtime();
    DirectStokes( Stokes, PardisoStokes, Stokes->b, Stokes->x );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDefect( SparseMat *Stokes, DirectSolver *PardisoStokes, Nparams* Nmodel, grid* mesh, params* model, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, scale scaling ) {

    printf("---- Solve Stokes in defect correction mode ----\n");

    double alpha = -1.0;
    double *dx;
    clock_t t_omp;

    dx = DoodzCalloc( Stokes->neq, sizeof(double));

    // Call direct solver
    t_omp = (double)omp_get_wtime();
    DirectStokes( Stokes, PardisoStokes, Stokes->F, dx );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

    // Line search
    if ( model->line_search == 1 ) {
        alpha = LineSearch( Stokes, dx, mesh, model, Nmodel, particles, topo_chain, topo, materials, scaling  );
    }

    // Return updated solution vector
    if ( Nmodel->stagnated == 0) {
        printf( "Correct x\n");
        ArrayPlusScalarArray( Stokes->x, alpha, dx, Stokes->neq );
    }

    DoodzFree(dx);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDecoupled( SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, SparseMat *Stokes, DirectSolver* PardisoStokes, params model, grid *mesh, scale scaling ) {

    printf("---- Solve Stokes in a decoupled/segregated fashion, lsolver = %d ----\n", model.lsolver);

    clock_t t_omp;

    // Call direct solver
    t_omp = (double)omp_get_wtime();

    if (model.lsolver== 0) DirectStokesDecoupled    ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes );
    if (model.lsolver== 1) KSPStokesDecoupled       ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes, Stokes, StokesA, StokesB, StokesC );
    if (model.lsolver== 2) KillerSolver             ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes, Stokes, StokesA, StokesB, StokesC );
    if (model.lsolver==-1) DirectStokesDecoupledComp( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes );

    ScaleVelocitiesRHSBack(StokesA, StokesD, Stokes->x);

    printf("** Time for direct decoupled Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDefectDecoupled( SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, SparseMat *Stokes, DirectSolver *PardisoStokes, Nparams* Nmodel, grid* mesh, params* model, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, scale scaling, SparseMat *JacobA, SparseMat *JacobB, SparseMat *JacobC ) {

    printf("---- Solve Stokes in a decoupled/segregated fashion and defect correction mode, lsolver = %d ----\n", model->lsolver);

    double alpha = -1.0;
    double *dx;//, *fu, *fp;
    clock_t t_omp;

    dx = DoodzCalloc( Stokes->neq,      sizeof(double));
//    fu = DoodzCalloc( Stokes->neq_mom,  sizeof(double));
//    fp = DoodzCalloc( Stokes->neq_cont, sizeof(double));

//    printf("NEQ %d \n", StokesA->neq );
//    for (int i =0;i<StokesA->neq;i++) {
//        if (fabs(StokesA->F[i])>1e5) printf("%2.2e \n", StokesA->F[i]);
////        printf("%2.2e \n", StokesA->b[i]);
//    }

    // Call direct solver
    t_omp = (double)omp_get_wtime();
    if ( model->lsolver == 0 ) DirectStokesDecoupled    ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes );
    if ( model->lsolver == 1 ) KSPStokesDecoupled       ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes, Stokes, JacobA, JacobB, JacobC );
    if ( model->lsolver == 2 ) KillerSolver             ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes, Stokes, JacobA, JacobB, JacobC );
    if ( model->lsolver ==-1 ) DirectStokesDecoupledComp( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

    if ( model->diag_scaling == 1 ) {
        ScaleVelocitiesRHSBack(StokesA, StokesD, dx);
    }

    // Line search
    if ( model->line_search == 1 ) {
        alpha = LineSearchDecoupled( Stokes, StokesA, StokesB, StokesC, StokesD, dx, mesh, model, Nmodel, particles, topo_chain, topo, materials, scaling  );
    }

    // Same solution extraction than in line search - consistent
    ExtractSolutions2( Stokes, mesh, model, dx, alpha );

//    InitialiseSolutionVector( mesh, Stokes, model );


//    // Return updated solution vector
//    if ( Nmodel->stagnated == 0) ArrayPlusScalarArray( Stokes->x, alpha, dx, Stokes->neq );
//
//    //    MinMaxArray( mesh->u_in, 1, mesh->Nx*(mesh->Nz+1), "u in");
//    SumArray( mesh->u_in, 1, mesh->Nx*(mesh->Nz+1), "u in");
//
//    ExtractSolutions( Stokes, mesh, model );

    //    MinMaxArray( mesh->u_in, 1, mesh->Nx*(mesh->Nz+1), "u in");
//    SumArray( mesh->u_in, 1, mesh->Nx*(mesh->Nz+1), "u in");

    //EvaluateStokesResidualDecoupled( Stokes, StokesA, StokesB, StokesC, StokesD, Nmodel, mesh, *model, scaling, 0 );


//    exit(1);

    DoodzFree(dx);
//    DoodzFree(fu);
//    DoodzFree(fp);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateStokesResidual( SparseMat *Stokes, Nparams *Nmodel, grid *mesh, params model, scale scaling, int quiet ) {

    int cc, nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    double celvol, Area=0.0;
    int ndofx=0, ndofz=0, ndofp=0;

    celvol = model.dx*model.dz;

    // Residuals
    double resx = 0.0;
    double resz = 0.0;
    double resp = 0.0;

    // Function evaluation
    BuildStokesOperator( mesh, model, 0, mesh->p_in,  mesh->u_in,  mesh->v_in, Stokes, 0 );

    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate(nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12) {
            ndofx++;
            resx        += Stokes->F[Stokes->eqn_u[cc]]*Stokes->F[Stokes->eqn_u[cc]];
            mesh->ru[cc] = Stokes->F[Stokes->eqn_u[cc]];
        }
    }
    Nmodel->resx = resx;

#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12) {
            ndofz++;
            resz        += Stokes->F[Stokes->eqn_v[cc]]*Stokes->F[Stokes->eqn_v[cc]];
            mesh->rv[cc] = Stokes->F[Stokes->eqn_v[cc]];
        }
    }
    Nmodel->resz = resz;

#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area        += celvol;
            resp        += Stokes->F[Stokes->eqn_p[cc]]*Stokes->F[Stokes->eqn_p[cc]];
            mesh->rp[cc] = Stokes->F[Stokes->eqn_p[cc]];
        }
    }
    Nmodel->resp = resp;

    // Sqrt
    Nmodel->resx =  sqrt(Nmodel->resx/ndofx);
    Nmodel->resz =  sqrt(Nmodel->resz/ndofz);
    Nmodel->resp =  sqrt(Nmodel->resp/ndofp);


    if ( quiet == 0 ) {
        printf("Fu = %2.6e\n", Nmodel->resx ); // Units of momentum
        printf("Fv = %2.6e\n", Nmodel->resz ); // Units of momentum
        printf("Fp = %2.6e\n", Nmodel->resp ); // Units of velocity gradient
    }

    if ( isnan(Nmodel->resx) || isnan(Nmodel->resz) || isnan(Nmodel->resp) ) {
        printf("Fu = %2.6e\n", Nmodel->resx ); // Units of momentum
        printf("Fv = %2.6e\n", Nmodel->resz ); // Units of momentum
        printf("Fp = %2.6e\n", Nmodel->resp );// Units of velocity gradient
        printf("Solve went wrong - Nan residuals...\nExiting...\n");
        exit(122);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateStokesResidualDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, Nparams *Nmodel, grid *mesh, params model, scale scaling, int quiet ) {

    int cc, nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    int ndofx=0, ndofz=0, ndofp=0;

    // Residuals
    double resx = 0.0;
    double resz = 0.0;
    double resp = 0.0;

    // Function evaluation
    if ( model.aniso == 0 ) BuildStokesOperatorDecoupled( mesh, model, 0, mesh->p_corr, mesh->p_in,  mesh->u_in,  mesh->v_in, Stokes, StokesA, StokesB, StokesC, StokesD, 0 );
    if ( model.aniso == 1 ) BuildJacobianOperatorDecoupled( mesh, model, 0, mesh->p_corr, mesh->p_in,  mesh->u_in,  mesh->v_in, Stokes, StokesA, StokesB, StokesC, StokesD, 0 );


    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes, StokesA ) private( cc ) firstprivate( nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12) {
            ndofx++;
            resx                          += StokesA->F[Stokes->eqn_u[cc]]*StokesA->F[Stokes->eqn_u[cc]];
            mesh->ru[cc]                   = StokesA->F[Stokes->eqn_u[cc]];
            StokesA->F[Stokes->eqn_u[cc]] *= StokesA->d[Stokes->eqn_u[cc]]; // Need to scale the residual here for Defect Correction formulation (F is the RHS)
//            if ( fabs(StokesA->F[Stokes->eqn_u[cc]]) > 1e5 ) printf( "%2.2e %2.2e %2.2e\n", StokesA->F[Stokes->eqn_u[cc]], StokesA->d[Stokes->eqn_u[cc]], mesh->ru[cc] );
        }
    }
    Nmodel->resx = resx;

#pragma omp parallel for shared( mesh, Stokes, StokesA ) private( cc ) firstprivate( nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12) {
            ndofz++;
            resz                          += StokesA->F[Stokes->eqn_v[cc]]*StokesA->F[Stokes->eqn_v[cc]];
            mesh->rv[cc]                   = StokesA->F[Stokes->eqn_v[cc]];
            StokesA->F[Stokes->eqn_v[cc]] *= StokesA->d[Stokes->eqn_v[cc]]; // Need to scale the residual here for Defect Correction formulation (F is the RHS)
//            printf("%2.2e %d\n", StokesA->F[Stokes->eqn_v[cc]], mesh->BCv.type[cc]);
        }
    }
    Nmodel->resz = resz;

#pragma omp parallel for shared( mesh, Stokes, StokesD ) private( cc ) firstprivate( ncz, ncx ) reduction(+:resp,ndofp)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            resp                                           += StokesD->F[Stokes->eqn_p[cc]- Stokes->neq_mom]*StokesD->F[Stokes->eqn_p[cc]- Stokes->neq_mom];
            mesh->rp[cc]                                    = StokesD->F[Stokes->eqn_p[cc]- Stokes->neq_mom];
            StokesD->F[Stokes->eqn_p[cc]- Stokes->neq_mom] *= StokesD->d[Stokes->eqn_p[cc]- Stokes->neq_mom]; // Need to scale the residual here for Defect Correction formulation (F is the RHS)
        }
    }
    Nmodel->resp = resp;

    // Sqrt
    Nmodel->resx =  sqrt(Nmodel->resx/ndofx);
    Nmodel->resz =  sqrt(Nmodel->resz/ndofz);
    Nmodel->resp =  sqrt(Nmodel->resp/ndofp);

    // Store initial residuals
    if ( Nmodel->nit == 0 ) {
        Nmodel->resx0 = Nmodel->resx;
        Nmodel->resz0 = Nmodel->resz;
        Nmodel->resp0 = Nmodel->resp;
    }

    if ( quiet == 0 ) {
        printf("Fu abs. = %2.6e --- Fu rel. = %2.6e\n", Nmodel->resx, Nmodel->resx/Nmodel->resx0 ); // Units of momentum
        printf("Fv abs. = %2.6e --- Fv rel. = %2.6e\n", Nmodel->resz, Nmodel->resz/Nmodel->resz0 ); // Units of momentum
        printf("Fp abs. = %2.6e --- Fp rel. = %2.6e\n", Nmodel->resp, Nmodel->resp/Nmodel->resp0 ); // Units of velocity gradient
    }

    if ( isnan(Nmodel->resx) || isnan(Nmodel->resz) || isnan(Nmodel->resp) ) {
        printf("Fu = %2.6e\n", Nmodel->resx * scaling.S * scaling.L ); // Units of momentum
        printf("Fv = %2.6e\n", Nmodel->resz * scaling.S * scaling.L ); // Units of momentum
        printf("Fp = %2.6e\n", Nmodel->resp * scaling.E * scaling.L * scaling.L ); // Units of velocity gradient
        printf("Solve went wrong - Nan residuals...\nExiting...\n");
        exit(122);
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateRHS( grid* mesh, params model, scale scaling, double RHO_REF ) {

    int c, k, l, c1, c2;
    int Nx, Nz, NzVx, NxVz;
    int Ncx, Ncz;
    double rhoVx, rhoVz, gx, gz, g, tet, x, z;
    double dx, dz;
    double inW, inE, inS, inN;

    Nx   = mesh->Nx;
    NxVz = mesh->Nx+1;
    Nz   = mesh->Nz;
    NzVx = mesh->Nz+1;
    Ncx=Nx-1; Ncz=Nz-1;
    dx   = mesh->dx;
    dz   = mesh->dz;

    int iPrW, iPrE, ixyN, ixyS, iVxS, iVxN;
    int iPrS, iPrN, ixyE, ixyW, iVyW, iVyE;
    
    int k1, c0;

    /* --------------------------------------------------------------------*/
    /* Here we calculate the forcing term -rho*gx on the finest grid level */
    /* --------------------------------------------------------------------*/

    // U POINTS
    for (l=0; l<NzVx; l++) {
        for (k=0; k<Nx; k++) {

            c  = k + l*Nx;
            c2 = k + (l-1)*Ncx;

            mesh->roger_x[c] = 0.0;
            gx  = model.gx;

            if (model.polar == 1 ) {
                x   = mesh->xg_coord[k];
                z   = mesh->zvx_coord[l];
                tet = atan(z/x);
                if (tet>0.0) gx     = model.gz*cos(tet);
                if (tet<0.0) gx     =-model.gz*cos(tet);
            }
            mesh->gx[c]         = gx;

            if (l>0 && l<NzVx-1) {

                if ( mesh->BCu.type[c] == -1 || mesh->BCu.type[c] == 2 || mesh->BCu.type[c] == -2 ) {

                    iPrW   = c2-1;
                    iPrE   = c2;
                    ixyN   = c;
                    ixyS   = c-Nx;
                    iVxS   = c-Nx;
                    iVxN   = c+Nx;
                    // Periodic BC !!!!!!!!!!!
                    if ( mesh->BCu.type[c] == -2 ) iPrW  = c2+Ncx-1;

                    // Free surface
                    inW=0.0, inE = 0.0, inS=0.0, inN = 0.0;
                    if (mesh->BCp.type[iPrW] == -1) inW = 1.0;
                    if (mesh->BCp.type[iPrE] == -1) inE = 1.0;
                    if (mesh->BCg.type[ixyS] != 30 && mesh->BCu.type[iVxS] != 13) inS = 1.0;
                    if (mesh->BCg.type[ixyN] != 30 && mesh->BCu.type[iVxN] != 13) inN = 1.0;

                    // Gravity force: take apparent viscosity (free surface correction)
                    if (model.polar==0) rhoVx             = 0.5*(mesh->rho_s[ixyS] + mesh->rho_s[ixyN]);
                    if (model.polar==1) rhoVx             = 0.5*(mesh->rho_n[iPrW] + mesh->rho_n[iPrE]);
                    mesh->roger_x[c]  = - gx * rhoVx;

                    // // Elastic force
                    // if ( model.iselastic == 1 && model.residual_form==0 ) {
                        // Inner nodes
                        // if (k>0 && k<Nx-1) {
                        //     if ( inE ) mesh->roger_x[c]  -= 1.0/dx * ( mesh->eta_n[iPrE] / (mesh->mu_n[iPrE]*model.dt)  * mesh->sxxd0[iPrE] );
                        //     if ( inW ) mesh->roger_x[c]  -= 1.0/dx * (-mesh->eta_n[iPrW] / (mesh->mu_n[iPrW]*model.dt)  * mesh->sxxd0[iPrW] );
                        //     if ( inN ) mesh->roger_x[c]  -= 1.0/dz * ( mesh->eta_s[ixyN] / (mesh->mu_s[ixyN]*model.dt)  * mesh->sxz0[ixyN] );
                        //     if ( inS ) mesh->roger_x[c]  -= 1.0/dz * (-mesh->eta_s[ixyS] / (mesh->mu_s[ixyS]*model.dt)  * mesh->sxz0[ixyS] );
                        // }
                    // }
                }
            }
            mesh->roger_x[c] = -mesh->roger_x[c];
        }
    }


    /* --------------------------------------------------------------------*/
    /* Here we calculate the forcing term -rho*gz on the finest grid level */
    /* --------------------------------------------------------------------*/
    
    double drhodz, om = model.free_surf_stab;
    double Wvalley = model.surf_Winc;
    double Vinc    = -model.surf_Vinc;


    // V POINTS
    for (l=0; l<Nz; l++) {
        for (k=0; k<NxVz; k++) {

            c  = k + l*(NxVz);
            c1 = k + l*(Nx)-1;
            c2 = k-1 + (l-1)*(Nx-1);

            mesh->roger_z[c]  = 0.0;
            gz  = model.gz;

            if (model.polar == 1 ) {
                x   = mesh->xvz_coord[k];
                z   = mesh->zg_coord[l];
                tet = atan(z/x);
                if (tet>0.0) gz     = model.gz*sin(tet);
                if (tet<0.0) gz     =-model.gz*sin(tet);
            }
            mesh->gz[c]         = gz;

            if (k>0 && k<NxVz-1) {

                if ( mesh->BCv.type[c] == -1  ) {

                    iPrS   = c2;
                    iPrN   = c2+Ncx;
                    ixyE   = c1+1;
                    ixyW   = c1;
                    iVyW   = c-1;
                    iVyE   = c+1;

                    // Free surface
                    inS=0.0, inN = 0.0, inW=0.0, inE = 0.0;
                    if (mesh->BCp.type[iPrS] == -1) inS = 1.0;
                    if (mesh->BCp.type[iPrN] == -1) inN = 1.0;
                    if (mesh->BCg.type[ixyW] != 30 && mesh->BCv.type[iVyW] != 13 ) inW = 1.0;
                    if (mesh->BCg.type[ixyE] != 30 && mesh->BCv.type[iVyE] != 13 ) inE = 1.0;

                    // Gravity force: use apparent density (free surface correction)
//                    if (model.polar==0) rhoVz             = 0.5 * (mesh->rho_s[ixyW] + mesh->rho_s[ixyE]);  // USE THIS ALWAYS but not for polar!
//                    if (model.polar==1) rhoVz             = 0.5 * (mesh->rho_n[iPrS] + mesh->rho_n[iPrN]);
                    rhoVz             = 0.5 * (mesh->rho_n[iPrS] + mesh->rho_n[iPrN]);
                    
                    mesh->roger_z[c]  = - gz * rhoVz;

//                    // Additional stabilisation term from surface processes
//                    if (model.surf_processes >= 1 && ( mesh->BCp.type[iPrS] == 30 || mesh->BCp.type[iPrS] == 31 || mesh->BCp.type[iPrN] == 30 || mesh->BCp.type[iPrN] == 31 ) ) {
//                        drhodz = ( mesh->rho_n[iPrN] - mesh->rho_n[iPrS] ) / dz;
//                        if ( fabs(mesh->xvz_coord[k]) <= 0.5*Wvalley ) {
//                            mesh->roger_z[c]  += -1.0*om*Vinc*model.dt*gz*drhodz;
//                        }
//                    }
                    
                    // // Elastic force
                    // if  (model.iselastic == 1 && model.residual_form==0 ) {
                        // if ( l>0 && l<Nz-1 ) {
                        //     if ( inN ) mesh->roger_z[c]  -= 1.0/dz * (  mesh->eta_n[iPrN] / (mesh->mu_n[iPrN] *model.dt ) * (mesh->szzd0[iPrN]) );
                        //     if ( inS ) mesh->roger_z[c]  -= 1.0/dz * ( -mesh->eta_n[iPrS] / (mesh->mu_n[iPrS] *model.dt ) * (mesh->szzd0[iPrS]) );
                        //     if ( inE ) mesh->roger_z[c]  -= 1.0/dx * (  mesh->eta_s[ixyE] / (mesh->mu_s[ixyE] *model.dt ) *  mesh->sxz0[ixyE]) ;
                        //     if ( inW ) mesh->roger_z[c]  -= 1.0/dx * ( -mesh->eta_s[ixyW] / (mesh->mu_s[ixyW] *model.dt ) *  mesh->sxz0[ixyW]) ;
                        // }
                    // }
                }
            }
            mesh->roger_z[c] = -mesh->roger_z[c];
        }
    }


    /* ------------------------------------------------------------------------*/
    /* Here we calculate the forcing term of the continuity and heat equation  */
    /* ------------------------------------------------------------------------*/

    // P-T RHS
    for (l=0; l<Nz-1; l++) {
        for (k=0; k<Nx-1; k++) {

            c  = k + l*(Nx-1);

            if ( mesh->BCp.type[c] != 30 && mesh->BCp.type[c] != 31 ) {

                // Heat equation
                mesh->rhs_t[c] = mesh->T[c];
                mesh->rhs_p[c] = 0.0;
                mesh->Qrho[c]  = 0.0;

                // Continuity equation
                mesh->rhs_p[c] = 0.0;

                if (model.compressible == 1 ) {
                    if (mesh->comp_cells[c] == 1) {
                        if ( model.VolChangeReac == 0 ) mesh->rhs_p[c] += mesh->p0_n[c]*mesh->bet_n[c]/model.dt;
                        if ( model.VolChangeReac == 1 ) mesh->rhs_p[c] += log(mesh->rho0_n[c])/model.dt;
                        if (model.adiab_heat > 0 ) {
                            mesh->rhs_p[c] += mesh->divth0_n[c];
                        }
                    }
                }
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvalNumberOfEquations( grid* mesh, SparseMat *Stokes ) {
    int inc = 0, incP=0;
    int k,l,kk;

    // Allocate
    Stokes->eqn_u = DoodzMalloc((mesh->Nz+1)*mesh->Nx * sizeof(int));
    Stokes->eqn_v = DoodzMalloc(mesh->Nz*(mesh->Nx+1) * sizeof(int));
    Stokes->eqn_p = DoodzMalloc((mesh->Nz-1)*(mesh->Nx-1) * sizeof(int));

    // Number U dofs
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            Stokes->eqn_u[kk] = -1;
            if ( mesh->BCu.type[kk] ==-1 || mesh->BCu.type[kk] ==2 || mesh->BCu.type[kk] ==-2  ) {
                Stokes->eqn_u[kk] = inc;
                inc++;
            }
            // Copy equation numbers from the West side
            if ( mesh->BCu.type[kk] ==-12) {
                Stokes->eqn_u[kk] = Stokes->eqn_u[kk-mesh->Nx+1];
            }
        }
    }

    // Number V dofs
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            Stokes->eqn_v[kk] = -1;
            if ( mesh->BCv.type[kk] ==-1 || mesh->BCv.type[kk] ==2 ) {
                Stokes->eqn_v[kk] = inc;
                inc++;
            }
        }
    }

    Stokes->neq_mom = inc;

    // Number P dofs
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            Stokes->eqn_p[kk] = -1;
            if ( mesh->BCp.type[kk] ==-1 ) {
                Stokes->eqn_p[kk] = inc;
                inc++;
                incP++;
            }
        }
    }
    Stokes->neq = inc;
    Stokes->neq_cont = incP;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseSolutionVector( grid* mesh, SparseMat *Stokes, params* model ) {
    int k,l,kk;

    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            if ( mesh->BCu.type[kk] ==-1 || mesh->BCu.type[kk] ==2 || mesh->BCu.type[kk] ==-2  ) {
                Stokes->x[ Stokes->eqn_u[kk] ] = mesh->u_in[kk];
            }

        }
    }

    // Number V dofs
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            if ( mesh->BCv.type[kk] ==-1 || mesh->BCv.type[kk] ==2 ) {
                Stokes->x[ Stokes->eqn_v[kk] ] = mesh->v_in[kk];
            }
        }
    }

    // Number P dofs
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            if ( mesh->BCp.type[kk] ==-1 ) {
                Stokes->x[ Stokes->eqn_p[kk] ] = mesh->p_in[kk];
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ExtractDiagonalScale(SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD ) {

// Extract Diagonal of A - Viscous block
int i, j, locNNz;
int I1, J1;

#pragma omp parallel for shared(StokesA, StokesB) private(I1, J1, i, j, locNNz )
    for (i=0;i<StokesA->neq; i++) {
        I1     = StokesA->Ic[i];
        locNNz = StokesA->Ic[i+1] - StokesA->Ic[i];
        for (J1=0;J1<locNNz; J1++) {
            j = StokesA->J[I1 + J1];
            if (i==j) StokesA->d[i] = StokesA->A[I1 + J1];
        }
        StokesB->d[i] = StokesA->d[i] ;
    }

//for (i=0;i<StokesA->neq; i++)
//    if (fabs(StokesA->d[i])>1e5) printf("WTF %2.2e\n", StokesA->d[i] );
//    // Extract Diagonal of D - Pressure block
//#pragma omp parallel for shared(StokesD) private(i)
//    for (i=0;i<StokesD->neq; i++) {
//        StokesC->d[i] = 1.0;
//    }

//    MinMaxArray(StokesA->d, 1, StokesA->neq, "diag. A" );

#pragma omp parallel for shared(StokesC,StokesD) private(I1,J1, i, j, locNNz )
    for (i=0;i<StokesD->neq; i++) {
        I1     = StokesD->Ic[i];
        locNNz = StokesD->Ic[i+1] - StokesD->Ic[i];
        for (J1=0;J1<locNNz; J1++) {
            j = StokesD->J[I1 + J1];
            if (i==j) StokesD->d[i] = StokesD->A[I1 + J1];
        }
        StokesC->d[i] = StokesD->d[i];
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ScaleMatrix(SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD ) {

    int i, j, locNNz;
    int I1, J1;

// Scale A
#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNz )
for (i=0;i<StokesA->neq; i++) {

    StokesA->b[i] *= StokesA->d[i];  // scale RHS
    StokesA->F[i] *= StokesA->d[i];  // scale RHS

    I1     = StokesA->Ic[i];
    locNNz = StokesA->Ic[i+1] - StokesA->Ic[i];
    for (J1=0;J1<locNNz; J1++) {
        j = StokesA->J[I1 + J1];
        StokesA->A[I1 + J1] *= StokesA->d[i]*StokesA->d[j];
    }
}

// Scale B
#pragma omp parallel for shared(StokesB,StokesA) private(I1,J1, i, j, locNNz )
for (i=0;i<StokesB->neq; i++) {

    StokesB->b[i] *= StokesB->d[i]; // scale RHS
    StokesB->F[i] *= StokesB->d[i]; // scale RHS

    I1     = StokesB->Ic[i];
    locNNz = StokesB->Ic[i+1] - StokesB->Ic[i];
    for (J1=0;J1<locNNz; J1++) {
        j = StokesB->J[I1 + J1];
        StokesB->A[I1 + J1] *= StokesA->d[i]*StokesD->d[j];
    }
}

// Scale C
#pragma omp parallel for shared(StokesC,StokesA) private(I1,J1, i, j, locNNz )
for (i=0;i<StokesC->neq; i++) {

    StokesC->b[i] *= StokesC->d[i]; // scale RHS
    StokesC->F[i] *= StokesC->d[i]; // scale RHS

    I1     = StokesC->Ic[i];
    locNNz = StokesC->Ic[i+1] - StokesC->Ic[i];
    for (J1=0;J1<locNNz; J1++) {
        j = StokesC->J[I1 + J1];
        StokesC->A[I1 + J1] *= StokesD->d[i]*StokesA->d[j];
    }
}

    // Scale D
#pragma omp parallel for shared(StokesD) private(I1,J1, i, j, locNNz )
    for (i=0;i<StokesC->neq; i++) {

        StokesD->b[i] *= StokesD->d[i]; // scale RHS
        StokesD->F[i] *= StokesD->d[i]; // scale RHS

        I1     = StokesD->Ic[i];
        locNNz = StokesD->Ic[i+1] - StokesD->Ic[i];
        for (J1=0;J1<locNNz; J1++) {
            j = StokesD->J[I1 + J1];
            StokesD->A[I1 + J1] *= StokesD->d[i]*StokesD->d[j];
        }
    }


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ScaleVelocitiesRHSBack(SparseMat * StokesA, SparseMat * StokesD, double* x) {

    int k;
#pragma omp parallel for shared(StokesA, x)
    for (k=0; k<StokesA->neq;k++) {
        x[k] *= StokesA->d[k];
        //        Stokes->b[k] /= StokesA->d[k];
    }

#pragma omp parallel for shared(StokesD, x)
    for (k=0; k<StokesD->neq;k++) {
        x[StokesA->neq + k] *= StokesD->d[k];
        //        Stokes->b[k] /= StokesA->d[k];
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
