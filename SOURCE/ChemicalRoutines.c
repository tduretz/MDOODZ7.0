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
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "header_MDOODZ.h"

#ifdef _OMP_
#include "omp.h"
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int EvalNumEqChem( grid *mesh, int *eqn_t ) {
    
    int inc = 0;
    int k,l,kk;
    
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            eqn_t[kk] = -1;
            if ( mesh->BCc.type[kk] != 30 ) {
                eqn_t[kk] = inc;
                inc++;
            }
        }
    }
    return inc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void AddCoeffChem( int* J, double*A, int neq, int jeq, int *nnzc, double coeff ) {
    J[*nnzc]         = jeq;
    A[*nnzc]         = coeff;
    *nnzc           += 1;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ChemicalDirectSolve( grid *mesh, params model, markers *particles, mat_prop* materials, double dt, scale scaling ) {
    
    // Energy matrix
    double  *A;
    double  *b, *bbc, *x;      // Coefs, rhs, sol
    int *Ic, *J;      // lines, colums
    int k, l, c1, c2, c3, nx, nz, nxvz, nzvx, ncx, ncz;
    double AW, AE, AS, AN;
    int *eqn_t, nnzc, eqn = 0;
    double val, theta=1.0;
    double transient=1.0;
    double aW=0.0, aE=0.0;
    double yW, yE, xn, zn, tet;
    int neq, nnz, p, cond;
    
    double Xeq, Xeq_phase;
    double *tau_kin, *prod, tau_kin_phase;
    double Pr,  Pr_phase;
    double dPr, dPr_phase;
//    double tau_kin = model.tau_kin;
//    double k_chem  = model.k_chem;
//    double dPr     = model.dPr;
//    double Pr      = model.Pr;
    
    double ks;//   = k_chem;
    double mink;// = k_chem;
    
    // Pre-calculate FD coefs
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    
    nx   = mesh->Nx;
    nz   = mesh->Nz;
    ncx  = nx-1;
    ncz  = nz-1;
    nxvz = nx+1;
    nzvx = nz+1;
    
    // Solve system using CHOLMOD
    cholmod_common c;
    cholmod_start( &c );
    cholmod_factor *Afact;
    cs_di *At;
    
    // Number of equations
    eqn_t = DoodzCalloc(ncx*ncz, sizeof(int));
    neq   = EvalNumEqChem( mesh, eqn_t );
    bbc   = DoodzCalloc(neq    , sizeof(double));
    
    // Guess number of non-zeros
    nnz   = 5*neq;
    
    // Set nnz count to 0
    nnzc = 0;
    
    // Allocate arrays
    b  = DoodzCalloc(neq     , sizeof(double));
    x  = DoodzCalloc(neq     , sizeof(double));
    Ic = DoodzCalloc((neq+1) , sizeof(int));
    J  = DoodzMalloc(nnz     * sizeof(int));
    A  = DoodzMalloc(nnz     * sizeof(double));
    tau_kin = DoodzCalloc(ncx*ncz, sizeof(double));
    prod    = DoodzCalloc(ncx*ncz, sizeof(double));
    
    //----------------------------------------------------//
    
    // Build right-hand side
    
    for( c2=0; c2<ncz*ncx; c2++) {
        k   = mesh->kp[c2];
        l   = mesh->lp[c2];
        c1  = k + l*nx;
        c3  = k + l*nxvz;
        eqn = eqn_t[c2];
        
        if ( mesh->BCc.type[c2] != 30 ) {
            
            // Contribution from transient term
            if ( model.UnsplitDiffReac == 0 ) b[eqn]  = transient * mesh->X_n[c2]/ dt;
            
            if ( model.UnsplitDiffReac == 1 ) {
                b[eqn]  = transient * mesh->X0_n[c2]/ dt;

                // Pr, dPr, tau_kin depends on phases
                Pr          = 0.0;
                dPr         = 0.0;
                tau_kin[c2] = 0.0;

                for ( p=0; p<model.Nb_phases; p++ ) {
                    cond =  fabs(mesh->phase_perc_n[p][c2])>1.0e-13;

                    if ( cond == 1 ) Pr          += mesh->phase_perc_n[p][c2] * materials->Pr[p];
                    if ( cond == 1 ) dPr         += mesh->phase_perc_n[p][c2] * materials->dPr[p];
                    if ( cond == 1 ) tau_kin[c2] += mesh->phase_perc_n[p][c2] * materials->tau_kin[p];
                }

                // Contribution for source
                Xeq     = 0.5*erfc((-mesh->p_in[c2] + Pr)/dPr) ;
                
                if (Xeq >= mesh->X0_n[c2]                      ) prod[c2] = 1.0;
                if (Xeq <  mesh->X0_n[c2] && model.NoReturn==1 ) prod[c2] = 0.0;

    //           printf("Xeq = %2.2e p_in = %2.2e Pr = %2.2e dPr = %2.2e\n", Xeq, mesh->p_in[c2]*scaling.S, Pr*scaling.S, dPr*scaling.S);

                b[eqn] += prod[c2] * Xeq / tau_kin[c2];
            }

        }
        
    }
    
    MinMaxArray(b, scaling.E, ncz*ncx, "b chemical");
    
    //----------------------------------------------------//

    printf("Assembling Energy matrix... with %d discrete equations\n", neq);
    
    // LOOP ON THE GRID TO CALCULATE FD COEFFICIENTS
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            
            c1  = k   + (l+1)*nx;
            c2  = k   + l*ncx;
            c3  = k   + l*nxvz+1;
            
            if (mesh->BCc.type[c2] != 30 ) {
                
                // Central coefficient
                eqn     = eqn_t[c2];
                Ic[eqn] = nnzc;
                
                AW    = mesh->kc_x[c1]       ;
                AE    = mesh->kc_x[c1+1]     ;
                AS    = mesh->kc_z[c3]       ;
                AN    = mesh->kc_z[c3+nxvz]  ;
                
                // Contribution to S dof
                if ( l>0 ) {
                    if (mesh->BCc.type[c2-ncx] != 30 ) {
                        val   = -theta*AS*one_dz_dz;
                        AddCoeffChem( J, A, eqn, eqn_t[c2-ncx],   &nnzc,  val );
                    }
                }
                
                // Periodic on the right side
                if (mesh->BCc.type[c2]==-2 && k==ncx-1 ) {
                    val = -theta*AE*one_dx_dx;
                    AddCoeffChem( J, A, eqn, eqn_t[c2-(ncx-1)],   &nnzc,  val );
                }
                
                // Contribution to W dof
                if (k>0) {
                    if (mesh->BCc.type[c2-1] != 30 ) {
                        val   = -theta*AW*one_dx_dx;
                        AddCoeffChem( J, A, eqn, eqn_t[c2-1],   &nnzc,  val );
                    }
                }
                
                // ----------------- Center point ----------------- //
                if ( model.UnsplitDiffReac == 1 ) val = transient/dt + prod[c2]/tau_kin[c2];
                if ( model.UnsplitDiffReac == 0 ) val = transient/dt;
                
                // Average conductivity for surface values (avoid zero conductivity)
                ks = 0.25*(AE+AW+AN+AS);
                if ( ks < mink ) {
                    printf("ACHTUNG: interloplated surface conductivity is set lower cutoff value!!\n");
                    ks = mink;
                }
                
                // Flux from sides: SOUTH
                if ( l>0 ) {
                    
                    if (mesh->BCc.type[c2-ncx] != 30 ) {
                        val +=  AS*one_dz_dz;
                    }
                    // Dirichlet contribution (free surface)
                    else {
                        if (mesh->BCc.type[c2] == 1 ) {
                            val      += 2.0*ks*one_dz_dz;
                            bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCc.val[c2];
                        }
                    }
                }
                
                // Flux from sides: WEST
                if ( (k>0) || (k==0 && mesh->BCc.type[c2]==-2) ) {
                    if (mesh->BCc.type[c2-1] != 30 ) {
                        val +=  AW*one_dx_dx;
                    }
                    // Dirichlet contribution (free surface)
                    else {
                        if (mesh->BCc.type[c2] == 1 ) {
                            val      += 2.0*ks*one_dx_dx;
                            bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCc.val[c2];
                        }
                    }
                }
                
                // Flux from sides: EAST
                if ( (k<ncx-1) || (k==ncx-1 && mesh->BCc.type[c2]==-2) ) {
                    if (mesh->BCc.type[c2+1] != 30 ) {
                        val +=  AE*one_dx_dx;
                    }
                    // Dirichlet contribution (free surface)
                    else {
                        if (mesh->BCc.type[c2] == 1 ) {
                            val      += 2.0*ks*one_dx_dx;
                            bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCc.val[c2];
                        }
                    }
                }
                
                // Flux from sides: NORTH
                if (l<ncz-1)  {
                    if ( mesh->BCc.type[c2+ncx] != 30 ) {
                        val +=  AN*one_dz_dz;
                    }
                    // Dirichlet contribution (free surface)
                    else {
                        if (mesh->BCc.type[c2] == 1 ) {
                            val      += 2.0*ks*one_dz_dz;
                            bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCc.val[c2];
                        }
                    }
                }
                
                // -------------- Dirichlet contibutions -------------- //
                
                if (l==0 && mesh->BCc.typS[k]==1) {
                    val      += 2.0*AS*one_dz_dz;
                    bbc[eqn] += 2.0*AS*one_dz_dz * mesh->BCc.valS[k];
                }
                
                if (l==0 && mesh->BCc.typS[k]==2) {
                    bbc[eqn] += -1.0/mesh->dx * mesh->BCc.valS[k];
                }
                
                if (k==0 && mesh->BCc.typW[l]==1) {
                    val      += 2.0*AW*one_dx_dx;
                    bbc[eqn] += 2.0*AW*one_dx_dx * mesh->BCc.valW[l];
                }
                
                if (k==ncx-1 && mesh->BCc.typE[l]==1) {
                    val      += 2.0*AE*one_dx_dx;
                    bbc[eqn] += 2.0*AE*one_dx_dx * mesh->BCc.valE[l];
                }
                
                if (l==ncz-1 && mesh->BCc.typN[k]==1)  {
                    val      += 2.0*AN*one_dz_dz;
                    bbc[eqn] += 2.0*AN*one_dz_dz * mesh->BCc.valN[k];
                }
                
                // Zero flux in radial coordinates
                if (k==0 && model.polar==1 && mesh->BCc.type[c2+ncx] != 30 ) {
                    
                    xn   = mesh->xg_coord[k]-model.dx/2.0;
                    zn   = mesh->zc_coord[l];
                    tet  = atan(zn/xn);
                    if (tet<0.0) tet  *= -1;
                    yW   = atan(M_PI/2.0-tet) * model.dx;
                    aW   = yW/model.dz;
                    val +=  (1.0-aW)*AW*one_dx_dx;
                }
                
                if (k==ncx-1 && model.polar==1) {
                    
                    xn   = mesh->xg_coord[k+1]+model.dx/2.0;
                    zn   = mesh->zc_coord[l];
                    tet  = atan(zn/xn);
                    if (tet<0.0) tet  *= -1;
                    yE   = atan(M_PI/2.0-tet) * model.dx;
                    aE   = yE/model.dz;
                    val +=  (1.0-aE)*AE*one_dx_dx;
                }
                
                AddCoeffChem( J, A, eqn, eqn_t[c2],   &nnzc, val );
                
                // -------------- Dirichlet contibutions -------------- //
                
                // Contribution to E dof
                if (k<ncx-1) {
                    if (mesh->BCc.type[c2+1] != 30 ) {
                        val   = -theta*AE*one_dx_dx;
                        AddCoeffChem( J, A, eqn, eqn_t[c2+1],   &nnzc,  val );
                    }
                }
                
                // Periodic on the left side
                if (mesh->BCc.type[c2]==-2 && k==0 ) {
                    val = -theta*AW*one_dx_dx;
                    AddCoeffChem( J, A, eqn, eqn_t[c2+(ncx-1)],   &nnzc,  val );
                }
                
                // Contribution to N dof
                if (l<ncz-1)  {
                    if (mesh->BCc.type[c2+ncx] != 30 ) {
                        val   = -theta*AN*one_dz_dz;
                        
                        // Zero flux in radial coordinates
                        if (k==0 && model.polar==1) {
                            val +=  (aW-1.0)*AW*one_dx_dx;
                        }
                        
                        if (k==ncx-1 && model.polar==1) {
                            val +=  (aE-1.0)*AE*one_dx_dx;
                        }
                        
                        AddCoeffChem( J, A, eqn, eqn_t[c2+ncx],   &nnzc,  val );
                    }
                }
            }
        }
    }
    
    // Final index
    Ic[neq] = nnzc;
    
    // Resize arrays to proper number of non-zeros
    double *bufd;
    int *bufi;
    bufi = DoodzRealloc(J, nnzc*sizeof(int));
    bufd = DoodzRealloc(A, nnzc*sizeof(double));
    J    = bufi;
    A    = bufd;
    printf("Chemical diffusion --> System size: ndof = %d, nnz = %d\n", neq, nnzc);
    
    // ------------------------------------------- SOLVER ------------------------------------------- //
    
    // Compute transpose of A matrix (used in polar mode only)
    At = TransposeA( &c, A, Ic, J, neq, nnzc );
    
    // Factor matrix only at the first step
    Afact = FactorEnergyCHOLMOD( &c, At, A, Ic, J, neq, nnzc, model.polar );
    
    // Solve
    ArrayPlusScalarArray( b, 1.0, bbc, neq );
    SolveEnergyCHOLMOD( &c, At, Afact, x, b, neq, nnzc, model.polar );
    
    // ------------------------------------------- SOLVER ------------------------------------------- //
    
    // Extract temperature T from solution vector x, compute temperature increments dT and integrate heat increments dUt
    //#pragma omp parallel for shared( mesh, scaling, x ) private( c2, eqn ) firstprivate( ncx, ncz, zero_celsius, model ) reduction (+:dUt)
    for( c2=0; c2<ncx*ncz; c2++) {
        eqn = eqn_t[c2];
        if ( mesh->BCc.type[c2] != 30 ) {
            //                mesh->dT[c2] = 1.0*(x[eqn] -  mesh->T[c2]);
            mesh->X_n[c2]  = x[eqn];
            if (mesh->X_n[c2] < 0.0) {
                printf("Negative X --- Are you crazy! (ChemicalDirectSolve)\n");
                exit(1);
            }
        }
        else {
            // Outside the free surface
            mesh->X_n[c2]  = 0.0;
        }
    }
    MinMaxArrayTag( mesh->X_n, 1.0, (mesh->Nx-1)*(mesh->Nz-1), "X_n", mesh->BCc.type );
    
    // Free arrays
    DoodzFree(x);
    DoodzFree(b);
    DoodzFree(Ic);
    DoodzFree(J);
    DoodzFree(A);
    DoodzFree(tau_kin);
    DoodzFree(prod);
    
    // Free factors
    cholmod_free_factor ( &Afact, &c) ;
    cs_spfree(At);
    cholmod_finish( &c ) ;
    DoodzFree(bbc);
    DoodzFree(eqn_t);
}
