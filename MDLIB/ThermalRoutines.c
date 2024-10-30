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
#include "mdoodz-private.h"

#ifdef _OMP_
#include "omp.h"
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int EvalNumEqTemp( grid *mesh, int *eqn_t ) {

    int inc = 0;
    int k,l,kk;

    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            eqn_t[kk] = -1;
            if ( mesh->BCt.type[kk] != 30 ) {
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


void AddCoeffThermal( int* J, double*A, int neq, int jeq, int *nnzc, double coeff ) {
    J[*nnzc]         = jeq;
    A[*nnzc]         = coeff;
    *nnzc           += 1;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EnergyDirectSolve( grid *mesh, params model, double *rhs_t, markers *particles, double dt, int shear_heating, int adiab_heating, scale scaling, int nit ) {

    // Energy matrix
    double  *A;
    double  *b, *bbc, *x;      // Coefs, rhs, sol
    int *Ic, *J;      // lines, colums
    int k, l, c0, c1, c2, c3, nx, nz, nxvz, nzvx, ncx, ncz;
    double AW, AE, AS, AN, rhoCp;
    int *eqn_t, nnzc, eqn;
    double val, theta=1.0;
    double dW = 0.0, dUt = 0.0, dUe = 0.0;
    double Wtot, Wth, Wel;
    eqn = 0;
    double transient=1.0;
    double Hr = 1.0, zero_celsius, ks, mink, maxk;
    double aW=0.0, aE=0.0;
    double yW, yE, xn, zn, tet;
    double syyd, eyyd, eyyd_el, eyyd_diss;

    // Pre-calculate FD coefs
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;

    nx   = mesh->Nx;
    nz   = mesh->Nz;
    ncx  = nx-1;
    ncz  = nz-1;
    nxvz = nx+1;
    nzvx = nz+1;

    // Get minimum value of k - safety: will be used a a lower limit if interpolated surface conductivity is too low
    MinMaxArrayVal( mesh->kx, nx*(nz+1), &mink, &maxk );

    // Solve system using CHOLMOD
    cholmod_common c;
    cholmod_start( &c );
    cholmod_factor *Afact;
    cs_di *At;

    double *Hs, *Ha, dexz_th, dexz_el, dexz_tot, dexx_th, dexx_el, dexx_tot, diss_limit=1.0e-8/(scaling.S*scaling.E);
    double dezz_el, dezz_th, dezz_tot, deyy_el, deyy_th, deyy_tot;
    int neq, nnz, it;

    // Number of equations
    eqn_t = DoodzCalloc(ncx*ncz, sizeof(int));
    neq   = EvalNumEqTemp( mesh, eqn_t );
    bbc   = DoodzCalloc(neq     , sizeof(double));

    // Guess number of non-zeros
    nnz   = 5*neq;

    for ( it=0; it<nit; it++ ) {

        // Set nnz count to 0
        nnzc = 0;

        // Allocate arrays
        b  = DoodzCalloc(neq     , sizeof(double));
        x  = DoodzCalloc(neq     , sizeof(double));
        Ic = DoodzCalloc((neq+1) , sizeof(int));
        J  = DoodzMalloc(nnz     * sizeof(int));
        A  = DoodzMalloc(nnz     * sizeof(double));
        Hs = DoodzCalloc(ncx*ncz, sizeof(double));
        Ha = DoodzCalloc(ncx*ncz, sizeof(double));

        //----------------------------------------------------//

        // Build right-hand side
//#pragma omp parallel for shared ( b, Hs, Ha, eqn_t, mesh ) private ( l, k, c0, c1, c2, c3, eqn, rhoCp, dexx_el, dexx_th, dexx_tot, dezz_el, dezz_th, dezz_tot, deyy_el, deyy_th, deyy_tot, dexz_el, dexz_th, dexz_tot, Wth, Wel, Wtot, syyd, eyyd, eyyd_el, eyyd_diss  ) firstprivate ( nx, ncx, nxvz, shear_heating, adiab_heating, model, transient, dt, Hr, diss_limit, scaling ) reduction( +:dUe, dW )
        for( c2=0; c2<ncz*ncx; c2++) {
            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c0  = (k+1) + (l+1)*(ncx+2); // refers to extended arrays
            c1  = k + l*nx;
            c3  = k + l*nxvz;
            eqn = eqn_t[c2];

            if ( mesh->BCT_exp.type[c0] != 30 ) {

                // Contribution from transient
                // printf("mesh->Cp = %2.2e\n", mesh->Cp[c2]*scaling.Cp);
                rhoCp   = mesh->rho_n[c2]*mesh->Cp[c2];
                b[eqn]  = transient * rhoCp * mesh->T[c2]/ dt;

                // Contribution from radiogenic sources
                b[eqn] += Hr*mesh->Qr[c2];

                // Contribution from dissipation
                if ( shear_heating == 1 ) {
                    b[eqn] += mesh->Wdiss[c2];
                    Hs[c2]  = mesh->Wdiss[c2];
                }

                // Contribution from adiabatic heat
                if ( adiab_heating == 1 ) Ha[c2]  = mesh->T[c2]*mesh->alp[c2]*0.5*(mesh->v_in[c3+1] + mesh->v_in[c3+1+nxvz])*model.gz*mesh->rho_n[c2];
                if ( adiab_heating == 2 ) Ha[c2]  = mesh->T[c2]*mesh->alp[c2]*(mesh->p_in[c2] - mesh->p0_n[c2])/model.dt; //mesh->T[c2]*mesh->alp[c2]*(mesh->dp[c2])/model.dt;
                if ( adiab_heating  > 0 ) b[eqn] += Ha[c2];
            }
        }

        MinMaxArrayTag( mesh->T, scaling.T, ncx*ncz, "T0", mesh->BCt.type );
        MinMaxArrayTag( Ha, scaling.S*scaling.E, ncx*ncz, "Ha", mesh->BCt.type );
        MinMaxArrayTag( Hs, scaling.S*scaling.E, ncx*ncz, "Hs", mesh->BCt.type );
        mesh->Uelastic += dUe;
        mesh->Work  += dW;

        //----------------------------------------------------//
        if ( it == 0 ) {

            printf("Assembling Energy matrix... with %d discrete equations\n", neq);

            // LOOP ON THE GRID TO CALCULATE FD COEFFICIENTS
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {

                    c0  = (k+1) + (l+1)*(ncx+2);
                    c1  = k   + (l+1)*nx;
                    c2  = k   + l*ncx;
                    c3  = k   + l*nxvz+1;

                    if (mesh->BCT_exp.type[c0] != 30 ) {

                        // Central coefficient
                        eqn     = eqn_t[c2];
                        Ic[eqn] = nnzc;

                        rhoCp = mesh->rho_n[c2]*mesh->Cp[c2];
                        AW    = mesh->kx[c1]       ;
                        AE    = mesh->kx[c1+1]     ;
                        AS    = mesh->kz[c3]       ;
                        AN    = mesh->kz[c3+nxvz]  ;

                        // Contribution to S dof
                        if ( l>0 ) {
                            if (mesh->BCT_exp.type[c0-(ncx+2)] != 30 ) {
                                val   = -theta*AS*one_dz_dz;
                                AddCoeffThermal( J, A, eqn, eqn_t[c2-ncx],   &nnzc,  val );
                            }
                        }

                        // Periodic on the E side
                        if (mesh->BCT_exp.type[c0+1]==-2 && k==ncx-1 ) {
                            val = -theta*AE*one_dx_dx;
                            AddCoeffThermal( J, A, eqn, eqn_t[c2-(ncx-1)],   &nnzc,  val );
                        }

                        // Contribution to W dof
                        if (k>0) {
                            if (mesh->BCT_exp.type[c0-1] != 30 ) {
                                val   = -theta*AW*one_dx_dx;
                                AddCoeffThermal( J, A, eqn, eqn_t[c2-1],   &nnzc,  val );
                            }
                        }

                        // ----------------- Center point ----------------- //
                        val = transient*rhoCp/dt;

                        // Average conductivity for surface values (avoid zero conductivity)
                        ks = 0.25*(AE+AW+AN+AS);
                        if ( ks < mink ) {
                            printf("ACHTUNG: interloplated surface conductivity is set lower cutoff value!!\n");
                            ks = mink;
                        }

                        // Flux from sides: SOUTH
                        if ( l>0 ) {

                            if (mesh->BCT_exp.type[c0-(ncx+2)] != 30 ) {
                                val +=  AS*one_dz_dz;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCT_exp.type[c0] == 1 ) {
                                    val      += 2.0*ks*one_dz_dz;
                                    bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCT_exp.val[c0];
                                }
                            }
                        }

                        // Flux from sides: WEST
                        if ( (k>0) || (k==0 && mesh->BCT_exp.type[c0-1]==-2) ) {
                            if (mesh->BCT_exp.type[c0-1] != 30 ) {
                                val +=  AW*one_dx_dx;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCT_exp.type[c0] == 1 ) {
                                    val      += 2.0*ks*one_dx_dx;
                                    bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCT_exp.val[c0];
                                }
                            }
                        }

                        // Flux from sides: EAST
                        if ( (k<ncx-1) || (k==ncx-1 && mesh->BCT_exp.type[c0+1]==-2) ) {
                            if (mesh->BCT_exp.type[c0+1] != 30 ) {
                                val +=  AE*one_dx_dx;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCT_exp.type[c0] == 1 ) {
                                    val      += 2.0*ks*one_dx_dx;
                                    bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCT_exp.val[c0];
                                }
                            }
                        }

                        // Flux from sides: NORTH
                        if (l<ncz-1)  {
                            if ( mesh->BCT_exp.type[c0+(ncx+2)] != 30 ) {
                                val +=  AN*one_dz_dz;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCT_exp.type[c0] == 1 ) {
                                    val      += 2.0*ks*one_dz_dz;
                                    bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCT_exp.val[c0];
                                }
                            }
                        }

                        // -------------- Dirichlet contibutions -------------- //

                        if (l==0 && mesh->BCT_exp.type[c0-(ncx+2)]==1) { // SOUTH
                            val      += 2.0*AS*one_dz_dz;
                            bbc[eqn] += 2.0*AS*one_dz_dz * mesh->BCT_exp.val[c0-(ncx+2)];
                        }

                        // Bottom flux
                        if (l==0 && mesh->BCT_exp.type[c0-(ncx+2)]==0) { // SOUTH
                            bbc[eqn] += -1.0/mesh->dx * mesh->BCT_exp.val[c0-(ncx+2)];
                        }

                        if (k==0 && mesh->BCT_exp.type[c0-1]==1) {  // WEST
                            val      += 2.0*AW*one_dx_dx;
                            bbc[eqn] += 2.0*AW*one_dx_dx * mesh->BCT_exp.val[c0-1];
                        }

                        if (k==ncx-1 && mesh->BCT_exp.type[c0+1]==1) { // EAST
                            val      += 2.0*AE*one_dx_dx;
                            bbc[eqn] += 2.0*AE*one_dx_dx * mesh->BCT_exp.val[c0+1];
                        }

                        if (l==ncz-1 && mesh->BCT_exp.type[c0+(ncx+2)]==1)  { // NORTH
                            val      += 2.0*AN*one_dz_dz;
                            bbc[eqn] += 2.0*AN*one_dz_dz * mesh->BCT_exp.val[c0+(ncx+2)];
                        }

                        // Zero flux in radial coordinates
                        if (k==0 && model.polar==1 && mesh->BCT_exp.type[c0+(ncx+2)]!=30 ) {
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

                        AddCoeffThermal( J, A, eqn, eqn_t[c2],   &nnzc, val );

                        // -------------- Dirichlet contibutions -------------- //

                        // Contribution to E dof
                        if (k<ncx-1) {
                            if (mesh->BCT_exp.type[c0+1] != 30 ) {
                                val   = -theta*AE*one_dx_dx;
                                AddCoeffThermal( J, A, eqn, eqn_t[c2+1],   &nnzc,  val );
                            }
                        }

                        // Periodic on W side
                        if (mesh->BCT_exp.type[c0-1]==-2 && k==0 ) {
                            val = -theta*AW*one_dx_dx;
                            AddCoeffThermal( J, A, eqn, eqn_t[c2+(ncx-1)],   &nnzc,  val );
                        }

                        // Contribution to N dof
                        if (l<ncz-1)  {
                            if (mesh->BCT_exp.type[c0+(ncx+2)] != 30 ) {
                                val   = -theta*AN*one_dz_dz;

                                // Zero flux in radial coordinates
                                if (k==0 && model.polar==1) {
                                    val +=  (aW-1.0)*AW*one_dx_dx;
                                }

                                if (k==ncx-1 && model.polar==1) {
                                    val +=  (aE-1.0)*AE*one_dx_dx;
                                }

                                AddCoeffThermal( J, A, eqn, eqn_t[c2+ncx],   &nnzc,  val );
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
            printf("Energy balance --> System size: ndof = %d, nnz = %d\n", neq, nnzc);
        }

        // ------------------------------------------- SOLVER ------------------------------------------- //

        // Compute transpose of A matrix (used in polar mode only)
        if ( it == 0 ) At = TransposeA( &c, A, Ic, J, neq, nnzc );

        // Factor matrix only at the first step
        if ( it == 0 ) Afact = FactorEnergyCHOLMOD( &c, At, A, Ic, J, neq, nnzc, model.polar );

        // Solve
        ArrayPlusScalarArray( b, 1.0, bbc, neq );
        SolveEnergyCHOLMOD( &c, At, Afact, x, b, neq, nnzc, model.polar );

        // ------------------------------------------- SOLVER ------------------------------------------- //

        // Extract temperature T from solution vector x, compute temperature increments dT and integrate heat increments dUt
        MinMaxArrayTag( mesh->T, scaling.T, (mesh->Nx-1)*(mesh->Nz-1), "T", mesh->BCt.type );
        dUt          = 0.0;
        double dT    = 0.0;
#pragma omp parallel for shared( mesh ) private( c2, c0, k, l, eqn, dT ) firstprivate( ncx, ncz, model, scaling, x, eqn_t ) reduction (+:dUt)
        // LOOP ON THE GRID TO CALCULATE FD COEFFICIENTS
        for( l=0; l<ncz; l++) {
            for( k=0; k<ncx; k++) {
                c0  = (k+1) + (l+1)*(ncx+2);
                c2  = k   + l*ncx;
                eqn = eqn_t[c2];
                if ( mesh->BCT_exp.type[c0] != 30 ) {
                    dT           = x[eqn] - mesh->T[c2];
                    mesh->T[c2]  = x[eqn];
                    dUt         += mesh->rho_n[c2]*mesh->Cp[c2]*dT;
                    if (mesh->T[c2] < 0.0) {
                        printf("Negative temperature --- Are you crazy! (EnergyDirectSolve)\n");
                        exit(1);
                    }
                }
                else {
                    // Outside the free surface
                    mesh->T[c2]  = zeroC/scaling.T;
                }
            }
        }
        mesh->Uthermal += dUt*model.dx*model.dz;
        MinMaxArrayTag( mesh->T, scaling.T, (mesh->Nx-1)*(mesh->Nz-1), "T", mesh->BCt.type );
#ifndef _VG_
        // Write matrix in a HDF5 file for debugging purposes
        if ( model.writer_debug == 1 ) {

            char *filename;
            asprintf( &filename, "EnergyMatrix.gzip.h5" );

            // Fill in DD data structure
            int cc;
            OutputSparseMatrix OutputDD;
            OutputDD.V         = A;
            OutputDD.Ic        = Ic;
            OutputDD.J         = J;
            OutputDD.b         = b;
            OutputDD.params[0] = nx;
            OutputDD.params[1] = nz;
            OutputDD.params[2] = mesh->dx;
            OutputDD.params[3] = mesh->dz;
            OutputDD.eqn_p     = DoodzMalloc(ncx*ncz*sizeof(int));

            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    OutputDD.eqn_p[cc]=eqn_t[cc];
                }
            }

            // Send data to file
            CreateOutputHDF5( filename );
            AddGroupToHDF5( filename, "model" );
            AddGroupToHDF5( filename, "matrix" );
            AddGroupToHDF5( filename, "numbering" );
            AddGroupToHDF5( filename, "fields" );
            AddFieldToGroup( filename, "model", "params" , 'd', 4, OutputDD.params,  1 );
            AddFieldToGroup( filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDD.eqn_p,  1 );
            AddFieldToGroup( filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup( filename, "matrix", "I" , 'i', neq+1, OutputDD.Ic,  1 );
            AddFieldToGroup( filename, "matrix", "J" , 'i', nnzc,  OutputDD.J,  1 );
            AddFieldToGroup( filename, "matrix", "V" , 'd', nnzc,  OutputDD.V,  1 );
            AddFieldToGroup( filename, "matrix", "rhs" , 'd', neq, OutputDD.b,  1 );
            DoodzFree( OutputDD.eqn_p );
            free(filename);
        }
#endif
        // Free arrays
        DoodzFree(x);
        DoodzFree(b);
        DoodzFree(Ic);
        DoodzFree(J);
        DoodzFree(A);
        DoodzFree(Hs);
        DoodzFree(Ha);
    }
    // Free factors
    cholmod_free_factor ( &Afact, &c) ;
    cs_spfree(At);
    cholmod_finish( &c ) ;
    DoodzFree(bbc);
    DoodzFree(eqn_t);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ThermalSteps( grid *mesh, params model, double *rhs_t, markers *particles, double total_time, scale scaling ) {

    double dt = total_time/10.0;
    int nit = floor(total_time/dt);

    printf( "Total thermal time = %2.2e s\n", total_time*scaling.t );
    printf( "Number of thermal steps = %d\n", nit );

    // Run thermal diffusion timesteps to generate initial temperature field
    EnergyDirectSolve( mesh, model, rhs_t, particles, dt, 0, 0, scaling, nit );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetThermalPert( grid* mesh, params model, scale scaling ) {

    double x, z;
    int    k, l, c2, nx, nz, ncx, ncz;
    double pert = model.therm_perturb_dT, x0 = model.therm_perturb_x0, z0 = model.therm_perturb_z0;
    double rad_x = model.therm_perturb_rad_x, rad_z = model.therm_perturb_rad_z;

    nx   = mesh->Nx;
    nz   = mesh->Nz;
    ncx  = nx-1;
    ncz  = nz-1;

    printf("Setting thermal perturbation of %3.1lf C...\n", pert*scaling.T);

    //#pragma omp parallel for shared( mesh ) private( x, z, l, k, c2 ) firstprivate( rad, pert, x0, z0, ncx)
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            c2 = k + l*ncx;
            x  = mesh->xc_coord[k];
            z  = mesh->zc_coord[l];
            if ( pow((x-x0)/rad_x,2.0) + pow((z-z0)/rad_z,2.0) < 1.0 ) {
                mesh->T[c2]  += pert;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/