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

void EnergyDirectSolve( grid *mesh, params model, double *rhoE, double *drhoE, double *rhs_t, double *rhoE0, markers *particles, double dt, int shear_heating, int adiab_heating, scale scaling, int nit ) {

    // Energy matrix
    double  *A;
    double  *b, *bbc, *x;      // Coefs, rhs, sol
    int *Ic, *J;      // lines, colums
    int k, l, c1, c2, c3, nx, nz, nxvz, nzvx, ncx, ncz;
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

        // Increment time
//        *time += dt;

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
#pragma omp parallel for shared ( b, Hs, Ha, eqn_t, mesh ) private ( l, k, c1, c2, c3, eqn, rhoCp, dexx_el, dexx_th, dexx_tot, dezz_el, dezz_th, dezz_tot, deyy_el, deyy_th, deyy_tot, dexz_el, dexz_th, dexz_tot, Wth, Wel, Wtot, syyd, eyyd, eyyd_el, eyyd_diss  ) firstprivate ( nx, ncx, nxvz, shear_heating, adiab_heating, model, transient, dt, Hr, diss_limit, scaling ) reduction( +:dUe, dW )
        for( c2=0; c2<ncz*ncx; c2++) {
                k   = mesh->kp[c2];
                l   = mesh->lp[c2];
                c1  = k + l*nx;
                c3  = k + l*nxvz;
                eqn = eqn_t[c2];

                if ( mesh->BCt.type[c2] != 30 ) {

                    // Contribution from transient
                    rhoCp   = mesh->rho_n[c2]*mesh->Cv[c2];
                    b[eqn]  = transient * rhoCp * mesh->T[c2]/ dt;

                    // Contribution from radiogenic sources
                    b[eqn] += Hr*mesh->Qr[c2];

                    // Contribution from dissipation
                    if ( shear_heating == 1 ) {
                        b[eqn] += mesh->Wdiss[c2];
                        Hs[c2]  = mesh->Wdiss[c2];
                    }

                    // Contribution from adiabatic heat
                    if ( adiab_heating == 1 ) Ha[c2]  = mesh->T[c2]*mesh->alp[c2]*0.5*(mesh->v_in[c3+1]+mesh->v_in[c3+1+nxvz])*model.gz*mesh->rho_n[c2];
                    if ( adiab_heating == 2 ) Ha[c2]  = mesh->T[c2]*mesh->alp[c2]*(mesh->p_in[c2] - mesh->p0_n[c2])/model.dt; //mesh->T[c2]*mesh->alp[c2]*(mesh->dp[c2])/model.dt;
                    if ( adiab_heating  > 0 ) b[eqn] += Ha[c2];

                }

        }

        MinMaxArrayTag( mesh->T, scaling.T, ncx*ncz, "T0", mesh->BCt.type );
        MinMaxArrayTag( Ha, scaling.S*scaling.E, ncx*ncz, "Ha", mesh->BCt.type );
        MinMaxArrayTag( Hs, scaling.S*scaling.E, ncx*ncz, "Hs", mesh->BCt.type );
        mesh->Ue += dUe;
        mesh->W  += dW;

        //----------------------------------------------------//
        if ( it == 0 ) {

            printf("Assembling Energy matrix... with %d discrete equations\n", neq);

            // LOOP ON THE GRID TO CALCULATE FD COEFFICIENTS
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {

                    c1  = k   + (l+1)*nx;
                    c2  = k   + l*ncx;
                    c3  = k   + l*nxvz+1;

                    if (mesh->BCt.type[c2] != 30 ) {

                        // Central coefficient
                        eqn     = eqn_t[c2];
                        Ic[eqn] = nnzc;

                        rhoCp = mesh->rho_n[c2]*mesh->Cv[c2];
                        AW    = mesh->kx[c1]       ;
                        AE    = mesh->kx[c1+1]     ;
                        AS    = mesh->kz[c3]       ;
                        AN    = mesh->kz[c3+nxvz]  ;

                        // Contribution to S dof
                        if ( l>0 ) {
                            if (mesh->BCt.type[c2-ncx] != 30 ) {
                                val   = -theta*AS*one_dz_dz;
                                AddCoeffThermal( J, A, eqn, eqn_t[c2-ncx],   &nnzc,  val );
                            }
                        }

                        // Periodic on the right side
                        if (mesh->BCt.type[c2]==-2 && k==ncx-1 ) {
                            val = -theta*AE*one_dx_dx;
                            AddCoeffThermal( J, A, eqn, eqn_t[c2-(ncx-1)],   &nnzc,  val );
                        }

                        // Contribution to W dof
                        if (k>0) {
                            if (mesh->BCt.type[c2-1] != 30 ) {
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

                            if (mesh->BCt.type[c2-ncx] != 30 ) {
                                val +=  AS*one_dz_dz;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCt.type[c2] == 1 ) {
                                    val      += 2.0*ks*one_dz_dz;
                                    bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCt.val[c2];
                                }
                            }
                        }

                        // Flux from sides: WEST
                        if ( (k>0) || (k==0 && mesh->BCt.type[c2]==-2) ) {
                            if (mesh->BCt.type[c2-1] != 30 ) {
                                val +=  AW*one_dx_dx;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCt.type[c2] == 1 ) {
                                    val      += 2.0*ks*one_dx_dx;
                                    bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCt.val[c2];
                                }
                            }
                        }

                        // Flux from sides: EAST
                        if ( (k<ncx-1) || (k==ncx-1 && mesh->BCt.type[c2]==-2) ) {
                            if (mesh->BCt.type[c2+1] != 30 ) {
                                val +=  AE*one_dx_dx;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCt.type[c2] == 1 ) {
                                    val      += 2.0*ks*one_dx_dx;
                                    bbc[eqn] += 2.0*ks*one_dx_dx * mesh->BCt.val[c2];
                                }
                            }
                        }

                        // Flux from sides: NORTH
                        if (l<ncz-1)  {
                            if ( mesh->BCt.type[c2+ncx] != 30 ) {
                                val +=  AN*one_dz_dz;
                            }
                            // Dirichlet contribution (free surface)
                            else {
                                if (mesh->BCt.type[c2] == 1 ) {
                                    val      += 2.0*ks*one_dz_dz;
                                    bbc[eqn] += 2.0*ks*one_dz_dz * mesh->BCt.val[c2];
                                }
                            }
                        }

                        // -------------- Dirichlet contibutions -------------- //

                        if (l==0 && mesh->BCt.typS[k]==1) {
                            val      += 2.0*AS*one_dz_dz;
                            bbc[eqn] += 2.0*AS*one_dz_dz * mesh->BCt.valS[k];
                        }

                        if (l==0 && mesh->BCt.typS[k]==2) {
                            bbc[eqn] += -1.0/mesh->dx * mesh->BCt.valS[k];
                        }

                        if (k==0 && mesh->BCt.typW[l]==1) {
                            val      += 2.0*AW*one_dx_dx;
                            bbc[eqn] += 2.0*AW*one_dx_dx * mesh->BCt.valW[l];
                        }

                        if (k==ncx-1 && mesh->BCt.typE[l]==1) {
                            val      += 2.0*AE*one_dx_dx;
                            bbc[eqn] += 2.0*AE*one_dx_dx * mesh->BCt.valE[l];
                        }

                        if (l==ncz-1 && mesh->BCt.typN[k]==1)  {
                            val      += 2.0*AN*one_dz_dz;
                            bbc[eqn] += 2.0*AN*one_dz_dz * mesh->BCt.valN[k];
                        }

                        // Zero flux in radial coordinates
                        if (k==0 && model.polar==1 && mesh->BCt.type[c2+ncx] != 30 ) {
//                            xn   = mesh->xg_coord[k];
//                            zn   = mesh->zc_coord[l];
//                            tet  = -atan(zn/xn) / (2.0*model.dx/model.dz);
//                            if (tet>0.0) yW  = atan(tet) * model.dx;
//                            if (tet<0.0) yW  =-atan(tet) * model.dx;
//                            aW   = yW/model.dz;


                            xn   = mesh->xg_coord[k]-model.dx/2.0;
                            zn   = mesh->zc_coord[l];
                            tet  = atan(zn/xn);
                            if (tet<0.0) tet  *= -1;
                            yW   = atan(M_PI/2.0-tet) * model.dx;
                            aW   = yW/model.dz;
                            val +=  (1.0-aW)*AW*one_dx_dx;
                        }

                        if (k==ncx-1 && model.polar==1) {
//                            xn   = mesh->xg_coord[k+1];
//                            zn   = mesh->zc_coord[l];
//                            tet  = -atan(zn/xn) / (2.0*model.dx/model.dz);
//                            if (tet>0.0) yE  = atan(tet) * model.dx;
//                            if (tet<0.0) yE  =-atan(tet) * model.dx;
//                            aE   = yE/model.dz;

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
                            if (mesh->BCt.type[c2+1] != 30 ) {
                                val   = -theta*AE*one_dx_dx;
                                AddCoeffThermal( J, A, eqn, eqn_t[c2+1],   &nnzc,  val );
                            }
                        }

                        // Periodic on the left side
                        if (mesh->BCt.type[c2]==-2 && k==0 ) {
                            val = -theta*AW*one_dx_dx;
                            AddCoeffThermal( J, A, eqn, eqn_t[c2+(ncx-1)],   &nnzc,  val );
                        }

                        // Contribution to N dof
                        if (l<ncz-1)  {
                            if (mesh->BCt.type[c2+ncx] != 30 ) {
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
        zero_celsius = zeroC;
#pragma omp parallel for shared( mesh, scaling, x ) private( c2, eqn ) firstprivate( ncx, ncz, zero_celsius, model ) reduction (+:dUt)
        for( c2=0; c2<ncx*ncz; c2++) {
            eqn = eqn_t[c2];
            if ( mesh->BCt.type[c2] != 30 ) {
                mesh->dT[c2] = 1.0*(x[eqn] -  mesh->T[c2]);
                mesh->T[c2]  = mesh->T[c2] + mesh->dT[c2];
                dUt         += mesh->rho_n[c2]*mesh->Cv[c2]*mesh->dT[c2];
                if (mesh->T[c2] < 0.0) {
                    printf("Negative temperature --- Are you crazy! (EnergyDirectSolve)\n");
                    exit(1);
                }
            }
            else {
                // Outside the free surface
                mesh->dT[c2] = 0.0;
                mesh->T[c2]  = zero_celsius/scaling.T;
            }
        }
        mesh->Ut += dUt*model.dx*model.dz;
        MinMaxArrayTag( mesh->T, scaling.T, (mesh->Nx-1)*(mesh->Nz-1), "T", mesh->BCt.type );
#ifndef _VG_
        // Write matrix in a HDF5 file for debugging purposes
        if ( model.write_debug == 1 ) {

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
            create_output_hdf5( filename );
            AddGroup_to_hdf5( filename, "model" );
            AddGroup_to_hdf5( filename, "matrix" );
            AddGroup_to_hdf5( filename, "numbering" );
            AddGroup_to_hdf5( filename, "fields" );
            AddFieldToGroup_generic( _TRUE_, filename, "model", "params" , 'd', 4, OutputDD.params,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDD.eqn_p,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "I" , 'i', neq+1, OutputDD.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "J" , 'i', nnzc,  OutputDD.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "V" , 'd', nnzc,  OutputDD.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs" , 'd', neq, OutputDD.b,  1 );

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

void ThermalSteps( grid *mesh, params model, double *rhoE, double *drhoE, double *rhs_t, double *rhoE0, markers *particles, double total_time, scale scaling ) {

    double dt = total_time/10.0;
    int nit = floor(total_time/dt);

    printf( "Total thermal time = %2.2e s\n", total_time*scaling.t );
    printf( "Number of thermal steps = %d\n", nit );

    // Run thermal diffusion timesteps to generate initial temperature field
    EnergyDirectSolve( mesh, model, mesh->T, mesh->dT, rhs_t, mesh->T, particles, dt, 0, 0, scaling, nit );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void AdvectHeat ( grid* mesh, params *model, scale scaling ) {
//
//    DoodzFP *VxC, *VzC;
//    VxC = DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1), sizeof(DoodzFP) );
//    VzC = DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1), sizeof(DoodzFP) );
//
//    VelocitiesOnCenters(  mesh->u_in, mesh->v_in, VxC, VzC, mesh->Nx, mesh->Nz, scaling );
//    FirstOrderUpwindAdvection( VxC, VzC, mesh->T, mesh->T, mesh, mesh->Nx-1, mesh->Nz-1, *model, scaling, 1 );
//
//    DoodzFree( VxC );
//    DoodzFree( VzC );
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Energies( grid* mesh, params model, scale scaling ) {

    int noisy=1;

    if ( noisy==1 ) {
        printf("Thermal energy = %2.8e\n",  mesh->Ut*scaling.S*scaling.L*scaling.L );
        printf("Mechnical work = %2.8e\n",  mesh->W *scaling.S*scaling.L*scaling.L );
        if (model.iselastic == 1) printf("Elastic energy = %2.8e\n", mesh->Ue*scaling.S*scaling.L*scaling.L );
        if (model.iselastic == 1) printf("Remainder      = %2.8e\n", (mesh->W - mesh->Ut - mesh->Ue)/mesh->W );
        else {
            if (noisy==1) printf("Difference     = %2.8e\n", (mesh->W - mesh->Ut)/mesh->W );
        }
    }

    mesh->Uthermal[model.step-1] = mesh->Ut*scaling.S*scaling.L*scaling.L;
    mesh->Uelastic[model.step-1] = mesh->Ue*scaling.S*scaling.L*scaling.L;
    mesh->Work[model.step-1]     = mesh->W*scaling.S*scaling.L*scaling.L;
    mesh->Time[model.step-1]     = model.time*scaling.t;
    mesh->Short[model.step-1]    = (model.L0-(model.xmax - model.xmin))/model.L0 * 100.0;

    if ( model.step == model.Nt ) {

        create_output_hdf5( "Energies.gzip.h5" );
        AddGroup_to_hdf5( "Energies.gzip.h5", "Fields" );
        AddFieldToGroup_generic( _TRUE_,  "Energies.gzip.h5", "Fields", "Time" ,     'd', model.Nt,  mesh->Time,      1 );
        AddFieldToGroup_generic( _TRUE_,  "Energies.gzip.h5", "Fields", "Short" ,    'd', model.Nt,  mesh->Short,     1 );
        AddFieldToGroup_generic( _TRUE_,  "Energies.gzip.h5", "Fields", "Work" ,     'd', model.Nt,  mesh->Work,      1 );
        AddFieldToGroup_generic( _TRUE_,  "Energies.gzip.h5", "Fields", "Uthermal" , 'd', model.Nt,  mesh->Uthermal,  1 );
        AddFieldToGroup_generic( _TRUE_,  "Energies.gzip.h5", "Fields", "Uelastic" , 'd', model.Nt,  mesh->Uelastic,  1 );

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetThermalPert( grid* mesh, params model, scale scaling ) {

    double x, z;
    int    k, l, c2, nx, nz, ncx, ncz;
    double pert = model.therm_pert_dT, x0 = model.therm_pert_x0, z0 = model.therm_pert_z0, rad = model.therm_pert_rad;

    nx   = mesh->Nx;
    nz   = mesh->Nz;
    ncx  = nx-1;
    ncz  = nz-1;

    printf("Setting thermal perturbation of %3.1lf C...\n", pert*scaling.T);

    #pragma omp parallel for shared( mesh ) private( x, z, l, k, c2 ) firstprivate( rad, pert, x0, z0, ncx)
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            c2 = k + l*ncx;
            x  = mesh->xc_coord[k];
            z  = mesh->zc_coord[l];
            if ( pow((x-x0),2.0) + pow((z-z0),2.0) < rad*rad ) {
                mesh->T[c2]  += pert;
            }
        }
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
