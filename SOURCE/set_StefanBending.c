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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {

    int k;
    double TopoLevel = 2.e-1/scaling.L; // sets zero initial topography
    double a = 100e3/scaling.L, b = 10e3/scaling.L, Rad=6370e3/scaling.L;
    double maxAngle = 18.0/2.0*M_PI/180.0;
    double maxX     = model.xmax, X, tet;
    double amp      = Rad*sin(M_PI/2.0) - Rad*sin(maxAngle+M_PI/2);

    for ( k=0; k<topo_chain->Nb_part; k++ ) {

        if (model.polar==0) {
            topo_chain->z[k]     = Rad;//b / (1 + pow(topo_chain->x[k]/a,2.0));
        }
        if (model.polar==1) {
            // see PolarCoordinatesStuff.py
            topo_chain->z[k]     =  sqrt((Rad - topo_chain->x[k])*(Rad + topo_chain->x[k]));
        }
//        printf("topo = %2.2e tet = %2.2e\n", topo_chain->z[k]*scaling.L/1e3, tet*180/M_PI);
        topo_chain->phase[k] = 0;
    }

    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    int np;
    FILE *read;
    int s1, s2;

    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    double T_init = (model.user2 + zeroC)/scaling.T;
    int    cyl = (int)model.user1;
    double X, Z;

//    printf("%2.2e %2.2e %2.2e\n", Ttop, Tbot, Tgrad );
    double a = 100e3/scaling.L, b = 10e3/scaling.L, Rad=6370e3/scaling.L;
    double maxAngle = 18.0/2.0*M_PI/180.0;
    double maxX     = model.xmax;
    double amp      = Rad*sin(M_PI/2.0) - Rad*sin(maxAngle+M_PI/2);
    double zMoho = Rad + model.user0/scaling.L;
    double zLAB  = Rad + model.user1/scaling.L;

    double Ttop  = 293.0/(scaling.T);
    double Tbot  = (model.user3 + zeroC)/scaling.T;
    double Tgrad = (Ttop-Tbot)/ (Lz - (model.zmax-Rad));
    double x0, z0, x_ell, z_ell, angle = 35.*M_PI/180;
    double a_ell = 2.0*model.user5/scaling.L, b_ell = 0.5*model.user5/scaling.L;

    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {

        // Standard initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;                   // set initial particle velocity (unused)
        particles->phase[np] = 0;                                               // same phase number everywhere
        particles->d[np]     = 0;                                               // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->rho[np]   = 0;
//        particles->T[np]     = T_init;
        particles->T[np]     = Ttop + Tgrad*(particles->z[np] - Rad);


        // ------------------------- //

//        // Mantle
        if ( model.polar==0 ) {
            if (particles->z[np] < zMoho) {// + b / (1.0 + pow(particles->x[np]/a, 2.0)) ) {
                particles->phase[np] = 1;
            }
            if (particles->z[np] < zLAB) {// + b / (1.0 + pow(particles->x[np]/a, 2.0)) ) {
                particles->phase[np] = 2;
            }
        }
        if ( model.polar==1 ) {
            Z = sqrt((zMoho - particles->x[np])*(zMoho + particles->x[np]));
            if ( particles->z[np] <  Z  ) particles->phase[np] = 1;
            Z = sqrt((zLAB  - particles->x[np])*(zLAB  + particles->x[np]));
            if ( particles->z[np] <  Z  ) particles->phase[np] = 2;

        }

        Z = Rad - 55e3/scaling.L;
        X = 20e3/scaling.L;
//        double rad = model.user5/scaling.L;
//        if ( ( pow(particles->x[np]-X,2) + pow(particles->z[np]-Z,2) ) < rad*rad  ) {
//            particles->phase[np] = 0;
//        }
        x0 = X;
        z0 = Z;
        x_ell = (particles->x[np]-x0)*cos(angle) + (particles->z[np]-z0)*sin(angle);
        z_ell =-(particles->x[np]-x0)*sin(angle) + (particles->z[np]-z0)*cos(angle);
        if (pow(x_ell/a_ell,2.0) + pow(z_ell/b_ell,2.0) < 1.0) particles->phase[np] = 0;

        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }

        //--------------------------//
        // DENSITY
        if ( model.eqn_state > 0 ) {
            particles->rho[np] = materials->rho[particles->phase[np]] * (1 -  materials->alp[particles->phase[np]] * (T_init - materials->T0[particles->phase[np]]) );
        }
        else {
            particles->rho[np] = materials->rho[particles->phase[np]];
        }
        //--------------------------//
    }
    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T,  scaling.T, particles->Nb_part, "Tp init"  );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {

    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->EpsBG;
    double Inflow = 0.0, VzOutflow = 0.0, x, z, V, tet, r, Vx, Vz, Vt, Vr;
    int    OutflowOnSides = (int)model->user4;

    // Inflow/Outflow velocities
    double Vx_tot =  (model->xmax-model->xmin) * model->EpsBG;
    double VxOut_W = 0.0, VzOut_W = 0.0, VxIn_W = 0.5*Vx_tot, VzIn_W = 0.0;
    double VxOut_E = 0.0, VzOut_E = 0.0, VxIn_E =-0.5*Vx_tot, VzIn_E = 0.0;
    double VzOut_S = 0.0;
    double maxAngle = 18.0/2.0*M_PI/180.0, tet_W, alp_W, tet_E, alp_E;

    if (model->polar == 1) {
        tet_W  = -maxAngle;
        alp_W  = M_PI/2.0 - tet_W;
        VxIn_W = (0.5*Vx_tot) * sin(alp_W);
        VzIn_W =-(0.5*Vx_tot) * cos(alp_W);
        printf("VxIn_W = %2.2e m/s VzIn_W = %2.2e m/s \n", VxIn_W*scaling.V, VzIn_W*scaling.V);

        tet_E  = maxAngle;
        alp_E  = M_PI/2.0 - tet_E;
        VxIn_E =-(0.5*Vx_tot) * sin(alp_E);
        VzIn_E = (0.5*Vx_tot) * cos(alp_E);
        printf("VxIn_E = %2.2e m/s VzIn_E = %2.2e m/s \n", VxIn_E*scaling.V, VzIn_E*scaling.V);
    }



    printf("Vx_tot = %2.2e m/s --- VxW = %2.2e m/s --- VxE = %2.2e m/s --- EpsBG = %2.2e 1/s \n", Vx_tot*scaling.V, VxIn_W*scaling.V, VxIn_E*scaling.V, model->EpsBG*scaling.E);

    //-----------------------------------------------------//

    // Declare and allocate boundary velocity arrays
    double HLit   = -(double)model->user2/scaling.L; // scale values of lithosphere thickness
    double *VxBC_W, *VxBC_E, *VzBC_W, *VzBC_E, top_W =  model->zmin, top_E = model->zmin;
    int c_E, c_W;
    double dx = model->dx, dz = model->dz;
    VxBC_W = DoodzCalloc(mesh->Nz+1, sizeof(double));
    VxBC_E = DoodzCalloc(mesh->Nz+1, sizeof(double));
    VzBC_W = DoodzCalloc(mesh->Nz+0, sizeof(double));
    VzBC_E = DoodzCalloc(mesh->Nz+0, sizeof(double));

    // Find how many cells have inflow/outflow on each side of the model
    int nc_in_W=0, nc_out_W=0, nc_in_E=0, nc_out_E=0;

    // Find top W/E
    for (l=0; l<mesh->Nz+1; l++) {

        // Vx node index
        c_W =           0  + l*(mesh->Nx);
        c_E = (mesh->Nx-1) + l*(mesh->Nx);

        if (mesh->BCu.type[c_W] != 30  && mesh->zvx_coord[l]>top_W) top_W = mesh->zvx_coord[l];
        if (mesh->BCu.type[c_E] != 30  && mesh->zvx_coord[l]>top_E) top_E = mesh->zvx_coord[l];
    }

    printf("Top W = %2.2e --- Top E = %2.2e\n", top_W*scaling.L, top_E*scaling.L);

    // Loop through Vx nodes and count which cells are inside the mode and have inflow/outflow
    for (l=0; l<mesh->Nz+1; l++) {

            // Vx node index
            c_W =           0  + l*(mesh->Nx);
            c_E = (mesh->Nx-1) + l*(mesh->Nx);

            // West
            if (l>0 && l<mesh->Nz) {
                if ( mesh->BCu.type[c_W] != 30 ) {
                    if (mesh->zvx_coord[l] > top_W - HLit) {
                        nc_in_W++;
                    }
                    else {
                        nc_out_W++;
                    }
                }
            }

            // East
            if (l>0 && l<mesh->Nz) {
                if ( mesh->BCu.type[c_E] != 30 ) {
                    if (mesh->zvx_coord[l] > top_E - HLit) {
                        nc_in_E++;
                    }
                    else {
                        nc_out_E++;
                    }
                }
            }
    }

    // Compute outflow velocities
    if (OutflowOnSides==1) {
        VxOut_W = -VxIn_W*nc_in_W/nc_out_W;
        VxOut_E = -VxIn_E*nc_in_E/nc_out_E;
    }
    else {
        VzOut_S = -( fabs(VxIn_W*nc_in_W*dz) + fabs(VxIn_E*nc_in_E*dz) ) / (mesh->Nx-1)/dx;
        printf( "VzOut_S = %2.2e m/s\n", VzOut_S*scaling.V);
    }

    if (model->polar == 1) {
        VzOut_W =-VxOut_W * cos(alp_W);
        VxOut_W = VxOut_W * sin(alp_W);
        VzOut_E =-VxOut_E * cos(alp_E);
        VxOut_E = VxOut_E * sin(alp_E);
    }

    printf("Vxin_W = %2.2e Vzin_W = %2.2e VxOut_W = %2.2e VzOut_W = %2.2e\n", VxIn_W*scaling.V, VzIn_W*scaling.V, VxOut_W*scaling.V, VzOut_W*scaling.V);
    printf("Vxin_E = %2.2e Vzin_E = %2.2e VxOut_E = %2.2e VzOut_E = %2.2e\n", VxIn_E*scaling.V, VzIn_E*scaling.V, VxOut_E*scaling.V, VzOut_E*scaling.V);

    // Loop through Vx nodes and assert BC values
    for (l=0; l<mesh->Nz+1; l++) {

            // Vx node index
            c_W = 0 + l*(mesh->Nx);
            c_E = (mesh->Nx-1) + l*(mesh->Nx);

            // West
            if  (l>0 && l<mesh->Nz) {
                if ( mesh->BCu.type[c_W] != 30 ) {
                    if (mesh->zvx_coord[l] > top_W - HLit) {
                        VxBC_W[l] = VxIn_W;
                    }
                    else {
                        VxBC_W[l] = VxOut_W;
                    }
                }
            }

            // East
            if  (l>0 && l<mesh->Nz) {
                if ( mesh->BCu.type[c_E] != 30 ) {
                    if (mesh->zvx_coord[l] > top_E - HLit) {
                        VxBC_E[l] = VxIn_E;
                    }
                    else {
                        VxBC_E[l] = VxOut_E;
                    }
                }
            }
    }

    // Loop through Vz nodes and assert BC values
    for (l=0; l<mesh->Nz; l++) {

        // Vz node index
        c_W = 0 + l*(mesh->Nx+1);
        c_E = (mesh->Nx) + l*(mesh->Nx+1);

        // West
            if ( mesh->BCv.type[c_W] != 30 ) {
                if (mesh->zg_coord[l] > top_W - HLit) {
                    VzBC_W[l] = VzIn_W;
                }
                else {
                    VzBC_W[l] = VzOut_W;
                }
            }

        // East
            if ( mesh->BCv.type[c_E] != 30 ) {
                if (mesh->zg_coord[l] > top_E - HLit) {
                    VzBC_E[l] = VzIn_E;
                }
                else {
                    VzBC_E[l] = VzOut_E;
                }
            }
    }

    // Smooth using explicit diffusion solver
    double qzS, qzN, knum=1, dtnum=1/knum*model->dz*model->dz/2, tnum=0.001, VxN, VxS, VzN, VzS;
    int    it, nit=(int)(tnum/(dtnum));

    // Loop through Vx nodes and assert BC values
    for (it=0; it<nit; it++) {
        for (l=0; l<mesh->Nz+1; l++) {

                // Vx node index
                c_W = 0 + l*(mesh->Nx);
                c_E = (mesh->Nx-1) + l*(mesh->Nx);

                // West
                    if ( mesh->BCu.type[c_W] != 30 ) {
                        if (l==0          ) VxS = VxOut_W;
                        else VxS = VxBC_W[l-1];
                        if (l==mesh->Nz || mesh->BCu.type[c_W+mesh->Nx] == 30 ) VxN = VxIn_W;
                        else VxN = VxBC_W[l+1];
                        qzN = (VxN - VxBC_W[l  ])/model->dz;
                        qzS = (VxBC_W[l  ] - VxS)/model->dz;
                        VxBC_W[l] = VxBC_W[l] + knum*dtnum/model->dz*(qzN-qzS);
                    }

                // East
            if ( mesh->BCu.type[c_E] != 30 ) {
                        if (l==0          ) VxS = VxOut_E;
                        else VxS = VxBC_E[l-1];
                        if (l==mesh->Nz || mesh->BCu.type[c_E+mesh->Nx] == 30 ) VxN = VxIn_E;
                        else VxN = VxBC_E[l+1];
                        qzN = (VxN - VxBC_E[l  ])/model->dz;
                        qzS = (VxBC_E[l  ] - VxS)/model->dz;
                        VxBC_E[l] = VxBC_E[l] + knum*dtnum/model->dz*(qzN-qzS);
                    }
            }

        for (l=0; l<mesh->Nz; l++) {

            // Vx node index
            c_W = 0 + l*(mesh->Nx+1);
            c_E = (mesh->Nx) + l*(mesh->Nx+1);

            // West
            if ( mesh->BCv.type[c_W] != 30 ) {
                if (l==0          ) VzS = VzOut_W;
                else VzS = VzBC_W[l-1];
                if (l==mesh->Nz-1 || mesh->BCv.type[c_W+mesh->Nx+1] == 30 ) VzN = VzIn_W;
                else VzN = VzBC_W[l+1];
                qzN = (VzN - VzBC_W[l  ])/model->dz;
                qzS = (VzBC_W[l  ] - VzS)/model->dz;
                VzBC_W[l] += knum*dtnum/model->dz*(qzN-qzS);
            }

            // East
            if ( mesh->BCv.type[c_E] != 30 ) {
                if (l==0          ) VzS = VzOut_E;
                else VzS = VzBC_E[l-1];
                if (l==mesh->Nz-1 || mesh->BCv.type[c_E+mesh->Nx+1] == 30 ) VzN = VzIn_E;
                else VzN = VzBC_E[l+1];
                qzN = (VzN - VzBC_E[l  ])/model->dz;
                qzS = (VzBC_E[l  ] - VzS)/model->dz;
                VzBC_E[l] += knum*dtnum/model->dz*(qzN-qzS);
            }

//            printf("%2.2e\n", VzBC_E[l]-VzBC_W[l]);

        }

    }


    //-----------------------------------------------------//


    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;

    /* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vx on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -2: periodic in the x direction (matches the physical boundary)                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
	/* --------------------------------------------------------------------------------------------------------*/

		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;

		for (l=0; l<mesh->Nz+1; l++) {
			for (k=0; k<mesh->Nx; k++) {

				c = k + l*(mesh->Nx);

                x  = mesh->xg_coord[k];
                z  = mesh->zvx_coord[l];

                if ( model->ismechanical ==0 && model->polar==1 ) {
                    r           = sqrt(x*x + z*z);
                    tet         = atan(z/x);
                    Vt          = 2.0*model->EpsBG*r*sin(tet);
                    mesh->u_in[c] = Vt*sin(tet);
                }
                
                if ( mesh->BCu.type[c] != 30 ) {

                    // Internal points:  -1
                    mesh->BCu.type[c] = -1;
                    mesh->BCu.val[c]  =  0;

                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = VxBC_W[l];
                    }

                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = VxBC_E[l];
                    }

//                    // Matching BC nodes WEST
//                    if (k==0 ) {
//                        mesh->BCu.type[c] = 0;
//                        mesh->BCu.val[c]  = Vx;
//                        Inflow           += fabs(mesh->BCu.val[c])*model->dz;
//                    }
//
//                    // Matching BC nodes EAST
//                    if (k==mesh->Nx-1 ) {
//                        mesh->BCu.type[c] = 0;
//                        mesh->BCu.val[c]  = Vx;
//                        Inflow           += fabs(mesh->BCu.val[c])*model->dz;
//                    }

                    // Free slip SOUTH
                    if (l==0  ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }

                    // Free slip NORTH
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                }

			}
		}

	/* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vz on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -2: periodic in the x direction (does not match the physical boundary)                             */
	/* Type-10: useless point (set to zero)                                                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
	/* --------------------------------------------------------------------------------------------------------*/

		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;

		for (l=0; l<mesh->Nz; l++) {
			for (k=0; k<mesh->Nx+1; k++) {

				c  = k + l*(mesh->Nx+1);
                
                x  = mesh->xvz_coord[k];
                z  = mesh->zg_coord[l];
                
                 if ( model->ismechanical ==0 && model->polar==1 ) {
                     r           = sqrt(x*x + z*z);
                     tet         = atan(z/x);
                     Vt          = 2.0*model->EpsBG*r*sin(tet);
                     mesh->v_in[c] =-Vt*cos(tet);
                 }

                if ( mesh->BCv.type[c] != 30 ) {

                    // Internal points:  -1
                    mesh->BCv.type[c] = -1;
                    mesh->BCv.val[c]  =  0;

                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = VzOut_S;//VzOutflow;//mesh->zg_coord[l] * model->EpsBG;
                    }

                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;//mesh->zg_coord[l] * model->EpsBG;
                    }

                    // Non-matching boundary WEST
                    if ( (k==0) ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   VzBC_W[l];
                    }

                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx) ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   VzBC_E[l];
                    }
                }

			}
		}


    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* Type 31: surface pressure (Dirichlet)                                                                   */
    /* --------------------------------------------------------------------------------------------------------*/


        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;

        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {

                c  = k + l*(NCX);

                if (mesh->BCt.type[c] != 30) {

                    // Internal points:  -1
                    mesh->BCp.type[c] = -1;
                    mesh->BCp.val[c]  =  0;
                }
            }
        }

	/* -------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for T on all grid levels                                                                   */
	/* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)      */
	/* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)             */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                   */
	/* Type -1: not a BC point (tag for inner points)                                                         */
    /* Type 30: not calculated (part of the "air")                                                            */
	/* -------------------------------------------------------------------------------------------------------*/

    double Ttop = 293.15/scaling.T;
    double Tbot= (model->user3 + zeroC)/scaling.T;
    double Tleft, Tright;

    double k_crazy = 1000*materials->k[2];

    if (model->step == 0) {
        materials->k_eff[2] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }

    else {
        materials->k_eff[2] = materials->k[2];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }

    printf("Ttop=%2.2e Tbot=%2.2e\n", Ttop*scaling.T, Tbot*scaling.T);


		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;

		for (l=0; l<mesh->Nz-1; l++) {
			for (k=0; k<mesh->Nx-1; k++) {

				c = k + l*(NCX);

                if ( mesh->BCt.type[c] != 30 ) {

                    // WEST
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typW[l] = 0;
                        mesh->BCt.valW[l] = 0;
                    }

                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = 0;
                    }

                    // SOUTH
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typS[k] = 1;
                        mesh->BCt.valS[k] = Tbot;
                    }

                    // NORTH
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typN[k] = 1;
                        mesh->BCt.valN[k] = Ttop;
                    }

                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = Ttop;
                        }
                    }


                }

			}
		}

    DoodzFree( VxBC_W );
    DoodzFree( VxBC_E );
    DoodzFree( VzBC_W );
    DoodzFree( VzBC_E );

	printf("Velocity and pressure were initialised\n");
	printf("Boundary conditions were set up\n");

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
