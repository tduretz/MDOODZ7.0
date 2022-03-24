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
    double TopoLevel = 0.0e3/scaling.L; // sets zero initial topography

    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {    int np;
    FILE *read;
    int s1, s2;
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    double T_init = (model.user0 + zeroC)/scaling.T;
    double radius = 0.1/scaling.L;
    double X, Z, xc = 0.0, zc = 0.0;
    int pos = 1;
    int il;
    double spacing = 0.1/scaling.L;
    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;                   // set initial particle velocity (unused)
        particles->phase[np] = 0;                                               // same phase number everywhere
        particles->d[np]     = 0;                                            // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->rho[np]   = 0;
        particles->T[np]     = T_init;

    
        // ------------------------- //
//        if (pos==0) {
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // DRAW INCLUSION
//        if (X*X + Z*Z < radius*radius) {
//            particles->phase[np] = 1;
//        }
//        }
//        else {
//        X = particles->x[np]-model.xmax;
//        Z = particles->z[np]-zc;
//        // DRAW INCLUSION
//        if (X*X + Z*Z < radius*radius) {
//            particles->phase[np] = 1;
//        }
//
//        // DRAW INCLUSION
//        X = particles->x[np]-model.xmin;
//        Z = particles->z[np]-zc;
//        if (X*X + Z*Z < radius*radius) {
//            particles->phase[np] = 1;
//        }
//        }
        //=================================================================
        // Layering (passive markers) - Matrice - phase 0 et 3 => Granulite
        if (particles->phase[np] == 0) {

            for( il=0; il<= 100; il=il+2 ) {
                if ((particles->z[np]-Lz/2.0-spacing/2.0)>=-((il+1)*spacing) && (particles->z[np]-Lz/2.0-spacing/2.0)<-(il*spacing)) particles->phase[np] = 4;
            }
            for( il=0; il<= 100; il=il+2 ) {
                if (particles->phase[np]==0 && (particles->x[np]-Lx/2.0-spacing/2.0)>=-((il+1)*spacing) && (particles->x[np]-Lx/2.0-spacing/2.0)<-(il*spacing)) particles->phase[np] = 5;
                if (particles->phase[np]==4 && (particles->x[np]-Lx/2.0-spacing/2.0)>=-((il+1)*spacing) && (particles->x[np]-Lx/2.0-spacing/2.0)<-(il*spacing)) particles->phase[np] = 0;
            }
            if (particles->phase[np]==5 || particles->phase[np]==4) particles->phase[np] = 3;
        }
        //==================================================================
                // 1ere elliptical inclusion
        double rad=0.1/scaling.L;
        double X,Xn,Z,Zn, xc=0.0, zc=0.0, la= 1.00*rad, sa = 1.00*rad, theta=(30.0)*M_PI/180.0;
        xc=-0.0*scaling.L;
        zc=-0.0*scaling.L;
        
        X = particles->x[np]-xc;
        Z = particles->z[np]-zc;
        
        Xn = X*cos(theta) - Z*sin(theta);
        Zn = X*sin(theta) + Z*cos(theta);
        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
        
        // ------------------------------------------------
//        // 1ere elliptical inclusion
//        double rad=0.25/scaling.L;
//        double X,Xn,Z,Zn, xc=model.xmax, zc=model.zmin, la= 1.00*rad, sa = 1.00*rad, theta=(90.0)*M_PI/180.0;
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // elliptical inclusion
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//
//        // 2eme elliptical inclusion
//        xc=model.xmin;
//        zc=model.zmin;
//        la= 1.0*rad;
//        sa = 1.0*rad;
//        //theta=(90.0)*M_PI/180.0;
//
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//
//        // 3eme elliptical inclusion
//        xc=+0.4/scaling.L;
//        zc=+0.2/scaling.L;
//        la= 1.0*rad;
//        sa = 9.0*rad;
//        //theta=(0.0)*M_PI/180.0;
//        
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // elliptical inclusion
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//        
//        // 4eme elliptical inclusion
//        xc=-0.4/scaling.L;
//        zc=-0.2/scaling.L;
//        la= 1.0*rad;
//        sa = 9.0*rad;
//        //theta=(0.0)*M_PI/180.0;
//        //
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // elliptical inclusion
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//        //
//        // 5eme elliptical inclusion
//        xc=-0.01/scaling.L;
//        zc=+0.35/scaling.L;
//        la= 1.0*rad;
//        sa = 15.0*rad;
//        //theta=(0.0)*M_PI/180.0;
//        
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // elliptical inclusion
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//        
//        // 6eme elliptical inclusion
//        xc=+0.01/scaling.L;
//        zc=-0.35/scaling.L;
//        la= 1.0*rad;
//        sa = 15.0*rad;
//        //theta=(0.0)*M_PI/180.0;
//        
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
//        // elliptical inclusion
//        Xn = X*cos(theta) - Z*sin(theta);
//        Zn = X*sin(theta) + Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
//        //--------------------------------------------------------------------------------
        
        
//        if (pos==0) {
//            X = particles->x[np]-xc;
//            Z = particles->z[np]-zc;
//            // DRAW INCLUSION
//            if (X*X + Z*Z < radius*radius) {
//                particles->phase[np] = 1;
//            }
//        }
//        else {
//            X = particles->x[np]-model.xmax;
//            Z = particles->z[np]-zc;
//            // DRAW INCLUSION
//            if (X*X + Z*Z < radius*radius) {
//                particles->phase[np] = 1;
//            }
//            
//            // DRAW INCLUSION
//            X = particles->x[np]-model.xmin;
//            Z = particles->z[np]-zc;
//            if (X*X + Z*Z < radius*radius) {
//                particles->phase[np] = 1;
//            }
//        }
        
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
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "Tp init" );
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
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
    
    // ---- T-Dependent marker types
    // -------------------- SPECIFIC TO YOANN's SETUP -------------------- //
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    X  = malloc (NX*sizeof(double));
    Z  = malloc (NZ*sizeof(double));
    XC = malloc (NCX*sizeof(double));
    ZC = malloc (NCZ*sizeof(double));
    
    for (k=0; k<NX; k++) {
        X[k] = mesh->xg_coord[k];
    }
    for (k=0; k<NCX; k++) {
        XC[k] = mesh->xc_coord[k];
    }
    for (l=0; l<NZ; l++) {
        Z[l] = mesh->zg_coord[l];
    }
    for (l=0; l<NCZ; l++) {
        ZC[l] = mesh->zc_coord[l];
    }
    
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
                
                if ( mesh->BCu.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[c] = -1;
                    mesh->BCu.val[c]  =  0;
                    
                    if (model->shear_style== 0 ) {
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
                    }
                    
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
                    if (model->shear_style== 1 ) {
                        
                        // Matching BC nodes WEST
                        if (k==0 ) {
                            mesh->BCu.type[c] = -2;
                            mesh->BCu.val[c]  = 0.0*model->EpsBG*Lx;
                        }
                        
                        // Matching BC nodes EAST
                        if (k==mesh->Nx-1 ) {
                            mesh->BCu.type[c] =  -12;
                            mesh->BCu.val[c]  = -0.0*model->EpsBG*Lx;
                        }
                        
                        // Free slip S
                        if (l==0 ) { //&& (k>0 && k<NX-1) ) {
                            mesh->BCu.type[c] =  11;
//                            mesh->BCu.val[c]  = -model->EpsBG*Lz/1.0;
                            mesh->BCu.val[c]  = 2.0*(mesh->zvx_coord[l]-0.0)*model->EpsBG;
                        }
                        
                        // Free slip N
                        if ( l==mesh->Nz) {// && (k>0 && k<NX-1)) {
                            mesh->BCu.type[c] =  11;
//                            mesh->BCu.val[c]  =  model->EpsBG*Lz/1.0;
                            mesh->BCu.val[c]  = 2.0*(mesh->zvx_coord[l]-0.0)*model->EpsBG;
                        }
                        
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
                
                if ( mesh->BCv.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCv.type[c] = -1;
                    mesh->BCv.val[c]  =  0;
                    
                    if (model->shear_style== 0 ) {
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Non-matching boundary WEST
                    if ( k==0 ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary EAST
                    if ( k==mesh->Nx ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                    }
                    
                    
                    if (model->shear_style== 1 ) {
                        // Matching BC nodes SOUTH
                        if (l==0 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = -0.0*model->EpsBG*Lz;
                        }
                        
                        // Matching BC nodes NORTH
                        if (l==mesh->Nz-1 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = 0.0*model->EpsBG*Lz;
                        }
                        
                        // Non-matching boundary points
                        if ( k==0   ) {    //&& (l>0 && l<NZ-1)
                            mesh->BCv.type[c] =  -12;
                            mesh->BCv.val[c]  =   0;
                        }
                        
                        // Non-matching boundary points
                        if ( k==mesh->Nx  ) { // && (l>0 && l<NZ-1)
                            mesh->BCv.type[c] =  -12;
                            mesh->BCv.val[c]  =   0;
                        }
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
    
    double Temperature = (model->user0+zeroC)/scaling.T;
    double scaled_Ma = 1e6*365.25*24*3600 / scaling.t;
    double TN = Temperature, TS = Temperature;
    double TW = Temperature, TE = Temperature;
    double Tbot, Tleft, Tright;
    
    printf("****************************************************************\n");
    printf("** Model time  = %2.10e Ma\n", model->time*scaling.t/(365.25*24.0*3600.0*1e6) );
    printf("** Temperature = %2.10e K\n", (Temperature * scaling.T));
    printf("****************************************************************\n");
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    for (l=0; l<NZ-1; l++) {
        for (k=0; k<NX-1; k++) {
            
            c = k + l*(NCX);
            
            //mesh->T[c] = Temperature ;
            
            if ( mesh->BCt.type[c] != 30 ) {
                
                if( model->isperiodic_x == 1 ){
                    
                    // WEST
                    
                    if ( k==0 ) {
                        mesh->BCt.type[c] = -2;
                        mesh->BCt.typW[l] = -2;
                        mesh->BCt.valW[l] = TW;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = -2;
                        mesh->BCt.typE[l] = -2;
                        mesh->BCt.valE[l] = TE;
                    }
                }
                else {
                    // WEST
                    
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typW[l] = 0;
                        mesh->BCt.valW[l] = TW;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = TE;
                    }
                }
                
                // SOUTH
                if ( l==0 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typS[k] = 1;
                    mesh->BCt.valS[k] = TS;
                }
                
                // NORTH
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typN[k] = 1;
                    mesh->BCt.valN[k] = TN;
                }
                
                // FREE SURFACE
                else {
                    if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.val[c]  = TN;
                    }
                }
                
                
            }
            
        }
    }
    
    
    //--------------------------------- CHEMICAL
    
    for (l=0; l<NZ-1; l++) {
           for (k=0; k<NX-1; k++) {
               
               c = k + l*(NCX);
                              
               if ( mesh->BCc.type[c] != 30 ) {
                   
                   if( model->isperiodic_x == 1 ){
                       
                       // WEST
                       
                       if ( k==0 ) {
                           mesh->BCc.type[c] = -2;
                           mesh->BCc.typW[l] = -2;
                           mesh->BCc.valW[l] = 0.0;
                       }
                       
                       // EAST
                       if ( k==NCX-1 ) {
                           mesh->BCc.type[c] = -2;
                           mesh->BCc.typE[l] = -2;
                           mesh->BCc.valE[l] = 0.0;
                       }
                   }
                   else {
                       // WEST
                       if ( k==0 ) {
                           mesh->BCc.type[c] = 0;
                           mesh->BCc.typW[l] = 0;
                           mesh->BCc.valW[l] = 0.0;
                       }
                       
                       // EAST
                       if ( k==NCX-1 ) {
                           mesh->BCc.type[c] = 0;
                           mesh->BCc.typE[l] = 0;
                           mesh->BCc.valE[l] = 0.0;
                       }
                   }
                   
                   // SOUTH
                   if ( l==0 ) {
                       mesh->BCc.type[c] = 0;
                       mesh->BCc.typS[k] = 0;
                       mesh->BCc.valS[k] = 0.0;
                   }
                   
                   // NORTH
                   if ( l==NCZ-1 ) {
                       mesh->BCc.type[c] = 0;
                       mesh->BCc.typN[k] = 0;
                       mesh->BCc.valN[k] = 0.0;
                   }
                   
                   // FREE SURFACE
                   else {
                       if ((mesh->BCc.type[c] == -1 || mesh->BCc.type[c] == 1 || mesh->BCc.type[c] == 0) && mesh->BCc.type[c+NCX] == 30) {
                           mesh->BCc.type[c] = 0;
                           mesh->BCc.val[c]  = 0.0;
                       }
                   }
                   
                   
               }
               
           }
       }


		
    
	free(X);
	free(Z);
	free(XC);
	free(ZC);
	printf("Velocity and pressure were initialised\n");
	printf("Boundary conditions were set up\n");
	
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
