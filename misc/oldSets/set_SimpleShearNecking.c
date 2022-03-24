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
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {
    
    int k;
    double TopoLevel = 0.0e3/scaling.L;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
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
	
	s1 = sizeof(int);
	s2 = sizeof(double);
    
    double H, xl, lx;
    int    k, nb_layers = 20;

    H = (model.zmax - model.zmin);
    lx = (model.xmax - model.xmin);
    
    double gsbg = 2e-3/scaling.L;
    double Temperature = (model.user2+zeroC)/scaling.T;
    double rad=model.user1/scaling.L, xc=-0e3/scaling.L, zc=-0e3/scaling.L;
    double spacing = 5.0e-3/scaling.L;
    
    double alpha = (model.user3)/180*M_PI; // angle radians
    double h     = model.user1/scaling.L;  // layer thickness
    double a     = tan(alpha);
    double b     = (h/2.0)/cos(alpha);    

    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->d[np]     = gsbg;                                            // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->T[np]     = Temperature;
        particles->phase[np] = 2;

        //----------------------------------
        // INCLUSIONS INCLUSIONS INCLUSIONS
        //----------------------------------
        
        
        
        // Dipping layer with angle alpha
        if ( (particles->z[np] < (a*particles->x[np] + b) ) && (particles->z[np] > (a*particles->x[np] - b) ) ) {
            particles->phase[np] = 0;
        }
        
        // Central inclusion
        if ( pow(particles->x[np],2) + pow(particles->z[np],2) < pow(rad/10,2) ) {
            particles->phase[np] = 2;
        }
        
        
        
//        // SE inclusion
//        if ( pow(particles->x[np]-0.5/scaling.L,2) + pow(particles->z[np]+0.5/scaling.L,2) < pow(rad,2) ) {
//            particles->phase[np] = 0;
//        }
//        // NW inclusion
//        if ( pow(particles->x[np]+0.5/scaling.L,2) + pow(particles->z[np]-0.5/scaling.L,2) < pow(rad,2) ) {
//            particles->phase[np] = 0;
//        }
//        // SW inclusion
//        if ( pow(particles->x[np]+0.5/scaling.L,2) + pow(particles->z[np]+0.5/scaling.L,2) < pow(rad,2) ) {
//            particles->phase[np] = 0;
//        }
//        // NE inclusion
//        if ( pow(particles->x[np]-0.5/scaling.L,2) + pow(particles->z[np]-0.5/scaling.L,2) < pow(rad,2) ) {
//            particles->phase[np] = 0;
//        }
//
//        // phase 1 for visualisation:
//        if (particles->phase[np] == 0 && particles->x[np]>0 && particles->z[np]>0) {
//            particles->phase[np] = 1;
//        }
//        if (particles->phase[np] == 0 && particles->x[np]<0 && particles->z[np]<0) {
//            particles->phase[np] = 1;
//        }
//
//        //----------------------------------
//        // STRIPES STRIPES STRIPES STRIPES
//        //----------------------------------
//        // Horizontal Layering
//        if (particles->phase[np]==2 && particles->z[np]> 1.0*H/6.0 && particles->z[np]< 2.0*H/6.0){
//            //if (particles->x[np]> 1.0*H/6.0 && particles->x[np]< 2.0*H/6.0)
//            particles->phase[np] = 3;
//        }
//        if (particles->phase[np]==2 && particles->z[np]< 0.0*H/6.0 && particles->z[np]> -1.0*H/6.0){
//            particles->phase[np] = 3;
//        }
//        if (particles->phase[np]==2 && particles->z[np]< -2.0*H/6.0 && particles->z[np]> -3.0*H/6.0){
//            particles->phase[np] = 3;
//        }
////        // Vertical alternance:
//        if (particles->x[np]> 1.0*H/6.0 && particles->x[np]< 2.0*H/6.0){
//            if (particles->phase[np] == 3) particles->phase[np] = 30;
//            if (particles->phase[np] == 2) particles->phase[np] = 20;
//        }
//        if (particles->x[np]> -1.0*H/6.0 && particles->x[np]< 0.0*H/6.0){
//            if (particles->phase[np] == 3) particles->phase[np] = 30;
//            if (particles->phase[np] == 2) particles->phase[np] = 20;
//        }
//
//        if (particles->x[np]> -3.0*H/6.0 && particles->x[np]< -2.0*H/6.0){
//            if (particles->phase[np] == 3) particles->phase[np] = 30;
//            if (particles->phase[np] == 2) particles->phase[np] = 20;
//        }
//
//        if (particles->phase[np] == 20) particles->phase[np] = 3;
//        if (particles->phase[np] == 30) particles->phase[np] = 2;
//
//
//        // here we define the initial particules velocity that should comply with the background homogeneous shear conditions
//        particles->Vx[np]    = -0.5*model.EpsBG*(particles->z[np] + (model.zmax-model.zmin)/2);
//        particles->Vz[np]    = 0;
//        particles->rho[np]   = materials->rho[particles->phase[np]];

    
        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
        
    }
    
    MinMaxArray(particles->Vx,  scaling.V, particles->Nb_part,   "Vxp init" );
    MinMaxArray(particles->Vz,  scaling.V, particles->Nb_part,   "Vzp init" );
    MinMaxArray(particles->rho, scaling.rho, particles->Nb_part, "Rho init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "T init" );
        
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
                        mesh->BCu.val[c]  = -1*model->EpsBG*Lz;
                    }
                    
                    // Free slip N
                    if ( l==mesh->Nz) {// && (k>0 && k<NX-1)) {
                        mesh->BCu.type[c] =  11;
                        mesh->BCu.val[c]  =  1*model->EpsBG*Lz;
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
                    if ( (k==0) ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx) ) {
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
                    if ( (k==0)   ) {    //&& (l>0 && l<NZ-1)
                        mesh->BCv.type[c] =  -12;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary points
                    if ( (k==mesh->Nx)  ) { // && (l>0 && l<NZ-1)
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
    
    double Ttop = 273.15/scaling.T;
    double Tbot, Tleft, Tright, Tbg=(model->user2+zeroC)/scaling.T;
    
    
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
                
                // LEFT
                if ( k==0 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typW[l] = 1;
                    mesh->BCt.valW[l] = Tbg;
                }
                
                // RIGHT
                if ( k==NCX-1 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typE[l] = 1;
                    mesh->BCt.valE[l] = Tbg;
                }
                
                // BOT
                if ( l==0 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typS[k] = 1;
                    mesh->BCt.valS[k] = Tbg;
                }
                
                // TOP
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 1;
                    mesh->BCt.typN[k] = 1;
                    mesh->BCt.valN[k] = Tbg;
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
