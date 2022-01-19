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

void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
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
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
	
    double xc=0.0, zc=-0.0, Tc=673/scaling.T, Tbg=673/scaling.T;
    
    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->Vx[np]    = 0;
        particles->Vz[np]    = 0;
        particles->T[np]     = 0;
        
        particles->phase[np] = 0;
        if (particles->z[np]>200e3/scaling.L && particles->z[np]<800e3/scaling.L && particles->x[np]<800e3/scaling.L) {
            particles->phase[np] = 1;
        }
    }
    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
    
    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta, Ttop=273.15/scaling.T;//(scaling.T+273)/scaling.T;
    double Lx, Lz;
    
    //    // For Piston Setup
    //	VzBC = -model->Pdriv*width*width/24/eta;
    //	printf("Pdriv = %2.2e VzBC = %2.2e m/s %2.2e %2.2e\n", model->Pdriv*(scaling.S/scaling.L), VzBC*scaling.V, scaling.S, scaling.V);
    
    
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
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
    
    // Switch off gravity
    if (model->time > 20e3*(3600*24*365.25)/scaling.t) {
        printf("%%%%%\n%%%%%%\nACHTUNG !!!!!!!!!!!!!!!!!!\ GRAVITY WAS SWITCHED OFF !!!!!!\n%%%%%%%\n%%%%%%\n");
        model->gx = 0;
        model->gz = 0;
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
                    mesh->BCv.type[c] =   11;
                    mesh->BCv.val[c]  =   0;
                }
                
                // Non-matching boundary EAST
                if ( k==mesh->Nx ) {
                    mesh->BCv.type[c] =   13;
                    mesh->BCv.val[c]  =   0;
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
    
    double Tbot, Tleft, Tright;
    
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
                    mesh->BCt.val[c]  = mesh->T[c];
                    mesh->BCt.typW[l] = 0;
                    mesh->BCt.valW[l] = mesh->T[c];
                }
                
                // EAST
                if ( k==NCX-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.val[c]  = mesh->T[c];
                    mesh->BCt.typE[l] = 0;
                    mesh->BCt.valE[l] = mesh->T[c];
                }
                
                // SOUTH
                if ( l==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.val[c]  = mesh->T[c];
                    mesh->BCt.typS[k] = 0;
                    mesh->BCt.valS[k] = mesh->T[c];
                }
                
                // NORTH
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.val[c]  = mesh->T[c];
                    mesh->BCt.typN[k] = 0;
                    mesh->BCt.valN[k] = mesh->T[c];
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

