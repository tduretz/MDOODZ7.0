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
    
    int k;
    double TopoLevel = 11.0e3/scaling.L;
    
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
    // Define dimensions;
    double rho_crust  = materials->rho[0];
    double top = model.user0/scaling.L;
    double slope = (top-model.zmin)/(model.xmax-model.xmin);
    double A     = 1.0/10*model.dx;
    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standard initialisation of particles
        particles->x[np]    += A * (((double)rand() / RAND_MAX) - 0.5);
        particles->x[np]    += A * (((double)rand() / RAND_MAX) - 0.5);
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;
        particles->phase[np] = 0;
        particles->T[np]     = zeroC;
        particles->d[np]     = 0.0;
        particles->phi[np]   = 0.0;
        particles->rho[np]   = rho_crust;
        
        // Draw inclined interface
        if (particles->z[np] < slope*(particles->x[np]-model.xmin) +  model.zmin) particles->phase[np] = 1;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {

	int   kk, k, l, c;
	double *X, *Z, *XC, *ZC;
	int   NX, NZ, NCX, NCZ;
	double VxBC, frac= 0.05;
    // Define dimensions;
    double Lx = (double) (model->xmax - model->xmin) ;
//    double Lz = (double) (model->zmax - model->zmin);
		
    NX  = mesh->Nx;
	NZ  = mesh->Nz;
	NCX = NX-1;
	NCZ = NZ-1;
	
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

    VxBC = -model->EpsBG*mesh->xg_coord[0];
	
	/* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vx on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type  3: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -1: not a BC point (tag fo inner points)                                                           */
	/* Type -2: periodic in the x direction (matches the physical boundary)                                    */
	/* --------------------------------------------------------------------------------------------------------*/
        NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
	        
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
                        mesh->BCu.val[c]  = VxBC;                // printf("uW_BC=%lf\n", mesh->BCu.val[c]);
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;               // printf("uE_BC=%lf\n", mesh->BCu.val[c]);
                    }
                    
                    // Free slip
                    if ( l==0 && (k>0 && k<mesh->Nx-1 ) ) {
                        mesh->BCu.type[c] =  11;
                        mesh->BCu.val[c]  =  VxBC;
                        if (mesh->xg_coord[k]> model->xmax - frac*Lx) {
                            mesh->BCu.val[c]  =  -(VxBC + VxBC/(frac*Lx) * (mesh->xg_coord[k] + (model->xmax - frac*Lx)));
                        }
//                        printf("VxBC = %2.2e\n", mesh->BCu.val[c]*scaling.V);
                    }
                    
                    if ( l==mesh->Nz && (k>0 && k<mesh->Nx-1 ) ) {
                        mesh->BCu.type[c] =  13;
                        mesh->BCu.val[c]  =  0;
                    }
                }
                
			}
		}
			
	/* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vz on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type  3: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -1: not a BC point (tag fo inner points)                                                           */
	/* Type -2: periodic in the x direction (does not match the physical boundary)                             */
	/* Type -10: useless point (set to zero)                                                                   */
	/* --------------------------------------------------------------------------------------------------------*/
			
		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		
		for (l=0; l<mesh->Nz; l++) {
			for (k=0; k<mesh->Nx+1; k++) {
				
				c  = k + l*(mesh->Nx+1);
                
                if ( mesh->BCv.type[c] != 30 ) {
                    
                    
                    // Boundary Conditions
                    
                    // Internal points:  -1
                    mesh->BCv.type[c] = -1;
                    mesh->BCv.val[c]  =  0;
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = model->EpsBG*mesh->zg_coord[l];
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = model->EpsBG*mesh->zg_coord[l];
                    }
                    
                    // Non-matching boundary points
                    if ( (k==0) && (l>0 && l<mesh->Nz-1 ) ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary points
                    if ( (k==mesh->Nx) && (l>0 && l<mesh->Nz-1 )) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                }
				
			}
		}
		    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
    /* --------------------------------------------------------------------------------------------------------*/
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        
        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {
                
                c  = k + l*(NCX);
                
                if (mesh->BCt.type[c] != 30) {
                    
                    // Boundary Conditions
                    
                    // Internal points:  -1
                    mesh->BCp.type[c] = -1;
                    mesh->BCp.val[c]  =  0;

                }
                
            }
        }
    
	/* -------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for T on all grid levels                                                                   */
	/* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)      */
	/* Type  3: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)             */
	/* Type -1: not a BC point (tag fo inner points)                                                          */
	/* Type -2: periodic in the x direction (matches the physical boundary)                                   */
	/* -------------------------------------------------------------------------------------------------------*/
    
    double Ttop = 273.15/scaling.T;    
	
	for(kk=0; kk<1; kk++) {
		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		
		for (l=0; l<mesh->Nz-1; l++) {
			for (k=0; k<mesh->Nx-1; k++) {
				
				c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
                    
                    //				// Internal points:  -1
                    //				mesh->BCt.type[c] = -1;
                    //				mesh->BCt.val[c]  =  0;
                    
                    // LEFT
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  =  mesh->T[c];
                    }
                    
                    // RIGHT
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // BOT
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // TOP
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = (Ttop*mesh->rho_n[c] * mesh->Cv[c]) / (mesh->rho_n[c] * mesh->Cv[c]);
                        }
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
