// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2022  MDOODZ Developper team
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
#include "mdoodz-private.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    for ( int k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = 0.0;
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {    int np;
    
    const double T_init = (model.user0 + zeroC)/scaling.T;
    const double P_init = model.PrBG;
    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standard initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;                   // set initial particle velocity (unused)
        particles->phase[np] = 0;                                               // same phase number everywhere
        particles->d[np]     = 0;                                               // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->T[np]     = T_init;
        particles->P[np]     = P_init;
        particles->noise[np] = ((double)rand() / (double)RAND_MAX) - 0.5;
                
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases-1) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }

        // SET DEFAULT DUAL PHASE NUMBER
        particles->dual[np] = particles->phase[np];

        //--------------------------//
    }

    const double Lx = (double) (model.xmax - model.xmin) ;
    const double Lz = (double) (model.zmax - model.zmin) ;
    double Ax, Az;

    // Generate checkerboard
    for ( np=0; np<particles->Nb_part; np++ ) {
            Ax = cos( 6.0*2.0*M_PI*particles->x[np] / Lx  );
            Az = sin( 6.0*2.0*M_PI*particles->z[np] / Lz  );

            if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && particles->dual[np]==0 ) {
                particles->dual[np] += model.Nb_phases; 
            }
    }

    // Generate checkerboard
    for ( np=0; np<particles->Nb_part; np++ ) {
            Ax = cos( 12.0*2.0*M_PI*particles->x[np] / Lx  );
            Az = sin( 12.0*2.0*M_PI*particles->z[np] / Lz  );

            if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && particles->dual[np]==1 ) {
                particles->dual[np] += model.Nb_phases; 
            }
    }

    // Generate checkerboard
    for ( np=0; np<particles->Nb_part; np++ ) {
            Ax = cos( 24.0*2.0*M_PI*particles->x[np] / Lx  );
            Az = sin( 24.0*2.0*M_PI*particles->z[np] / Lz  );

            if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && particles->dual[np]==2 ) {
                particles->dual[np] += model.Nb_phases; 
            }
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {

    int   kk, k, l, c, c1;
    const double a = model->user3;
    const double b = model->user4;
    
    // Define dimensions;
    const double  Lx = (double) (model->xmax - model->xmin) ;
    const double  Lz = (double) (model->zmax - model->zmin) ;
        
    // Sizes
    const int NX  = mesh->Nx; 
    const int NCX = NX-1;
    const int NZ  = mesh->Nz;  
    const int NCZ = NZ-1;

    // Define initial T
    const double Tfix = (pow(mesh->p_in[0]*scaling.S/1e9/b, 1/a) + zeroC)/scaling.T;
    printf("Pressure: %2.2e , Temperature: %2.2e\n", mesh->p_in[0]*scaling.S,  pow(mesh->p_in[0]*scaling.S/1e9/b, 1/a));
    
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
        
    for (l=0; l<mesh->Nz+1; l++) {
        for (k=0; k<mesh->Nx; k++) {
            
            c = k + l*(mesh->Nx);
            
            if ( mesh->BCu.type[c] != 30 ) {
                
                // Internal points:  -1
                mesh->BCu.type[c] = -1;
                mesh->BCu.val[c]  =  0;
                
                // ------------------- PURE SHEAR ------------------- // 
                if (model->shear_style== 0 ) {
                
                    // Matching BC nodes WEST
                    if ( k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * (model->EpsBG - model->DivBG/3.0);
                    }
                    
                    // Matching BC nodes EAST
                    if ( k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * (model->EpsBG - model->DivBG/3.0);
                    }
                    
                    // Free slip SOUTH
                    if ( l==0 ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                    
                }

                // ------------------- PERIODIC SIMPLE SHEAR ------------------- // 
                if (model->shear_style== 1 ) {
                    
                    // Matching BC nodes WEST
                    if ( k==0 ) {
                        mesh->BCu.type[c] = -2;
                        mesh->BCu.val[c]  = 0.0*model->EpsBG*Lx;
                    }
                    
                    // Matching BC nodes EAST
                    if ( k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] =  -12;
                        mesh->BCu.val[c]  = -0.0*model->EpsBG*Lx;
                    }
                    
                    // Free slip S
                    if ( l==0 ) { //&& (k>0 && k<NX-1) ) {
                        mesh->BCu.type[c] =  11;
                        mesh->BCu.val[c]  = -1*model->EpsBG*Lz;
                    }
                    
                    // Free slip N
                    if ( l==mesh->Nz ) {// && (k>0 && k<NX-1)) {
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
        
    for (l=0; l<mesh->Nz; l++) {
        for (k=0; k<mesh->Nx+1; k++) {
            
            c  = k + l*(mesh->Nx+1);
            
            if ( mesh->BCv.type[c] != 30 ) {
                
                // Internal points:  -1
                mesh->BCv.type[c] = -1;
                mesh->BCv.val[c]  =  0;
                
                // ------------------- PURE SHEAR ------------------- // 
                if ( model->shear_style== 0 ) {
                
                    // Matching BC nodes SOUTH
                    if ( l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * (model->EpsBG + model->DivBG/3.0);
                    }
                    
                    // Matching BC nodes NORTH
                    if ( l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * (model->EpsBG + model->DivBG/3.0);
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
                
                // ------------------- PERIODIC SIMPLE SHEAR ------------------- // 
                if (model->shear_style== 1 ) {
                    // Matching BC nodes SOUTH
                    if ( l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = -0.0*model->EpsBG*Lz + model->DivBG/3.0;
                    }
                    
                    // Matching BC nodes NORTH
                    if ( l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0*model->EpsBG*Lz + model->DivBG/3.0;
                    }
                    
                    // Non-matching boundary points
                    if ( k==0 ) {    //&& (l>0 && l<NZ-1)
                        mesh->BCv.type[c] =  -12;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary points
                    if ( k==mesh->Nx ) { // && (l>0 && l<NZ-1)
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
        
    for (l=0; l<NCZ; l++) {
        for (k=0; k<NCX; k++) {
            
            c  = k + l*(NCX);
            
            if ( mesh->BCt.type[c] != 30 ) {
                        
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
        
        for (l=0; l<mesh->Nz-1; l++) {
            for (k=0; k<mesh->Nx-1; k++) {
                
                c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {

                    // !!!!!!!!!! FORCE TEMPERATURE
                    mesh->T[c] = Tfix;
                    
                    // LEFT
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // RIGHT
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // BOT
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // TOP
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
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
        
    printf("Boundary conditions were set up\n");
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
