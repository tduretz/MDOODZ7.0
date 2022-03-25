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
#include "mdoodz-private.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = 0.0e3/scaling.L; // sets zero initial topography
    double A = model.zmax/2.0;
    double s = model.zmax/4.0;

    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = A*exp(-pow(topo_chain->x[k],2) / 2.0/s/s  );// TopoLevel;
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
    double P_init = model.PrBG;
    double setup  = (int)model.user1;
    double radius = model.user2/scaling.L;
    double X, Z, Xn, Zn, xc = 0.0, zc = 0.0, sa=radius/2.0, la=radius*2.0, theta=-30.0*M_PI/180.0;
    double Ax, Az;

    // Fixed random seed
    srand(69);

    char *ph_hr;
    int nb_elems;

    if (setup==2) {
        // Read input file
        FILE *fid;
        nb_elems = 1921*1921;
        ph_hr = malloc((nb_elems)*sizeof(char));
        fid = fopen(model.input_file, "rb"); // Open file
        if (!fid){
            fprintf(stderr, "\nUnable to open file %s. I will exit here ... \n", model.input_file); exit(2);
        }
        fread( ph_hr, sizeof(char), nb_elems, fid);
        fclose(fid);
    }
    // image properties
    int nx_hr = 1922, nz_hr = 1922, ix, iz;
    double dx_hr = (model.xmax - model.xmin)/(nx_hr-1.0);
    double dz_hr = (model.zmax - model.zmin)/(nz_hr-1.0);
    double dstx, dstz;
    
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
        
        if (setup==2) {
            // ------------------------- //
            // Locate markers in the image files
            // Find index of minimum/west temperature node
            dstx = ( particles->x[np] - model.xmin );
            ix   = ceil( dstx/dx_hr - 0.0 ) - 1;
            // Find index of minimum/west pressure node
            dstz = ( particles->z[np] - model.zmin );
            iz   = ceil( dstz/dz_hr  - 0.0 ) - 1;
            // Attribute phase
            if (ix<0) {printf("sauceisse!!!\n"); exit(1);}
            if (iz<0) {printf("puréee!!!\n"); exit(1);}
            if (ix + iz*(nx_hr-1) > nb_elems) {printf("puréee!!!\n"); exit(1);}
    //        printf("%d\n", (int)ph_hr[ix + iz*(nx_hr-1)]);
            particles->phase[np] = (int)ph_hr[ix + iz*(nx_hr-1)];
        }
        if (setup==1) {
            // DRAW INCLUSION
            X  = particles->x[np]-xc;
            Z  = particles->z[np]-zc;
            Xn = X*cos(theta) - Z*sin(theta);
            Zn = X*sin(theta) + Z*cos(theta);
            if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 1;
        }
        
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases-1) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }

        // SET DEFAULT DUAL PHASE NUMBER
        particles->dual[np] = particles->phase[np];

        //--------------------------//
    }

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

    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "Tp init" );
    if (setup==1) DoodzFree(ph_hr);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {

    int    kk, k, l, c, c1;
    int    NX, NZ, NCX, NCZ;
    double Lx, Lz, a = model->user3, b = model->user4;
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
        
    // Sizes
    NX  = mesh->Nx; NCX = NX-1;
    NZ  = mesh->Nz; NCZ = NZ-1;

    double Tfix = (pow(mesh->p_in[0]*scaling.S/1e9/b, 1/a) + zeroC)/scaling.T;
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
    double Tbot, Tleft, Tright;
        
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
