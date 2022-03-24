
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "header_MDOODZ.h"
#include "time.h"

#define USE_BC_USER 0

// Function declarations
void SetBCs_freeSlipBox( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials );


// Set physical properties on the grid and boundary conditions
void SetBCs_new( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials ) {
    switch (model->BC_setup_type) {
        case 0:
            //            // now points to the function SetBCs in old fashion input files
            //            // eventually it would point towards the function SetBCs_user in User.c
            //            SetBCs( mesh, model, scaling, particles, materials );
            //            break;
        case 1:
            SetBCs_freeSlipBox( mesh, model, scaling, particles, materials );
            break;
        case 2:
#if USE_BC_USER
            SetBCs_user( mesh, model, scaling, particles, materials );
            break;
#else
            printf("error: you requested the user defined boundary conditions but you did not compile with this option. Recompile with the option USE_BC_USER=1 and try again.\n");
#endif
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


// Set physical properties on the grid and boundary conditions
void SetBCs_freeSlipBox( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials ) {
    
    int   k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ;
    
    int StressBC_W=0, StressBC_E=0;
    //if ( model->user0 == 1 ) StressBC_E = 1;
    double TN = 273.15/scaling.T, TS = 1330/scaling.T;
    double TW = 273.15/scaling.T, TE = 1330/scaling.T;
    
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
    
    //--------------------------------------------------------------------------------
    
    for (l=0; l<mesh->Nz+1; l++) {
        for (k=0; k<mesh->Nx; k++) {
            
            c = k + l*(mesh->Nx);
            
            if ( mesh->BCu.type[c] != 30 ) {
                
                // Internal points:  -1
                mesh->BCu.type[c] = -1;
                mesh->BCu.val[c]  =  0;
                
                // Matching BC nodes WEST
                if (k==0 && l>0 && l<mesh->Nz ) { //
                    
                    if ( StressBC_W==1 ) {
                        mesh->BCu.type[c] = 2;
                        c1 = k + (l-1)*(mesh->Nx-1) ;
                        mesh->BCu.val[c]  = -mesh->p_lith[c1];
                    }
                    else {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -X[k]*model->EpsBG;//0.0;
                    }
                }
                
                // Matching BC nodes EAST
                if (k==mesh->Nx-1 && l>0 && l<mesh->Nz) { //
                    if ( StressBC_E==1 ) {
                        mesh->BCu.type[c] = 2;
                        c1 = k + (l-1)*(mesh->Nx-1) -1;
                        mesh->BCu.val[c]  = -mesh->p_lith[c1];
                    }
                    else {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -X[k]*model->EpsBG;//0.0;
                    }
                }
                
                
                // Free slip
                if (l==0  ) {
                    mesh->BCu.type[c] = 13;
                    mesh->BCu.val[c]  = 0.0;
                }
                
                if ( l==mesh->Nz) {
                    mesh->BCu.type[c] = 13;
                    mesh->BCu.val[c]  = 0.0;
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
                
                // Matching BC nodes SOUTH
                if (l==0  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = Z[l]*model->EpsBG;//0.0;
                }
                
                // Matching BC nodes NORTH
                if (l==mesh->Nz-1  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = Z[l]*model->EpsBG;//0.0;
                }
                
                // Non-matching boundary points
                if ( (k==0) ) {
                    if ( StressBC_W==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0.0;
                    }
                    else {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;
                    }
                }
                
                // Non-matching boundary points
                if ( (k==mesh->Nx) ) {
                    if ( StressBC_E==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0.0;
                    }
                    else{
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;
                    }
                }
                
                //                    // Normal stress applied to the free surface
                //                    if (l<mesh->Nz-1 && k>0 && k<mesh->Nx) {
                //                        c1 = k + (l)*(mesh->Nx-1)-1;
                //                        if (mesh->BCp.type[0][c1]==31) mesh->BCv.val[c]  =  -0*mesh->p_lith[c1-(mesh->Nx-1)];
                //                    }
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
            
            if (mesh->BCt.type[c] != 30) {
                
                // Internal points:  -1
                mesh->BCp.type[c] = -1;
                mesh->BCp.val[c]  =  0.0;
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
    
    for (l=0; l<mesh->Nz-1; l++) {
        for (k=0; k<mesh->Nx-1; k++) {
            
            c = k + l*(NCX);
            
            if ( mesh->BCt.type[c] != 30 ) {
                
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
                
                // SOUTH
                if ( l==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typS[k] = 1;
                    mesh->BCt.valS[k] = TS;
                }
                
                // NORTH
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 0;
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
