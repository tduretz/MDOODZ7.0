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

double zConvertToPolar( double x, double y_limit, scale* scaling ) {
    double r = Rad_Earth/scaling->L - y_limit;
    return sqrt((r - x)*(r + x));
}

double xConvertToPolar( double r, double x_limit, scale* scaling ) {
    double   tet = acos( x_limit / (Rad_Earth/scaling->L) );
    return   r*cos(tet);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = 0.0e3/scaling.L; // sets zero initial topography

    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        if (model.polar==0) {
            topo_chain->z[k]     = Rad_Earth/scaling.L;
        }
        if (model.polar==1) {
            // see PolarCoordinatesStuff.py
            topo_chain->z[k]     =  zConvertToPolar(topo_chain->x[k], 0.0, &scaling);
        }
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np;
    double Lx = (double) (model.xmax - model.xmin);
    double Lz = (double) (model.zmax - model.zmin);
    double Tbot = 773.15/scaling.T;
    double gsbg   = 2e-3/scaling.L;                           // reference grain size
    double Rad=Rad_Earth/scaling.L;
    
    double  l2 = model.user3/scaling.L;                    // longueur du biseau
    
    double  H = 61e3/scaling.L;								// lower plate thickness
    double C = 10e3/scaling.L;								//crust thickness
    //double W = 0e3/scaling.L;								//channel width
    double L = 2500e3/scaling.L;								//length upper plate
    double T = 330e3/scaling.L;									//length mantle wedge
    double  bis = model.user3/scaling.L;						//length biseau
    double	bis_fin = 0e3/scaling.L;							// épaisseur fin de biseau
    
    //lower plate horizontal
    
    double x0 = -150e-3/scaling.L, y0 = -H;
    
    //slab
    
    double x1 = 850e3/scaling.L, y1 =-C/scaling.L; //trench
    double x4 = x1 - 400e3/scaling.L, y4 = y1 - 215e3/scaling.L; //tip of the slab
    double x2 = x1 + 40e3/scaling.L, y2 = y1 - H; //courbure du slab
    double x3 = x4 + 40e3/scaling.L, y3 = y4 - H; //bottom of the slab
    double x5 = x1 , y5  = 0e3/scaling.L; //debut plaque sup coté trench
    double x6 = x5 - L, y6 = y5; // fin de la plaque sup
    double x7 = x5 - T, y7 = 87e3/scaling.L; //bas du wedge plaque sup
    double x8 =  2750e3/scaling.L,  y8 = y2;              // début biseau
    double x9 = x8 + bis, y9 = -bis_fin; // fin biseau
    double x10 = x6+100e3/scaling.L, y10 = y7; // biseau plaque sup
    
    //Slab limits
    
    //Upper slab limit
    double a0 = (y1-(-y4))/(x1-x4);
    double b0 = y1 - (x1*a0);
    
    //Lower slab limit
    double a2 = (y2-(-y3))/(x2-x3);
    double b2 = y2 - (x2*a2);
    
    //Left slab limit
    double a1 = (y3-(-y4))/(x3-x4);
    double b1 = y3 - (x3*a1);
    
    //Right slab limit
    double a3 = (y2-(-y1))/(x2-x1);
    double b3 = y2 -(x2*a3);
    
    //Upper plate limit right
    double a4 = (y5-y7)/(x5-x7);
    double b4 = y5 -(x5*a4);
    
    //Upper plate limit left
    double a5 = (y10-y6)/(x10-x6);
    double b5 = y10 -(x10*a5);
    
    //Biseau
    double  abiseau = (y9-(-y2)) / (x9-x8);                   // pente du biseau
    double  bbiseau = y9 - (x9*abiseau);
    
    double zS, zN, xW, xE, zW, zE, tet, r;
    
    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;
        particles->d[np]     = gsbg;                          // same grain size everywhere
        particles->phi[np]   = 0.0;                           // zero porosity everywhere
        particles->X[np]     = 0.0;                           // X set to 0
        particles->T[np]     = Tbot;                          // same temperature everywhere
        
        particles->phase[np] = 1;
        
        // Upper plate 1
        
        // Regarde là Ben:
        // chaque boite est limitées par les coordoneés, xW, xE, yS, yN
        // Coordonnée polaire (r) du marker:
        r   = sqrt( pow(particles->x[np],2) + pow(particles->z[np],2) );
        
        
//        double angle = 0;
//        double x0 = 0.0;
//        double z0 = Rad-200e3/scaling.L;
//        double a_ell = 50e3/scaling.L;
//        double b_ell = 50e3/scaling.L;
//        double x_ell = (particles->x[np]-x0)*cos(angle) + (particles->z[np]-z0)*sin(angle);
//        double z_ell =-(particles->x[np]-x0)*sin(angle) + (particles->z[np]-z0)*cos(angle);
//        if (pow(x_ell/a_ell,2) + pow(z_ell/b_ell,2) < 1) particles->phase[np] = 0;
        
        
        // Upper plate biseau West
        xW  = xConvertToPolar( r, x6 , &scaling );
        xE  = xConvertToPolar( r, x10, &scaling );
        zS  = zConvertToPolar( particles->x[np], a5*particles->x[np]+b5, &scaling );
        if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1) {
            particles->phase[np] = 2;
        }

        // Upper plate central

        xW  = xConvertToPolar( r, x10, &scaling );
        xE  = xConvertToPolar( r, x7 , &scaling );
        zS  = zConvertToPolar( particles->x[np], y7, &scaling );
                if (particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS && model.user1 == 1 && model.user2 == 1 ) {
            particles->phase[np] = 2;
        }

        // Upper plate biseau Est
        xW  = xConvertToPolar( r, x7 , &scaling );
        zS  = zConvertToPolar( particles->x[np], a4*particles->x[np]+b4, &scaling );
        if ( particles->x[np] > xW && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1) {
            particles->phase[np] = 2;
        }


        // Weak layer dipping slab

        xW  = xConvertToPolar( r, x4, &scaling );
        xE  = xConvertToPolar( r, x1, &scaling );
        zN  = zConvertToPolar( particles->x[np], a0*particles->x[np]+b0, &scaling );
        zS  = zConvertToPolar( particles->x[np], a0*particles->x[np]+b0+C, &scaling );


        if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS && particles->z[np] < zN && model.user2 == 1) {
                    particles->phase[np] = 3;
        }

        // Dipping slab

        xW  = xConvertToPolar( r, x4, &scaling );
        xE  = xConvertToPolar( r, x1, &scaling );
        zN  = zConvertToPolar( particles->x[np], a0*particles->x[np]+b0+C, &scaling );
        zS  = zConvertToPolar( particles->x[np], a0*particles->x[np]+b0+C+H, &scaling );


        if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS && particles->z[np] < zN && model.user2 == 1) {
                    particles->phase[np] = 0;
        }



        // Weak layer lower plate central

        xW  = xConvertToPolar( r, x1, &scaling );
        xE  = xConvertToPolar( r, x8 , &scaling );
        zS  = zConvertToPolar( particles->x[np], C, &scaling );
                if (particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS && model.user1 == 1 && model.user2 == 1 ) {
            particles->phase[np] = 3;
        }

        // Lower plate central

        xW  = xConvertToPolar( r, x1, &scaling );
        xE  = xConvertToPolar( r, x8 , &scaling );
        zN  = zConvertToPolar( particles->x[np], C, &scaling );
        zS  = zConvertToPolar( particles->x[np], C+H, &scaling );
                if (particles->x[np] > xW && particles->x[np] < xE && particles->z[np] < zN && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1 ) {
            particles->phase[np] = 0;
        }



        // Lower plate biseau est
        xW  = xConvertToPolar( r, x8 , &scaling );
        zS  = zConvertToPolar( particles->x[np], abiseau*particles->x[np]+bbiseau, &scaling );
        if ( particles->x[np] > xW && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1) {
            particles->phase[np] = 0;
        }


        // Lower plate biseau est
        xW  = xConvertToPolar( r, x8 , &scaling );
        zS  = zConvertToPolar( particles->x[np], abiseau*particles->x[np]+bbiseau, &scaling );
        if ( particles->x[np] > xW && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1) {
            particles->phase[np] = 0;
        }

        // Weak layer lower plate biseau est
        xE  = xConvertToPolar( r, x9 , &scaling );
        xW  = xConvertToPolar( r, x8 , &scaling );
        zS  = zConvertToPolar( particles->x[np], C, &scaling );
        if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS  && model.user1 == 1 && model.user2 == 1) {
            particles->phase[np] = 3;
        }
         
        // lower mantle
        xW  = xConvertToPolar( r, -3500e3/scaling.L , &scaling );
        xE  = xConvertToPolar( r, 3500e3/scaling.L , &scaling );
        zN  = zConvertToPolar( particles->x[np], 660e3/scaling.L, &scaling );
        if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] < zN && model.user1 == 1 && model.user2 == 1) {particles->phase[np] = 4;
        }
        
        // lower mantle
               xW  = xConvertToPolar( r, -3500e3/scaling.L , &scaling );
               xE  = xConvertToPolar( r, 3500e3/scaling.L , &scaling );
               zN  = zConvertToPolar( particles->x[np], Rad, &scaling );
               if ( particles->x[np] > xW && particles->x[np] < xE && particles->z[np] > zS && model.user1 == 1 && model.user2 == 1) {particles->phase[np] = 5;
               }

        
        //
        //      // biseau
        //        if (particles->x[np] > x8 &&  particles->z[np] > abiseau*particles->x[np]+bbiseau-C+Rad) {
        //            particles->phase[np] = 0;
        //        }
        //
        //              // biseau weak layer
        //
        //        if (particles->x[np] > x8 && particles->x[np] < x9 && particles->z[np] > Rad-C && model.user2 == 1) {
        //            particles->phase[np] = 3;
        //        }
        //
        //		// end right model
        //        if (particles->x[np] > x9 && model.user2 == 1) {
        //            particles->phase[np] = 1;
        //        }
        
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
    
    int   k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ;
    
    int StressBC_W=0;
    int StressBC_E=0;
    if ( model->user0 == 1 ) StressBC_E = 1;
    
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
                        mesh->BCu.val[c]  = 0.0;
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
                        mesh->BCu.val[c]  = 0.0;
                    }
                }
                
                
                // Free slip
                if (l==0  ) {
                    mesh->BCu.type[c] = 13;
                    mesh->BCu.val[c]  =  0;
                }
                
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
    int numVz = 0, numVz0 = 0;
    
    for (l=0; l<mesh->Nz; l++) {
        for (k=0; k<mesh->Nx+1; k++) {
            
            c  = k + l*(mesh->Nx+1);
            
            if ( mesh->BCv.type[c] != 30 ) {
                
                numVz++;
                
                // Internal points:  -1
                mesh->BCv.type[c] = -1;
                mesh->BCv.val[c]  =  0;
                
                // Matching BC nodes SOUTH
                if (l==0  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = 0.0;
                }
                
                // Matching BC nodes NORTH
                if (l==mesh->Nz-1  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = 0.0;
                }
                
                // Non-matching boundary points
                if ( (k==0) ) {
                    if ( StressBC_W==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0;
                    }
                    else {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                }
                
                // Non-matching boundary points
                if ( (k==mesh->Nx) ) {
                    if ( StressBC_E==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0;
                    }
                    else{
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                }
                
                //                    // Normal stress applied to the free surface
                //                    if (l<mesh->Nz-1 && k>0 && k<mesh->Nx) {
                //                        c1 = k + (l)*(mesh->Nx-1)-1;
                //                        if (mesh->BCp.type[0][c1]==31) mesh->BCv.val[c]  =  -0*mesh->p_lith[c1-(mesh->Nx-1)];
                //                    }
            }
            if ( mesh->BCv.type[c] == 0) numVz0 ++;
            
            //            if ( mesh->BCv.type[c] != 30 ) {
            //
            //                numVz++;
            //
            ////                                // Internal points:  -1
            ////                                mesh->BCv.type[c] = -1;
            ////                                mesh->BCv.val[c]  =  0;
            //
            //                // Non-matching boundary points
            //                if ( (k==0) ) {
            //                    mesh->BCv.type[c] =   13;
            //                    mesh->BCv.val[c]  =   0;
            //                }
            //
            //                // Non-matching boundary points
            //                if ( (k==mesh->Nx) ) {
            //                    mesh->BCv.type[c] =   13;
            //                    mesh->BCv.val[c]  =   0;
            //                }
            //
            //                // Matching BC nodes SOUTH
            //                if (l==0 ) {
            //                    mesh->BCv.type[c] = 0;
            //                    mesh->BCv.val[c]  = 0.0;
            //                }
            //
            //                // Matching BC nodes NORTH
            //                if (l==mesh->Nz-1 ) {
            //                    mesh->BCv.type[c] = 0;
            //                    mesh->BCv.val[c]  = 0.0;
            //                }
            //
            //            }
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
    
    double TN = 273.15/scaling.T, TS = 1330/scaling.T;
    double TW = 273.15/scaling.T, TE = 1330/scaling.T;
    
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
