#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "header_MDOODZ.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SET INITIAL TOPOGRPAHIC MARKER CHAIN
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = -0.0e3/scaling.L; // sets zero initial topography
    double h_pert = model.user3/scaling.L;   
    double xmin   = model.xmin;                                     // xmin
    double xmax   = model.xmax;                                     // xmax
    
    //for ( k=0; k<topo_chain->Nb_part; k++ ) {
    //    topo_chain->z[k]     = TopoLevel;
    //    topo_chain->phase[k] = 0;
    //}
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel + h_pert*(3330.0-2800.0)/2800.0*cos(2*M_PI*topo_chain->x[k]/(xmax-xmin)); //Hplateau = Hroot*(RHOmantle-RHOcrust)/RHOcrust
        topo_chain->phase[k] = 0;
    }
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SET INITIAL MARKER FIELD
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np; 
    double HLit   = model.user1/scaling.L;                          // lithosphere thickness
    double HCrust = model.user2/scaling.L;                          // crust thickness
    double Rad    = model.user0/scaling.L;                          // seed radius
    double h_pert = model.user3/scaling.L;                          // perturbation amplitude
    double Tsurf  = 273.15/scaling.T, Tpart;                        // surface temperature
    double Tmant  = (1330.0+273.15)/scaling.T;                        // adiabatic mantle temperature
    double xmin   = model.xmin;                                     // xmin
    double xmax   = model.xmax;                                     // xmax
    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -particles->x[np]*model.EpsBG;           // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;           // set initial particle velocity (unused)
        particles->phase[np] = 3;                                       // same phase number everywhere
        particles->d[np]     = materials->gs_ref[particles->phase[np]]; // Grain size
        particles->phi[np]   = 0.0;                                     // zero porosity everywhere
        
        //--------------------------//
        // TEMPERATURE - linear geotherm --- Will be overriden
        particles->T[np]       = 0.0;
        Tpart = ((Tmant-Tsurf)/HLit)*(-particles->z[np]) + Tsurf;
        if (Tpart>Tmant) Tpart = Tmant;
        particles->T[np]       = Tpart;
        
        //--------------------------//
        // Phases - lithosphere-asthenosphere
        particles->phase[np] = 2;   //lithos. M. everywhere
        if ( particles->z[np] < -HLit)  particles->phase[np] = 3; // astheno. M. below base lithos
        // Moho_depth equation :  -h_pert*cos(2*pi*x_current/(xmax-xmin));
        if ( particles->z[np] > -HCrust-h_pert*cos(2*M_PI*particles->x[np]/(xmax-xmin))) particles->phase[np] = 1; // crust above the Moho
     
        
        if ( model.user4 == 1 ) {
            if ( particles->z[np] > -HCrust+2.5e3/scaling.L-h_pert*cos(2*M_PI*particles->x[np]/(xmax-xmin)) && particles->z[np] < -HCrust+4.5e3/scaling.L-h_pert*cos(2*M_PI*particles->x[np]/(xmax-xmin))) particles->phase[np] = 7; // marker between 2.5 and 3 km above the Moho
        
            if ( particles->z[np] > -HCrust+7.0e3/scaling.L-h_pert*cos(2*M_PI*particles->x[np]/(xmax-xmin)) && particles->z[np] < -HCrust+9.0e3/scaling.L-h_pert*cos(2*M_PI*particles->x[np]/(xmax-xmin))) particles->phase[np] = 8; // marker between 4.5 and 5 km above the Moho
        }
        
        //if ( particles->z[np] > -HCrust                             )  particles->phase[np] = 1;
        //if ( particles->z[np] < -HCrust && particles->z[np] > -HLit )  particles->phase[np] = 2;
        //if ( particles->z[np] < -HLit                               )  particles->phase[np] = 3;
        //if ( pow( particles->z[np] - z_seed, 2 ) + pow( particles->x[np] - x_seed, 2 ) < Rad*Rad )  { particles->phase[np] = 0; };
       
        
        
        //--------------------------//
        // DENSITY --- Will be overriden
        particles->rho[np] = materials->rho[particles->phase[np]];
        
        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
    }

    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part,  "Tp init" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
    
    int   kk, k, l, c, c1, np;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double TN =  273.15/scaling.T, TS = (1330.+273.15)/scaling.T;
    double TW = (1330.+273.15)/scaling.T, TE = (1330.+273.15)/scaling.T;
    double Tbot, Tleft, Tright;

    // Set unresaonable conductivity in the mantle to generate and adiabatic mantle as initial condition
    int    phase_ast = 3;                     // This has to be the ASTHENOSPHERE phase number
    double k_crazy = 1000*materials->k[phase_ast];

    if (model->step == 0) {
        materials->k_eff[phase_ast] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }
    else {
        materials->k_eff[phase_ast] = materials->k[phase_ast];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }

	// Mantle Impregnation => mantle phase change if Pressure < 1GPa & Temperature > 1050°C 
	      // Loop on particles
//               for( np=0; np<particles->Nb_part; np++ ) {
//                                if (particles->z[np]>-50.0e3/scaling.L && particles->P[np]<1.0e9/scaling.S && particles->T[np]> 1323.0/scaling.T && particles->phase[np]==2) { particles->phase[np] = 6; }
//                                else if (particles->z[np]>-50.0e3/scaling.L && particles->P[np]<1.0e9/scaling.S && particles->T[np]> 1323.0/scaling.T && particles->phase[np]==3) { particles->phase[np] = 6; }
//                                else { ; }
//               }


 
    // Some stuff...
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
                    
                    if ( (k==0 || k==NCX-1) && l==NCZ-1 ) {
                        mesh->BCp.type[c] =  0;
                        mesh->BCp.val[c]  =  0;
                    }
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
                        mesh->BCt.valS[k] = mesh->T[c];
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
