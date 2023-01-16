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
#include"head.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, Mgrid Mmesh, scale scaling ) {
    
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
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
   
    int np;
    FILE *read;
    int s1=sizeof(int), s2=sizeof(double);
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    double T_init = (model.user0 + zeroC)/scaling.T;
    
    // Amplitude of random noise
    double **pert_lo0, **pert_up0, **pert_lo, **pert_up;
    double A = 100.0e-2/scaling.L;
    int    ncx = model.Nx-1, ic, it;
    double H = model.user2/scaling.L;
    double z_up0= H/2, z_lo0=-H/2;
    
    // Layers
    int nlayers = 9, il;
    double layers[9] = { -6.5, -5.0, -3.5, -2.0, -0.5, 1.0, 2.5, 4.0, 5.5};
    double layer_bottom, spacing;
    spacing = (Lz - nlayers*H - 2*6/scaling.L)/8;
 
    // Allocate cellwise contours perturbations
    pert_lo0 = DoodzCalloc( nlayers, sizeof(double*) );
    pert_up0 = DoodzCalloc( nlayers, sizeof(double*) );
    pert_lo  = DoodzCalloc( nlayers, sizeof(double*) );
    pert_up  = DoodzCalloc( nlayers, sizeof(double*) );
    for ( il=0; il<nlayers; il++ ) {
        pert_lo0[il] = DoodzCalloc( ncx, sizeof(double) );
        pert_up0[il] = DoodzCalloc( ncx, sizeof(double) );
        pert_lo[il]  = DoodzCalloc( ncx, sizeof(double) );
        pert_up[il]  = DoodzCalloc( ncx, sizeof(double) );
    }

    
    // Generation of perturbed layer contours on the grid
    if ((int)model.user1 == 0) {
        
        printf("Layer thickness = %2.2e m \n", H*scaling.L);
        
        
        for ( il=0; il<nlayers; il++ ) {
            
            // Generate perturbed contours
            for( ic=30; ic<ncx-30; ic++ ) {
                pert_lo[il][ic] = 1*A*(double)rand() / RAND_MAX - A/2;
                pert_up[il][ic] = 1*A*(double)rand() / RAND_MAX - A/2;
            }
            
            // Smooth contours with diffusion
            double D = 0.4;
            for( it=0; it<ncx; it++ ) {
                for( ic=0; ic<ncx-0; ic++ ) {
                    pert_lo0[il][ic] = pert_lo[il][ic];
                    pert_up0[il][ic] = pert_up[il][ic];
                }
                for( ic=1; ic<ncx-1; ic++ ) {
                    pert_lo[il][ic] = pert_lo0[il][ic] + D*(pert_lo0[il][ic-1] + pert_lo0[il][ic+1] - 2*pert_lo0[il][ic]);
                    pert_up[il][ic] = pert_up0[il][ic] + D*(pert_up0[il][ic-1] + pert_up0[il][ic+1] - 2*pert_up0[il][ic]);
                }
            }
        }
        
    }

    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;                   // set initial particle velocity (unused)
//        particles->phase[np] = 0;                                               // same phase number everywhere
        particles->d[np]     = 0;                                            // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->rho[np]   = 0;
        particles->T[np]     = T_init;
    
        // ------------------------- //
        // DRAW LAYER
        if ((int)model.user1 == 0) {
            
            for ( il=0; il<nlayers; il++ ) {
                
                layer_bottom = model.zmin + 6/scaling.L + il*(H+spacing);
                
                // Get the column:
                ic   = ceil( (( particles->x[np] - (model.xmin+model.dx/2) )/model.dx) + 0.5) - 1;
                if (ic<0)     ic = 0;
                if (ic>ncx-1) ic = ncx-1;
                // Corresponding layer contour perturbation
                if (particles->z[np]>layer_bottom+pert_lo[il][ic] && particles->z[np]<layer_bottom+ H+pert_up[il][ic] ) particles->phase[np] = il+1;
                
            }
        }

        
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
    
    if ((int)model.user1 == 0) {
        printf("WRITE\n");
        // Write the Doudou file
        if (fopen(model.input_file, "wb")!=NULL){
            read = fopen(model.input_file, "wb");
            printf("Writing %d particles from file %s...\n", particles->Nb_part, model.input_file );
        }
        else {
            printf("Cannot open file %s, check if the file exists in the current location !\n Exiting", model.input_file);
            exit(1);
        }
        fwrite( particles->x, s2, particles->Nb_part, read);
        fwrite( particles->z, s2, particles->Nb_part, read);
        fwrite( particles->phase, s1, particles->Nb_part, read);
        fclose(read);
    }
    
    for ( il=0; il<nlayers; il++ ) {
        DoodzFree(pert_lo0[il]);
        DoodzFree(pert_up0[il]);
        DoodzFree(pert_lo[il] );
        DoodzFree(pert_up[il] );
    }
    DoodzFree(pert_lo0);
    DoodzFree(pert_up0);
    DoodzFree(pert_lo );
    DoodzFree(pert_up );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void MSetFields( Mgrid *mesh, params *model, Mparams *Mmodel, scale scaling, markers* particles, mat_prop *materials ) {
    
    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    
    if (model->step >= 1){
        materials->k[4] = materials->k[3];
        printf("Running with normal conductivity in the asthenosphere!\n");
    }
    
    
    // ---- T-Dependent marker types
    // -------------------- SPECIFIC TO YOANN's SETUP -------------------- //
    
    NX  = mesh->Nx[0];
    NZ  = mesh->Nz[0];
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    X  = malloc (NX*sizeof(double));
    Z  = malloc (NZ*sizeof(double));
    XC = malloc (NCX*sizeof(double));
    ZC = malloc (NCZ*sizeof(double));
    
    for (k=0; k<NX; k++) {
        X[k] = mesh->xg_coord[0][k];
    }
    for (k=0; k<NCX; k++) {
        XC[k] = mesh->xc_coord[0][k];
    }
    for (l=0; l<NZ; l++) {
        Z[l] = mesh->zg_coord[0][l];
    }
    for (l=0; l<NCZ; l++) {
        ZC[l] = mesh->zc_coord[0][l];
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
	
	for(kk=0; kk<Mmodel->n_level; kk++) {
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
        
		for (l=0; l<mesh->Nz[kk]+1; l++) {
			for (k=0; k<mesh->Nx[kk]; k++) {
				
				c = k + l*(mesh->Nx[kk]);
                
                if ( mesh->BCu.type[kk][c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[kk][c] = -1;
                    mesh->BCu.val[kk][c]  =  0;
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[kk][c] = 0;
                        mesh->BCu.val[kk][c]  = -mesh->xg_coord[0][k] * model->EpsBG;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx[kk]-1 ) {
                        mesh->BCu.type[kk][c] = 0;
                        mesh->BCu.val[kk][c]  = -mesh->xg_coord[0][k] * model->EpsBG;
                    }
                    
                    // Free slip SOUTH
                    if (l==0  ) {
                        mesh->BCu.type[kk][c] = 13;
                        mesh->BCu.val[kk][c]  =  0;
                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz[kk] ) {
                        mesh->BCu.type[kk][c] = 13;
                        mesh->BCu.val[kk][c]  =  0;
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
	
	for(kk=0; kk<Mmodel->n_level; kk++) {
		
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz[kk]; l++) {
			for (k=0; k<mesh->Nx[kk]+1; k++) {
				
				c  = k + l*(mesh->Nx[kk]+1);
                
                if ( mesh->BCv.type[kk][c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCv.type[kk][c] = -1;
                    mesh->BCv.val[kk][c]  =  0;
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[kk][c] = 0;
                        mesh->BCv.val[kk][c]  = mesh->zg_coord[0][l] * model->EpsBG;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz[kk]-1 ) {
                        mesh->BCv.type[kk][c] = 0;
                        mesh->BCv.val[kk][c]  = mesh->zg_coord[0][l] * model->EpsBG;
                    }
                    
                    // Non-matching boundary WEST
                    if ( (k==0) ) {
                        mesh->BCv.type[kk][c] =   13;
                        mesh->BCv.val[kk][c]  =   0;
                    }
                    
                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx[kk]) ) {
                        mesh->BCv.type[kk][c] =   13;
                        mesh->BCv.val[kk][c]  =   0;
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
    
    for(kk=0; kk<Mmodel->n_level; kk++) {
        
        NX  = mesh->Nx[kk];
        NZ  = mesh->Nz[kk];
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {
                
                c  = k + l*(NCX);
                
                if (mesh->BCt.type[c] != 30) {
                            
                    // Internal points:  -1
                    mesh->BCp.type[kk][c] = -1;
                    mesh->BCp.val[kk][c]  =  0;
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
    
    double Ttop = 273.15/scaling.T;
    double Tbot, Tleft, Tright;
    
	
	for(kk=0; kk<1; kk++) {
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz[kk]-1; l++) {
			for (k=0; k<mesh->Nx[kk]-1; k++) {
				
				c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
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
