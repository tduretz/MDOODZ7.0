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
    double TopoLevel = 0.75; // sets zero initial topography
    
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
	
    int k;

    // Define dimensions of box;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    double x0 =0.0, z0 =-0*0.0*Lz ;
    double Smajor = model.user1 / scaling.L, Sminor = model.user2 / scaling.L, theta = (360-model.user3) * PI / 180;
    double xp, zp;
    double a , b, crit;
    
    // Loop over particles
    for( k=0; k<particles->Nb_part; k++ ) {
        
        particles->phase[k] = 0;   // integer --- phase type
        particles->Vx[k]    = 0.0;
        particles->Vz[k]    = 0.0;
        particles->T[k]     = 0.0;
        particles->phi[k]   = 0.0;
        particles->X[k]     = 0.0;
        
        xp = particles->x[k];
        zp = particles->z[k];

        // Check if particle lies within the ellipse. General eq.: a^2/Smajor^2 + b^2/Sminor^2 <= 1
        a = (cos(theta)*(xp - x0) + sin(theta)*(zp - z0));
        b = (sin(theta)*(xp - x0) - cos(theta)*(zp - z0));
        crit = (a*a)/(Smajor*Smajor) + (b*b)/(Sminor*Sminor);
        if (crit <= 1) {
            if ( model.diffuse_X == 0 ) particles->phase[k] = 1;   // Sharp boundary
            if ( model.diffuse_X == 1 ) particles->X[k]     = 1.0; // Smooth boundary
        }
        
        // DO THIS LAST
        particles->d[k]     = materials->gs_ref[particles->phase[k]];
        
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
	
	int   kk, k, l, c, c1;
	double *X, *Z, *XC, *ZC;
	int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
	double dmin, VzBC, eta = 1e4 / scaling.eta, Ttop=273.15/scaling.T;//(scaling.T+273)/scaling.T;
    double Lx, Lz;
    double Stress = -7e6/scaling.S;
    double vx, vz, p, m, sxx, syy, dx=model->dx, dz=model->dz, f=0.0;
    
    double Tfix = (model->user4+273.15) / scaling.T;
    
    int StressBC_W = 0;
    int StressBC_E = 0;
    int StressBC_S = 0;
    int StressBC_N = 0;
    
    if ( model->user0 == 1 ) {
        StressBC_W = 1;
        StressBC_E = 1;
        StressBC_S = 1;
        StressBC_N = 1;
    }
	
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
    
    // Fix temperature
    for (l=0; l<mesh->Nz-1; l++) {
        for (k=0; k<mesh->Nx-1; k++) {
            c = k + l*(NCX);
            mesh->T[c]     = Tfix;
        }
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
                
                c  = k + l*(mesh->Nx);
                
                
                if ( mesh->BCu.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[c] = -1;
                    mesh->BCu.val[c]  =  0;
                    
                    // Matching BC nodes WEST
                    if (k==0 && l>0 && l<mesh->Nz) {
                        
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k], mesh->zvx_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
//                        mesh->BCu.val[c]  = vx;
                        
//                        if ( StressBC_W==1 ) {
//                            eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k]+f*dx, mesh->zvx_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                            mesh->BCu.type[c] = 2;
//                            mesh->BCu.val[c]  = sxx;
//                        }
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 && l>0 && l<mesh->Nz) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k], mesh->zvx_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1]);
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
//                        mesh->BCu.val[c]  = vx;
                        
//                        if ( StressBC_E==1 ) {
//                            eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k]-f*dx, mesh->zvx_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                            mesh->BCu.type[c] = 2;
//                            mesh->BCu.val[c]  = sxx;
//                        }
                    }
                
//                    // Free slip SOUTH
//                    if ( l==0  ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k], model->zmin, 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                        mesh->BCu.type[c] =  11;
//                        mesh->BCu.val[c]  =  vx;
//
//                    }
//
//                    // Free slip NORTH
//                    if ( l==mesh->Nz ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xg_coord[0][k],  model->zmax, 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                        mesh->BCu.type[c] =  11;
//                        mesh->BCu.val[c]  =  vx;
//
//                    }
                    // Free slip SOUTH
                    if ( l==0  ) {
                        mesh->BCu.type[c] =  13;
                        mesh->BCu.val[c]  =  0.0;

                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] =  13;
                        mesh->BCu.val[c]  =  0.0;

                    }
                
                }
                
            }
        }
        
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    

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
                    if (l==0 && k>0 && k<mesh->Nx) {
//                       eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xvz_coord[0][k], mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
//                        mesh->BCv.val[c]  = vz;
                        
//                        if ( StressBC_S==1 ) {
//                            eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xvz_coord[0][k], mesh->zg_coord[0][l]+f*dz, 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                            mesh->BCv.type[c] = 2;
//                            mesh->BCv.val[c]  = syy;
//                        }


                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 && k>0 && k<mesh->Nx) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xvz_coord[0][k], mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1]);
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
//                        mesh->BCv.val[c]  = vz;
                        
//                        if ( StressBC_N==1 ) {
//                            eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xvz_coord[0][k], mesh->zg_coord[0][l]-f*dz, 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                            mesh->BCv.type[c] = 2;
//                            mesh->BCv.val[c]  = syy;
//                        };
                    }
                    
//                    // Non-matching boundary WEST
//                    if ( k==0 ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, model->xmin, mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                        mesh->BCv.type[c] =   11;
//                        mesh->BCv.val[c]  =   vz;
//                    }
//
//                    // Non-matching boundary EAST
//                    if ( k==mesh->Nx ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, model->xmax, mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                        mesh->BCv.type[c] =   11;
//                        mesh->BCv.val[c]  =   vz;
//
//                    }
                    
                    // Non-matching boundary WEST
                    if ( k==0 ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, model->xmin, mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;
                    }
                    
                    // Non-matching boundary EAST
                    if ( k==mesh->Nx ) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, model->xmax, mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;

                    }
                    
//                    // Normal stress applied to the free surface
//                    if (l<mesh->Nz-1 && k>0 && k<mesh->Nx) {
//                        eval_anal_Dani( &vx, &vz, &p, &m, &sxx, &syy, mesh->xvz_coord[0][k], mesh->zg_coord[0][l], 1, model->user1, materials->eta0[0], materials->eta0[1] );
//                        c1 = k + (l)*(mesh->Nx-1)-1;
//                       if (mesh->BCp.type[0][c1]==31) mesh->BCv.val[c]  = syy;
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
    
    for(kk=0; kk<1; kk++) {
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

