// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2019  MDOODZ Developper team
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

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ExpandCentroidArray( double* CentroidArray, double* temp, grid* mesh, params *model ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    int per = model->isperiodic_x;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Fill interior points
    for (k=0; k<Ncx; k++) {
        for (l=0; l<Ncz; l++) {
            c0 = k + l*(Ncx);
            c1 = k + (l+1)*(Ncx+2) + 1;
//            temp[c1] = CentroidArray[c0];
            if (mesh->BCp.type[c0] == -1) temp[c1] = CentroidArray[c0];
            if (mesh->BCp.type[c0] == 31) temp[c1] = CentroidArray[c0-Ncx];
        }
    }
    
    // Fill sides - avoid corners - assume zero flux
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (0)*(Ncx+2);       // South
        c1 = k + (1)*(Ncx+2);       // up neighbour
        temp[c0] = temp[c1];
    }
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (Ncz+1)*(Ncx+2);   // North
        c1 = k + (Ncz  )*(Ncx+2);   // down neighbour
        temp[c0] = temp[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = 0 + (l)*(Ncx+2);       // West
        if (per == 0) c1 = 1           + (l)*(Ncx+2);       // right neighbour
        if (per == 1) c1 = (Ncx+2-1-1) + (l)*(Ncx+2);       // right neighbour
        temp[c0] = temp[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = (Ncx+1) + (l)*(Ncx+2); // East
        if (per==0) c1 = (Ncx  ) + (l)*(Ncx+2); // left neighbour
        if (per==1) c1 = 1       + (l)*(Ncx+2); // left neighbour
        temp[c0] = temp[c1];
    }
    
    // Corners - assume zero flux
    c0 = (0) + (0)*(Ncx+2);         // South-West
    if (per==0) c1 = (1) + (1)*(Ncx+2);         // up-right neighbour
    if (per==1) c1 = (0) + (1)*(Ncx+2);         // up       neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (0)*(Ncx+2);     // South-East
    if (per==0) c1 = (Ncx  ) + (1)*(Ncx+2);     // up-left neighbour
    if (per==1) c1 = (Ncx+1) + (1)*(Ncx+2);     // up      neighbour
    temp[c0] = temp[c1];
    c0 = (0) + (Ncz+1)*(Ncx+2);     // North-West
    if (per==0) c1 = (1) + (Ncz  )*(Ncx+2);     // down-right neighbour
    if (per==1) c1 = (0) + (Ncz  )*(Ncx+2);     // down       neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (Ncz+1)*(Ncx+2); // North-West
    if (per==0) c1 = (Ncx  ) + (Ncz  )*(Ncx+2); // down-left neighbour
    if (per==1) c1 = (Ncx+1) + (Ncz  )*(Ncx+2); // down      neighbour
    temp[c0] = temp[c1];
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InterpCentroidsToVerticesDouble( double* CentroidArray, double* VertexArray, grid* mesh, params *model ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;    
    int per = model->isperiodic_x;
    
    Nx  = mesh->Nx;
    Nz  = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Allocate temporary swelled centroid of size = (ncx+2) * (ncz+2)
    DoodzFP* temp      = DoodzCalloc((Ncx+2)*(Ncz+2), sizeof(DoodzFP));
    int*     temp_flag = DoodzCalloc((Ncx+2)*(Ncz+2), sizeof(int));
    
    // Fill interior points
    for (k=0; k<Ncx; k++) {
        for (l=0; l<Ncz; l++) {
            c0 = k + l*(Ncx);
            c1 = k + (l+1)*(Ncx+2) + 1;
            temp[c1]      = CentroidArray[c0];
            temp_flag[c1] = mesh->BCp.type[c0];
        }
    }
    
    // Fill sides - avoid corners - assume zero flux
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (0)*(Ncx+2);       // South
        c1 = k + (1)*(Ncx+2);       // up neighbour
        temp[c0]      = temp[c1];
        temp_flag[c0] = temp_flag[c1];
    }
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (Ncz+1)*(Ncx+2);   // North
        c1 = k + (Ncz  )*(Ncx+2);   // down neighbour
        temp[c0]      = temp[c1];
        temp_flag[c0] = temp_flag[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = 0 + (l)*(Ncx+2);       // West
        if (per == 0) c1 = 1           + (l)*(Ncx+2);       // right neighbour
        if (per == 1) c1 = (Ncx+2-1-1) + (l)*(Ncx+2);       // right neighbour
        temp[c0]      = temp[c1];
        temp_flag[c0] = temp_flag[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = (Ncx+1) + (l)*(Ncx+2); // East
        if (per==0) c1 = (Ncx  ) + (l)*(Ncx+2); // left neighbour
        if (per==1) c1 = 1       + (l)*(Ncx+2); // left neighbour
        temp[c0]      = temp[c1];
        temp_flag[c0] = temp_flag[c1];
    }
    
    // Corners - assume zero flux
    c0 = (0) + (0)*(Ncx+2);         // South-West
    if (per==0) c1 = (1) + (1)*(Ncx+2);         // up-right neighbour
    if (per==1) c1 = (0) + (1)*(Ncx+2);         // up       neighbour
    temp[c0]      = temp[c1];
    temp_flag[c0] = temp_flag[c1];
    c0 = (Ncx+1) + (0)*(Ncx+2);     // South-East
    if (per==0) c1 = (Ncx  ) + (1)*(Ncx+2);     // up-left neighbour
    if (per==1) c1 = (Ncx+1) + (1)*(Ncx+2);     // up      neighbour
    temp[c0]      = temp[c1];
    temp_flag[c0] = temp_flag[c1];
    c0 = (0) + (Ncz+1)*(Ncx+2);     // North-West
    if (per==0) c1 = (1) + (Ncz  )*(Ncx+2);     // down-right neighbour
    if (per==1) c1 = (0) + (Ncz  )*(Ncx+2);     // down       neighbour
    temp[c0]      = temp[c1];
    temp_flag[c0] = temp_flag[c1];
    c0 = (Ncx+1) + (Ncz+1)*(Ncx+2); // North-West
    if (per==0) c1 = (Ncx  ) + (Ncz  )*(Ncx+2); // down-left neighbour
    if (per==1) c1 = (Ncx+1) + (Ncz  )*(Ncx+2); // down      neighbour
    temp[c0]      = temp[c1];
    temp_flag[c0] = temp_flag[c1];
    
#pragma omp parallel for shared( temp, temp_flag, VertexArray, mesh ) private( k,l,c1, c0 )  firstprivate( Ncx,Ncz, Nx, Nz )
    // interpolate from temp array to actual vertices array
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            
            c1 = k + (l)*(Nx);
            c0 = k + (l+1)*(Ncx+2) + 1;
            
            // Default 0 - above free surface
            VertexArray[c1] = 0.0;
            
            // Else interpolate
            if ( mesh->BCg.type[c1] != 30 ) {
                int inSW = 1, inSE = 1, inNW = 1, inNE = 1;
                if (temp_flag[c0-1-(Ncx+2)] == 30 || temp_flag[c0-1-(Ncx+2)] == 31) inSW = 0; 
                if (temp_flag[c0-0-(Ncx+2)] == 30 || temp_flag[c0-0-(Ncx+2)] == 31) inSE = 0;
                if (temp_flag[c0-1]         == 30 || temp_flag[c0-1]         == 31) inNW = 0;
                if (temp_flag[c0-0]         == 30 || temp_flag[c0-0]         == 31) inNE = 0;
                VertexArray[c1] = 0.25*( inSW*temp[c0-1-(Ncx+2)] + inSE*temp[c0-0-(Ncx+2)] + inNW*temp[c0-1] + inNE*temp[c0] );
  
            }
        }
    }
    
    // Free temporary array
    DoodzFree(temp);
    DoodzFree(temp_flag);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InterpVerticesToCentroidsDouble( double* CentroidArray, double* VertexArray, grid* mesh, params *model ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    double *temp;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
#pragma omp parallel for shared( CentroidArray, VertexArray, mesh ) private( k,l,c1, c0 )  firstprivate( Ncx,Ncz, Nx )
    // Fill interior points
    for (k=0; k<Ncx; k++) {
        for (l=0; l<Ncz; l++) {
            c0 = k + l*(Ncx);
            c1 = k + l*(Nx);
            
            if ( mesh->BCp.type[c0] != 30 &&  mesh->BCp.type[c0] != 31 ) {
                CentroidArray[c0] = 0.25*( VertexArray[c1] + VertexArray[c1+1] + VertexArray[c1+Nx] + VertexArray[c1+1+Nx] );
            }
        }
    }
    
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Initialise grid
void SetGridCoordinates( grid *mesh, params *model, int nx, int nz ) {
    
    // Set coordinates of mesh vertices and cell centroids, Vx and Vz nodes
	
    int i;
    double dx0, dz0;
    
    // Set spacing
    mesh->dx = (model->xmax-model->xmin) / (nx-1);
    mesh->dz = (model->zmax-model->zmin) / (nz-1);
    
    model->dx = mesh->dx;
    model->dz = mesh->dz;
	
    // Gridlines positions
    mesh->xg_coord[0] = model->xmin;
    mesh->zg_coord[0] = model->zmin;
	
    for (i = 1; i<nx; i++) {
        mesh->xg_coord[i] = mesh->xg_coord[i-1] + mesh->dx;
    }
	
    for (i = 1; i<nz; i++) {
        mesh->zg_coord[i] = mesh->zg_coord[i-1] + mesh->dz;
    }
    
    //---------------
    
    // Initial mesh
    // Set spacing
    dx0 = (model->xmax0-model->xmin0) / (nx-1);
    dz0 = (model->zmax0-model->zmin0) / (nz-1);
    
    // Gridlines positions
    mesh->xg_coord0[0] = model->xmin0;
    mesh->zg_coord0[0] = model->zmin0;
    
    for (i = 1; i<nx; i++) {
        mesh->xg_coord0[i] = mesh->xg_coord0[i-1] + dx0;
    }
    
    for (i = 1; i<nz; i++) {
        mesh->zg_coord0[i] = mesh->zg_coord0[i-1] + dz0;
    }
    
    //---------------
    
    // Cell centers positions
    mesh->xc_coord[0] = model->xmin + mesh->dx/2.0;
    mesh->zc_coord[0] = model->zmin + mesh->dz/2.0;
    
    for (i = 1; i<nx-1; i++) {
        mesh->xc_coord[i]  = mesh->xc_coord[i-1]  + mesh->dx;
    }
    for (i = 1; i<nz-1; i++) {
        mesh->zc_coord[i]  = mesh->zc_coord[i-1]  + mesh->dz;
    }
    
    // Ghosted nodes
    mesh->xvz_coord[0] = model->xmin - mesh->dx/2.0;
    mesh->zvx_coord[0] = model->zmin - mesh->dz/2.0;
    
    for (i = 1; i<nx+1; i++) {
        mesh->xvz_coord[i] = mesh->xvz_coord[i-1] + mesh->dx;
    }
    for (i = 1; i<nz+1; i++) {
        mesh->zvx_coord[i] = mesh->zvx_coord[i-1] + mesh->dz;
    }
    
    // Ghosted nodes
    mesh->xg_coord_ext[0] = model->xmin - mesh->dx;
    mesh->zg_coord_ext[0] = model->zmin - mesh->dz;
    
    for (i = 1; i<nx+2; i++) {
        mesh->xg_coord_ext[i] = mesh->xg_coord_ext[i-1] + mesh->dx;
    }
    for (i = 1; i<nz+2; i++) {
        mesh->zg_coord_ext[i] = mesh->zg_coord_ext[i-1] + mesh->dz;
    }
	
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseSolutionFields( grid *mesh, params *model ) {
    
    // Set initial velocity and pressure fields
	
    int nx, nz, nxvz, nzvx, ncx, ncz;
    int k, l, c;
    double eps = 1e-13; // perturbation to avoid zero pressure that results in Nan d(eta)dP in numerical differentiation
	
    nx   = mesh->Nx;
    nz   = mesh->Nz;
    nxvz = nx+1;
    nzvx = nz+1;
    ncx  = nx-1;
    ncz  = nz-1;

    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            
            c = k + l*nx;
            
//            // This works with cylindrical
//            mesh->u_in[c]  = 0.0;
//            if (mesh->BCu.type[c] == 0) mesh->u_in[c]  = mesh->BCu.val[c];
            
            if ( mesh->BCu.type[c] == 30 )  mesh->u_in[c]  = 0.0;
            else {        
                if (model->step==0) {
                    // Initial velocity field (zero or pure shear)
                    if (model->EpsBG == 0) mesh->u_in[c]  = 0.0;
                    // Pure shear
                    else mesh->u_in[c]  = -mesh->xg_coord[k]*(model->EpsBG - model->DivBG/3.0);
                    if (model->isperiodic_x == 1) mesh->u_in[c] = 2.0*(mesh->zvx_coord[l])*model->EpsBG + mesh->xg_coord[k]*model->DivBG/3.0; // Simple shear
                }
                // Force Dirichlets
                if (mesh->BCu.type[c] == 0) mesh->u_in[c]  = mesh->BCu.val[c]; 
            }
    
        }
    }
	
    // Initial velocity field Vz
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            
            c = k + l*nxvz;
            
            if ( mesh->BCv.type[c] == 30 )  mesh->v_in[c]  = 0.0;
            else {
                if (model->step==0) {
                    // Initial velocity field (zero or pure shear)
                    if (model->EpsBG == 0) mesh->v_in[c]  = 0.0;
                    else mesh->v_in[c]  = mesh->zg_coord[l]*(model->EpsBG + model->DivBG/3.0);
                    if (model->isperiodic_x == 1) mesh->v_in[c]  = 0.0 + mesh->zg_coord[l]*model->DivBG/3.0;
                }
                // Force Dirichlets
                if (mesh->BCv.type[c] == 0) mesh->v_in[c]  = mesh->BCv.val[c];
            }
            
        }
    }
    
    // Initial pressure field
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            c = k + l*ncx;               
            if ( mesh->BCp.type[c] == 30 ||  mesh->BCp.type[c] == 31 ) mesh->p_in[c]  = 0.0;
            if ( mesh->BCp.type[c] != 30 ||  mesh->BCp.type[c] != 31 ) {
                if (model->step==0) {
                    mesh->p_in[c]  = 0.0 + model->PrBG;
                }
            }
        }
    }

    // Very important step: Set BC's here!!!!!!!
    ApplyBC( mesh, model ); 
    printf("Velocity field was set to background pure shear\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ComputeLithostaticPressure( grid *mesh, params *model, double RHO_REF, scale scaling, int mode ) {
	
    // Compute lithostatic pressure by cumulative sum of rho*g across model thickness
    
    int nx, nz, ncx, ncz;
    int k, l, c, c1;
    double rho_eff;
    double eps = 1e-13; // perturbation to avoid zero pressure that results in Nan d(eta)dP in numerical differentiation

    
    nx  = mesh->Nx;
    nz  = mesh->Nz;
    ncx = nx-1;
    ncz = nz-1;
        
    Initialise1DArrayDouble( mesh->p_lith,  (mesh->Nx-1)*(mesh->Nz-1), 0.0 );
    
    // Cell center arrays
    for( l=ncz-2; l>=0; l--) {
        for( k=0; k<ncx; k++) {

            // Initialise vertices variables
            c  = k + l*ncx;

            mesh->p_lith[c]  = 0.0;

            // density
            if ( mode == 0 ) rho_eff = RHO_REF;
            if ( mode == 1 ) rho_eff = mesh->rho_n[c];

            // Initialise pressure variables : Compute lithostatic pressure
            if ( mesh->BCp.type[c] != 30 && mesh->BCp.type[c] != 31 ) { // First row (surface)
                mesh->p_lith[c]  = mesh->p_lith[c+ncx] -  model->gz * mesh->dz * rho_eff;

            }
        }
    }

    // Add confining pressure
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            c  = k + l*ncx;
            // density
            if ( mode == 0 ) rho_eff = RHO_REF;
            if ( mode == 1 ) rho_eff = mesh->rho_n[c];
            
            if ( mesh->BCp.type[c] != 30 && mesh->BCp.type[c] != 31 ) {
                mesh->p_lith[c] += model->PrBG;
                mesh->p_lith[c]  += 0.5*model->gz * mesh->dz * rho_eff;
            }
        }
    }
    
    
//    // + eps
//    // Cell center arrays
//    for( l=0; l<ncz; l++) {
//        for( k=0; k<ncx; k++) {            
//            c  = k + l*ncx;
//            mesh->p_lith[c] += eps;
//            mesh->p[c]       = mesh->p_in[c];
//        }
//    }
    

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GridIndices( grid* mesh) {
    
// Compute i, j indices of flattened arrays (vertices, centroids, Vx, Vz)
    
    int k,l,kk;
    
    // Cell vertices
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*(mesh->Nx);
            mesh->kn[kk] = k;
            mesh->ln[kk] = l;
        }
    }
    
    // Vx nodes
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            mesh->kvx[kk] = k;
            mesh->lvx[kk] = l;
        }
    }
    
    // Vz nodes
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            
            mesh->kvz[kk] = k;
            mesh->lvz[kk] = l;
        }
    }
    
    // Cell Centroids
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            mesh->kp[kk] = k;
            mesh->lp[kk] = l;
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetUpModel_NoMarkers ( grid* mesh, params *model, scale *scaling ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, k1, c0, c1;
    double x, z;
    double radius = model->user1/scaling->L;
    double z0     = 0.5*(model->zmax + model->zmin );
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    printf("Setting up mode without using markers --- DEBUG !!!!\n");
    
    // Vertices
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k + l*(Ncx);
        
        x = mesh->xc_coord[k];
        z = mesh->zc_coord[l] - z0;
        
        mesh->T[k1] = 0.05;
        
        mesh->phase_perc_n[0][k1] = 1.0;
        mesh->phase_perc_n[1][k1] = 0.0;
        
        if (x*x + z*z < radius*radius) {
            mesh->phase_perc_n[0][k1] = 0.0;
            mesh->phase_perc_n[1][k1] = 1.0;
        }
    }
    
    
    
    // Vertices
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        
        x = mesh->xg_coord[k];
        z = mesh->zg_coord[l];
        
        mesh->phase_perc_s[0][k1] = 1.0;
        mesh->phase_perc_s[1][k1] = 0.0;
        
        if (x*x + z*z < radius*radius) {
            mesh->phase_perc_s[0][k1] = 0.0;
            mesh->phase_perc_s[1][k1] = 1.0;
        }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Diffuse_X( grid* mesh, params* model, scale* scaling ) {
    
    double K = 1.0e-6 / (pow(scaling->L,2.0)/scaling->t);
    double dx = mesh->dx, dz = mesh->dz;
    double L_diff = model->diffusion_length;
    double dt = pow(MAXV(dx,dz), 2.0) / K / 4.1;
    double t_diff = pow(L_diff,2.0)/K;
    int    nsteps = (int) (t_diff/dt);
    int    k, l, it, c;
    int    ncx = mesh->Nx-1;
    int    ncz = mesh->Nz-1;
    double *X0, qE, qW, qN, qS;
    int * flag;
    X0   = DoodzCalloc(sizeof(DoodzFP), ncx*ncz);
    flag = DoodzCalloc(sizeof(int), ncx*ncz);
    printf("Diffusion length = %2.2e  Diffusion time = %2.2e dt = %2.2e nsteps = %03d\n", L_diff*scaling->L, t_diff*scaling->t, dt*scaling->t, nsteps);
    
    
    for (l=1;l<ncz-1;l++) {
        
        for (k=1;k<ncx-1;k++) {
            c = k + l*ncx;
            
            if ( mesh->X_n[c]>0.99 ) flag[c] = 1;
        }
    }
    
    
    for (it=0;it<nsteps;it++) {
        ArrayEqualArray( X0, mesh->X0_n, ncx*ncz );
        
        // !!!! BCs are not included!
        for (l=1;l<ncz-1;l++) {
            
            for (k=1;k<ncx-1;k++) {
                c = k + l*ncx;
                
                if (flag[c]!=1) {
                    
                    qW = - K*(X0[c]-X0[c-1])/dx;
                    qE = - K*(X0[c+1]-X0[c])/dx;
                    qS = - K*(X0[c]-X0[c-ncz])/dz;
                    qN = - K*(X0[c+ncz]-X0[c])/dz;
                    mesh->X_n[c] = X0[c] - dt* ( 1.0/dx*(qE-qW) + 1.0/dz*(qN-qS));
                }
            }
            
        }
    }
    MinMaxArrayTag( mesh->X_n,   1.0,    (mesh->Nx-1)*(mesh->Nz-1),   "Xreac_n",   mesh->BCp.type );
    
    DoodzFree(X0);
    DoodzFree(flag);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
