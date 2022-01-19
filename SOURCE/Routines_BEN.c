#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
//---- M-Doodz header file
#include "header_MDOODZ.h"
#include "umfpack.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateCourantCriterion_BEN( double* Vx, double* Vz, params *model, scale scaling,  grid* mesh, int quiet ) {
    
    int k, l, c;
    double minVx=0.0, minVz=0.0, maxVx=0.0, maxVz=0.0, dmin, dtc=0, vmax;
    double C = model->Courant;
    
    for (k=0; k<model->Nx; k++) {
        for (l=0; l<model->Nz+1; l++) {
            c = k + l*model->Nx;
            maxVx = MAXV(maxVx, (Vx[c]));
            minVx = MINV(minVx, (Vx[c]));
        }
    }
    
    for (k=0; k<model->Nx+1; k++) {
        for (l=0; l<model->Nz; l++) {
            c = k + l*(model->Nx+1);
            maxVz = MAXV(maxVz, (Vz[c]));
            minVz = MINV(minVz, (Vz[c]));
        }
    }
    if (quiet==0) printf("Min Vxm = %2.2e m/s / Max Vxm = %2.2e m/s\n", minVx * scaling.V, maxVx * scaling.V);
    if (quiet==0) printf("Min Vzm = %2.2e m/s / Max Vzm = %2.2e m/s\n", minVz * scaling.V, maxVz * scaling.V);
    
    dmin = MINV(model->dx, model->dz);
    vmax = MAXV(maxVx, maxVz);
    
    if (model->dt_constant == 0) {
        
        double fact;
        
        // Timestep increase factor
        if ( model->iselastic == 1 ) {
            fact = 1.05;
        }
        else {
            fact = 2.0;
        }
        
        // Courant dt
        dtc = C * dmin / fabs(vmax);
        
        // Timestep cutoff : Do not allow for very large timestep increase
        if (dtc > fact*model->dt0 ) {
            dtc = fact*model->dt0;
        }
        
        // If timestep is adaptive
        if ( model->dt_constant != 1 ) {
            model->dt = dtc;
        }
        
        // If there is no motion, then the timestep becomes huge: cut off the motion.
        if( model->dt>1e30 ) {
            dtc = 0;
            model->dt = model->dt_start;
        }
        
        // Cutoff infinitely slow motion
        if (vmax<1e-30) {
            dtc = 0;
            model->dt = model->dt_start;
        }
        
       
        if (quiet==0) printf("Current dt = %2.2e s / Courant dt = %2.2e s\n", model->dt * scaling.t, dtc * scaling.t );
    }
    else {
        if (quiet==0) printf("Fixed timestep dt = %2.2e s\n", model->dt * scaling.t );
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/



void FindClosestPhase_BEN( markers* particles, int ic, int jc, grid mesh, int *ind_list, int new_ind ) {
    
    // Find phase of a new particle according to its closest neighbour. This function works at the scale of a cell.
    // The array 'ind_list' is list of indices of particle already owing to that cell.
    
    double dst, min_dst;
    int ind, icell, min_index=ind_list[0];
    
    min_dst = ( pow( mesh.dx, 2) + pow( mesh.dz, 2) );
    icell = ic + jc * (mesh.Nx-1);
    
    // Loop on the particles inside the current cell, find which one is the closest to the newly created one.
    for (ind=0; ind<mesh.nb_part_cell[icell]; ind++) { //!
        dst = ( pow (particles->x[new_ind] - particles->x[ind_list[ind]], 2) + pow(particles->z[new_ind] - particles->z[ind_list[ind]], 2) );
        if (dst<min_dst) {
            min_dst   = dst;
            min_index = ind_list[ind];
        }
    }
    
    // Set new marker point properties
//    particles->x0[new_ind]            = particles->x[new_ind];
//    particles->z0[new_ind]            = particles->z[new_ind];
    particles->phase[new_ind]         = particles->phase[min_index];
    particles->sxxd[new_ind]          = particles->sxxd[min_index];
//    particles->szzd[new_ind]          = particles->szzd[min_index];
    particles->sxz[new_ind]           = particles->sxz[min_index];
//    particles->sxxd0[new_ind]         = particles->sxxd0[min_index];
//    particles->szzd0[new_ind]         = particles->szzd0[min_index];
//    particles->sxz0[new_ind]          = particles->sxz0[min_index];
    particles->Vx[new_ind]            = particles->Vx[min_index];
    particles->Vz[new_ind]            = particles->Vz[min_index];
//    particles->Vx0[new_ind]           = particles->Vx0[min_index];
//    particles->Vz0[new_ind]           = particles->Vz0[min_index];
    particles->strain[new_ind]        = particles->strain[min_index];
//    particles->rhoE[new_ind]          = particles->rhoE[min_index];
//    particles->rhoE0[new_ind]         = particles->rhoE0[min_index];
//    particles->eII[new_ind]           = particles->eII[min_index];
    particles->P[new_ind]             = particles->P[min_index];
//    particles->eta[new_ind]           = particles->eta[min_index];
    particles->generation[new_ind]    = particles->generation[min_index];
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FindClosestPhaseVertex_BEN( markers* particles, int ic, int jc, grid mesh, int *ind_list, int new_ind, int nb_neigh ) {
    
    // Find phase of a new particle according to its closest neighbour.
    // This function works at the scale of a cell.
    // The array 'ind_list' is list of indices of particle already owing to that cell.
    
    double dst, min_dst;
    int ind, icell, min_index=ind_list[0];
    
    min_dst = sqrt( pow( 4*mesh.dx, 2) + pow( 4*mesh.dz, 2) );
    icell = ic + jc * (mesh.Nx);
    
    for (ind=0; ind<nb_neigh; ind++) {
        
        dst = sqrt( pow (particles->x[new_ind] - particles->x[ind_list[ind]], 2) + pow(particles->z[new_ind] - particles->z[ind_list[ind]], 2) );
        
        if (dst<min_dst) {
            min_dst = dst;
            min_index = ind_list[ind];
        }
    }
    
    if (particles->phase[min_index] <-1 || particles->phase[min_index]>10 ) {
        printf("Could not find matching phase for newly created particle, this is a bug (attribute phase id : %d)\n Error at node ic = %d jc = %d \nExiting...", particles->phase[min_index], ic, jc);
        exit(50);
    }
    
    // Set new marker point properties
//    particles->x0[particles->Nb_part] = particles->x[particles->Nb_part];
//    particles->z0[particles->Nb_part] = particles->z[particles->Nb_part];
    particles->phase[new_ind]         = particles->phase[min_index];
    particles->sxxd[new_ind]          = particles->sxxd[min_index];
//    particles->szzd[new_ind]          = particles->szzd[min_index];
    particles->sxz[new_ind]           = particles->sxz[min_index];
//    particles->sxxd0[new_ind]         = particles->sxxd0[min_index];
//    particles->szzd0[new_ind]         = particles->szzd0[min_index];
//    particles->sxz0[new_ind]          = particles->sxz0[min_index];
    particles->Vx[new_ind]            = particles->Vx[min_index];
    particles->Vz[new_ind]            = particles->Vz[min_index];
//    particles->Vx0[new_ind]           = particles->Vx0[min_index];
//    particles->Vz0[new_ind]           = particles->Vz0[min_index];
    particles->strain[new_ind]        = particles->strain[min_index];
//    particles->rhoE[new_ind]          = particles->rhoE[min_index];
//    particles->rhoE0[new_ind]         = particles->rhoE0[min_index];
//    particles->eII[new_ind]           = particles->eII[min_index];
    particles->P[new_ind]             = particles->P[min_index];
//    particles->eta[new_ind]           = particles->eta[min_index];
    particles->generation[new_ind]    = particles->generation[min_index];
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartCell_BEN( markers *particles, grid mesh, int ic, int jc, int* ind_list, params model, int *ind_part_reuse, int *inc_reuse, int nb_part_reuse, surface topo  ) {
    
    // This functions add particles to cells that are tag as 'particle deficients'.
    // So far, it adds 4 new particles at different locations within the cell.
    // The new particles get the phase and stored stress of the closest particle.
    
    int Nb_add = 4;
    double new_x, new_z, h=model.zmax;
    int new_ind;
    double xmin = model.xmin + model.dx/2, distance;
    int cell;
    
    
    // Reallocate particle arrays if no more space in arrays (TO BE IMPLEMENTED)
    if (particles->Nb_part + Nb_add > particles->Nb_part_max && Nb_add > nb_part_reuse) {
        printf("You have reached the maximum number of particles currently available (%d), please increase it...\n", particles->Nb_part_max );
        printf("Exiting...\n");
        exit(1);
    }
    
    if ( Nb_add == 4 ) {
        
        //-------------------------------------------------------------------------------------------------------------------------//
        
        new_x = mesh.xg_coord[ic] + mesh.dx/3 - mesh.dx/10;
        new_z = mesh.zg_coord[jc] + mesh.dz/3 - mesh.dz/10;
        distance=fabs(new_x - xmin);
        cell = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x > model.xmin && new_z > model.zmin && new_z<h  ) {
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
            }
            
            // Add 4 particules
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhase_BEN( particles, ic, jc, mesh, ind_list, new_ind  );
        }
        
        //-------------------------------------------------------------------------------------------------------------------------//
        
        new_x = mesh.xg_coord[ic] + 2*mesh.dx/3 + mesh.dx/10;
        new_z = mesh.zg_coord[jc] + mesh.dz/3   - mesh.dz/10;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x > model.xmin && new_z > model.zmin && new_z<h  ) {
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhase_BEN( particles, ic, jc, mesh, ind_list, new_ind );
        }
        
        //-------------------------------------------------------------------------------------------------------------------------//
        
        new_x = mesh.xg_coord[ic] + mesh.dx/3   - mesh.dx/10;
        new_z = mesh.zg_coord[jc] + 2*mesh.dz/3 + mesh.dz/10;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x > model.xmin && new_z > model.zmin && new_z<h  ) {
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhase_BEN( particles, ic, jc, mesh, ind_list, new_ind );
        }
        
        //-------------------------------------------------------------------------------------------------------------------------//
        
        new_x = mesh.xg_coord[ic] + 2*mesh.dx/3  + mesh.dx/10;
        new_z = mesh.zg_coord[jc] + 2*mesh.dz/3  + mesh.dz/10;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x > model.xmin && new_z > model.zmin && new_z<h  ) {
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
            }
            else {
                new_ind = particles->Nb_part;
                //                printf( "generated new index : %d\n", new_ind);
                particles->Nb_part++;
            }
            
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhase_BEN( particles, ic, jc, mesh, ind_list, new_ind );
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartVert_BEN( markers *particles, grid mesh, int ic, int jc, int* ind_list, params model, int nb_neigh, int *ind_part_reuse, int *inc_reuse, int nb_part_reuse, surface topo ) {
    
    // This functions add particles to cells that are tag as 'particle deficients'.
    // So far, it adds 4 new particles at different locations within the cell.
    // The new particles get the phase and stored stress of the closest particle.
    
    int Nb_add = 4;
    double new_x, new_z, h=model.zmax;
    int new_ind;
    double xmin = model.xmin + model.dx/2, distance;
    int cell;
    
    // Reallocate particle arrays if no more space in arrays (TO BE IMPLEMENTED)
    if (particles->Nb_part + Nb_add > particles->Nb_part_max && Nb_add > nb_part_reuse ) {
        printf("You have reached the maximum number of particles currently available (%d), please increase it...\n", particles->Nb_part_max );
        printf("Exiting...\n");
        exit(1);
    }
    
    // Add 4 particules
    if ( Nb_add == 4 ) {
        
        //-------------------------------------------------------------------------------------------------------------------------//
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/2 + mesh.dx/3;
        new_z = mesh.zg_coord[jc] - mesh.dz/2 + mesh.dz/3;
        distance=fabs(new_x - xmin);
        
        cell = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x > model.xmin && new_z > model.zmin && new_z<h   ) {
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
                //printf("reusing p #%d\n", new_ind);
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
                //printf("new p #%d\n", new_ind);
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex_BEN( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh );
        }
        
        //-------------------------------------------------------------------------------------------------------------------------//
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/2 + 2*mesh.dx/3;
        new_z = mesh.zg_coord[jc] - mesh.dz/2 + mesh.dz/3;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_x < model.xmax && new_z > model.zmin && new_z<h   ) {
            
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
                //printf("reusing p #%d\n", new_ind);
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
                //printf("new p #%d\n", new_ind);
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex_BEN( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh );
        }
        
        //-------------------------------------------------------------------------------------------------------------------------//
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/2 + mesh.dx/3;
        new_z = mesh.zg_coord[jc] - mesh.dz/2 + 2*mesh.dz/3;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_z < model.zmax && new_x > model.xmin && new_z<h  ) {
            
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[(*inc_reuse)];
                (*inc_reuse)++;
                //printf("reusing p #%d inc_reuse=%d\n", new_ind, *inc_reuse);
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
                //printf("new p #%d\n", new_ind);
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex_BEN( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh );
        }
        
        //------------------------------------------------------------------------------------------------------------------------/
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/2 + 2*mesh.dx/3;
        new_z = mesh.zg_coord[jc] - mesh.dz/2 + 2*mesh.dz/3;
        distance=fabs(new_x - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        if ( model.free_surf==1 ) h = topo.a[cell]*new_x + topo.b[cell];
        
        if ( new_z < model.zmax && new_x < model.xmax && new_z<h ) {
            
            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                new_ind = ind_part_reuse[*inc_reuse];
                (*inc_reuse)++;
            }
            else {
                new_ind = particles->Nb_part;
                particles->Nb_part++;
            }
            
            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex_BEN( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh );
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CountPartCell_BEN ( markers* particles, grid *mesh, params model, surface topo, int reseed, scale scaling  ) {
    
    // This function counts the number of particle that are currently in each cell of the domain.
    // The function detects cells that are lacking of particle and call the particle re-seeding routine.
    
    int k, l, ic, jc, kc, nthreads, thread_num, **NPPC, **NPPV;
    int ncx=mesh->Nx-1, ncz=mesh->Nz-1, Nx=mesh->Nx, Nz=mesh->Nz, Nb_part=particles->Nb_part;
    double distance;
    int nb_part_reuse=0, *ind_part_reuse, inc_reuse=0;
    double dx=mesh->dx, dz=mesh->dz;
    
    int **ind_per_cell, *ind_part;
    int kk = 0;
    int jj, n_neigh=0, n_neigh0, ic2, *neigh_part, ii;
    int p, count, sum;
    
    // Initialise
    Initialise2DArrayInt( mesh->nb_part_cell, ncx, ncz, 0 );
    Initialise2DArrayInt( mesh->nb_part_vert, Nx, Nz, 0 );
    
    //------------------------- NODES -----------------------//
    
    
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    
    // Allocate
    NPPV = DoodzCalloc( nthreads, sizeof(int*));
    for (k=0; k<nthreads; k++) {
        NPPV[k] = DoodzCalloc( Nx*Nz, sizeof(int));
    }
    
        double ***PHASE_s;
    
        // Allocate
        PHASE_s = DoodzCalloc( nthreads, sizeof(double**));
        for (k=0; k<nthreads; k++) {
            PHASE_s[k] = DoodzCalloc( model.Nb_phases, sizeof(double*));
            for (l=0; l<model.Nb_phases; l++) {
                PHASE_s[k][l] = DoodzCalloc( (Nx)*(Nz), sizeof(double));
            }
        }
    
    
    // Count number of particles per nodes
    //#pragma omp parallel for shared ( particles, model, NPPV, Nb_part, PHASE_s  )    reduction(+:nb_part_reuse)        \
    private ( k, distance, ic, jc, thread_num )                    \
    firstprivate( ncx, mesh ) schedule( static )
    
#pragma omp parallel for shared ( particles, model, NPPV, Nb_part  )    reduction(+:nb_part_reuse)        \
private ( k, distance, ic, jc, thread_num )                    \
firstprivate( ncx, mesh ) schedule( static )
    
    
    for (k=0; k<Nb_part; k++) {
        
        thread_num = omp_get_thread_num();
        
        if ( particles->phase[k] != -1 ) {
            
            // -------- Check NODES ----------
            
            // Get column/line indices
            distance=fabs(particles->x[k]-mesh->xg_coord[0]);
            ic=ceil((distance/dx)+0.5) -1;
            distance=fabs(particles->z[k]-mesh->zg_coord[0]);
            jc=ceil((distance/dz)+0.5) -1;
            
            // If particles are around the nodes: count them
            if ( particles -> phase[k]>=0 ) {
                NPPV[thread_num][ic + jc * Nx]++;
                 PHASE_s[thread_num][particles->phase[k]][ic + jc * Nx] ++;
            }
        }
        else {
            // -------- Check What's OUT  ----------
            // If particles are taged as 'out' (-1): count them for re-use
            (nb_part_reuse)++;
        }
    }
    
    printf ("Before NODE re-seeding, nb_part_reuse: %d\n ", nb_part_reuse);
    
        // Initialise Phase arrays
    #pragma omp parallel for shared ( mesh, Nx, Nz, model ) private( p, l ) schedule( static )
        for ( l=0; l<Nx*Nz; l++ ) {
            for ( p=0; p<model.Nb_phases; p++ ) {
                mesh->phase_perc_s[p][l] = 0;
            }
        }
    
    
    // Final reduction for NODES
    //#pragma omp parallel for shared ( mesh, NPPV, Nx, Nz, nthreads, PHASE_s ) private( l, k ) schedule( static )
#pragma omp parallel for shared ( mesh, NPPV, Nx, Nz, nthreads ) private( l, k ) schedule( static )
    for ( l=0; l<Nx*Nz; l++ ) {
        for ( k=0; k<nthreads; k++ ) {
                        for ( p=0; p<model.Nb_phases; p++ ) {
                            mesh->phase_perc_s[p][l] += PHASE_s[k][p][l];
                        }
            mesh->nb_part_vert[l]   += NPPV[k][l];
        }
    }
    
    // Normalise to phase percentages
#pragma omp parallel for shared ( mesh, Nx, Nz, model ) private( p, l ) schedule( static )
    for ( l=0; l<Nx*Nz; l++ ) {
        //        if ( mesh->BCg.type[l] != 30 ) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            if (mesh->nb_part_vert[l]>0) mesh->phase_perc_s[p][l] /= mesh->nb_part_vert[l];
        }
        
        //        }
    }
    
//    for ( p=0; p<model.Nb_phases; p++ ) {
//        MinMaxArrayTag( mesh->phase_perc_s[p], 1.0, Nx*Nz,  "perc s", mesh->BCg.type);
//    }
    //-------------------------------------------------------------------------------------------------------------------------//
    
    // NODES CYCLES
    // Loop over cells and check for the minimum amount of particles
    
    ind_per_cell = DoodzMalloc(sizeof(int*)*(mesh->Nx) * (mesh->Nz));
    ind_part     = DoodzCalloc((mesh->Nx) * (mesh->Nz), sizeof(int));
    
    // Loop on cells and allocate the required number of particles
    for (k=0; k<(mesh->Nx) * (mesh->Nz); k++) {
        ind_per_cell[k] = DoodzMalloc(sizeof(int)*(mesh->nb_part_vert[k]));
    }
    
    // Initialise the array containing particles indices to reuse
    ind_part_reuse = DoodzMalloc(sizeof(int)*(nb_part_reuse));
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    
    // Build list of particles to reuse
    
    //#pragma omp parallel for shared ( particles, ind_part_reuse ) reduction ( +:kk )
    for (k=0; k<Nb_part; k++) {
        if ( particles->phase[k] == -1 ) {
            ind_part_reuse[kk] = k;
            kk++;
        }
    }
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    
    // LOOP ON PARTICLES, GIVE THEIR ID'S TO THE CELL THAT OWNS THEM
    //#pragma omp parallel for shared ( particles, mesh, ind_part, ind_per_cell )    \
    private ( k, kc, distance, ic, jc )                    \
    firstprivate ( Nb_part, dx, dz ) \
    schedule ( static )
    
    for (k=0; k<Nb_part; k++) {
        
        if ( particles->phase[k] !=-1 ) {
            
            // Get the column:
            distance=fabs(particles->x[k]-mesh->xg_coord[0]);
            ic=ceil((distance/dx)+0.5) -1;
            
            // Get the line:
            distance=fabs(particles->z[k]-mesh->zg_coord[0]);
            jc=ceil((distance/dz)+0.5) -1;
            
            // Global index
            kc = ic + jc * (mesh->Nx);
            ind_per_cell[kc][ind_part[kc]] = k;
            
            //#pragma omp atomic
            (ind_part[kc])++;
        }
        
    }
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    if (reseed == 1) {
    
    // LOOP ON NODES - RESEED PARTICLES IF NEEDED
    for (k=0; k<mesh->Nx; k++) {
        for (l=0; l<mesh->Nz; l++) {
            ic = k + l * (mesh->Nx);
            
            // Detect unsufficient number of points
            if (mesh->nb_part_vert[ic] < particles -> min_part_cell && mesh->BCg.type[ic] != 30 ) {
                
                // IF there are no particles in the current cell, look in neighbouring cells
                if ( mesh->nb_part_vert[ic] < 1 ) {
                    
                    n_neigh = 0;
                    
                    // Determine number of particles in neighbouring cells, take care of the boundaries
                    if (k>0) {
                        ic2 = k + l * (mesh->Nx) - 1;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if (k<mesh->Nx-1) {
                        ic2 = k + l * (mesh->Nx) + 1;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if (l>0) {
                        ic2 = k + l * (mesh->Nx) - mesh->Nx;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if (l<mesh->Nz-1) {
                        ic2 = k + l * (mesh->Nx) + mesh->Nx;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    
                    // Diagonal neighbours
                    if ( k>0 && l<mesh->Nz-1 ) { // NW
                        ic2 = k + l * (mesh->Nx) - 1 + mesh->Nx;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if ( k>0 && l>0 ) {  // SW
                        ic2 = k + l * (mesh->Nx) - 1 - mesh->Nx;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if ( k<mesh->Nx-1 && l<mesh->Nz-1 ) { // NE
                        ic2 = k + l * (mesh->Nx) + 1 + mesh->Nx;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    if ( k<mesh->Nx-1 && l>0 ) { // SE
                        ic2 = k + l * (mesh->Nx) - mesh->Nx + 1;
                        n_neigh += mesh->nb_part_vert[ic2];
                    }
                    
                    
                    // Allocate array of neighbouring cell particles
                    neigh_part = DoodzMalloc( sizeof(int)*n_neigh);
                    
                    // Reset n-neigh, re-evaluate number of neighbourgs and save their indices
                    n_neigh = 0;
                    
                    // Determine the list of particle indices that are lying in the neighbouring cells, take care of the boundaries
                    if (k>0) {
                        ic2 = k + l * (mesh->Nx) - 1;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    if (k<mesh->Nx-1) {
                        ic2 = k + l * (mesh->Nx) + 1;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    if (l>0) {
                        ic2 = k + l * (mesh->Nx) - mesh->Nx;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    if (l<mesh->Nz-1) {
                        ic2 = k + l * (mesh->Nx) + mesh->Nx;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    
                    // Diagonal neighbours
                    if ( k>0 && l<mesh->Nz-1 ) { // NW
                        ic2 = k + l * (mesh->Nx) - 1 + mesh->Nx;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    
                    if ( k>0 && l>0 ) {  // SW
                        ic2 = k + l * (mesh->Nx) - 1 - mesh->Nx;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    
                    if ( k<mesh->Nx-1 && l<mesh->Nz-1 ) { // NE
                        ic2 = k + l * (mesh->Nx) + 1 + mesh->Nx;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    
                    if ( k<mesh->Nx-1 && l>0 ) { // SE
                        ic2 = k + l * (mesh->Nx) - mesh->Nx + 1;
                        n_neigh0 = n_neigh;
                        n_neigh += mesh->nb_part_vert[ic2];
                        ii = 0;
                        for (jj=n_neigh0; jj<n_neigh; jj++) {
                            neigh_part[jj] = ind_per_cell[ic2][ii];
                            ii++;
                        }
                    }
                    
                    if (n_neigh==0) {
                        printf("no neighbouring particles for node k = %d l = %d tag=%d z=%lf... The code is gonna breakdown...\n", k,l, mesh->BCg.type[ic], mesh->zg_coord[l]*1000);
                        exit(1);
                    }
                    else {
                        AddPartVert_BEN( particles, *mesh, k, l, neigh_part, model, n_neigh, ind_part_reuse, &inc_reuse, nb_part_reuse, topo );
                    }
                    DoodzFree(neigh_part);
                }
                else {
                    //AddPartVert( particles, *mesh, k, l, ind_per_cell[ic], model,  mesh->nb_part_vert[ic], ind_part_reuse, &inc_reuse, nb_part_reuse );
                }
            }
        }
    }
        
    }
    
    printf("%d re-used particles, %d particles left out\n", inc_reuse, nb_part_reuse-inc_reuse);
    
    /* free table-table */
    for (k=0; k<(mesh->Nx) * (mesh->Nz); k++) {
        DoodzFree(ind_per_cell[k]);
    }
    DoodzFree(ind_per_cell);
    DoodzFree(ind_part);
    DoodzFree(ind_part_reuse);
    
        // Free PHASE
        for (k=0; k<nthreads; k++) {
            for (l=0; l<model.Nb_phases; l++) {
                DoodzFree( PHASE_s[k][l] );
            }
            DoodzFree( PHASE_s[k] );
        }
        DoodzFree( PHASE_s );
    
    
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    
    
    // CELL CYCLES
    
    nb_part_reuse = 0;
    int *ind_part_reuse2;
    double ***PHASE_n;
    
    // Allocate
    NPPC = DoodzCalloc( nthreads, sizeof(int*));
    for (k=0; k<nthreads; k++) {
        NPPC[k] = DoodzCalloc( Nx*Nz, sizeof(int));
    }
    
    // Allocate
    PHASE_n = DoodzCalloc( nthreads, sizeof(double**));
    for (k=0; k<nthreads; k++) {
        PHASE_n[k] = DoodzCalloc( model.Nb_phases, sizeof(double*));
        for (l=0; l<model.Nb_phases; l++) {
            PHASE_n[k][l] = DoodzCalloc( (Nx)*(Nz), sizeof(double));
        }
    }
    
    
    // Count number of particles per cell
#pragma omp parallel for shared ( particles, model, NPPC, Nb_part, PHASE_n )    reduction(+:nb_part_reuse )        \
private ( k, distance, ic, jc, thread_num )                    \
firstprivate( ncx, mesh ) schedule( static )
    
    for (k=0; k<Nb_part; k++) {
        
        thread_num = omp_get_thread_num();
        
        if ( particles->phase[k] != -1 ) {
            
            // -------- Check CELLS ----------
            
            // Get column/line indices
            distance=fabs(particles->x[k] - mesh->xc_coord[0]);
            ic = ceil((distance/mesh->dx)+0.5) - 1;
            distance=fabs(particles->z[k] - mesh->zc_coord[0]);
            jc = ceil((distance/mesh->dz)+0.5) - 1;
            
            // If particles within cells then count them
            if ( particles -> phase[k]>=0 ) {
                NPPC[thread_num][ic + jc * ncx] ++;
                PHASE_n[thread_num][particles->phase[k]][ic + jc * ncx] ++;
            }
        }
        else {
            // -------- Check What's OUT  ----------
            // If particles are taged as 'out' (-1): count them for re-use
            (nb_part_reuse)++;
        }
    }
    
    printf ("Before CELL re-seeding, nb_part_reuse: %d\n", nb_part_reuse);
    
    // Initialise Phase arrays
#pragma omp parallel for shared ( mesh, ncx, ncz, model ) private( p, l ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            mesh->phase_perc_n[p][l] = 0;
        }
    }
    
    // Final reduction for CELLS
#pragma omp parallel for shared ( mesh, NPPC, ncx, ncz, nthreads, PHASE_n, model ) private( p, l, k ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        
        for ( k=0; k<nthreads; k++ ) {
            for ( p=0; p<model.Nb_phases; p++ ) {
                mesh->phase_perc_n[p][l] += PHASE_n[k][p][l];
            }
            mesh->nb_part_cell[l]   += NPPC[k][l];
        }
    }
    
    //    // Find number of different phases per cell
    //#pragma omp parallel for shared ( mesh, ncx, ncz, model ) private( count, p, l ) schedule( static )
    //    for ( l=0; l<ncx*ncz; l++ ) {
    //        count = 0;
    //            for ( p=0; p<model.Nb_phases; p++ ) {
    //                //if (mesh->phase_perc_n[p][l]>0) count++;
    //            }
    //        //mesh->num_phases_n[l] = count;
    //    }
    
    // Normalise to phase percentages
#pragma omp parallel for shared ( mesh, ncx, ncz, model ) private( p, l ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        //        if ( mesh->BCp.type[0][l] != 30 && mesh->BCp.type[0][l] != 31) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            if (mesh->nb_part_cell[l]>0) mesh->phase_perc_n[p][l] /= mesh->nb_part_cell[l];
        }
        //        }
    }
    //    int cc;
    //
    //    for( l=0; l<ncz; l++) {
    //        for( k=0; k<ncx; k++) {
    //            cc = k + l*ncx;
    //            printf("%d ", mesh->nb_part_cell[cc]);
    //        }
    //        printf("\n");
    //    }
    //
    //     printf("\n");
    //
    //    for( l=0; l<ncz; l++) {
    //        for( k=0; k<ncx; k++) {
    //            cc = k + l*(ncx);
    //            printf("%d ", mesh->BCp.type[cc]);
    //        }
    //        printf("\n");
    //    }
    //
    //
    //    printf("\n");
    //
    //    int cc;
    //    for( l=0; l<ncz; l++) {
    //        for( k=0; k<ncx; k++) {
    //            cc = k + l*ncx;
    //            printf("%lf ", mesh->phase_perc_n[cc]);
    //        }
    //        printf("\n");
    //    }
    //
    //    printf("\n");
    //
    //    for( l=0; l<ncz; l++) {
    //        for( k=0; k<ncx; k++) {
    //            cc = k + l*ncx;
    //            printf("%lf ", mesh->phase_perc_n[1][cc]);
    //        }
    //        printf("\n");
    //    }
    
    
    //-------------------------------------------------------------------------------------------//
    
    // Initialise the array containing particles indices to reuse
    ind_part_reuse2 = DoodzMalloc(sizeof(int)*(nb_part_reuse));
    
    // Loop over cells and check for the minimum amount of particles
    ind_per_cell =  DoodzMalloc(sizeof(int*)*ncx*ncz);
    ind_part     =  DoodzCalloc(sizeof(int), ncx*ncz);
    
    // Loop on cells and allocate the required number of particles
    for (k=0; k<ncx*ncz; k++) {
        ind_per_cell[k] =  DoodzMalloc(sizeof(int)*(mesh->nb_part_cell[k]));
    }
    
    // Build list of particles to reuse
    kk = 0;
    //#pragma omp parallel for shared ( particles, ind_part_reuse ) reduction ( +:kk )
    for (k=0; k<Nb_part; k++) {
        if ( particles->phase[k] == -1 ) {
            ind_part_reuse2[kk] = k;
            kk++;
        }
    }
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    
    // Loop on particles, give their ID's to the cell that owns them
    //#pragma omp parallel for shared ( particles, mesh, ind_part, ind_per_cell )    \
    private ( k, kc, distance, ic, jc )                    \
    firstprivate( Nb_part )
    for (k=0; k<Nb_part; k++) {
        
        if ( particles->phase[k] != -1 ) {
            distance=fabs(particles->x[k] - mesh->xc_coord[0]);
            ic = ceil((distance/mesh->dx)+0.5) - 1;
            distance=fabs(particles->z[k] - mesh->zc_coord[0]);
            jc = ceil((distance/mesh->dz)+0.5) - 1;
            
            kc = ic + jc * (mesh->Nx-1);
            
            ind_per_cell[kc][ind_part[kc]] = k;
            //#pragma omp atomic
            ind_part[kc] += 1;
            
        }
    }
    
    //-------------------------------------------------------------------------------------------------------------------------//
    
    
    int *ind_list, neighs, oo, nb;
    
     if (reseed == 1) {
    
    // Loop on cells and add particles
    for (k=0; k<ncx; k++) {
        for (l=0; l<ncz; l++) {
            ic = k + l * ncx;
            
            if (mesh->BCt.type[ic] != 30) {
                
                if ( mesh->nb_part_cell[ic] < particles -> min_part_cell ) {
                    
                    if ( mesh->nb_part_cell[ic]>0 ) {
                        
                        // search closest particles in the current cell
                        AddPartCell_BEN( particles, *mesh, k, l, ind_per_cell[ic], model, ind_part_reuse2, &inc_reuse, nb_part_reuse, topo   );
                    }
                    else {
                        
                        neighs = 0;
                        oo     = 0;
                        
                        //                    printf("No particle in the CELL i=%d j=%d, checking the neigbours...\n",  k, l);
                        // list the number of particles in the neighbouring cells
                        if (k>0)     neighs+= mesh->nb_part_cell[ic-1];
                        if (k<ncx-1) neighs+= mesh->nb_part_cell[ic+1];
                        if (l>0)     neighs+= mesh->nb_part_cell[ic-ncx];
                        if (l<ncz-1) neighs+= mesh->nb_part_cell[ic+ncx];
                        
                        if (neighs == 0) {
                            printf("All the neighbouring CELLS of ix = %d iz = %d are empty, simulation will stop\n", k, l);
                            //                        exit(200);
                        }
                        
                        ind_list = DoodzMalloc( neighs * sizeof(int));
                        if (k>0) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic-1]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic-1][nb];
                                oo++;
                            }
                        }
                        if (k<ncx-1) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic+1]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic+1][nb];
                                oo++;
                            }
                        }
                        if (l>0) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic-ncx]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic-ncx][nb];
                                oo++;
                            }
                        }
                        if (l<ncz-1) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic+ncx]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic+ncx][nb];
                                oo++;
                            }
                        }
                        
                        //                    printf("CELL i=%d j=%d, Xc=%lf, Zc=%lf\n", k, l, mesh->xc_coord[k]*1000, mesh->zc_coord[l]*1000);
                        //                    for ( nb=0; nb<neighs; nb++ ) {
                        //                        printf( "%d ", ind_list[nb] );
                        //                    }
                        //                    printf( "\n" );
                        
                        AddPartCell_BEN( particles, *mesh, k, l, ind_list, model, ind_part_reuse2, &inc_reuse, nb_part_reuse, topo );
                        DoodzFree(ind_list);
                    }
                }
            }
        }
    }
     }
    printf("%d re-used particles, %d particles left out\n", inc_reuse, nb_part_reuse-inc_reuse);
    
    
    /* free table-table */
    for (k=0; k<ncx*ncz; k++) {
        DoodzFree(ind_per_cell[k]);
    }
    DoodzFree(ind_per_cell);
    DoodzFree(ind_part);
    DoodzFree(ind_part_reuse2);
    
    for (l=0; l<nthreads; l++) {
        DoodzFree(NPPC[l]);
    }
    DoodzFree(NPPC);
    
    for ( k=0; k<nthreads; k++ ) {
        DoodzFree(NPPV[k]);
    }
    DoodzFree(NPPV);
    
    // Free PHASE
    for (k=0; k<nthreads; k++) {
        for (l=0; l<model.Nb_phases; l++) {
            DoodzFree( PHASE_n[k][l] );
        }
        DoodzFree( PHASE_n[k] );
    }
    DoodzFree( PHASE_n );
    
    
    
    
}
/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PutPartInBox_BEN( markers *particles, grid *mesh, params model, surface topo, scale scaling )

// Set the original particle layout throughout the mesh.
// The particles are set with regular spacing.

{
    int i, j, ki, kj, np;
    double dx_particles, dz_particles;
    double H, W;
    
    //int add_noise, double noise, double* random_x, double* random_z
    
    // Model size
    W = mesh->xg_coord[ mesh->Nx-1 ] - mesh->xg_coord[0];
    H = mesh->zg_coord[ mesh->Nz-1 ] - mesh->zg_coord[0];
    
    // Compute the spacing between particles:
    dx_particles = mesh->dx/(double)(particles->Nx_part+1);
    dz_particles = mesh->dz/(double)(particles->Nz_part+1);
    
    printf("dx = %2.6e dz = %2.6e Nx_part = %d\n", mesh->dx, mesh->dz, (particles->Nx_part));

    
    printf("Initial particle spacing : dxm = %lf dzm = %lf m\n", dx_particles*scaling.L, dz_particles*scaling.L);
    
    // Loop over loop over loop over loop (on s'en branle, c'est du C)
    np = 0;
    double h, a, b, coord_x, coord_z;
    
    for (j=0; j<mesh->Nx-1; j++) {
        for (i=0; i<mesh->Nz-1; i++) {
            
            if (model.free_surf == 1) {
                a = topo.a[j];
                b = topo.b[j];
            }
            for (kj=0; kj<particles->Nx_part; kj++) {
                for (ki=0; ki<particles->Nz_part; ki++) {
                    
                    if (model.free_surf == 1) {
                        
                        coord_x = mesh->xg_coord[j] + dx_particles*kj + dx_particles;
                        coord_z = mesh->zg_coord[i] + dz_particles*ki + dz_particles;
                        
                        h = b + a*coord_x;
                        
                        if ( coord_z < h ) {
                            particles->x[np]  = coord_x;
                            particles->z[np]  = coord_z;
                            np++;
                        }
                    }
                    
                    if (model.free_surf == 0) {
                        coord_x = mesh->xg_coord[j] + dx_particles*kj + dx_particles;
                        coord_z = mesh->zg_coord[i] + dz_particles*ki + dz_particles;
                        particles->x[np]  = coord_x;
                        particles->z[np]  = coord_z;
                        np++;
                    }
                }
            }
        }
    }
    particles->Nb_part = np;
    
    /* //------------------------------------ */
    /* // Adding noise in x and z direction: */
    /* //------------------------------------ */
    /* if (add_noise == 1) { */
    
    /*   for (np=0;np<particles.x.n;np++) { */
    /*     random_x[np]=(double)rand()/(double)RAND_MAX; */
    /*     random_z[np]=(double)rand()/(double)RAND_MAX; */
    /*     random_x[np]=((2*random_x[np])-1)*dx_particles*noise; */
    /*     random_z[np]=((2*random_z[np])-1)*dz_particles*noise; */
    
    /*     //      printf("random_x : %f\n",random_x[np]); */
    /*     //      printf("random_z : %f\n",random_z[np]); */
    
    /*     particles.x.val[np]=particles.x.val[np]+random_x[np]; */
    /*     particles.z.val[np]=particles.z.val[np]+random_z[np]; */
    /*       } */
    /*   //disp('   Check if particles are inside the model box'); */
    /*   eps=0.0001; */
    /*   for (np=0;np<particles.x.n;np++) { */
    /*     if(particles.x.val[np]<=eps) { */
    /*     //printf("Oooops... particle outside the box - left side!\n"); */
    /*         particles.x.val[np]=0+eps; */
    /*     } */
    
    /*     if(particles.x.val[np]>=W-eps) { */
    /*     //printf("Oooops... particle outside the box - right side!\n"); */
    /*         particles.x.val[np]=W-eps; */
    /*     } */
    
    /*     if(particles.z.val[np]>=-eps) { */
    /*     //printf("Oooops... particle outside the box - top side!\n"); */
    /*         particles.z.val[np]=-eps; */
    /*     } */
    
    /*     if(particles.z.val[np]<=-H+eps) { */
    /*     //printf("Oooops... particle outside the box - bottom side!\n"); */
    /*         particles.z.val[np]=-H+eps; */
    /*     } */
    
    /*     } */
    /* } */
    
    printf("Initial number of particles = %d\n", particles->Nb_part);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs_BEN( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials ) {
    
    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz;
    
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
    
    //    printf ("BEFORE FLaG\n");
    //    for (l=0; l<mesh->Nz-1; l++) {
    //        for (k=0; k<mesh->Nx-1; k++) {
    //            c = k + l*(mesh->Nx-1);
    //
    //            printf ("%04.0lf ", mesh->Cv[c]*scaling.Cv);
    //        }
    //        printf("\n");
    //    }
    
    
    
    //    printf ("AFTER FLaG\n");
    //    for (l=0; l<mesh->Nz-1; l++) {
    //        for (k=0; k<mesh->Nx-1; k++) {
    //            c = k + l*(mesh->Nx-1);
    //
    //            printf ("%04.0lf ", mesh->Cv[c]*scaling.Cv);
    //        }
    //        printf("\n");
    //    }
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
    
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
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<mesh->Nz+1; l++) {
            for (k=0; k<mesh->Nx; k++) {
                
                c = k + l*(mesh->Nx);
                
                if ( mesh->BCu.type[c] != 30 ) {
                    
                    //                // Internal points:  -1
                    //                mesh->BCu.type[c] = -1;
                    //                mesh->BCu.val[c]  =  0;
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -0.0;
                    }
                    
                    // Free slip
                    if (l==0  ) {
                        mesh->BCu.type[c] =  3;
                        mesh->BCu.val[c]  =  0;
                    }
                    
                    // N
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] =  3;
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
    int numVz = 0, numVz0 = 0;
    
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<mesh->Nz; l++) {
            for (k=0; k<mesh->Nx+1; k++) {
                
                c  = k + l*(mesh->Nx+1);
                
                // Boundary Conditions
                
                if ( mesh->BCv.type[c] != 30 ) {
                    
                    numVz++;
                    
                    //                // Internal points:  -1
                    //                mesh->BCv.type[c] = -1;
                    //                mesh->BCv.val[c]  =  0;
                    
                    // Non-matching boundary points
                    if ( k==0 ) {
                        mesh->BCv.type[c] =   3;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Non-matching boundary points
                    if ( k==mesh->Nx ) {
                        mesh->BCv.type[c] =   3;
                        mesh->BCv.val[c]  =   0;
                    }
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;
                    }
                    
                }
                if ( mesh->BCv.type[c] == 0) numVz0 ++;
            }
        }
        
    printf("%d %d\n",numVz, numVz0);
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
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
                
                // Boundary Conditions
                
                if ( mesh->BCp.type[c] != 30 ) {
                    
                    //                // Internal points:  -1
                    //                mesh->BCp.type[c] = -1;
                    //                mesh->BCp.val[c]  =  0;
                    
                    //                // Matching boundary points (4 corners)
                    //                if ( (l==0 && k==0) || (l==0 && k==NCX-1) || (l==NCZ-1 && k==NCX-1) || (l==NCZ-1 && k==0)  ) {
                    //                   mesh->BCp.type[c] = 0;
                    //                   mesh->BCp.val[c]  = 0;
                    //                }
                    
                    //                // Top corner pressures
                    //                if ((k==0 && l==NCZ-1) || (k==NCX-1 && l==NCZ-1) ) {
                    //                        mesh->BCp.type[c] = 0;
                    //                        mesh->BCp.val[c]  = 0;
                    //                }
                    
                    // Top left corner pressure
                    if ((k==1 && l==NCZ-1) ) {
                        mesh->BCp.type[c] = 0;
                        mesh->BCp.val[c]  = 0;//mesh->p_in[c+1];
                        //printf("BCP %lf\n",mesh->BCp.val[c]);
                    }
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
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {
                
                c = k + l*(NCX);
                
                //                // Internal points:  -1
                //                mesh->BCt.type[c] = -1;
                //                mesh->BCt.val[c]  =  0;
                
                if (mesh->BCt.type[c] != 30) {
                    
                    //
                    //                    if ( mesh->BCt.type[c+NCX] == 30 ) {
                    //                        mesh->BCt.type[c+NCX] = 1;
                    //                        mesh->BCt.val[c+NCX]  = 273.15/scaling.T;
                    //                    }
                    
//                    // LEFT
//                    if ( k==0 ) {
//                        mesh->BCt.type[c] = 1;
//                        mesh->BCt.val[c]  = mesh->rhoE[c] / (mesh->rho_n[c] * mesh->Cv[c]);
//                    }
//
//                    // RIGHT
//                    if ( k==NCX-1 ) {
//                        mesh->BCt.type[c] = 1;
//                        mesh->BCt.val[c]  = mesh->rhoE[c] / (mesh->rho_n[c] * mesh->Cv[c]);
//                    }
//
//                    // BOT
//                    if ( l==0 ) {
//                        mesh->BCt.type[c] = 1;
//                        mesh->BCt.val[c]  = 773.15/scaling.T;
//                    }
//
//                    // TOP
//                    if ( l==NCZ-1 ) {
//                        mesh->BCt.type[c] = 1;
//                        mesh->BCt.val[c]  = 273.15/scaling.T;
//                    }
//                    else {
//                        if (mesh->BCt.type[c] == -1 && mesh->BCt.type[c+NCX] == 30) {
//                            mesh->BCt.type[c] = 1;
//                            mesh->BCt.val[c]  = mesh->rhoE[c] / (mesh->rho_n[c] * mesh->Cv[c]);
//                        }
//                    }
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

void SetParticles_BEN( markers *particles, scale scaling, params model, mat_prop* materials  ) {
    
    int np;
    double Lx = (double) (model.xmax - model.xmin);
    double Lz = (double) (model.zmax - model.zmin);
    double Tbot = 773.15/scaling.T;
    double gsbg   = 2e-3/scaling.L;                           // reference grain size
    
    double  H = 1.2e-2/scaling.L;
    double x0 = 20.0e-2/scaling.L, y0 = -H;
    double x1 = x0 - 4.974e-2/scaling.L, y1 = -3.3552e-2/scaling.L-H;
    double L = sqrt( pow(x0-x1,2) + pow(y0-y1,2) );
    double dip = acos( (x0-x1)/L );
    double delta_z = H/ cos(dip);
    //    double x2 = x1, y2 = y1 + delta_z;
    //    double x4 = x0, y4 = 0.0;
    //    double a1 = (y0-y1)/(x0-x1);
    //    double a2 = (y4-y2)/(x4-x1);
    //    double b1 = y1  - (x1*a1);
    //    double b2 = y2  - (x2*a2);
    //    double width_channel = 30e3/scaling.L;
    //    double a3 = a2;
    //    double b3 = y2  - ((x2-width_channel)*a3);
    
    
    double ax0 = 20.0e-2/scaling.L, ay0 = 0e-2/scaling.L;
    double ax1 = ax0 - 4.974e-2/scaling.L, ay1 = ay0 - 3.3552e-2/scaling.L;
    double ax2 = ax1 + 0.0049/scaling.L, ay2 = ay1 - 0.0110/scaling.L;
    double ax3 = ax0 + 0.0049/scaling.L, ay3 = ay0 - 0.0110/scaling.L;
    
    // Upper slab limit
    double aa0 = (ay0-ay1)/(ax0-ax1);
    double bb0 = ay1  - (ax1*aa0);
    
    // Lower slab limit
    double aa2 = (ay3-ay2)/(ax3-ax2);
    double bb2 = ay2  - (ax2*aa2);
    
    // Left slab limit
    double aa1 = (ay2-ay1)/(ax2-ax1);
    double bb1 = ay2  - (ax2*aa1);
    
    // Right slab limit
    double aa3 = (ay3-ay0)/(ax3-ax0);
    double bb3 = ay0  - (ax0*aa3);
    
    
    // right end of plate
    double xright = 1000e3/scaling.L;
    
    
    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;
        //        particles->d[np]     = gsbg;                          // same grain size everywhere
        //        particles->phi[np]   = 0.0;                           // zero porosity everywhere
        //        particles->X[np]     = 0.0;                           // X set to 0
        //        particles->T[np]     = Tbot;                          // same temperature everywhere
        
        particles->phase[np] = 1;
        
        // Draw dipping slap
        if ( particles->z[np] < aa0*particles->x[np]+bb0 && particles->z[np] > aa2*particles->x[np]+bb2 && particles->z[np] > aa1*particles->x[np]+bb1 && particles->z[np] < aa3*particles->x[np]+bb3 ) {
            particles->phase[np] = 0;
        }
        
        // Lower plate
        if (particles->x[np] > 20.0e-1 && particles->z[np] > -H) {
            particles->phase[np] = 0;
        }
        
        // Sticky air
        if ( particles->z[np] > 0.0 ) {
            particles->phase[np] = 2;
        }
        
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokes_BEN( SparseMat *mat, DirectSolver *pardi, double *rhs, double *sol ) {
    
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
    int status, i;
    void *Symbolic = NULL;
    void *Numeric = NULL;
    int msglvl = 0;
    
//    clock_t t_omp;
    
//    t_omp = (double)omp_get_wtime();
    
    if (msglvl==1) printf ("\nUMFPACK V%d.%d (%s) demo: _di_ version\n",  UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;
    
    /* get the default control parameters */
    umfpack_di_defaults (Control) ;
    
    /* change the default print level for this demo */
    /* (otherwise, nothing will print) */
    if (msglvl==1) Control [UMFPACK_PRL] = 6 ;
    
    /* print the license agreement */
    if (msglvl==1) umfpack_di_report_status (Control, UMFPACK_OK) ;
    if (msglvl==1) Control [UMFPACK_PRL] = 5 ;
    
    Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    //    Control [UMFPACK_PIVOT_TOLERANCE] = 1;
    Control [UMFPACK_ORDERING] = 1; // not 2
    Control [UMFPACK_AGGRESSIVE] = 1;
    Control [UMFPACK_IRSTEP]   = 0;
    Control [UMFPACK_FIXQ] = 1;
    
    /* print the control parameters */
    if (msglvl==1) umfpack_di_report_control (Control) ;
    
    int n = mat->neq;
    
    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_symbolic (n, n, mat->Ic, mat->J, mat->A, &Symbolic,
                                  Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
//        error ("umfpack_di_symbolic failed") ;
    }
    
    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_numeric (mat->Ic, mat->J, mat->A, Symbolic, &Numeric,
                                 Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
//        error ("umfpack_di_numeric failed") ;
    }
    umfpack_di_free_symbolic(&Symbolic);
    
    /* ---------------------------------------------------------------------- */
    /* solve Ax=b */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_solve (UMFPACK_Aat, mat->Ic, mat->J, mat->A, sol, rhs,
                               Numeric, Control, Info) ;
    if (msglvl==1)umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    if (status < 0)
    {
//        error ("umfpack_di_solve failed") ;
    }
    umfpack_di_free_numeric(&Numeric);
    
//    printf("** Time for UMFPACK = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokes_BEN( SparseMat *Stokes, DirectSolver* PardisoStokes ) {
    
    // Call MKL Pardiso
    DirectStokes( Stokes, PardisoStokes, Stokes->b, Stokes->x );
}

