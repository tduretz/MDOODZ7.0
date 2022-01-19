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
#include "math.h"
#include "stdlib.h"
#include "header_MDOODZ.h"
#include "time.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CountPartCell2 ( markers* particles, grid *mesh, params model, surface topo, surface topo_ini, int reseed, scale scaling ) {

    // This function counts the number of particle that are currently in each cell of the domain.
    // The function detects cells that are lacking of particle and call the particle re-seeding routine.

    int k, l, ic, jc, kc, nthreads, thread_num, **NPPC, **NPPV;
    int ncx=mesh->Nx-1, ncz=mesh->Nz-1, Nx=mesh->Nx, Nz=mesh->Nz, Nb_part=particles->Nb_part;
    double distance;
    int nb_part_reuse=0, *ind_part_reuse, inc_reuse=0;
    double dx=mesh->dx, dz=mesh->dz;
    int **ind_per_cell, *ind_part;
    int kk = 0;
    int jj, n_neigh=0, n_neigh0, ic2, *neigh_part, ii, p;

    printf("USING CountPartCell2\n");

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
#pragma omp parallel for shared ( particles, model, NPPV, Nb_part, PHASE_s  )    reduction(+:nb_part_reuse)        \
private ( k, distance, ic, jc, thread_num )                    \
firstprivate( ncx, mesh ) schedule( static )
    for (k=0; k<Nb_part; k++) {

        thread_num = omp_get_thread_num();

        if ( particles->phase[k] != -1 ) {

            // -------- Check NODES ----------

            // Get column/line indices
            distance = (particles->x[k]-mesh->xg_coord[0]);
            ic       = ceil((distance/dx)+0.5) -1;
            distance = (particles->z[k]-mesh->zg_coord[0]);
            jc       = ceil((distance/dz)+0.5) -1;

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

    printf ("Before NODE re-seeding, nb_part_reuse: %d\n", nb_part_reuse);

    // Final reduction for NODES
    for ( l=0; l<Nx*Nz; l++ ) {
        mesh->nb_part_vert[l] = 0;
        //        if ( mesh->BCg.type[l] != 30 ) {
        for ( k=0; k<nthreads; k++ ) {
            //            for ( p=0; p<model.Nb_phases; p++ ) {
            //                mesh->phase_perc_s[p][l] += PHASE_s[k][p][l];
            //            }
            mesh->nb_part_vert[l]   += NPPV[k][l];
        }
        //        }
    }

    // Initialise Phase arrays
#pragma omp parallel for shared ( mesh, Nx, Nz, model ) private( p, l ) schedule( static )
    for ( l=0; l<Nx*Nz; l++ ) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            if (mesh->nb_part_vert[l]>0) mesh->phase_perc_s[p][l] = 0;
        }
    }


    // Final reduction for NODES
    for ( l=0; l<Nx*Nz; l++ ) {
        //        if ( mesh->BCg.type[l] != 30 ) {
        for ( k=0; k<nthreads; k++ ) {
            for ( p=0; p<model.Nb_phases; p++ ) {
                if (mesh->nb_part_vert[l]>0) mesh->phase_perc_s[p][l] += PHASE_s[k][p][l];
            }
            //            mesh->nb_part_vert[l]   += NPPV[k][l];
        }
        //        }
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
    for (k=0; k<Nb_part; k++) {
        if ( particles->phase[k] == -1 ) {
            ind_part_reuse[kk] = k;
            kk++;
        }
    }

    //-------------------------------------------------------------------------------------------------------------------------//

    // LOOP ON PARTICLES, GIVE THEIR ID'S TO THE CELL THAT OWNS THEM
    for (k=0; k<Nb_part; k++) {

        if ( particles->phase[k] !=-1 ) {

            // Get the column:
            distance=(particles->x[k]-mesh->xg_coord[0]);
            ic=ceil((distance/dx)+0.5) -1;

            // Get the line:
            distance=(particles->z[k]-mesh->zg_coord[0]);
            jc=ceil((distance/dz)+0.5) -1;

            // Global index
            kc = ic + jc * (mesh->Nx);
            ind_per_cell[kc][ind_part[kc]] = k;

            //#pragma omp atomic
            (ind_part[kc])++;
        }
    }

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

    // ----------------

    // ----------------

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
            distance = (particles->x[k] - mesh->xc_coord[0]);
            ic       = ceil((distance/mesh->dx)+0.5) - 1;
            distance = (particles->z[k] - mesh->zc_coord[0]);
            jc       = ceil((distance/mesh->dz)+0.5) - 1;

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


    // Final reduction for CELLS
#pragma omp parallel for shared ( mesh, NPPC, ncx, ncz, nthreads, PHASE_n, model ) private( p, l, k ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        mesh->nb_part_cell[l] = 0.0;
        for ( k=0; k<nthreads; k++ ) {
            mesh->nb_part_cell[l]   += NPPC[k][l];
        }
    }


    // Initialise Phase arrays
#pragma omp parallel for shared ( mesh, ncx, ncz, model ) private( p, l ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            if (mesh->nb_part_cell[l]>0) mesh->phase_perc_n[p][l] = 0;
        }
    }

    // Final reduction for CELLS
#pragma omp parallel for shared ( mesh, NPPC, ncx, ncz, nthreads, PHASE_n, model ) private( p, l, k ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {

        for ( k=0; k<nthreads; k++ ) {
            for ( p=0; p<model.Nb_phases; p++ ) {
                if (mesh->nb_part_cell[l]>0) mesh->phase_perc_n[p][l] += PHASE_n[k][p][l];
            }
        }
    }

    // Normalise to phase percentages
#pragma omp parallel for shared ( mesh, ncx, ncz, model ) private( p, l ) schedule( static )
    for ( l=0; l<ncx*ncz; l++ ) {
        //        if ( mesh->BCp.type[0][l] != 30 && mesh->BCp.type[0][l] != 31) {
        for ( p=0; p<model.Nb_phases; p++ ) {
            if (mesh->nb_part_cell[l]>0) mesh->phase_perc_n[p][l] /= mesh->nb_part_cell[l];
        }
    }

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
    for (k=0; k<Nb_part; k++) {

        if ( particles->phase[k] != -1 ) {
            distance=(particles->x[k] - mesh->xc_coord[0]);
            ic = ceil((distance/mesh->dx)+0.5) - 1;
            distance=(particles->z[k] - mesh->zc_coord[0]);
            jc = ceil((distance/mesh->dz)+0.5) - 1;

            kc = ic + jc * (mesh->Nx-1);

            ind_per_cell[kc][ind_part[kc]] = k;
            ind_part[kc] ++;
        }
    }

    double xC, zC, dst, dst_min;
    int    Np = particles->Nb_part;
    int    neigh, min_index, new_index;
    int *ind_list, neighs, oo, nb, deact;

    int Np_add  = 0;

    if (reseed == 1) {



        // Loop on cells and add particles
        for (k=0; k<ncx; k++) {
            for (l=0; l<ncz; l++) {
                ic = k + l * ncx;

                if ( mesh->BCt.type[ic] != 30 && mesh->nb_part_cell[ic] < particles -> min_part_cell ) {

                    //----------------

                    neighs = 0;
                    oo     = 0;

                    // list the number of particles in the neighbouring cells
                    neighs              += mesh->nb_part_cell[ic];

                    if (neighs == 0) {
                                           printf("All the neighbouring CELLS of ix = %d iz = %d are empty, simulation will stop\n", k, l);
                                           exit(1);
                                       }
                    if (k>0     && mesh->BCt.type[ic-1  ] != 30) neighs += mesh->nb_part_cell[ic-1];
                    if (k<ncx-1 && mesh->BCt.type[ic+1  ] != 30) neighs += mesh->nb_part_cell[ic+1];
                    if (l>0     && mesh->BCt.type[ic-ncx] != 30) neighs += mesh->nb_part_cell[ic-ncx];
                    if (l<ncz-1 && mesh->BCt.type[ic+ncx] != 30) neighs += mesh->nb_part_cell[ic+ncx];

                    // Diagonal neighbours
                    neighs += mesh->nb_part_cell[ic-1-ncx];
                    neighs += mesh->nb_part_cell[ic+1-ncx];
                    neighs += mesh->nb_part_cell[ic-1+ncx];
                    neighs += mesh->nb_part_cell[ic+1+ncx];
                    // Extended neighbours
                    neighs += mesh->nb_part_cell[ic-2];
                    neighs += mesh->nb_part_cell[ic+2];
                    neighs += mesh->nb_part_cell[ic-2*ncx];
                    neighs += mesh->nb_part_cell[ic+2*ncx];

                    if (neighs == 0) {
                        printf("All the neighbouring CELLS of ix = %d iz = %d are empty, simulation will stop\n", k, l);
                        exit(1);
                    }

                    if (neighs > 0) {

                        ind_list = DoodzMalloc( neighs * sizeof(int));


                        for ( nb=0; nb<mesh->nb_part_cell[ic]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic][nb];
                            oo++;
                        }

                        if (k>0 && mesh->BCt.type[ic-1] != 30) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic-1]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic-1][nb];
                                oo++;
                            }
                        }
                        if (k<ncx-1 && mesh->BCt.type[ic+1] != 30) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic+1]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic+1][nb];
                                oo++;
                            }
                        }
                        if (l>0 && mesh->BCt.type[ic-ncx] != 30) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic-ncx]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic-ncx][nb];
                                oo++;
                            }
                        }
                        if (l<ncz-1 && mesh->BCt.type[ic+ncx] != 30) {
                            for ( nb=0; nb<mesh->nb_part_cell[ic+ncx]; nb++ ) {
                                ind_list[oo] = ind_per_cell[ic+ncx][nb];
                                oo++;
                            }
                        }




                        // Diagonal cells: add particles from neigbours
                        for ( nb=0; nb<mesh->nb_part_cell[ic-1-ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic-1-ncx][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic+1-ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic+1-ncx][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic-1+ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic-1+ncx][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic+1+ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic+1+ncx][nb];
                            oo++;
                        }

                        // Extended neighbours
                        for ( nb=0; nb<mesh->nb_part_cell[ic-2]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic-2][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic+2]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic+2][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic-2*ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic-2*ncx][nb];
                            oo++;
                        }

                        for ( nb=0; nb<mesh->nb_part_cell[ic+2*ncx]; nb++ ) {
                            ind_list[oo] = ind_per_cell[ic+2*ncx][nb];
                            oo++;
                        }

                        //                        printf("%d %d %d\n",mesh->nb_part_cell[ic],mesh->nb_part_cell[ic-ncx], neighs);
//
                        //----------------

                        xC = mesh->xc_coord[k];
                        zC = mesh->zc_coord[l];

                        // Center
                        new_index = Np+Np_add;
                        particles->x[new_index] = xC;
                        particles->z[new_index] = zC;

//                        neigh = 0;
//                        dst_min = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
//                        dst_min = model.xmax-model.xmin;
//                        min_index = 0;
//
//                        for ( nb=0; nb<particles->Nb_part; nb++ ) {
//                            neigh = nb;
////                            printf("%d %d\n", nb,oo);
////                            if (nb>oo) exit(67);
//                            dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
//                            if ( dst < dst_min && particles->phase[neigh]!=-1)
//                                dst_min = dst;
//                            min_index = neigh;
//                        }

                        neigh = ind_list[0];
                        dst_min = model.xmax-model.xmin;//sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
                        min_index = ind_list[0];

                        for ( nb=0; nb<neighs; nb++ ) {

                            neigh = ind_list[nb];
//                            if (nb>oo) exit(67);
                            dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
//                            if (nb< mesh->nb_part_cell[ic]) printf("current cell dst %2.2f dst_min %2.2f\n", dst, dst_min);
//                            if (nb>=mesh->nb_part_cell[ic]) printf("below   cell dst %2.2f dst_min %2.2f\n", dst, dst_min);
                            if ( dst < dst_min ) {
                                dst_min = dst;
                                min_index = neigh;
                            }
                        }

                        AssignMarkerProperties ( particles, new_index, min_index, &model, mesh, model.DirectNeighbour );
                        Np_add++;
                        //
                        //                             // SW
                        //                             new_index = Np+Np_add;
                        //                             particles->x[new_index] = xC-dx/4.0;
                        //                             particles->z[new_index] = zC-dz/4.0;
                        //
                        //                             neigh = ind_list[0];
                        //                             dst_min = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
                        //                             min_index = ind_list[0];
                        //
                        ////                             if (dst>sqrt(dx*dx+dz*dz)) {
                        ////                                 printf("WOUHLALAH\n");
                        ////                                 exit(1);
                        ////                             }
                        //
                        //                             for ( nb=0; nb<neighs; nb++ ) {
                        //                                 neigh = ind_list[nb];
                        //                                 dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
                        //                                 if ( dst < dst_min )
                        //                                     dst_min = dst;
                        //                                     min_index = neigh;
                        //                             }
                        //
                        //                             AssignMarkerProperties ( particles, new_index, min_index, &model, mesh, model.DirectNeighbour );
                        //                             Np_add++;
                        //
                        //                             // SE
                        //                             new_index = Np+Np_add;
                        //                             particles->x[new_index] = xC+dx/4.0;
                        //                             particles->z[new_index] = zC-dz/4.0;
                        //
                        //                             neigh = ind_list[0];
                        //                             dst_min = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
                        //                             min_index = ind_list[0];
                        //
                        ////                             if (dst>sqrt(dx*dx+dz*dz)) {
                        ////                                 printf("WOUHLALAH\n");
                        ////                                 exit(1);
                        ////                             }
                        //
                        //                             for ( nb=0; nb<neighs; nb++ ) {
                        //                                 neigh = ind_list[nb];
                        //                                 dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
                        //                                 if ( dst < dst_min )
                        //                                     dst_min = dst;
                        //                                     min_index = neigh;
                        //                             }
                        //
                        //                             AssignMarkerProperties ( particles, new_index, min_index, &model, mesh, model.DirectNeighbour );
                        //                             Np_add++;
                        //
                        //                             // NW
                        //                             new_index = Np+Np_add;
                        //                             particles->x[new_index] = xC-dx/4.0;
                        //                             particles->z[new_index] = zC+dz/4.0;
                        //
                        //                             neigh = ind_list[0];
                        //                             dst_min = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
                        //                             min_index = ind_list[0];
                        //
                        ////                             if (dst>sqrt(dx*dx+dz*dz)) {
                        ////                                 printf("WOUHLALAH\n");
                        ////                                 exit(1);
                        ////                             }
                        //
                        //                             for ( nb=0; nb<neighs; nb++ ) {
                        //                                 neigh = ind_list[nb];
                        //                                 dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
                        //                                 if ( dst < dst_min )
                        //                                     dst_min = dst;
                        //                                     min_index = neigh;
                        //                             }
                        //
                        //                             AssignMarkerProperties ( particles, new_index, min_index, &model, mesh, model.DirectNeighbour );
                        //                             Np_add++;
                        //
                        //                            // NE
                        //                             new_index = Np+Np_add;
                        //                             particles->x[new_index] = xC+dx/4.0;
                        //                             particles->z[new_index] = zC+dz/4.0;
                        //
                        //                             neigh = ind_list[0];
                        //                             dst_min = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0) );
                        //                             min_index = ind_list[0];
                        //
                        ////                             if (dst>sqrt(dx*dx+dz*dz)) {
                        ////                                 printf("WOUHLALAH\n");
                        ////                                 exit(1);
                        ////                             }
                        //
                        //                             for ( nb=0; nb<neighs; nb++ ) {
                        //                                 neigh = ind_list[nb];
                        //                                 dst   = sqrt( pow(particles->x[new_index] - particles->x[neigh], 2.0) + pow(particles->z[new_index] - particles->z[neigh], 2.0));
                        //                                 if ( dst < dst_min )
                        //                                     dst_min = dst;
                        //                                     min_index = neigh;
                        //                             }
                        //
                        //                             AssignMarkerProperties ( particles, new_index, min_index, &model, mesh, model.DirectNeighbour );
                        //                             Np_add++;

                        DoodzFree(ind_list);

                    }

                }

            }
        }


        deact = 0;

        for (k=0; k<ncx; k++) {
            for (l=0; l<ncz; l++) {
                ic = k + l * ncx;

                // Deactivate markers
                if ( mesh->BCt.type[ic] != 30 && mesh->nb_part_cell[ic] > particles -> min_part_cell + 4 ) {
                    // Loop through markers in the current cell and deactivate if more markers than wanted
                    for ( nb=0; nb<mesh->nb_part_cell[ic]; nb++ ) {
                        if (nb > particles -> min_part_cell + 4) {
                            //                             printf("Deleting marker\n");
                            new_index = ind_per_cell[ic][nb];
                            particles->phase[new_index] = -1;
                            deact++;
                        }

                    }
                }
            }
        }



    }

    particles->Nb_part = Np+Np_add;
    printf("Deactivated particles: %03d\n", deact);



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

    MinMaxArrayTagInt  ( mesh->nb_part_vert, 1.0, (mesh->Nx-0)*(mesh->Nz-0), "nb part vert", mesh->BCg.type );
    MinMaxArrayTagInt  ( mesh->nb_part_cell, 1.0, (mesh->Nx-1)*(mesh->Nz-1), "nb part cell", mesh->BCp.type );

}
