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
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {

    // Define dimensions;
    const double T_init = (model.user0 + zeroC)/scaling.T;
    const double radius = model.user1/scaling.L;

    // Loop on particles
    for ( int ip=0; ip<particles->Nb_part; ip++ ) {
        particles->phase[ip]    = 0;
        particles->T[ip]        = T_init;
        if ( ( pow(particles->x[ip],2) + pow(particles->z[ip],2) ) < pow(radius,2) ) {
            particles->phase[ip] = 1;
        }
        // SET DEFAULT DUAL PHASE NUMBER used for visualisation of phase on grid
        particles->dual[ip] = particles->phase[ip];
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {


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
        
    const int Nx  = mesh->Nx;
    const int Nz  = mesh->Nz;

    // Set the BCs for Vx on all grid levels
    for (int j=0; j < Nz+1; j++) {
        for (int i=0; i < Nx; i++) {
            // global index for Vx point
            int c = i + j*Nx;

            // Internal points:  -1
            mesh->BCu.type[c] = -1;
            mesh->BCu.val[c]  =  0;

            // West & East boundary
            if ( i == 0 || i == Nx-1) { 
                mesh->BCu.type[c] = 0;
                mesh->BCu.val[c]  = -mesh->xg_coord[i] * model->EpsBG;
            }

            // North & South boundary
            if ( j == 0 || j == Nz) { 
                mesh->BCu.type[c] = 13;
                mesh->BCu.val[c]  = 0.0;
            }

        }
    }

    // Set the BCs for Vz on all grid levels
    for (int j=0; j < Nz; j++) {
        for (int i=0; i < Nx+1; i++) {
        int c = (Nx + 1) * j + i;

        // Internal points:  -1
        mesh->BCv.type[c] = -1;
        mesh->BCv.val[c] = 0;

        // North & South boundary
        if (j == 0 || j == (Nz - 1)) {
          mesh->BCv.type[c] = 0;
          mesh->BCv.val[c] = mesh->zg_coord[j] * model->EpsBG; // corrected valgrind!
        }

        // West & East boundary
        if (i == 0 || i == Nx) {
          mesh->BCv.type[c] = 13;
          mesh->BCv.val[c] = 0.0;
        }
      }
    }

    // Set the BCs for P on all grid levels
    for (int j=0; j < Nz-1; j++) {
        for (int i=0; i < Nx-1; i++) {
        int c = (Nx - 1) * j + i;
        // Internal points:  -1
        mesh->BCp.type[c] = -1; // BCp is not used
        mesh->BCp.val[c] = 0;
      }
    }

    // Set the BCs for T on all grid levels
    for (int j=0; j < Nz-1; j++) {
        for (int i=0; i < Nx-1; i++) {
        int c = (Nx - 1) * j + i;

        if (i == 0) {
          mesh->BCt.type[c] = 0; // BCp is not used
          mesh->BCt.typW[j] = 0; // 0 corresponds to Neumann BC and 1 corresponds to Dirichlet
          mesh->BCt.valW[j] = 0.0;
        }

        if (i == Nx - 2) {
          mesh->BCt.type[c] = 0; // BCp is not used
          mesh->BCt.typE[j] = 0;
          mesh->BCt.valE[j] = 0.0;
        }

        if (j == 0) {
          mesh->BCt.type[c] = 0; // BCp is not used
          mesh->BCt.typS[i] = 0;
          mesh->BCt.valS[i] = 0.0;
        }

        if (j == Nz - 2) {
          mesh->BCt.type[c] = 0; // BCp is not used
          mesh->BCt.typN[i] = 0;
          mesh->BCt.valN[i] = 0.0;
        }
      }
    }
}
