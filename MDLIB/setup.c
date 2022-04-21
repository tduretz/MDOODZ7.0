#include "mdoodz-private.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"


void BuildInitialTopography(MdoodzInstance *instance, markers *topo_chain) {
  BuildInitialTopography_ff buildInitialTopography = instance->BuildInitialTopography;
  for (int k = 0; k < topo_chain->Nb_part; k++) {
    const double x_coord = topo_chain->x[k];
    topo_chain->z[k]     = buildInitialTopography.SetSurfaceZCoord(instance, x_coord);
    topo_chain->phase[k] = buildInitialTopography.SetSurfacePhase(instance, x_coord);
  }
  printf("Topographic chain initialised with %d markers\n",
         topo_chain->Nb_part);
}

void SetParticles(MdoodzInstance *instance, markers *particles) {
  SetParticles_ff setParticles = instance->SetParticles;
  for (int np = 0; np < particles->Nb_part; np++) {
    Coordinates coordinates = {
            .x = particles->x[np],
            .z = particles->z[np]};
    particles->Vx[np]    = setParticles.SetHorizontalStrainRate(instance, coordinates);
    particles->Vz[np]    = setParticles.SetVerticalStrainRate(instance, coordinates);
    particles->phase[np] = setParticles.SetPhase(instance, coordinates);
    particles->d[np]     = setParticles.SetGrainSize(instance, coordinates);
    particles->phi[np]   = setParticles.SetPorosity(instance, coordinates);
    particles->rho[np]   = setParticles.SetDensity(instance, coordinates);
    particles->T[np]     = setParticles.SetTemperature(instance, coordinates);
    particles->X[np]     = setParticles.SetXComponent(instance, coordinates);

    if (particles->phase[np] > instance->model.Nb_phases) {
      printf("Lazy bastard! Fix your particle phase ID! \n");
      exit(144);
    }
  }
}

void SetBCs(MdoodzInstance *instance, grid *mesh) {
  SetBCs_ff setBCs = instance->SetBCs;
  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for Vx on all grid levels */
  /* Type  0: Dirichlet point that matches the physical boundary (Vx:
   * left/right, Vz: bottom/top)            */
  /* Type 11: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)       */
  /* Type  2: Neumann point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)         */
  /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)              */
  /* Type -2: periodic in the x direction (matches the physical boundary) */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < mesh->Nz + 1; l++) {
    for (int k = 0; k < mesh->Nx; k++) {
      const int c = k + l * (mesh->Nx);
      if (mesh->BCu.type[c] != 30) {
        POSITION position;
        if (k == 0) {
          position = LEFT;
        } else if (k == mesh->Nx - 1) {
          position = RIGHT;
        } else if (l == 0) {
          position = BOTTOM;
        } else if (l == mesh->Nz) {
          position = TOP;
        } else {
          position = INTERNAL;
        }
        mesh->BCu.type[c]       = setBCs.SetBCVxType(instance, position);
        Coordinates coordinates = {
                .x = mesh->xg_coord[k],
                .z = mesh->zg_coord[l]};
        mesh->BCu.val[c] = setBCs.SetBCVxValue(instance, position, coordinates);
      }
    }
  }

  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for Vz on all grid levels */
  /* Type  0: Dirichlet point that matches the physical boundary (Vx:
   * left/right, Vz: bottom/top)            */
  /* Type 11: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)       */
  /* Type  2: Neumann point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)         */
  /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)              */
  /* Type -2: periodic in the x direction (does not match the physical boundary)
   */
  /* Type-10: useless point (set to zero) */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < mesh->Nz; l++) {
    for (int k = 0; k < mesh->Nx + 1; k++) {
      const int c = k + l * (mesh->Nx + 1);
      if (mesh->BCv.type[c] != 30) {
        POSITION position;
        if (l == 0) {
          position = BOTTOM;
        } else if (l == mesh->Nz - 1) {
          position = TOP;
        } else if (k == 0) {
          position = LEFT;
        } else if (k == mesh->Nx) {
          position = RIGHT;
        } else {
          position = INTERNAL;
        }
        mesh->BCv.type[c]       = setBCs.SetBCVzType(instance, position);
        Coordinates coordinates = {
                .x = mesh->xg_coord[k],
                .z = mesh->zg_coord[l]};
        mesh->BCv.val[c] = setBCs.SetBCVzValue(instance, position, coordinates);
      }
    }
  }

  const int NCX = mesh->Nx - 1;
  const int NCZ = mesh->Nz - 1;

  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for P on all grid levels */
  /* Type  0: Dirichlet within the grid */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* Type 31: surface pressure (Dirichlet) */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < NCZ; l++) {
    for (int k = 0; k < NCX; k++) {
      const int c = k + l * (NCX);
      POSITION  position;
      if (k == 0) {
        position = LEFT;
      } else if (k == NCX - 1) {
        position = RIGHT;
      } else if (l == 0) {
        position = BOTTOM;
      } else if (l == NCZ - 1) {
        position = TOP;
      } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c + NCX] == 30) {
        position = FREE_SURFACE;
      } else {
        position = INTERNAL;
      }
      if (mesh->BCt.type[c] != 30) {
        mesh->BCp.type[c] = setBCs.SetBCPType(instance, position);
        mesh->BCp.val[c]  = 0;
      }
    }
  }

  /* -------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for T on all grid levels */
  /* Type  1: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)      */
  /* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)             */
  /* Type -2: periodic in the x direction (matches the physical boundary) */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* -------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < NCZ; l++) {
    for (int k = 0; k < NCX; k++) {
      const int c = k + l * (NCX);
      POSITION  position;
      if (k == 0) {
        position = LEFT;
      } else if (k == NCX - 1) {
        position = RIGHT;
      } else if (l == 0) {
        position = BOTTOM;
      } else if (l == NCZ - 1) {
        position = TOP;
      } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c + NCX] == 30) {
        position = FREE_SURFACE;
      } else {
        position = INTERNAL;
      }
      if (mesh->BCt.type[c] != 30) {
        mesh->BCt.type[c] = setBCs.SetBCTType(instance, position);
        mesh->BCt.val[c]  = setBCs.SetBCTValue(instance, position, mesh->T[c]);
      }
    }
  }

  printf("Velocity and pressure were initialised\n");
  printf("Boundary conditions were set up\n");
}
