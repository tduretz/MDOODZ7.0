#include "header_MDOODZ.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

void BuildInitialTopography(markers *topo_chain, params model, scale scaling) {

  double Amplitude = 7e3 / scaling.L;
  double Wavelength = 2800e3 / scaling.L;

  for (int k = 0; k < topo_chain->Nb_part; k++) {
    topo_chain->z[k] =
        -Amplitude * cos(2.0 * M_PI * topo_chain->x[k] / Wavelength);
    topo_chain->phase[k] = 0;
  }

  printf("Topographic chain initialised with %d markers\n",
         topo_chain->Nb_part);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz
 * -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles(markers *particles, scale scaling, params model,
                  mat_prop *materials) {

  double z_LAB = -100.0e3 / scaling.L; // lithosphere thickness

  // Loop on particles
  for (int np = 0; np < particles->Nb_part; np++) {

    // Standart initialisation of particles
    particles->phase[np] = 0; // same phase number everywhere
    if (particles->z[np] > z_LAB)
      particles->phase[np] = 1;

    particles->dual[np] = particles->phase[np];

    //--------------------------//
    // SANITY CHECK
    if (particles->phase[np] > model.Nb_phases - 1) {
      printf("Lazy bastard! Fix your particle phase ID! \n");
      exit(144);
    }
    //--------------------------//
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz
 * -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs(grid *mesh, params *model, scale scaling, markers *particles,
            mat_prop *materials, surface *topo) {

  int k, l, c;
  double *X, *Z, *XC, *ZC;
  int NX, NZ, NCX, NCZ;
  double TN = 273.15 / scaling.T;
  double TW = (1330. + 273.15) / scaling.T, TE = (1330. + 273.15) / scaling.T;

  // Set unresaonable conductivity in the mantle to generate and adiabatic
  // mantle as initial condition
  int phase_ast1 = 2; // This has to be the ASTHENOSPHERE phase number
  int phase_ast2 = 4;
  double k_crazy = 1000 * materials->k[phase_ast1];

  if (model->step == 0) {
    materials->k_eff[phase_ast1] = k_crazy;
    materials->k_eff[phase_ast2] = k_crazy;
    printf("Running with crazy conductivity for the asthenosphere!!\n");
  } else {
    materials->k_eff[phase_ast1] = materials->k[phase_ast1];
    materials->k_eff[phase_ast2] = materials->k[phase_ast2];
    printf("Running with normal conductivity for the asthenosphere...\n");
  }

  // Some stuff...
  NX = mesh->Nx;
  NZ = mesh->Nz;
  NCX = NX - 1;
  NCZ = NZ - 1;

  X = malloc(NX * sizeof(double));
  Z = malloc(NZ * sizeof(double));
  XC = malloc(NCX * sizeof(double));
  ZC = malloc(NCZ * sizeof(double));

  for (k = 0; k < NX; k++) {
    X[k] = mesh->xg_coord[k];
  }
  for (k = 0; k < NCX; k++) {
    XC[k] = mesh->xc_coord[k];
  }
  for (l = 0; l < NZ; l++) {
    Z[l] = mesh->zg_coord[l];
  }
  for (l = 0; l < NCZ; l++) {
    ZC[l] = mesh->zc_coord[l];
  }

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

  for (l = 0; l < mesh->Nz + 1; l++) {
    for (k = 0; k < mesh->Nx; k++) {

      c = k + l * (mesh->Nx);

      if (mesh->BCu.type[c] != 30) {

        // Internal points:  -1
        mesh->BCu.type[c] = -1;
        mesh->BCu.val[c] = 0;

        // Matching BC nodes WEST
        if (k == 0) {
          mesh->BCu.type[c] = 0;
          mesh->BCu.val[c] = 0.0;
        }

        // Matching BC nodes EAST
        if (k == mesh->Nx - 1) {
          mesh->BCu.type[c] = 0;
          mesh->BCu.val[c] = 0.0;
        }

        // No slip SOUTH
        if (l == 0) {
          mesh->BCu.type[c] = 11;
          mesh->BCu.val[c] = 0;
        }

        // Free slip NORTH
        if (l == mesh->Nz) {
          mesh->BCu.type[c] = 13;
          mesh->BCu.val[c] = 0;
        }
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

  for (l = 0; l < mesh->Nz; l++) {
    for (k = 0; k < mesh->Nx + 1; k++) {

      c = k + l * (mesh->Nx + 1);

      if (mesh->BCv.type[c] != 30) {

        // Internal points:  -1
        mesh->BCv.type[c] = -1;
        mesh->BCv.val[c] = 0;

        // Matching BC nodes SOUTH
        if (l == 0) {
          mesh->BCv.type[c] = 0;
          mesh->BCv.val[c] = 0.0;
        }

        // Matching BC nodes NORTH
        if (l == mesh->Nz - 1) {
          mesh->BCv.type[c] = 0;
          mesh->BCv.val[c] = 0.0;
        }

        // Non-matching boundary WEST
        if ((k == 0)) {
          mesh->BCv.type[c] = 13;
          mesh->BCv.val[c] = 0.0;
        }

        // Non-matching boundary EAST
        if ((k == mesh->Nx)) {
          mesh->BCv.type[c] = 13;
          mesh->BCv.val[c] = 0.0;
        }
      }
    }
  }

  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for P on all grid levels */
  /* Type  0: Dirichlet within the grid */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* Type 31: surface pressure (Dirichlet) */
  /* --------------------------------------------------------------------------------------------------------*/

  NX = mesh->Nx;
  NZ = mesh->Nz;
  NCX = NX - 1;
  NCZ = NZ - 1;

  for (l = 0; l < NCZ; l++) {
    for (k = 0; k < NCX; k++) {

      c = k + l * (NCX);

      if (mesh->BCt.type[c] != 30) {

        // Internal points:  -1
        mesh->BCp.type[c] = -1;
        mesh->BCp.val[c] = 0;

        if ((k == 0 || k == NCX - 1) && l == NCZ - 1) {
          mesh->BCp.type[c] = 0;
          mesh->BCp.val[c] = 0;
        }
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

  NX = mesh->Nx;
  NZ = mesh->Nz;
  NCX = NX - 1;
  NCZ = NZ - 1;

  for (l = 0; l < mesh->Nz - 1; l++) {
    for (k = 0; k < mesh->Nx - 1; k++) {

      c = k + l * (NCX);

      if (mesh->BCt.type[c] != 30) {

        // WEST
        if (k == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.typW[l] = 0;
          mesh->BCt.valW[l] = TW;
        }

        // EAST
        if (k == NCX - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.typE[l] = 0;
          mesh->BCt.valE[l] = TE;
        }

        // SOUTH
        if (l == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.typS[k] = 1;
          mesh->BCt.valS[k] = mesh->T[c];
        }

        // NORTH
        if (l == NCZ - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.typN[k] = 1;
          mesh->BCt.valN[k] = TN;
        }

        // FREE SURFACE
        else {
          if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 ||
               mesh->BCt.type[c] == 0) &&
              mesh->BCt.type[c + NCX] == 30) {
            mesh->BCt.type[c] = 1;
            mesh->BCt.val[c] = TN;
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

int main(int nargs, char *args[]) {
  char *setupFileName = GetSetupFileName(nargs, args);
  RunMDOODZ(setupFileName, BuildInitialTopography, SetParticles, SetBCs);
}