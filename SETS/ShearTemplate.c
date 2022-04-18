#include "math.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"

void BuildInitialTopography(markers *topo_chain, params model, scale scaling) {

  int k;
  double A = model.zmax / 2.0;
  double s = model.zmax / 4.0;

  for (k = 0; k < topo_chain->Nb_part; k++) {
    topo_chain->z[k] =
        A * exp(-pow(topo_chain->x[k], 2) / 2.0 / s / s); // TopoLevel;
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
  int np;

  double T_init = (model.user0 + zeroC) / scaling.T;
  double radius = model.user1 / scaling.L;
  double X, Z, xc = 0.0, zc = 0.0;

  // Loop on particles
  for (np = 0; np < particles->Nb_part; np++) {

    // Standart initialisation of particles
    particles->Vx[np] = -1.0 * particles->x[np] *
                        model.EpsBG; // set initial particle velocity (unused)
    particles->Vz[np] = particles->z[np] *
                        model.EpsBG; // set initial particle velocity (unused)
    particles->phase[np] = 0;        // same phase number everywhere
    particles->d[np] = 0;            // same grain size everywhere
    particles->phi[np] = 0.0;        // zero porosity everywhere
    particles->rho[np] = 0;
    particles->T[np] = T_init;
    X = particles->x[np] - xc;
    Z = particles->z[np] - zc;

    // ------------------------- //
    // DRAW INCLUSION
    if (X * X + Z * Z < radius * radius) {
      particles->phase[np] = 1;
    }

    // SANITY CHECK
    if (particles->phase[np] > model.Nb_phases) {
      printf("Lazy bastard! Fix your particle phase ID! \n");
      exit(144);
    }

    //--------------------------//
    // DENSITY
    if (model.eqn_state > 0) {
      particles->rho[np] =
          materials->rho[particles->phase[np]] *
          (1 - materials->alp[particles->phase[np]] *
                   (T_init - materials->T0[particles->phase[np]]));
    } else {
      particles->rho[np] = materials->rho[particles->phase[np]];
    }
    //--------------------------//
  }
  MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init");
  MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init");
  MinMaxArray(particles->T, scaling.T, particles->Nb_part, "Tp init");
}

void SetBCs(grid *mesh, params *model, scale scaling, markers *particles,
            mat_prop *materials, surface *topo) {

  int k, l, c;
  double Lx = (double)(model->xmax - model->xmin);
  double Lz = (double)(model->zmax - model->zmin);

  int NX = mesh->Nx;
  int NZ = mesh->Nz;
  int NCX = NX - 1;
  int NCZ = NZ - 1;

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

        if (model->shear_style == 0) {

          // Matching BC nodes WEST
          if (k == 0) {
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
          }

          // Matching BC nodes EAST
          if (k == mesh->Nx - 1) {
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
          }

          // Free slip SOUTH
          if (l == 0) {
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c] = 0;
          }

          // Free slip NORTH
          if (l == mesh->Nz) {
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c] = 0;
          }
        }
        if (model->shear_style == 1) {

          // Matching BC nodes WEST
          if (k == 0) {
            mesh->BCu.type[c] = -2;
            mesh->BCu.val[c] = 0.0 * model->EpsBG * Lx;
          }

          // Matching BC nodes EAST
          if (k == mesh->Nx - 1) {
            mesh->BCu.type[c] = -12;
            mesh->BCu.val[c] = -0.0 * model->EpsBG * Lx;
          }

          // Free slip S
          if (l == 0) { //&& (k>0 && k<NX-1) ) {
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c] = -1 * model->EpsBG * Lz;
          }

          // Free slip N
          if (l == mesh->Nz) { // && (k>0 && k<NX-1)) {
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c] = 1 * model->EpsBG * Lz;
          }
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

        if (model->shear_style == 0) {

          // Matching BC nodes SOUTH
          if (l == 0) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
          }

          // Matching BC nodes NORTH
          if (l == mesh->Nz - 1) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
          }

          // Non-matching boundary WEST
          if (k == 0) {
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c] = 0;
          }

          // Non-matching boundary EAST
          if (k == mesh->Nx) {
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c] = 0;
          }
        }

        if (model->shear_style == 1) {
          // Matching BC nodes SOUTH
          if (l == 0) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = -0.0 * model->EpsBG * Lz;
          }

          // Matching BC nodes NORTH
          if (l == mesh->Nz - 1) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = 0.0 * model->EpsBG * Lz;
          }

          // Non-matching boundary points
          if (k == 0) {
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c] = 0;
          }

          // Non-matching boundary points
          if (k == mesh->Nx) {
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c] = 0;
          }
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

  for (l = 0; l < NCZ; l++) {
    for (k = 0; k < NCX; k++) {
      c = k + l * (NCX);
      if (mesh->BCt.type[c] != 30) {
        // Internal points:  -1
        mesh->BCp.type[c] = -1;
        mesh->BCp.val[c] = 0;
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

  double Ttop = 273.15 / scaling.T;

  NCX = NX - 1;
  NCZ = NZ - 1;

  for (l = 0; l < mesh->Nz - 1; l++) {
    for (k = 0; k < mesh->Nx - 1; k++) {

      c = k + l * (NCX);

      if (mesh->BCt.type[c] != 30) {

        // LEFT
        if (k == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }

        // RIGHT
        if (k == NCX - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }

        // BOT
        if (l == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }

        // TOP
        if (l == NCZ - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }
        // FREE SURFACE
        else {
          if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 ||
               mesh->BCt.type[c] == 0) &&
              mesh->BCt.type[c + NCX] == 30) {
            mesh->BCt.type[c] = 1;
            mesh->BCt.val[c] = Ttop;
          }
        }
      }
    }
  }

  printf("Velocity and pressure were initialised\n");
  printf("Boundary conditions were set up\n");
}

int main(int nargs, char *args[]) {
  char *setupFileName = GetSetupFileName(nargs, args);
  RunMDOODZ(setupFileName, BuildInitialTopography, SetParticles, SetBCs);
}
