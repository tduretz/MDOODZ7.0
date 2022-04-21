#include "math.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"

double SetZCoord(MdoodzInstance *instance, double x_coord) {
  const double A = instance->model.zmax / 2.0;
  const double s = instance->model.zmax / 4.0;
  return A * exp(-pow(x_coord, 2) / 2.0 / s / s);
}

int SetPhase(MdoodzInstance *instance, double x_coord) {
  return 0;
}

void BuildInitialTopography(MdoodzInstance *instance, markers *topo_chain) {
  for (int k = 0; k < topo_chain->Nb_part; k++) {
    const x_coord        = topo_chain->x[k];
    topo_chain->z[k]     = SetZCoord(instance, x_coord);
    topo_chain->phase[k] = SetPhase(instance, x_coord);
  }
  printf("Topographic chain initialised with %d markers\n",
         topo_chain->Nb_part);
}

double SetHorizontalStrainRate(MdoodzInstance *instance, double x, double z) {
  return -1.0 * x * instance->model.EpsBG;
}

double SetVerticalStrainRate(MdoodzInstance *instance, double x, double z) {
  return z * instance->model.EpsBG;
}

int P_SetPhase(MdoodzInstance *instance, double x, double z) {
  double xc = 0.0, zc = 0.0;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double X      = x - xc;
  const double Z      = z - zc;
  if (X * X + Z * Z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetGrainSize(MdoodzInstance *instance, double x, double z) {
  return 0;
}

double SetPorosity(MdoodzInstance *instance, double x, double z) {
  return 0;
}

double SetDensity(MdoodzInstance *instance, double x, double z) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  const int phase     = P_SetPhase(instance, x, z);
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] *
           (1 - instance->materials.alp[phase] *
                        (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

double SetTemperature(MdoodzInstance *instance, double x, double z) {
  return (instance->model.user0 + zeroC) / instance->scaling.T;
}

void SetParticles(MdoodzInstance *instance, markers *particles) {
  for (int np = 0; np < particles->Nb_part; np++) {
    particles->Vx[np]    = SetHorizontalStrainRate(instance, particles->x[np], particles->z[np]);
    particles->Vz[np]    = SetVerticalStrainRate(instance, particles->x[np], particles->z[np]);
    particles->phase[np] = P_SetPhase(instance, particles->x[np], particles->z[np]);
    particles->d[np]     = SetGrainSize(instance, particles->x[np], particles->z[np]);
    particles->phi[np]   = SetPorosity(instance, particles->x[np], particles->z[np]);
    particles->rho[np]   = SetDensity(instance, particles->x[np], particles->z[np]);
    particles->T[np]     = SetTemperature(instance, particles->x[np], particles->z[np]);

    if (particles->phase[np] > instance->model.Nb_phases) {
      printf("Lazy bastard! Fix your particle phase ID! \n");
      exit(144);
    }
  }
}

void SetBCs(MdoodzInstance *instance, grid *mesh) {
  const double Lx = (double) (instance->model.xmax - instance->model.xmin);
  const double Lz = (double) (instance->model.zmax - instance->model.zmin);

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
        if (instance->model.shear_style == 0) {
          if (k == 0) {
            // Matching BC nodes WEST
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c]  = -mesh->xg_coord[k] * instance->model.EpsBG;
          } else if (k == mesh->Nx - 1) {
            // Matching BC nodes EAST
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c]  = -mesh->xg_coord[k] * instance->model.EpsBG;
          } else if (l == 0) {
            // Free slip SOUTH
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c]  = 0;
          } else if (l == mesh->Nz) {
            // Free slip NORTH
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c]  = 0;
          } else {
            // Internal points:  -1
            mesh->BCu.type[c] = -1;
            mesh->BCu.val[c]  = 0;
          }
        } else if (instance->model.shear_style == 1) {
          if (k == 0) {
            // Matching BC nodes WEST
            mesh->BCu.type[c] = -2;
            mesh->BCu.val[c]  = 0.0 * instance->model.EpsBG * Lx;
          } else if (k == mesh->Nx - 1) {
            // Matching BC nodes EAST
            mesh->BCu.type[c] = -12;
            mesh->BCu.val[c]  = -0.0 * instance->model.EpsBG * Lx;
          } else if (l == 0) {
            // Free slip S
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c]  = -1 * instance->model.EpsBG * Lz;
          } else if (l == mesh->Nz) {
            // Free slip N
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c]  = 1 * instance->model.EpsBG * Lz;
          } else {
            // Internal points:  -1
            mesh->BCu.type[c] = -1;
            mesh->BCu.val[c]  = 0;
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

  for (int l = 0; l < mesh->Nz; l++) {
    for (int k = 0; k < mesh->Nx + 1; k++) {
      const int c = k + l * (mesh->Nx + 1);
      if (mesh->BCv.type[c] != 30) {
        if (instance->model.shear_style == 0) {
          if (l == 0) {
            // Matching BC nodes SOUTH
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c]  = mesh->zg_coord[l] * instance->model.EpsBG;
          } else if (l == mesh->Nz - 1) {
            // Matching BC nodes NORTH
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c]  = mesh->zg_coord[l] * instance->model.EpsBG;
          } else if (k == 0) {
            // Non-matching boundary WEST
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c]  = 0;
          } else if (k == mesh->Nx) {
            // Non-matching boundary EAST
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c]  = 0;
          } else {
            // Internal points:  -1
            mesh->BCv.type[c] = -1;
            mesh->BCv.val[c]  = 0;
          }
        } else if (instance->model.shear_style == 1) {
          if (l == 0) {
            // Matching BC nodes SOUTH
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c]  = -0.0 * instance->model.EpsBG * Lz;
          } else if (l == mesh->Nz - 1) {
            // Matching BC nodes NORTH
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c]  = 0.0 * instance->model.EpsBG * Lz;
          } else if (k == 0) {
            // Non-matching boundary points
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c]  = 0;
          } else if (k == mesh->Nx) {
            // Non-matching boundary points
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c]  = 0;
          } else {
            // Internal points:  -1
            mesh->BCv.type[c] = -1;
            mesh->BCv.val[c]  = 0;
          }
        }
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
      if (mesh->BCt.type[c] != 30) {
        // Internal points:  -1
        mesh->BCp.type[c] = -1;
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

  const double Ttop = 273.15 / instance->scaling.T;

  for (int l = 0; l < mesh->Nz - 1; l++) {
    for (int k = 0; k < mesh->Nx - 1; k++) {
      const int c = k + l * (NCX);
      if (mesh->BCt.type[c] != 30) {
        if (k == 0) {
          // LEFT
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c]  = mesh->T[c];
        } else if (k == NCX - 1) {
          // RIGHT
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c]  = mesh->T[c];
        } else if (l == 0) {
          // BOT
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c]  = mesh->T[c];
        } else if (l == NCZ - 1) {
          // TOP
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c]  = mesh->T[c];
        } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 ||
                    mesh->BCt.type[c] == 0) &&
                   mesh->BCt.type[c + NCX] == 30) {
          // FREE SURFACE
          mesh->BCt.type[c] = 1;
          mesh->BCt.val[c]  = Ttop;
        }
      }
    }
  }

  printf("Velocity and pressure were initialised\n");
  printf("Boundary conditions were set up\n");
}
