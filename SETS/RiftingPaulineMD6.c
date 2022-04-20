#include "math.h"
#include "mdoodz.h"
#include "stdio.h"

void BuildInitialTopography(markers *topo_chain, params model, scale scaling) {
  const double TopoLevel = -0.0e3 / scaling.L;
  const double h_pert = model.user3 / scaling.L; // perturbation amplitude
  for (int k = 0; k < topo_chain->Nb_part; k++) {
    topo_chain->z[k] = TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 *
                                       cos(2 * M_PI * topo_chain->x[k] /
                                           (model.xmax - model.xmin));
    topo_chain->phase[k] = 0;
  }
}

double GetTemperature(scale scaling, params model, double depth) {
  const double HLit = model.user1 / scaling.L; // lithosphere thickness
  const double TsurfScaled = zeroC / scaling.T;
  const double HLitScaled = HLit / scaling.L;
  const double TmantScaled = (1330.0 + zeroC) / scaling.T;
  const double Tpart =
      ((TmantScaled - TsurfScaled) / HLitScaled) * depth + TsurfScaled;
  if (Tpart > TmantScaled) {
    return TmantScaled;
  } else {
    return Tpart;
  }
}

int GetPhase(scale scaling, params model, double depth) {
  const double HLit = model.user1 / scaling.L;   // lithosphere thickness
  const double HCrust = model.user2 / scaling.L; // crust thickness
  const double h_pert = model.user3 / scaling.L; // perturbation amplitude
  const double xmin = model.xmin;                // xmin
  const double xmax = model.xmax;                // xmax
  const double HLitScaled = HLit / scaling.L;
  const double mohoDepth = h_pert * cos(2 * M_PI * depth / (xmax - xmin));

  if (depth < -HLitScaled) {
    // astheno. M. below base lithos
    return 3;
  }
  if (depth > -HCrust - mohoDepth) {
    // crust above the Moho
    return 1;
  }
  if (model.user4 == 1) {
    if (depth > -HCrust + 2.5e3 / scaling.L - mohoDepth &&
        depth < -HCrust + 4.5e3 / scaling.L - mohoDepth) {
      // marker between 2.5 and 3 km above the Moho
      return 7;
    } else if (depth > -HCrust + 7.0e3 / scaling.L - mohoDepth &&
               depth < -HCrust + 9.0e3 / scaling.L - mohoDepth) {
      // marker between 4.5 and 5 km above the Moho
      return 8;
    }
  }
  return 2;
}

double GetPorosity() { return 0.0; }

double GetGrainSize(mat_prop *materials, scale scaling, params model,
                    double depth) {
  return materials->gs_ref[GetPhase(scaling, model, depth)];
}

double GetHorizontalStrainRate(double x, double EpsBG) { return -x * EpsBG; }

double GetVerticalStrainRate(double z, double EpsBG) { return z * EpsBG; }

void SetParticles(markers *particles, scale scaling, params model,
                  mat_prop *materials) {
  for (int np = 0; np < particles->Nb_part; np++) {
    particles->Vx[np] = GetHorizontalStrainRate(particles->x[np], model.EpsBG);
    particles->Vz[np] = GetVerticalStrainRate(particles->z[np], model.EpsBG);
    particles->phi[np] = GetPorosity(); // zero porosity everywhere
    particles->T[np] = GetTemperature(scaling, model, -particles->z[np]);
    particles->phase[np] = GetPhase(scaling, model, -particles->z[np]);
    particles->d[np] =
        GetGrainSize(materials, scaling, model, -particles->z[np]);
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz
 * -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs(grid *mesh, params *model, scale scaling, markers *particles,
            mat_prop *materials, surface *topo) {
  // Set unresaonable conductivity in the mantle to generate and adiabatic
  // mantle as initial condition
  const int phase_ast = 3; // This has to be the ASTHENOSPHERE phase number
  const double k_crazy = 1000 * materials->k[phase_ast];
  if (model->step == 0) {
    materials->k_eff[phase_ast] = k_crazy;
    printf("Running with crazy conductivity for the asthenosphere!!\n");
  } else {
    materials->k_eff[phase_ast] = materials->k[phase_ast];
    printf("Running with normal conductivity for the asthenosphere...\n");
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

  for (int l = 0; l < mesh->Nz + 1; l++) {
    for (int k = 0; k < mesh->Nx; k++) {
      const int c = k + l * (mesh->Nx);
      if (mesh->BCu.type[c] != 30) {
        if (k == 0) {
          // Matching BC nodes WEST
          mesh->BCu.type[c] = 0;
          mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
        } else if (k == mesh->Nx - 1) {
          // Matching BC nodes EAST
          mesh->BCu.type[c] = 0;
          mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
        } else if (l == 0) {
          // Free slip SOUTH
          mesh->BCu.type[c] = 13;
          mesh->BCu.val[c] = 0;
        } else if (l == mesh->Nz) {
          // Free slip NORTH
          mesh->BCu.type[c] = 13;
          mesh->BCu.val[c] = 0;
        } else {
          // Internal points:  -1
          mesh->BCu.type[c] = -1;
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

  for (int l = 0; l < mesh->Nz; l++) {
    for (int k = 0; k < mesh->Nx + 1; k++) {
      const int c = k + l * (mesh->Nx + 1);
      if (mesh->BCv.type[c] != 30) {
        if (l == 0) {
          // Matching BC nodes SOUTH
          mesh->BCv.type[c] = 0;
          mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
        } else if (l == mesh->Nz - 1) {
          // Matching BC nodes NORTH
          mesh->BCv.type[c] = 0;
          mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
        } else if (k == 0) {
          // Non-matching boundary WEST
          mesh->BCv.type[c] = 13;
          mesh->BCv.val[c] = 0;
        } else if (k == mesh->Nx) {
          // Non-matching boundary EAST
          mesh->BCv.type[c] = 13;
          mesh->BCv.val[c] = 0;
        } else {
          // Internal points:  -1
          mesh->BCv.type[c] = -1;
          mesh->BCv.val[c] = 0;
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

  for (int l = 0; l < mesh->Nz - 1; l++) {
    for (int k = 0; k < mesh->Nx - 1; k++) {
      const int c = k + l * (mesh->Nx - 1);
      if (mesh->BCt.type[c] != 30) {
        if ((k == 0 || k == mesh->Nx - 1 - 1) && l == mesh->Nz - 1 - 1) {
          mesh->BCp.type[c] = 0;
          mesh->BCp.val[c] = 0;
        } else {
          // Internal points:  -1
          mesh->BCp.type[c] = -1;
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

  const double TN = 273.15 / scaling.T;
  const double TW = (1330. + 273.15) / scaling.T,
               TE = (1330. + 273.15) / scaling.T;

  for (int l = 0; l < mesh->Nz - 1; l++) {
    for (int k = 0; k < mesh->Nx - 1; k++) {
      const int c = k + l * (mesh->Nx - 1);
      if (mesh->BCt.type[c] != 30) {
        if (k == 0) {
          // WEST
          mesh->BCt.type[c] = 0;
          mesh->BCt.typW[l] = 0;
          mesh->BCt.valW[l] = TW;
        } else if (k == mesh->Nx - 1 - 1) {
          // EAST
          mesh->BCt.type[c] = 0;
          mesh->BCt.typE[l] = 0;
          mesh->BCt.valE[l] = TE;
        } else if (l == 0) {
          // SOUTH
          mesh->BCt.type[c] = 0;
          mesh->BCt.typS[k] = 1;
          mesh->BCt.valS[k] = mesh->T[c];
        } else if (l == mesh->Nz - 1 - 1) {
          // NORTH
          mesh->BCt.type[c] = 0;
          mesh->BCt.typN[k] = 1;
          mesh->BCt.valN[k] = TN;
        } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 ||
                    mesh->BCt.type[c] == 0) &&
                   mesh->BCt.type[c + mesh->Nx - 1] == 30) {
          // FREE SURFACE
          mesh->BCt.type[c] = 1;
          mesh->BCt.val[c] = TN;
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
