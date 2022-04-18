#include "math.h"
#include "mdoodz.h"

void BuildInitialTopography(markers *topo_chain, params model, scale scaling) {}

void SetParticles(markers *particles, scale scaling, params model,
                  mat_prop *materials) {
  const double T_init = (model.user0 + zeroC) / scaling.T;
  const double radius = model.user1 / scaling.L;
  for (int np = 0; np < particles->Nb_part; np++) {
    particles->T[np] = T_init;
    if ((pow(particles->x[np], 2) + pow(particles->z[np], 2)) <
        pow(radius, 2)) {
      particles->phase[np] = 1;
    } else {
      particles->phase[np] = 0;
    }
    particles->dual[np] = particles->phase[np];
  }
}

void SetBCs(grid *mesh, params *model, scale scaling, markers *particles,
            mat_prop *materials, surface *topo) {
  const int Nx = mesh->Nx;
  const int Nz = mesh->Nz;
  // Set the BCs for Vx on all grid levels
  for (int j = 0; j < Nz + 1; j++) {
    for (int i = 0; i < Nx; i++) {
      const int c = i + j * Nx;
      if (i == 0 || i == Nx - 1) {
        // West & East boundary
        mesh->BCu.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCu.val[c] = -mesh->xg_coord[i] * model->EpsBG;
      } else if (j == 0 || j == Nz) {
        // North & South boundary
        mesh->BCu.type[c] = BOUNDARY_NEUMANN;
        mesh->BCu.val[c] = 0.0;
      } else {
        mesh->BCu.type[c] = NOT_BC;
        mesh->BCu.val[c] = 0;
      }
    }
  }

  // Set the BCs for Vz on all grid levels
  for (int j = 0; j < Nz; j++) {
    for (int i = 0; i < Nx + 1; i++) {
      const int c = (Nx + 1) * j + i;
      if (j == 0 || j == (Nz - 1)) {
        // North & South boundary
        mesh->BCv.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCv.val[c] =
            mesh->zg_coord[j] * model->EpsBG; // corrected valgrind!
      } else if (i == 0 || i == Nx) {
        // West & East boundary
        mesh->BCv.type[c] = BOUNDARY_NEUMANN;
        mesh->BCv.val[c] = 0.0;
      } else {
        // Internal points:  -1
        mesh->BCv.type[c] = NOT_BC;
        mesh->BCv.val[c] = 0;
      }
    }
  }

  // Set the BCs for P on all grid levels
  for (int j = 0; j < Nz - 1; j++) {
    for (int i = 0; i < Nx - 1; i++) {
      const int c = (Nx - 1) * j + i;
      mesh->BCp.type[c] = NOT_BC;
      mesh->BCp.val[c] = 0;
    }
  }

  // Set the BCs for T on all grid levels
  for (int j = 0; j < Nz - 1; j++) {
    for (int i = 0; i < Nx - 1; i++) {
      const int c = (Nx - 1) * j + i;
      if (i == 0) {
        mesh->BCt.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCt.typW[j] = BOUNDARY_DIRICHLET;
        mesh->BCt.valW[j] = 0.0;
      } else if (i == Nx - 2) {
        mesh->BCt.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCt.typE[j] = BOUNDARY_DIRICHLET;
        mesh->BCt.valE[j] = 0.0;
      } else if (j == 0) {
        mesh->BCt.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCt.typS[i] = BOUNDARY_DIRICHLET;
        mesh->BCt.valS[i] = 0.0;
      } else if (j == Nz - 2) {
        mesh->BCt.type[c] = BOUNDARY_DIRICHLET;
        mesh->BCt.typN[i] = BOUNDARY_DIRICHLET;
        mesh->BCt.valN[i] = 0.0;
      }
    }
  }
}

int main(int nargs, char *args[]) {
  char *setupFileName = GetSetupFileName(nargs, args);
  RunMDOODZ(setupFileName, BuildInitialTopography, SetParticles, SetBCs);
}
