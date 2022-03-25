#include "hdf5.h"
#include "header_MDOODZ.h"
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#define FILENAME "Output00001.gzip.h5"

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

void SetParticles(markers *particles, scale scaling, params model,
                  mat_prop *materials) {
  int np;
  double T_init = (model.user0 + zeroC) / scaling.T;
  double radius = model.user1 / scaling.L;
  double X, Z, xc = 0.0, zc = 0.0;

  for (np = 0; np < particles->Nb_part; np++) {
    particles->Vx[np] = -1.0 * particles->x[np] *
                        model.EpsBG;
    particles->Vz[np] = particles->z[np] *
                        model.EpsBG;
    particles->phase[np] = 0;
    particles->d[np] = 0;
    particles->phi[np] = 0.0;
    particles->rho[np] = 0;
    particles->T[np] = T_init;
    X = particles->x[np] - xc;
    Z = particles->z[np] - zc;

    if (X * X + Z * Z < radius * radius) {
      particles->phase[np] = 1;
    }

    // SANITY CHECK
    if (particles->phase[np] > model.Nb_phases) {
      printf("Lazy bastard! Fix your particle phase ID! \n");
      exit(144);
    }

    if (model.eqn_state > 0) {
      particles->rho[np] =
          materials->rho[particles->phase[np]] *
          (1 - materials->alp[particles->phase[np]] *
                   (T_init - materials->T0[particles->phase[np]]));
    } else {
      particles->rho[np] = materials->rho[particles->phase[np]];
    }
  }
  MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init");
  MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init");
  MinMaxArray(particles->T, scaling.T, particles->Nb_part, "Tp init");
}

void SetBCs(grid *mesh, params *model, scale scaling, markers *particles,
            mat_prop *materials, surface *topo) {

  int k, l, c;
  double *X, *Z, *XC, *ZC;
  int NX, NZ, NCX, NCZ;
  double Lx, Lz;

  // Define dimensions;
  Lx = (double)(model->xmax - model->xmin);
  Lz = (double)(model->zmax - model->zmin);
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

  NX = mesh->Nx;
  NZ = mesh->Nz;

  for (l = 0; l < mesh->Nz + 1; l++) {
    for (k = 0; k < mesh->Nx; k++) {

      c = k + l * (mesh->Nx);

      if (mesh->BCu.type[c] != 30) {
        mesh->BCu.type[c] = -1;
        mesh->BCu.val[c] = 0;

        if (model->shear_style == 0) {
          if (k == 0) {
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
          }

          if (k == mesh->Nx - 1) {
            mesh->BCu.type[c] = 0;
            mesh->BCu.val[c] = -mesh->xg_coord[k] * model->EpsBG;
          }

          if (l == 0) {
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c] = 0;
          }

          if (l == mesh->Nz) {
            mesh->BCu.type[c] = 13;
            mesh->BCu.val[c] = 0;
          }
        }
        if (model->shear_style == 1) {
          if (k == 0) {
            mesh->BCu.type[c] = -2;
            mesh->BCu.val[c] = 0.0 * model->EpsBG * Lx;
          }

          if (k == mesh->Nx - 1) {
            mesh->BCu.type[c] = -12;
            mesh->BCu.val[c] = -0.0 * model->EpsBG * Lx;
          }

          if (l == 0) {
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c] = -1 * model->EpsBG * Lz;
          }

          if (l == mesh->Nz) {
            mesh->BCu.type[c] = 11;
            mesh->BCu.val[c] = 1 * model->EpsBG * Lz;
          }
        }
      }
    }
  }

  NX = mesh->Nx;
  NZ = mesh->Nz;

  for (l = 0; l < mesh->Nz; l++) {
    for (k = 0; k < mesh->Nx + 1; k++) {

      c = k + l * (mesh->Nx + 1);

      if (mesh->BCv.type[c] != 30) {

        // Internal points:  -1
        mesh->BCv.type[c] = -1;
        mesh->BCv.val[c] = 0;

        if (model->shear_style == 0) {

          if (l == 0) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
          }

          if (l == mesh->Nz - 1) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = mesh->zg_coord[l] * model->EpsBG;
          }

          if ((k == 0)) {
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c] = 0;
          }

          if ((k == mesh->Nx)) {
            mesh->BCv.type[c] = 13;
            mesh->BCv.val[c] = 0;
          }
        }

        if (model->shear_style == 1) {
          if (l == 0) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = -0.0 * model->EpsBG * Lz;
          }

          if (l == mesh->Nz - 1) {
            mesh->BCv.type[c] = 0;
            mesh->BCv.val[c] = 0.0 * model->EpsBG * Lz;
          }

          if ((k == 0)) {
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c] = 0;
          }

          if ((k == mesh->Nx)) {
            mesh->BCv.type[c] = -12;
            mesh->BCv.val[c] = 0;
          }
        }
      }
    }
  }

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
      }
    }
  }

  double Ttop = 273.15 / scaling.T;

  NX = mesh->Nx;
  NZ = mesh->Nz;
  NCX = NX - 1;
  NCZ = NZ - 1;

  for (l = 0; l < mesh->Nz - 1; l++) {
    for (k = 0; k < mesh->Nx - 1; k++) {

      c = k + l * (NCX);

      if (mesh->BCt.type[c] != 30) {
        if (k == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }
        if (k == NCX - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }
        if (l == 0) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        }
        if (l == NCZ - 1) {
          mesh->BCt.type[c] = 0;
          mesh->BCt.val[c] = mesh->T[c];
        } else {
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

  free(X);
  free(Z);
  free(XC);
  free(ZC);
}

int main() {
  RunMDOODZ("ShearTemplate.txt", BuildInitialTopography, SetParticles, SetBCs);
  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t IterationsGroup = H5Gopen(File, "Iterations", H5P_DEFAULT);
  hid_t NumberStepsDataset =
      H5Dopen(IterationsGroup, "NumberSteps", H5P_DEFAULT);
  int NumberStepsArray[1] = {0};
  H5Dread(NumberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          NumberStepsArray);
  int StepsCount = NumberStepsArray[0];

  assert(StepsCount == 1);
}