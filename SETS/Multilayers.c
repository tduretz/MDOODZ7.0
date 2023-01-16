#include "math.h"
#include "mdoodz.h"
#include "stdlib.h"

static bool isPerturbationInitialised = false;
static double **pert_lo, **pert_up;

void        InitialisePerturbation(MdoodzInput *input) {
  int      nlayers  = 9;
  int      ncx      = input->model.Nx - 1;
  double **pert_lo0 = calloc(nlayers, sizeof(double *));
  double **pert_up0 = calloc(nlayers, sizeof(double *));
  pert_lo    = calloc(nlayers, sizeof(double *));
  pert_up    = calloc(nlayers, sizeof(double *));

  for (int il = 0; il < nlayers; il++) {
    pert_lo0[il]       = calloc(ncx, sizeof(double));
    pert_up0[il]       = calloc(ncx, sizeof(double));
    pert_lo[il] = calloc(ncx, sizeof(double));
    pert_up[il] = calloc(ncx, sizeof(double));
  }

  double A = 0.5 / input->scaling.L;

  for (int il = 0; il < nlayers; il++) {
    for (int ic = 30; ic < ncx - 30; ic++) {
      pert_lo[il][ic] = 1 * A * (double) arc4random() / UINT32_MAX - A / 2;
      pert_up[il][ic] = 1 * A * (double) arc4random() / UINT32_MAX - A / 2;
    }

    double D = 0.4;
    for (int it = 0; it < ncx; it++) {
      for (int ic = 0; ic < ncx - 0; ic++) {
        pert_lo0[il][ic] = pert_lo[il][ic];
        pert_up0[il][ic] = pert_up[il][ic];
      }
      for (int ic = 1; ic < ncx - 1; ic++) {
        pert_lo[il][ic] = pert_lo0[il][ic] + D * (pert_lo0[il][ic - 1] + pert_lo0[il][ic + 1] - 2 * pert_lo0[il][ic]);
        pert_up[il][ic] = pert_up0[il][ic] + D * (pert_up0[il][ic - 1] + pert_up0[il][ic + 1] - 2 * pert_up0[il][ic]);
      }
    }
  }
  isPerturbationInitialised = true;
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  int    nlayers = (int) input->model.user2;
  double H       = input->model.user0 / input->scaling.L;
  int    ncx     = input->model.Nx - 1;

  if (!isPerturbationInitialised) {
    InitialisePerturbation(input);
  }

  double spacing = -input->model.user1 / input->scaling.L;
  for (int il = 0; il < nlayers; il++) {

    double layer_bottom = input->model.zmin + 6 / input->scaling.L + il * (H + spacing);

    // Get the column:
    int    ic           = ceil(((coordinates.x - (input->model.xmin + input->model.dx / 2)) / input->model.dx) + 0.5) - 1;
    if (ic < 0) ic = 0;
    if (ic > ncx - 1) ic = ncx - 1;
    // Corresponding layer contour perturbation


    double pert_lo_c = pert_lo[il][ic];
    double pert_up_c = pert_up[il][ic];

    if (coordinates.z > layer_bottom + pert_lo_c && coordinates.z < layer_bottom + H + pert_up_c) {
      return il + 1;
    }
  }
  return 0;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase   = SetPhase,
                  .SetDensity = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("Multilayers.txt", &setup);
}
