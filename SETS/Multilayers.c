#include "math.h"
#include "mdoodz.h"
#include "stdlib.h"

static bool isPerturbationInitialised = false;
static double **pert_lo, **pert_up;

void        InitialisePerturbation(MdoodzInput *input) {
  int      nlayers  = (int) input->model.user2;
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

  double A = 1 / input->scaling.L;

  for (int il = 0; il < nlayers; il++) {
    for (int ic = 30; ic < ncx - 30; ic++) {
#ifdef __APPLE__
      pert_lo[il][ic] = 1 * A * (double) arc4random() / UINT32_MAX - A / 2;
      pert_up[il][ic] = 1 * A * (double) arc4random() / UINT32_MAX - A / 2;
#else
      pert_lo[il][ic] = 1 * A * (double) rand() / RAND_MAX - A / 2;
      pert_up[il][ic] = 1 * A * (double) rand() / RAND_MAX - A / 2;
#endif
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
  if (!isPerturbationInitialised) {
    InitialisePerturbation(input);
  }

  const double H            = input->model.user0 / input->scaling.L;
  const double spacing      = input->model.user1 / input->scaling.L;
  const int    nlayers      = (int) input->model.user2;
  const int    ncx          = input->model.Nx - 1;
  const double layerBoxSize = (nlayers * (H + spacing)) - spacing;
  const double totalBoxSize = -input->model.zmin + input->model.zmax;

  const double offset = (totalBoxSize - layerBoxSize) / 2;

  for (int il = 0; il < nlayers; il++) {
    const double layer_bottom = input->model.zmin + offset + il * (H + spacing);
    int    ic           = ceil(((coordinates.x - (input->model.xmin + input->model.dx / 2)) / input->model.dx) + 0.5) - 1;
    if (ic < 0) ic = 0;
    if (ic > ncx - 1) ic = ncx - 1;
    const double pert_lo_c = pert_lo[il][ic];
    const double pert_up_c = pert_up[il][ic];
    if (coordinates.z > layer_bottom + pert_lo_c && coordinates.z < layer_bottom + H + pert_up_c) {
      return 1;
    }
  }
  return 0;
}


int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase   = SetPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("Multilayers.txt", &setup);
}
