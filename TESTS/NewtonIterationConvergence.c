#include "assert.h"
#include "hdf5.h"
#include "mdoodz.h"
#include "stdlib.h"
#define FILENAME "Output00001.gzip.h5"


int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

void MutateInput(MdoodzInput *input, MutateInputParams *mutateInputParams) {
  input->model.shear_style           = mutateInputParams->int1;
  input->model.aniso                 = mutateInputParams->int2;
  input->model.free_surf             = mutateInputParams->int3;
  input->model.Newton                = 1;
  const int matrixPhase              = 0;
  input->materials.npwl[matrixPhase] = mutateInputParams->double1;
  if (input->model.aniso) {
    input->materials.aniso_factor[matrixPhase] = 2.0;
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
          .MutateInput = MutateInput,
  };

  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  mutateInputParams->int1              = 0;  // shear_style
  mutateInputParams->int2              = 0;  // matrix aniso
  mutateInputParams->double1           = 1.0;// matrix nwpl
  setup.mutateInputParams              = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_pureshear_iso.h5");

  mutateInputParams->int1 = 1;// shear_style
  setup.mutateInputParams = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_simpleshear_iso.h5");

  mutateInputParams->int1 = 0;// shear_style
  mutateInputParams->int2 = 1;// matrix aniso
  setup.mutateInputParams = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_pureshear_aniso.h5");

  mutateInputParams->int1 = 1;// shear_style
  setup.mutateInputParams = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_simpleshear_aniso.h5");

  mutateInputParams->int1    = 0;  // shear_style
  mutateInputParams->int2    = 0;  // matrix aniso
  mutateInputParams->double1 = 3.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_pureshear_iso.h5");

  mutateInputParams->int1 = 1;// shear_style
  mutateInputParams->int2 = 0;// matrix aniso
  setup.mutateInputParams = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_simpleshear_iso.h5");

  mutateInputParams->int1    = 0;  // shear_style
  mutateInputParams->int2    = 1;  // matrix aniso
  mutateInputParams->double1 = 3.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_pureshear_aniso.h5");

  mutateInputParams->int1    = 1;  // shear_style
  mutateInputParams->double1 = 3.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_simpleshear_aniso.h5");

  mutateInputParams->int1      = 0;  // shear_style
  mutateInputParams->int2      = 0;  // matrix aniso
  mutateInputParams->int3      = 1;  // free_surface
  setup.BuildInitialTopography = (BuildInitialTopography_ff *) malloc(sizeof(BuildInitialTopography_ff));
  setup.SetParticles->SetPhase = NULL;
  setup.SetParticles->SetDensity = NULL;
  mutateInputParams->double1   = 1.0;// matrix nwpl
  setup.mutateInputParams      = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "freesurf_lin_pureshear_iso.h5");

  mutateInputParams->int2    = 1;  // matrix aniso
  mutateInputParams->double1 = 1.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "freesurf_lin_pureshear_iso.h5");

  mutateInputParams->int2    = 0;  // matrix aniso
  mutateInputParams->double1 = 3.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "freesurf_nonlin_pureshear_aniso.h5");

  mutateInputParams->int2    = 1;  // matrix aniso
  mutateInputParams->double1 = 3.0;// matrix nwpl
  setup.mutateInputParams    = mutateInputParams;
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "freesurf_nonlin_pureshear_aniso.h5");
}
