#include "assert.h"
#include "hdf5.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"


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
  input->model.shear_style = mutateInputParams->int1;
  input->model.free_surf   = 0;
  const int matrixPhase    = 0;
  if (mutateInputParams->int2) {
    input->materials.aniso_factor[matrixPhase] = 2.0;
    input->model.aniso                         = 1;
  } else {
    input->materials.aniso_factor[matrixPhase] = 1.0;
    input->model.aniso                         = 0;
  }
  if (mutateInputParams->int4) {
    input->materials.npwl[matrixPhase] = 3.0;
    input->materials.cstv[matrixPhase] = 0;
    input->materials.pwlv[matrixPhase] = 1;
  } else {
    input->materials.npwl[matrixPhase] = 1.0;
    input->materials.cstv[matrixPhase] = 1;
    input->materials.pwlv[matrixPhase] = 0;
  }
}

int getStepsCount(char *hdf5FileName) {
  hid_t file            = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t iterationsGroup = H5Gopen(file, "Iterations", H5P_DEFAULT);
  hid_t numberStepsDataset =
          H5Dopen(iterationsGroup, "NumberSteps", H5P_DEFAULT);
  int numberStepsArray[1] = {0};
  H5Dread(numberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          numberStepsArray);
  int stepsCount = numberStepsArray[0];
  printf("Number of iterations: %d\n", stepsCount);
  return stepsCount;
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
  setup.mutateInputParams    = mutateInputParams;

  printf("lin_pureshear_iso");
  mutateInputParams->int1    = 0;  // shear_style
  mutateInputParams->int2    = 0;  // matrix aniso
  mutateInputParams->int4    = 0;  // non-linear
  mutateInputParams->double1 = 1.0;// matrix nwpl
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_pureshear_iso.h5");
  int stepsCount = getStepsCount("lin_pureshear_iso.h5");
  assert(stepsCount == 1);

  printf("lin_simpleshear_iso");
  mutateInputParams->int1 = 1;// shear_style
  mutateInputParams->int4 = 0;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_simpleshear_iso.h5");
  stepsCount = getStepsCount("lin_simpleshear_iso.h5");
  assert(stepsCount == 1);

  printf("lin_pureshear_aniso");
  mutateInputParams->int1 = 0;// shear_style
  mutateInputParams->int2 = 1;// matrix aniso
  mutateInputParams->int4 = 0;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_pureshear_aniso.h5");
  stepsCount = getStepsCount("lin_pureshear_aniso.h5");
  printf("stepsCount: %d", stepsCount);

  printf("lin_simpleshear_aniso");
  mutateInputParams->int1 = 1;// shear_style
  mutateInputParams->int4 = 0;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "lin_simpleshear_aniso.h5");
  stepsCount = getStepsCount("lin_simpleshear_aniso.h5");
  assert(stepsCount == 1);

  printf("nonlin_pureshear_iso");
  mutateInputParams->int1 = 0;// shear_style
  mutateInputParams->int2 = 0;// matrix aniso
  mutateInputParams->int4 = 1;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_pureshear_iso.h5");
  stepsCount = getStepsCount("nonlin_pureshear_iso.h5");
  assert(stepsCount > 1);

  printf("nonlin_simpleshear_iso");
  mutateInputParams->int1 = 1;// shear_style
  mutateInputParams->int2 = 0;// matrix aniso
  mutateInputParams->int4 = 1;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_simpleshear_iso.h5");
  stepsCount = getStepsCount("nonlin_simpleshear_iso.h5");
  assert(stepsCount > 1);

  printf("nonlin_pureshear_aniso");
  mutateInputParams->int1 = 0;// shear_style
  mutateInputParams->int2 = 1;// matrix aniso
  mutateInputParams->int4 = 1;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_pureshear_aniso.h5");
  stepsCount = getStepsCount("nonlin_pureshear_aniso.h5");
  assert(stepsCount > 1);

  printf("nonlin_simpleshear_aniso");
  mutateInputParams->int1 = 1;// shear_style
  mutateInputParams->int4 = 1;// non-linear
  RunMDOODZ("ShearTemplate.txt", &setup);
  rename("Output00001.gzip.h5", "nonlin_simpleshear_aniso.h5");
  stepsCount = getStepsCount("nonlin_simpleshear_aniso.h5");
  assert(stepsCount > 1);
}
