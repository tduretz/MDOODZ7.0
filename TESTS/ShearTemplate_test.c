#include "assert.h"
#include "hdf5.h"
#include "mdoodz.h"
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

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == N || position == NE || position == NW || position == SE || position == SW) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = -coordinates.x * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == NE || position == NW || position == SE || position == SW) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = coordinates.z * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

int main() {
  MdoodzInput instance = {
          .inputFileName = "ShearTemplate.txt",
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ(&instance);

  hid_t File            = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t IterationsGroup = H5Gopen(File, "Iterations", H5P_DEFAULT);
  hid_t NumberStepsDataset =
          H5Dopen(IterationsGroup, "NumberSteps", H5P_DEFAULT);
  int NumberStepsArray[1] = {0};
  H5Dread(NumberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          NumberStepsArray);
  int StepsCount = NumberStepsArray[0];

  assert(StepsCount == 1);
}
