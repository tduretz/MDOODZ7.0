#include "mdoodz.h"
#include "math.h"

bool IsInCircle(Coordinates coordinates, double rad) {
  return pow(coordinates.x, 2) + pow(coordinates.z, 2) < pow(rad, 2);
}

bool IsInHorizontalZones(double zCoord, double H) {
  if (zCoord > 1.0 * H / 6.0 && zCoord < 2.0 * H / 6.0) {
    return true;
  } else if (zCoord < 0.0 * H / 6.0 && zCoord > -1.0 * H / 6.0) {
    return true;
  } else if (zCoord < -2.0 * H / 6.0 && zCoord > -3.0 * H / 6.0) {
    return true;
  } else {
    return false;
  }
}

bool IsInVerticalZones(double xCoord, double H) {
  if (xCoord > 1.0 * H / 6.0 && xCoord < 2.0 * H / 6.0) {
    return true;
  } else if (xCoord > -1.0 * H / 6.0 && xCoord < 0.0 * H / 6.0) {
    return true;
  } else if (xCoord > -3.0 * H / 6.0 && xCoord < -2.0 * H / 6.0) {
    return true;
  } else {
    return false;
  }
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  double H    = (input->model.zmax - input->model.zmin);
  double rad = 0.2337 / input->scaling.L;

  if (IsInCircle(coordinates, rad)) {
    if ((coordinates.x < 0 && coordinates.z < 0) || (coordinates.x > 0 && coordinates.z > 0)) {
      return 1;
    } else {
      return 0;
    }
  } else if (IsInVerticalZones(coordinates.x, H)) {
    if (IsInHorizontalZones(coordinates.z, H)) {
      return 2;
    } else {
      return 3;
    }
  } else {
    if (!IsInHorizontalZones(coordinates.z, H)) {
      return 2;
    } else {
      return 3;
    }
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double SetVerticalVelocity(MdoodzInput *input, Coordinates coordinates) {
  return 0.0;
}

double SetHorizontalVelocity(MdoodzInput *input, Coordinates coordinates) {
  return -0.5 * input->model.EpsBG * (coordinates.z + (input->model.zmax - input->model.zmin) / 2);
}

double SetGrainSize(MdoodzInput *input, Coordinates coordinates, int phase) {
  return 2e-3 / input->scaling.L;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetDensity            = SetDensity,
                  .SetVerticalVelocity   = SetVerticalVelocity,
                  .SetHorizontalVelocity = SetHorizontalVelocity,
                  .SetGrainSize          = SetGrainSize,

          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetSimpleShearBCVx,
                  .SetBCVz = SetSimpleShearBCVz,
          },
  };
  RunMDOODZ("SimpleShearAnisoHomo.txt", &setup);
}
