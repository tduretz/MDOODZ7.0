#include "math.h"
#include "mdoodz.h"

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
  const Ellipse circle = (Ellipse){
          .angle   = 0,
          .radiusX = input->model.user1,
          .radiusZ = input->model.user1,
          .centreX = 0.0,
          .centreZ = 0.0,
  };
  if (IsEllipseCoordinates(coordinates, circle, input->scaling.L)) {
    return 1;
  } else {
    return 0;
  }
}

int SetDual(MdoodzInput *input, Coordinates coordinates, int phase) {
  if (phase != 0) {
    return phase;
  }
  const double Lx = (double) (input->model.xmax - input->model.xmin);
  const double Lz = (double) (input->model.zmax - input->model.zmin);
  const double Ax = cos(6.0 * 2.0 * M_PI * coordinates.x / Lx);
  const double Az = sin(6.0 * 2.0 * M_PI * coordinates.z / Lz);
  if ((Az < 0.0 && Ax < 0.0) || (Az > 0.0 && Ax > 0.0)) {
    return input->model.Nb_phases;
  } else {
    return phase;
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase = SetPhase,
                  .SetDual  = SetDual,

          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("ShearInclusionAniso.txt", &setup);
}
