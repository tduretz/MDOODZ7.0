#include "mdoodz.h"
#include "math.h"

// bool IsInHorizontalZones(double zCoord, double H) {
//   if (zCoord > 1.0 * H / 6.0 && zCoord < 2.0 * H / 6.0) {
//     return true;
//   } else if (zCoord < 0.0 * H / 6.0 && zCoord > -1.0 * H / 6.0) {
//     return true;
//   } else if (zCoord < -2.0 * H / 6.0 && zCoord > -3.0 * H / 6.0) {
//     return true;
//   } else {
//     return false;
//   }
// }

// bool IsInVerticalZones(double xCoord, double H) {
//   if (xCoord > 1.0 * H / 6.0 && xCoord < 2.0 * H / 6.0) {
//     return true;
//   } else if (xCoord > -1.0 * H / 6.0 && xCoord < 0.0 * H / 6.0) {
//     return true;
//   } else if (xCoord > -3.0 * H / 6.0 && xCoord < -2.0 * H / 6.0) {
//     return true;
//   } else {
//     return false;
//   }
// }

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double  H      = (input->model.zmax - input->model.zmin);
  const double  rad    = 0.2337 * 2;
  const Ellipse circle = (Ellipse){
          .angle   = 0,
          .radiusX = rad,
          .radiusZ = rad,
          .centreX = 0.0,
          .centreZ = 0.0,
  };
  if (IsEllipseCoordinates(coordinates, circle, input->scaling.L)) {
     return 1;
    } else {
      return 0;
    }
}

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // This checker board looks weird at low resolution (?)

    // Set checkerboard for phase 0
    Ax = cos( 2*3.0*M_PI*coordinate.x / Lx  );
    Az = sin( 2*3.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Ax = cos( 2*3.0*M_PI*coordinate.x / Lx  );
    Az = sin( 2*3.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }

  return dual_phase;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double SetVerticalVelocity(MdoodzInput *input, Coordinates coordinates) {
  return 0.0;
}

double SetHorizontalVelocity(MdoodzInput *input, Coordinates coordinates) {
  return -0.5 * input->model.bkg_strain_rate * (coordinates.z + (input->model.zmax - input->model.zmin) / 2);
}

double SetGrainSize(MdoodzInput *input, Coordinates coordinates, int phase) {
  return 2e-3 / input->scaling.L;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetDualPhase          = SetDualPhase,
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
