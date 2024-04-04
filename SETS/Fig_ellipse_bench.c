#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {

  // Fig. 06a) Rectangle pure shear  
  double Xshift = 0.0;
  double Zshift = 0.0;
  double Majaxis = 0.35;
  double Minaxis = 0.15;
  double Angle = 45.0 / 180.0 * M_PI;

  // scaling geometry values
  Xshift  = Xshift  / input->scaling.L;
  Zshift  = Zshift  / input->scaling.L;
  Majaxis = Majaxis / input->scaling.L;
  Minaxis = Minaxis / input->scaling.L;
  // Angle doesn't need scaling (in radians btw)

  bool isInsideEllipse = false;
  isInsideEllipse = ((coordinates.x - Xshift)*cos(Angle) + (coordinates.z - Zshift)*sin(Angle))*((coordinates.x - Xshift)*cos(Angle) + (coordinates.z - Zshift)*sin(Angle))/(Majaxis * Majaxis) + (-(coordinates.z - Zshift)*cos(Angle) + (coordinates.x - Xshift)*sin(Angle))*(-(coordinates.z - Zshift)*cos(Angle) + (coordinates.x - Xshift)*sin(Angle))/ (Minaxis * Minaxis) < 1; // rotated ellipse
  // reminder ellipse function: inside = [(x-x0)*cos(phi)+(z-z0)*sin(phi)]^2/majaxis^2 + [(x-x0)*sin(phi)-(z-z0)*cos(phi)]^2/minaxis^2 < 1

  if (isInsideEllipse) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double SetAnisoAngle(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double radius = input->model.user1 / input->scaling.L;

  // Fig. 06a) Rectangle pure shear  
  double Xshift = 0.0;
  double Zshift = 0.0;
  double Majaxis = 0.35;
  double Minaxis = 0.15;
  double Angle = 45.0 / 180.0 * M_PI;

  // scaling geometry values
  Xshift  = Xshift  / input->scaling.L;
  Zshift  = Zshift  / input->scaling.L;
  Majaxis = Majaxis / input->scaling.L;
  Minaxis = Minaxis / input->scaling.L;
  // Angle doesn't need scaling (in radians btw)

  bool isInsideEllipse = false;
  isInsideEllipse = ((coordinates.x - Xshift)*cos(Angle) + (coordinates.z - Zshift)*sin(Angle))*((coordinates.x - Xshift)*cos(Angle) + (coordinates.z - Zshift)*sin(Angle))/(Majaxis * Majaxis) + (-(coordinates.z - Zshift)*cos(Angle) + (coordinates.x - Xshift)*sin(Angle))*(-(coordinates.z - Zshift)*cos(Angle) + (coordinates.x - Xshift)*sin(Angle))/ (Minaxis * Minaxis) < 1; // rotated ellipse
  // reminder ellipse function: inside = [(x-x0)*cos(phi)+(z-z0)*sin(phi)]^2/majaxis^2 + [(x-x0)*sin(phi)-(z-z0)*cos(phi)]^2/minaxis^2 < 1

  if (isInsideEllipse) {
    return Angle * 180.0 / M_PI + 90; // corresponding aniso_angle
  } else {
    return rand()*360; // random value in matrix
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetAnisoAngle         = SetAnisoAngle,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("Fig_ellipse_bench.txt", &setup);
}
