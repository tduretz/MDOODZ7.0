#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double f = input->model.user1;

  double Lx = 1.0;

  double radius = sqrt((Lx*Lx * f)/M_PI); // f = a/A with a = pi*r^2 and A = L^2
  // p.ex. for f = 0.3 and L = 1.0: radius = 0.30901936161855165

  // one single circle
  int nellipses = 1;
  double Xshift[1]  = { 0.0};
  double Zshift[1]  = { 0.0};
  double Majaxis[1] = { radius}; // semi-major axis (!)
  double Minaxis[1] = { radius};
  double Angle[1]   = { 0.7853981633974483 }; // pi/4 or 45°

  // scaling geometry values
  int i;
  for (i = 0; i < nellipses; ++i) {
    Xshift[i]  = Xshift[i]  / input->scaling.L;
    Zshift[i]  = Zshift[i]  / input->scaling.L;
    Majaxis[i] = Majaxis[i] / input->scaling.L;
    Minaxis[i] = Minaxis[i] / input->scaling.L;
    // Angle doesn't need scaling (in radians btw)
  }

  bool isInsideSomeEllipse = false;
  bool isInsideThisEllipse = false;
  for (i = 0; i < nellipses; ++i) {
    isInsideThisEllipse = ((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))*((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))/(Majaxis[i] * Majaxis[i]) + (-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))*(-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))/ (Minaxis[i] * Minaxis[i]) < 1; // rotated ellipse
    // reminder ellipse function: inside = [(x-x0)*cos(phi)+(z-z0)*sin(phi)]^2/majaxis^2 + [(x-x0)*sin(phi)-(z-z0)*cos(phi)]^2/minaxis^2 < 1
    if (isInsideThisEllipse) {
      isInsideSomeEllipse = true;
    }
  }
  if (isInsideSomeEllipse) {
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

  return 135.0;
  //const double f = input->model.user1;
  //double Lx = 1.0;
  //double radius = sqrt((Lx*Lx * f)/M_PI); // f = a/A with a = pi*r^2 and A = L^2

  //// one single circle
  //int nellipses = 1;
  //double Xshift[1]  = { 0.0};
  //double Zshift[1]  = { 0.0};
  //double Majaxis[1] = { radius}; // semi-major axis (!)
  //double Minaxis[1] = { radius};
  //double Angle[1]   = { 0.7853981633974483 }; // pi/4 or 45°

  //bool isInsideSomeEllipse = false;
  //bool isInsideThisEllipse = false;
  //for (int i = 0; i < nellipses; ++i) {
  //  isInsideThisEllipse = ((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))*((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))/(Majaxis[i] * Majaxis[i]) + (-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))*(-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))/ (Minaxis[i] * Minaxis[i]) < 1; // rotated ellipse
  //  // reminder ellipse function: inside = [(x-x0)*cos(phi)+(z-z0)*sin(phi)]^2/majaxis^2 + [(x-x0)*sin(phi)-(z-z0)*cos(phi)]^2/minaxis^2 < 1
  //  if (isInsideThisEllipse) {
  //    isInsideSomeEllipse = true;
  //    return Angle[i] * 180.0 / M_PI + 90; // corresponding aniso_angle
  //  }
  //}
  //if (!isInsideSomeEllipse) {
  //  //return 135; // fixed value in matrix
  //  return rand()*360; // random value in matrix
  //}
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
  RunMDOODZ("William_Dabrowski_1incl.txt", &setup);
}
