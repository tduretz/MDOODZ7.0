#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {

  // B6) large circles
  int nellipses = 19;
  double Xshift[19] = { -0.039839 , -0.039839 , -0.21756 , -0.40536 , -0.3594 , 0.082289 , 0.055881 , 0.39664 , -0.27811 , 0.35478 , 0.35478 , 0.49922 , -0.50078 , 0.23346 , -0.049993 , -0.33661 , -0.37744 , 0.085244 , 0.24523 };
  double Zshift[19] = { 0.47959 , -0.52041 , -0.25289 , 0.24778 , -0.17126 , -0.21801 , 0.11881 , -0.21808 , 0.3959 , 0.46174 , -0.53826 , -0.057673 , -0.057673 , -0.31448 , 0.27161 , 0.096272 , -0.38035 , -0.39411 , 0.10369 };
  double Majaxis[19] = { 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 };
  double Minaxis[19] = { 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 };
  double Angle[19] = { 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 };

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
  const double radius = input->model.user1 / input->scaling.L;

  // B6) large circles
  int nellipses = 19;
  double Xshift[19] = { -0.039839 , -0.039839 , -0.21756 , -0.40536 , -0.3594 , 0.082289 , 0.055881 , 0.39664 , -0.27811 , 0.35478 , 0.35478 , 0.49922 , -0.50078 , 0.23346 , -0.049993 , -0.33661 , -0.37744 , 0.085244 , 0.24523 };
  double Zshift[19] = { 0.47959 , -0.52041 , -0.25289 , 0.24778 , -0.17126 , -0.21801 , 0.11881 , -0.21808 , 0.3959 , 0.46174 , -0.53826 , -0.057673 , -0.057673 , -0.31448 , 0.27161 , 0.096272 , -0.38035 , -0.39411 , 0.10369 };
  double Majaxis[19] = { 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 };
  double Minaxis[19] = { 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 , 0.075 };
  double Angle[19] = { 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 , 0.7854 };

  bool isInsideSomeEllipse = false;
  bool isInsideThisEllipse = false;
  for (int i = 0; i < nellipses; ++i) {
    isInsideThisEllipse = ((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))*((coordinates.x - Xshift[i])*cos(Angle[i]) + (coordinates.z - Zshift[i])*sin(Angle[i]))/(Majaxis[i] * Majaxis[i]) + (-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))*(-(coordinates.z - Zshift[i])*cos(Angle[i]) + (coordinates.x - Xshift[i])*sin(Angle[i]))/ (Minaxis[i] * Minaxis[i]) < 1; // rotated ellipse
    // reminder ellipse function: inside = [(x-x0)*cos(phi)+(z-z0)*sin(phi)]^2/majaxis^2 + [(x-x0)*sin(phi)-(z-z0)*cos(phi)]^2/minaxis^2 < 1
    if (isInsideThisEllipse) {
      isInsideSomeEllipse = true;
      return Angle[i] * 180.0 / M_PI + 90; // corresponding aniso_angle
    }
  }
  if (!isInsideSomeEllipse) {
    //return 135; // fixed value in matrix
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
  RunMDOODZ("WilliamB6.txt", &setup);
}
