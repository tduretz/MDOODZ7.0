#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {

  // B5) circles, same size, 5x5 all aligned
  double d_Xshift[25] = { 0.9485484492371699 , 0.9333762393916095  , 0.3076373730588343  , 0.8635277035926671   , 0.9485484492371699 , 0.5890979596406679 , 0.6829027842212059 , 0.9942515475871702 , 0.9621266398494146    , 0.5890979596406679 , 0.908601239653215  , 0.8622361781698946   , 0.9925757704891546  , 0.8503822186310541  , 0.908601239653215 , 0.0013589198992898233 , 0.9334670539028086 , 0.07588996956074079 , 0.8657927415787832  , 0.0013589198992898233     , 0.9485484492371699 , 0.9333762393916095   , 0.3076373730588343  , 0.8635277035926671 , 0.9485484492371699 };
  double d_Zshift[25] = { 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.4273171754694547 , 0.5425866165938917 , 0.6201318911085705 , 0.90273532152766   , 0.0019465420865421024 , 0.5425866165938917 , 0.6583670587253173 , 0.052142170900339346 , 0.11051969187091559 , 0.12033430560852221 , 0.6583670587253173 , 0.7117852782330175    , 0.8241369329108513 , 0.8610212856101548  , 0.35619242168110343 , 0.7117852782330175       , 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.4273171754694547 };
  //double d_Xshift[25] = { 0.9485484492371699 , 0.9333762393916095 , 0.3076373730588343 , 0.8635277035926671 , 0.15496802379385022 , 0.5890979596406679 , 0.6829027842212059 , 0.9942515475871702 , 0.9621266398494146 , 0.4259237710938296 , 0.908601239653215 , 0.8622361781698946 , 0.9925757704891546 , 0.8503822186310541 , 0.6854127168974997 , 0.0013589198992898233 , 0.9334670539028086 , 0.07588996956074079 , 0.8657927415787832 , 0.367442674337337 , 0.2509770419357381 , 0.7662147966485621 , 0.49155362059076 , 0.34965262596847924 , 0.22472835082437814   };
  //double d_Zshift[25] = { 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.24442697436391425 , 0.5425866165938917 , 0.6201318911085705 , 0.90273532152766 , 0.0019465420865421024 , 0.8767279816239639 , 0.6583670587253173 , 0.052142170900339346 , 0.11051969187091559 , 0.12033430560852221 , 0.8371639193166704 , 0.7117852782330175 , 0.8241369329108513 , 0.8610212856101548 , 0.35619242168110343 , 0.0026585669431200554 , 0.098480050950839 , 0.3672892330306361 , 0.35165850677517885 , 0.40825308767366664 , 0.7289186882099035 };
  
  int nellipses = 25;
  double Xshift[25];
  double Zshift[25];
  double Majaxis[25];
  double Minaxis[25];
  double Angle[25];
  int kg;
  for (int i = 0; i < 5; ++i) {     // z idx
    for (int j = 0; j < 5; ++j) {   // x idx
      kg = 5 * i + j;               // global idx
      Xshift[kg]  = ( (1.0 / 4.0) * j - 0.5 ) ;// + (d_Xshift[kg]-0.5)*0.05; // including random delta
      Zshift[kg]  = ( (1.0 / 4.0) * i - 0.5 ) ;// + (d_Zshift[kg]-0.5)*0.05; // including random delta
      Majaxis[kg] = 0.075;
      Minaxis[kg] = 0.075;
      Angle[kg]   = 45.0 / 180.0 * M_PI;
      if (j == 1 || j == 3) {
        Zshift[kg] += 1.0 / 8.0 ;
      }
    }
  }
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

  // B5) circles, same size, 5x5 all aligned
  double d_Xshift[25] = { 0.9485484492371699 , 0.9333762393916095  , 0.3076373730588343  , 0.8635277035926671   , 0.9485484492371699 , 0.5890979596406679 , 0.6829027842212059 , 0.9942515475871702 , 0.9621266398494146    , 0.5890979596406679 , 0.908601239653215  , 0.8622361781698946   , 0.9925757704891546  , 0.8503822186310541  , 0.908601239653215 , 0.0013589198992898233 , 0.9334670539028086 , 0.07588996956074079 , 0.8657927415787832  , 0.0013589198992898233     , 0.9485484492371699 , 0.9333762393916095   , 0.3076373730588343  , 0.8635277035926671 , 0.9485484492371699 };
  double d_Zshift[25] = { 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.4273171754694547 , 0.5425866165938917 , 0.6201318911085705 , 0.90273532152766   , 0.0019465420865421024 , 0.5425866165938917 , 0.6583670587253173 , 0.052142170900339346 , 0.11051969187091559 , 0.12033430560852221 , 0.6583670587253173 , 0.7117852782330175    , 0.8241369329108513 , 0.8610212856101548  , 0.35619242168110343 , 0.7117852782330175       , 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.4273171754694547 };
  //double d_Xshift[25] = { 0.9485484492371699 , 0.9333762393916095 , 0.3076373730588343 , 0.8635277035926671 , 0.15496802379385022 , 0.5890979596406679 , 0.6829027842212059 , 0.9942515475871702 , 0.9621266398494146 , 0.4259237710938296 , 0.908601239653215 , 0.8622361781698946 , 0.9925757704891546 , 0.8503822186310541 , 0.6854127168974997 , 0.0013589198992898233 , 0.9334670539028086 , 0.07588996956074079 , 0.8657927415787832 , 0.367442674337337 , 0.2509770419357381 , 0.7662147966485621 , 0.49155362059076 , 0.34965262596847924 , 0.22472835082437814   };
  //double d_Zshift[25] = { 0.4273171754694547 , 0.10767404451273999 , 0.46873148134167775 , 0.012067110769540812 , 0.24442697436391425 , 0.5425866165938917 , 0.6201318911085705 , 0.90273532152766 , 0.0019465420865421024 , 0.8767279816239639 , 0.6583670587253173 , 0.052142170900339346 , 0.11051969187091559 , 0.12033430560852221 , 0.8371639193166704 , 0.7117852782330175 , 0.8241369329108513 , 0.8610212856101548 , 0.35619242168110343 , 0.0026585669431200554 , 0.098480050950839 , 0.3672892330306361 , 0.35165850677517885 , 0.40825308767366664 , 0.7289186882099035 };
  
  int nellipses = 25;
  double Xshift[25];
  double Zshift[25];
  double Majaxis[25];
  double Minaxis[25];
  double Angle[25];
  int kg;
  for (int i = 0; i < 5; ++i) {     // z idx
    for (int j = 0; j < 5; ++j) {   // x idx
      kg = 5 * i + j;               // global idx
      Xshift[kg]  = ( (1.0 / 4.0) * j - 0.5 ) ;// + (d_Xshift[kg]-0.5)*0.05; // including random delta
      Zshift[kg]  = ( (1.0 / 4.0) * i - 0.5 ) ;// + (d_Zshift[kg]-0.5)*0.05; // including random delta
      Majaxis[kg] = 0.075;
      Minaxis[kg] = 0.075;
      Angle[kg]   = 45.0 / 180.0 * M_PI;
      if (j == 1 || j == 3) {
        Zshift[kg] += 1.0 / 8.0 ;
      }
    }
  }

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
  RunMDOODZ("AnisoViscTest_evolv_multi_ellipses.txt", &setup);
}
