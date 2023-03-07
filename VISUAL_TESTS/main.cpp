extern "C" {
#include "mdoodz.h"
}
#include "visual-tests.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace std;
namespace fs = filesystem;

class RiftingChenin {
  static double RPSetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    const double TopoLevel = -0.0e3 / instance->scaling.L;
    const double h_pert    = instance->model.user3 / instance->scaling.L;
    return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
  }

  static int RPSetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double lithosphereThickness  = instance->model.user1 / instance->scaling.L;
    const double crustThickness        = instance->model.user2 / instance->scaling.L;
    const double perturbationAmplitude = instance->model.user3 / instance->scaling.L;
    const double mohoLevel             = -crustThickness - perturbationAmplitude * cos(2 * M_PI * coordinates.x / (instance->model.xmax - instance->model.xmin));
    const bool   isBelowLithosphere    = coordinates.z < -lithosphereThickness;
    const bool   isAboveMoho           = coordinates.z > mohoLevel;

    if (instance->model.user4 && isAboveMoho) {
      const bool is2500MAboveMoho = coordinates.z > mohoLevel + 2500 / instance->scaling.L;
      const bool is4500MAboveMoho = coordinates.z > mohoLevel + 4500 / instance->scaling.L;
      const bool is7000MAboveMoho = coordinates.z > mohoLevel + 7000 / instance->scaling.L;
      const bool is9000MAboveMoho = coordinates.z > mohoLevel + 9000 / instance->scaling.L;
      if (is2500MAboveMoho && !is4500MAboveMoho) {
        return 7;
      } else if (is7000MAboveMoho && !is9000MAboveMoho) {
        return 8;
      } else {
        return 1;
      }
    } else if (isAboveMoho) {
      return 1;
    } else if (isBelowLithosphere) {
      return 3;
    } else {
      return 2;
    }
  }

  static double RPSetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
    const double surfaceTemperature   = 273.15 / instance->scaling.T;
    const double mantleTemperature    = (1330.0 + 273.15) / instance->scaling.T;
    const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
    if (particleTemperature > mantleTemperature) {
      return mantleTemperature;
    } else {
      return particleTemperature;
    }
  }

  static double RPSetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const int astenospherePhase = 3;
    return instance->materials.gs_ref[astenospherePhase];
  }

  static double RPSetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
    return -coordinates.x * instance->model.bkg_strain_rate;
  }

  static double RPSetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
    return coordinates.z * instance->model.bkg_strain_rate;
  }

  static char RPSetBCPType(MdoodzInput *instance, POSITION position) {
    if (position == NE || position == NW) {
      return 0;
    } else {
      return -1;
    }
  }

  static SetBC RPSetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
    SetBC  bc;
    double surfaceTemperature = zeroC / instance->scaling.T;
    if (position == free_surface) {
      bc.value = surfaceTemperature;
      bc.type  = 1;
    } else {
      bc.value = 0.0;
      bc.type  = 0;
    }
    return bc;
  }

  static SetBC RPSetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
    SetBC  bc;
    double surfaceTemperature = zeroC / instance->scaling.T;
    double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
    if (position == S || position == SE || position == SW) {
      bc.value = particleTemperature;
      bc.type  = 1;
    } else if (position == N || position == NE || position == NW) {
      bc.value = surfaceTemperature;
      bc.type  = 1;
    } else if (position == W || position == E) {
      bc.value = mantleTemperature;
      bc.type  = 0;
    } else {
      bc.value = 0;
      bc.type  = 0;
    }
    return bc;
  }

  static void RPMutateInput(MdoodzInput *instance, MutateInputParams *params) {
    int *astenospherePhases     = (int *) malloc(sizeof(int));
    astenospherePhases[0]       = {3};
    instance->crazyConductivity = new CrazyConductivity{
            .multiplier = 1000,
            .phases     = astenospherePhases,
            .nPhases    = 1,
    };
  }

  public: void run() {
    MdoodzSetup riftingCheninModel = {
            .BuildInitialTopography = new BuildInitialTopography_ff{
                    .SetSurfaceZCoord = RPSetSurfaceZCoord,
            },
            .SetParticles = new SetParticles_ff{
                    .SetHorizontalVelocity = RPSetHorizontalVelocity,
                    .SetVerticalVelocity   = RPSetVerticalVelocity,
                    .SetPhase              = RPSetPhase,
                    .SetTemperature        = RPSetTemperature,
                    .SetGrainSize          = RPSetGrainSize,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx    = SetPureShearBCVx,
                    .SetBCVz    = SetPureShearBCVz,
                    .SetBCT     = RPSetBCT,
                    .SetBCPType = RPSetBCPType,
                    .SetBCTNew  = RPSetBCTNew,
            },
            .MutateInput = RPMutateInput};
    RunMDOODZ("RiftingChenin.txt", &riftingCheninModel);
  }
};

class TopoBenchCase1 {
  static double TBCSetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    double Amplitude  = 7e3 / instance->scaling.L;
    double Wavelength = 2800e3 / instance->scaling.L;
    return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
  }

  static int TBCSetPhase(MdoodzInput *instance, Coordinates coordinates) {
    double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
    if (coordinates.z > lithosphereBottomDepth) {
      return 1;
    } else {
      return 0;
    }
  }

  static SetBC TBCSetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == S || position == SW || position == SE) {
      bc.value = 0.0;
      bc.type  = 11;
    } else if (position == N || position == NW || position == NE) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == W || position == E) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static SetBC TBCSetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static char TBCSetBCPType(MdoodzInput *instance, POSITION position) {
    if (position == NE || position == NW) {
      return 0;
    } else {
      return -1;
    }
  }

  public:  void run() {
    MdoodzSetup topoBenchCaseModel = {
            .BuildInitialTopography = new BuildInitialTopography_ff{
                    .SetSurfaceZCoord = TBCSetSurfaceZCoord,
            },
            .SetParticles = new SetParticles_ff{
                    .SetPhase = TBCSetPhase,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx    = TBCSetBCVx,
                    .SetBCVz    = TBCSetBCVz,
                    .SetBCPType = TBCSetBCPType,
            },
    };
    RunMDOODZ("TopoBenchCase1.txt", &topoBenchCaseModel);
  }
};

class ShearTemplate {
  static int STSetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double STSetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (1 == 0 == 0) {
      return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
    } else {
      return instance->materials.rho[phase];
    }
  }

  public: void run() {
    MdoodzSetup shearTemplate = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = STSetPhase,
                    .SetDensity = STSetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureOrSimpleShearBCVx,
                    .SetBCVz = SetPureOrSimpleShearBCVz,
            },
    };
    RunMDOODZ("ShearTemplate.txt", &shearTemplate);
  }

  public: void run1() {
    MdoodzSetup shearTemplate = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = STSetPhase,
                    .SetDensity = STSetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureOrSimpleShearBCVx,
                    .SetBCVz = SetPureOrSimpleShearBCVz,
            },
    };
    RunMDOODZ("ShearTemplate1.txt", &shearTemplate);
  }

  public: void runAniso() {
    MdoodzSetup shearTemplate = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = STSetPhase,
                    .SetDensity = STSetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureOrSimpleShearBCVx,
                    .SetBCVz = SetPureOrSimpleShearBCVz,
            },
    };
    RunMDOODZ("ShearTemplateAniso.txt", &shearTemplate);
  }
};

void RenameTopoBenchCaseFiles() {
  rename("Output00000.gzip.h5", "TopoBenchCase00000.gzip.h5");
  rename("Output00010.gzip.h5", "TopoBenchCase00010.gzip.h5");
  rename("Output00020.gzip.h5", "TopoBenchCase00020.gzip.h5");
  rename("Output00030.gzip.h5", "TopoBenchCase00030.gzip.h5");
  rename("Output00040.gzip.h5", "TopoBenchCase00040.gzip.h5");
  rename("Output00050.gzip.h5", "TopoBenchCase00050.gzip.h5");
  rename("Output00060.gzip.h5", "TopoBenchCase00060.gzip.h5");
  rename("Output00070.gzip.h5", "TopoBenchCase00070.gzip.h5");
  rename("Output00080.gzip.h5", "TopoBenchCase00080.gzip.h5");
  rename("Output00090.gzip.h5", "TopoBenchCase00090.gzip.h5");
  rename("Output00100.gzip.h5", "TopoBenchCase00100.gzip.h5");
}

void RenameGSEFiles() {
  fs::create_directory("GSE");
  const auto copyOptions = fs::copy_options::overwrite_existing;
  fs::copy("Output00010.gzip.h5", "GSE/GSE00010.gzip.h5", copyOptions);
  fs::copy("Output00020.gzip.h5", "GSE/GSE00020.gzip.h5", copyOptions);
  fs::copy("Output00030.gzip.h5", "GSE/GSE00030.gzip.h5", copyOptions);
  fs::copy("Output00040.gzip.h5", "GSE/GSE00040.gzip.h5", copyOptions);
  fs::copy("Output00050.gzip.h5", "GSE/GSE00050.gzip.h5", copyOptions);
  fs::copy("Output00060.gzip.h5", "GSE/GSE00060.gzip.h5", copyOptions);
  fs::copy("Output00070.gzip.h5", "GSE/GSE00070.gzip.h5", copyOptions);
  fs::copy("Output00080.gzip.h5", "GSE/GSE00080.gzip.h5", copyOptions);
  fs::copy("Output00090.gzip.h5", "GSE/GSE00090.gzip.h5", copyOptions);
  fs::copy("Output00100.gzip.h5", "GSE/GSE00100.gzip.h5", copyOptions);
}

void RenameVEPFiles() {
  fs::create_directory("VEP");
  const auto copyOptions = fs::copy_options::overwrite_existing;
  fs::copy("Output00001.gzip.h5", "VEP/VEP00001.gzip.h5", copyOptions);
  fs::copy("Output00005.gzip.h5", "VEP/VEP00005.gzip.h5", copyOptions);
  fs::copy("Output00010.gzip.h5", "VEP/VEP00010.gzip.h5", copyOptions);
  fs::copy("Output00015.gzip.h5", "VEP/VEP00015.gzip.h5", copyOptions);
  fs::copy("Output00020.gzip.h5", "VEP/VEP00020.gzip.h5", copyOptions);
  fs::copy("Output00025.gzip.h5", "VEP/VEP00025.gzip.h5", copyOptions);
  fs::copy("Output00026.gzip.h5", "VEP/VEP00026.gzip.h5", copyOptions);
  fs::copy("Output00027.gzip.h5", "VEP/VEP00027.gzip.h5", copyOptions);
  fs::copy("Output00028.gzip.h5", "VEP/VEP00028.gzip.h5", copyOptions);
  fs::copy("Output00029.gzip.h5", "VEP/VEP00029.gzip.h5", copyOptions);
  fs::copy("Output00030.gzip.h5", "VEP/VEP00030.gzip.h5", copyOptions);
}

class ShearHeatingDuretz14 {
  static int SHD14SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SHD14SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SHD14SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (1 == 0 == 0) {
      return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
    } else {
      return instance->materials.rho[phase];
    }
  }

  static SetBC SHD14SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
    SetBC bc;
    if (position == W || position == E || position == S || position == N || position == SE || position == SW || position == NE || position == NW) {
      bc.type  = 0;
      bc.value = 0.0;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
    return bc;
  }


  static SetBC SHD14SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
    SetBC bc;
    if (position == W || position == E || position == S || position == N || position == SE || position == SW || position == NE || position == NW) {
      bc.type  = 0;
      bc.value = 0.0;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
    return bc;
  }

  public: void run() {
    MdoodzSetup shearHeatingDuretzModel = {
            .SetParticles = new (SetParticles_ff){
                    .SetPhase       = SHD14SetPhase,
                    .SetTemperature = SHD14SetTemperature,
                    .SetDensity     = SHD14SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx   = SetPureOrSimpleShearBCVx,
                    .SetBCVz   = SetPureOrSimpleShearBCVz,
                    .SetBCT    = SHD14SetBCT,
                    .SetBCTNew = SHD14SetBCTNew,
            },
    };
    RunMDOODZ("ShearHeatingDuretz14.txt", &shearHeatingDuretzModel);
  }
};

class GSE {
  static int GSESetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double A          = 2e-3 / instance->scaling.L;
    const double layer_bot0 = -5e-2 / instance->scaling.L;
    const double layer_top0 = 5e-2 / instance->scaling.L;
    const double Lx         = (instance->model.xmax - instance->model.xmin);
    const double layer_top  = layer_top0 - A * cos(coordinates.x * 2.0 * M_PI / Lx);
    const double layer_bot  = layer_bot0 + A * cos(coordinates.x * 2.0 * M_PI / Lx);
    if (coordinates.z > layer_bot && coordinates.z < layer_top) {
      return 1;
    } else {
      return 0;
    }
  }

  static double GSESetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.gs_ref[phase];
  }

  static double GSESetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {// phase
    return instance->materials.rho[phase];
  }

  static double GSESetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    const double T = (instance->model.user0 + zeroC) / instance->scaling.T;
    return T;
  }

  //----------------------------- THERMAL SetBC -----------------------------//


  static SetBC GSESetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  static SetBC GSESetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  public: void run() {
    MdoodzSetup gseModel = {
            .SetParticles = new (SetParticles_ff){
                    .SetPhase       = GSESetPhase,
                    .SetTemperature = GSESetTemperature,
                    .SetGrainSize   = GSESetGrainSize,
                    .SetDensity     = GSESetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx   = SetPureShearBCVx,
                    .SetBCVz   = SetPureShearBCVz,
                    .SetBCT    = GSESetBCT,
                    .SetBCTNew = GSESetBCTNew,
            },
    };
    RunMDOODZ("PinchSwellGSE.txt", &gseModel);
  }
};

class VEP_Duretz18 {
  static int VEPSetPhase(MdoodzInput *input, Coordinates coordinates) {
    const double radius = input->model.user1 / input->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double VEPSetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
    if (1 == 0) {
      return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
    } else {
      return input->materials.rho[phase];
    }
  }

  static SetBC GSESetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  static SetBC GSESetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  public: void run() {
    MdoodzSetup vepDuretzModel = {
            .SetParticles = new (SetParticles_ff){
                    .SetPhase       = VEPSetPhase,
                    .SetDensity = VEPSetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx   = SetPureShearBCVx,
                    .SetBCVz   = SetPureShearBCVz,
                    .SetBCT    = GSESetBCT,
                    .SetBCTNew = GSESetBCTNew,
            },
    };
    RunMDOODZ("VEP_Duretz18.txt", &vepDuretzModel);
  }
};

class Shrinking {
  static int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    const double radius = 0.25 / input->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {

    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Ax = cos( 6.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 6.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
      dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Ax = cos( 24.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 24.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
      dual_phase += input->model.Nb_phases;
    }

    return dual_phase;
  }

  static double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
    const double P_init = (input->model.bkg_pressure         ) / input->scaling.S;
    if (1 == 0) {
      return input->materials.rho[phase] * exp(input->materials.bet[phase]*P_init -  input->materials.alp[phase] * T_init);
    } else {
      return input->materials.rho[phase];
    }
  }

  public: int run() {
    MdoodzSetup setup = {
            .SetParticles  = new SetParticles_ff{
                    .SetPhase              = SetPhase,
                    .SetDualPhase          = SetDualPhase,
                    .SetDensity            = SetDensity,
            },
            .SetBCs = new SetBCs_ff {
                    .SetBCVx = SetPureOrSimpleShearBCVx,
                    .SetBCVz = SetPureOrSimpleShearBCVz,
            },
    };
    RunMDOODZ("Shrinking.txt", &setup);
  }
};

void RunTestCases() {
  (*new VEP_Duretz18).run();
  RenameVEPFiles();
  (*new GSE).run();
  RenameGSEFiles();
  (*new ShearHeatingDuretz14).run();
  rename("Output00005.gzip.h5", "ShearHeatingDuretz14.gzip.h5");
  (*new RiftingChenin).run();
  rename("Output00050.gzip.h5", "RiftingChenin50.gzip.h5");
  (*new ShearTemplate).run();
  rename("Output00005.gzip.h5", "ShearTemplate.gzip.h5");
  (*new ShearTemplate).run1();
  rename("Output00005.gzip.h5", "ShearTemplate1.gzip.h5");
  (*new ShearTemplate).runAniso();
  rename("Output00005.gzip.h5", "ShearTemplateAniso.gzip.h5");
  (*new TopoBenchCase1).run();
  RenameTopoBenchCaseFiles();
  (*new Shrinking).run();
}

string currentDateTime() {
  time_t    now = time(0);
  struct tm tstruct;
  char      buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

void UpdateReadmeTimestamp() {
  ifstream input("../VISUAL_TESTS/readme.md");
  ofstream myfile;
  myfile.open("readme.md");
  for (string line; getline(input, line);) {
    if (line.find("Last run date") != string::npos) {
      line = "Last run date: ";
      line.append(currentDateTime());
    }
    myfile << line << endl;
  }
  myfile.close();
  fs::copy("readme.md", "../VISUAL_TESTS/readme.md", fs::copy_options::update_existing);
}

int main() {
  RunTestCases();
  PlotShrinkingGif();
  PlotShrinkingGifRef();
  PlotGSE();
  PlotGSERef();
  PlotVEP();
  PlotVEPRef();
  PlotRiftingChenin();
  PlotRiftingCheninReference();
  PlotShearTemplate();
  PlotShearTemplateReference();
  PlotShearTemplate1();
  PlotShearTemplate1Reference();
  PlotShearTemplateAniso();
  PlotShearTemplateAnisoReference();
  PlotShearHeatingDuretz14();
  PlotShearHeatingDuretz14Reference();
  PlotTopoBenchCase1();
  UpdateReadmeTimestamp();
}