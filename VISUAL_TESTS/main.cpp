extern "C" {
#include "mdoodz.h"
}
#include "visual-tests.h"
#include <cmath>
#include <cstdio>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>

using namespace std;
namespace fs = std::experimental::filesystem;


double RPSetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int RPSetPhase(MdoodzInput *instance, Coordinates coordinates) {
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

double RPSetTemperature(MdoodzInput *instance, Coordinates coordinates) {
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

double RPSetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const int astenospherePhase = 3;
  return instance->materials.gs_ref[astenospherePhase];
}

double RPSetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double RPSetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

char RPSetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC RPSetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else {
    bc.value = 0.0;
    bc.type  = 0;
  }
  return bc;
}

SetBC RPSetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
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

void RPMutateInput(MdoodzInput *instance, MutateInputParams *params) {
  int *astenospherePhases     = (int *) malloc(sizeof(int));
  astenospherePhases[0]       = {3};
  instance->crazyConductivity = new CrazyConductivity{
          .multiplier = 1000,
          .phases     = astenospherePhases,
          .nPhases    = 1,
  };
}

MdoodzSetup CreateRiftingCheninInstance() {
  return (MdoodzSetup){
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
}

double TBCSetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  double Amplitude  = 7e3 / instance->scaling.L;
  double Wavelength = 2800e3 / instance->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
}

int TBCSetPhase(MdoodzInput *instance, Coordinates coordinates) {
  double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
  if (coordinates.z > lithosphereBottomDepth) {
    return 1;
  } else {
    return 0;
  }
}

SetBC TBCSetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
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

SetBC TBCSetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
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

char TBCSetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

MdoodzSetup CreateTopoBenchCase1Instance() {
  return (MdoodzSetup){
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
}

int STSetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double STSetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

MdoodzSetup CreateShearTemplateInstance() {
  return (MdoodzSetup){
          .SetParticles = new SetParticles_ff{
                  .SetPhase   = STSetPhase,
                  .SetDensity = STSetDensity,
          },
          .SetBCs = new SetBCs_ff{
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
}

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

int SHD14SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SHD14SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  return (instance->model.user0 + zeroC) / instance->scaling.T;
}

double SHD14SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

SetBC SHD14SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
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


SetBC SHD14SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
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

MdoodzSetup CreateShearHeatingDuretz14Instance() {
  return (MdoodzSetup){
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
}

int GSESetPhase(MdoodzInput *instance, Coordinates coordinates) {
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
double GSESetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.gs_ref[phase];
}

double GSESetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {// phase
  return instance->materials.rho[phase];
}

double GSESetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double T = (instance->model.user0 + zeroC) / instance->scaling.T;
  return T;
}

//----------------------------- THERMAL SetBC -----------------------------//


SetBC GSESetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
  return {
          .value = gridTemperature,
          .type  = 0,
  };
}

SetBC GSESetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
  return {
          .value = gridTemperature,
          .type  = 0,
  };
}

MdoodzSetup CreateGSEInstance() {
  return (MdoodzSetup){
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
}


int VEPSetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double VEPSetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}


MdoodzSetup CreateVEPInstance() {
  return (MdoodzSetup){
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
}

void RunTestCases() {
  MdoodzSetup vep = CreateVEPInstance();
  RunMDOODZ("VEP_Duretz18.txt", &vep);
  RenameVEPFiles();
  MdoodzSetup gse = CreateGSEInstance();
  RunMDOODZ("PinchSwellGSE.txt", &gse);
  RenameGSEFiles();
  MdoodzSetup shearHeatingDuretz14 = CreateShearHeatingDuretz14Instance();
  RunMDOODZ("ShearHeatingDuretz14.txt", &shearHeatingDuretz14);
  rename("Output00005.gzip.h5", "ShearHeatingDuretz14.gzip.h5");
  MdoodzSetup RiftingChenin = CreateRiftingCheninInstance();
  RunMDOODZ("RiftingChenin.txt", &RiftingChenin);
  rename("Output00050.gzip.h5", "RiftingChenin50.gzip.h5");
  MdoodzSetup shearTemplate = CreateShearTemplateInstance();
  RunMDOODZ("ShearTemplate.txt", &shearTemplate);
  rename("Output00005.gzip.h5", "ShearTemplate.gzip.h5");
  MdoodzSetup shearTemplate1 = CreateShearTemplateInstance();
  RunMDOODZ("ShearTemplate1.txt", &shearTemplate1);
  rename("Output00005.gzip.h5", "ShearTemplate1.gzip.h5");
  MdoodzSetup topoBenchCase1 = CreateTopoBenchCase1Instance();
  RunMDOODZ("TopoBenchCase1.txt", &topoBenchCase1);
  RenameTopoBenchCaseFiles();
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
  PlotGSE();
  PlotGSERef();
  PlotVEP();
  PlotVEPRef();
  PlotRiftingChenin();
  PlotShearHeatingDuretz14Reference();
  PlotRiftingCheninReference();
  PlotShearTemplate();
  PlotShearTemplateReference();
  PlotShearTemplate1();
  PlotShearTemplate1Reference();
  PlotShearHeatingDuretz14();
  PlotTopoBenchCase1();
  UpdateReadmeTimestamp();
}