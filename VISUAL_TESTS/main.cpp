extern "C" {
#include "mdoodz.h"
}
#include <cstdio>
#include <iostream>
#include <cmath>
#include "visual-tests.h"
#include <fstream>

namespace fs = std::filesystem;


double RPSetSurfaceZCoord(MdoodzInstance *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int RPSetPhase(MdoodzInstance *instance, Coordinates coordinates) {
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

double RPSetTemperature(MdoodzInstance *instance, Coordinates coordinates) {
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

double RPSetGrainSize(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  const int astenospherePhase = 3;
  return instance->materials.gs_ref[astenospherePhase];
}

double RPSetHorizontalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double RPSetVerticalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

SetBC RPSetBCVx(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = -coordinates.x * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}


SetBC RPSetBCVz(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = coordinates.z * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}

char RPSetBCPType(MdoodzInstance *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC RPSetBCT(MdoodzInstance *instance, POSITION position, double particleTemperature) {
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


SetBC RPSetBCTNew(MdoodzInstance *instance, POSITION position, double particleTemperature) {
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

MdoodzInstance CreateRiftingPaulineInstance() {
  int *astenospherePhases = (int *) malloc(sizeof(int));
  astenospherePhases[0]   = {3};
  return (MdoodzInstance){
          .inputFileName          = "RiftingPauline.txt",
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
                  .SetBCVx    = RPSetBCVx,
                  .SetBCVz    = RPSetBCVz,
                  .SetBCT     = RPSetBCT,
                  .SetBCPType = RPSetBCPType,
                  .SetBCTNew  = RPSetBCTNew,
          },
          .crazyConductivity = new CrazyConductivity{
                  .multiplier = 1000,
                  .phases     = astenospherePhases,
                  .nPhases    = 1,
          },
  };
}

double TBCSetSurfaceZCoord(MdoodzInstance *instance, double x_coord) {
  double Amplitude  = 7e3 / instance->scaling.L;
  double Wavelength = 2800e3 / instance->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
}

int TBCSetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
  if (coordinates.z > lithosphereBottomDepth) {
    return 1;
  } else {
    return 0;
  }
}

SetBC TBCSetBCVx(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
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

SetBC TBCSetBCVz(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
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

char TBCSetBCPType(MdoodzInstance *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

MdoodzInstance CreateTopoBenchCase1Instance() {
  int astenospherePhases[2] = {2, 4};
  return (MdoodzInstance){
          .inputFileName          = "TopoBenchCase1.txt",
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
          .crazyConductivity = new CrazyConductivity{
                  .multiplier = 1000,
                  .phases     = astenospherePhases,
                  .nPhases    = 2,
          }};
}

int STSetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double STSetDensity(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

SetBC STSetBCVx(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (instance->model.shear_style == 0) {
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
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == S || position == SE || position == SW) {
      bc.value = -instance->model.EpsBG * Lz;
      bc.type  = 11;
    } else if (position == N || position == NE || position == NW) {
      bc.value = instance->model.EpsBG * Lz;
      bc.type  = 11;
    } else if (position == E) {
      bc.value = 0.0;
      bc.type  = -12;
    } else if (position == W) {
      bc.value = 0.0;
      bc.type  = -2;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  }
  return bc;
}

SetBC STSetBCVz(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (instance->model.shear_style == 0) {
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
  } else {
    if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = -12;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  }
  return bc;
}

MdoodzInstance CreateShearTemplateInstance(int shear_style) {
  char *inputFileName;
  if (shear_style) {
    inputFileName = "ShearTemplate1.txt";
  } else {
    inputFileName = "ShearTemplate.txt";
  }
  return (MdoodzInstance){
          .inputFileName = inputFileName,
          .SetParticles  = new SetParticles_ff{
                   .SetPhase   = STSetPhase,
                   .SetDensity = STSetDensity,
          },
          .SetBCs = new SetBCs_ff{
                  .SetBCVx = STSetBCVx,
                  .SetBCVz = STSetBCVz,
          },
  };
}

void RenameTopoBenchCaseFiles() {
  std::rename("Output00000.gzip.h5", "TopoBenchCase00000.gzip.h5");
  std::rename("Output00010.gzip.h5", "TopoBenchCase00010.gzip.h5");
  std::rename("Output00020.gzip.h5", "TopoBenchCase00020.gzip.h5");
  std::rename("Output00030.gzip.h5", "TopoBenchCase00030.gzip.h5");
  std::rename("Output00040.gzip.h5", "TopoBenchCase00040.gzip.h5");
  std::rename("Output00050.gzip.h5", "TopoBenchCase00050.gzip.h5");
  std::rename("Output00060.gzip.h5", "TopoBenchCase00060.gzip.h5");
  std::rename("Output00070.gzip.h5", "TopoBenchCase00070.gzip.h5");
  std::rename("Output00080.gzip.h5", "TopoBenchCase00080.gzip.h5");
  std::rename("Output00090.gzip.h5", "TopoBenchCase00090.gzip.h5");
  std::rename("Output00100.gzip.h5", "TopoBenchCase00100.gzip.h5");
}

void RunTestCases() {
  MdoodzInstance riftingPauline = CreateRiftingPaulineInstance();
  RunMDOODZ(&riftingPauline);
  std::rename("Output00050.gzip.h5", "RiftingPauline50.gzip.h5");
  MdoodzInstance shearTemplate = CreateShearTemplateInstance(0);
  RunMDOODZ(&shearTemplate);
  std::rename("Output00005.gzip.h5", "ShearTemplate.gzip.h5");
  MdoodzInstance shearTemplate1 = CreateShearTemplateInstance(1);
  RunMDOODZ(&shearTemplate1);
  std::rename("Output00005.gzip.h5", "ShearTemplate1.gzip.h5");
  MdoodzInstance topoBenchCase1 = CreateTopoBenchCase1Instance();
  RunMDOODZ(&topoBenchCase1);
  RenameTopoBenchCaseFiles();
}

std::string currentDateTime() {
  time_t    now = time(0);
  struct tm tstruct;
  char      buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

void UpdateReadmeTimestamp() {
  std::ifstream input("../VISUAL_TESTS/readme.md");
  std::ofstream myfile;
  myfile.open("readme.md");
  for (std::string line; getline(input, line);) {
    if (line.find("Last update date") != std::string::npos) {
      line = "Last update date: ";
      line.append(currentDateTime());
    }
    myfile << line << std::endl;
  }
  myfile.close();
  fs::copy("readme.md", "../VISUAL_TESTS/readme.md", fs::copy_options::update_existing);
}

int main() {
  RunTestCases();
  PlotRiftingPauline();
  PlotRiftingPaulineReference();
  PlotTopoBenchCase1();
  PlotShearTemplateReference();
  PlotShearTemplate();
  PlotShearTemplate1Reference();
  PlotShearTemplate1();
  UpdateReadmeTimestamp();
}