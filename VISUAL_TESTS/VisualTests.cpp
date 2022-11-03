extern "C" {
#include "mdoodz.h"
}
#include "visual-tests.h"
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
namespace fs = filesystem;

class MinimalistPrinter : public testing::EmptyTestEventListener {
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

  void OnTestStart(const testing::TestInfo &test_info) override {
    printf("*** Test %s.%s starting.\n",
           test_info.test_suite_name(), test_info.name());
  }

  void OnTestPartResult(const testing::TestPartResult &test_part_result) override {
    printf("%s in %s:%d\n%s\n",
           test_part_result.failed() ? "*** Failure" : "Success",
           test_part_result.file_name(),
           test_part_result.line_number(),
           test_part_result.summary());
  }

  void OnTestEnd(const testing::TestInfo &test_info) override {
    printf("*** Test %s.%s ending.\n", test_info.test_suite_name(), test_info.name());
    UpdateReadmeTimestamp();
  }
};

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::TestEventListeners &listeners = testing::UnitTest::GetInstance()->listeners();
  listeners.Append(new MinimalistPrinter);
  return RUN_ALL_TESTS();
}


class RiftingChenin : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .BuildInitialTopography = new BuildInitialTopography_ff{
                           .SetSurfaceZCoord = RPSetSurfaceZCoord,
            },
                   .SetParticles = new SetParticles_ff{
                           .SetHorizontalVelocity = SetHorizontalVelocity,
                           .SetVerticalVelocity   = SetVerticalVelocity,
                           .SetPhase              = SetPhase,
                           .SetTemperature        = SetTemperature,
                           .SetGrainSize          = SetGrainSize,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx    = SetPureShearBCVx,
                           .SetBCVz    = SetPureShearBCVz,
                           .SetBCT     = SetBCT,
                           .SetBCPType = SetBCPType,
                           .SetBCTNew  = SetBCTNew,
            },
                   .MutateInput = MutateInput};
  }

  static double RPSetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    const double TopoLevel = -0.0e3 / instance->scaling.L;
    const double h_pert    = instance->model.user3 / instance->scaling.L;
    return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
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

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
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

  static double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const int astenospherePhase = 3;
    return instance->materials.gs_ref[astenospherePhase];
  }

  static double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
    return -coordinates.x * instance->model.EpsBG;
  }

  static double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
    return coordinates.z * instance->model.EpsBG;
  }

  static char SetBCPType(MdoodzInput *instance, POSITION position) {
    if (position == NE || position == NW) {
      return 0;
    } else {
      return -1;
    }
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
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

  static SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
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

  static void MutateInput(MdoodzInput *instance, MutateInputParams *params) {
    int *astenospherePhases     = (int *) malloc(sizeof(int));
    astenospherePhases[0]       = {3};
    instance->crazyConductivity = new CrazyConductivity{
            .multiplier = 1000,
            .phases     = astenospherePhases,
            .nPhases    = 1,
    };
  }
};


class TopoBenchCase1 : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .BuildInitialTopography = new BuildInitialTopography_ff{
                           .SetSurfaceZCoord = SetSurfaceZCoord,
            },
                   .SetParticles = new SetParticles_ff{
                           .SetPhase = SetPhase,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx    = SetBCVx,
                           .SetBCVz    = SetBCVz,
                           .SetBCPType = SetBCPType,
            },
    };
  }

  static double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    double Amplitude  = 7e3 / instance->scaling.L;
    double Wavelength = 2800e3 / instance->scaling.L;
    return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
    if (coordinates.z > lithosphereBottomDepth) {
      return 1;
    } else {
      return 0;
    }
  }

  static SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
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

  static SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
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

  static char SetBCPType(MdoodzInput *instance, POSITION position) {
    if (position == NE || position == NW) {
      return 0;
    } else {
      return -1;
    }
  }
};

class ShearTemplate : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .SetParticles = new SetParticles_ff{
                           .SetPhase   = SetPhase,
                           .SetDensity = SetDensity,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx = SetPureOrSimpleShearBCVx,
                           .SetBCVz = SetPureOrSimpleShearBCVz,
            },
                   .MutateInput = MutateInput,
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (instance->model.eqn_state > 0) {
      return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
    } else {
      return instance->materials.rho[phase];
    }
  }

  static void MutateInput(MdoodzInput *input, MutateInputParams *mutateInputParams) {
    input->model.shear_style     = mutateInputParams->int1;
    input->model.writerSubfolder = mutateInputParams->str1;
  }
};


class ShearHeatingDuretz14 : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .SetParticles = new (SetParticles_ff){
                           .SetPhase       = SetPhase,
                           .SetTemperature = SetTemperature,
                           .SetDensity     = SetDensity,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx   = SetPureOrSimpleShearBCVx,
                           .SetBCVz   = SetPureOrSimpleShearBCVz,
                           .SetBCT    = SetBCT,
                           .SetBCTNew = SetBCTNew,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (instance->model.eqn_state > 0) {
      return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
    } else {
      return instance->materials.rho[phase];
    }
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
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


  static SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
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
};


class PinchSwellGSE : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .SetParticles = new (SetParticles_ff){
                           .SetPhase       = SetPhase,
                           .SetTemperature = SetTemperature,
                           .SetGrainSize   = SetGrainSize,
                           .SetDensity     = SetDensity,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx   = SetPureShearBCVx,
                           .SetBCVz   = SetPureShearBCVz,
                           .SetBCT    = SetBCT,
                           .SetBCTNew = SetBCTNew,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
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

  static double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.gs_ref[phase];
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {// phase
    return instance->materials.rho[phase];
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    const double T = (instance->model.user0 + zeroC) / instance->scaling.T;
    return T;
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  static SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }
};

class VEPDuretz18 : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .SetParticles = new (SetParticles_ff){
                           .SetPhase   = SetPhase,
                           .SetDensity = SetDensity,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx   = SetPureShearBCVx,
                           .SetBCVz   = SetPureShearBCVz,
                           .SetBCT    = SetBCT,
                           .SetBCTNew = SetBCTNew,
            },
    };
  }

  static int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    const double radius = input->model.user1 / input->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
    if (input->model.eqn_state > 0) {
      return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
    } else {
      return input->materials.rho[phase];
    }
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }

  static SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
    return {
            .value = gridTemperature,
            .type  = 0,
    };
  }
};

TEST_F(ShearTemplate, PureShear) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 0;       // shear_style
  mutateInputParams->str1              = testName;// writerSubfolder
  RunMDOODZ("ShearTemplate.txt", &setup);
  PlotShearTemplateReference();
  PlotShearTemplate();
}

TEST_F(ShearTemplate, SimpleShear) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 1;       // shear_style
  mutateInputParams->str1              = testName;// writerSubfolder
  RunMDOODZ("ShearTemplate.txt", &setup);
  PlotShearTemplate1Reference();
  PlotShearTemplate1();
}

TEST_F(ShearHeatingDuretz14, test) {
  RunMDOODZ("ShearHeatingDuretz14.txt", &setup);
  PlotShearHeatingDuretz14Reference();
  PlotShearHeatingDuretz14();
}

TEST_F(RiftingChenin, test) {
  RunMDOODZ("RiftingChenin.txt", &setup);
  PlotRiftingChenin();
  PlotRiftingCheninReference();
}

TEST_F(TopoBenchCase1, test) {
  RunMDOODZ("TopoBenchCase1.txt", &setup);
  PlotTopoBenchCase1();
}

TEST_F(VEPDuretz18, test) {
  RunMDOODZ("VEP_Duretz18.txt", &setup);
  PlotVEP();
  PlotVEPRef();
}

TEST_F(PinchSwellGSE, test) {
  RunMDOODZ("PinchSwellGSE.txt", &setup);
  PlotGSE();
  PlotGSERef();
}