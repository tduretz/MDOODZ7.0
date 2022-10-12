extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include <hdf5.h>
#include <math.h>

using ::testing::Gt;
using ::testing::Lt;

class MinimalistPrinter : public testing::EmptyTestEventListener {
  void RemoveDir(const char *dirName) {
    char command[500];
    snprintf(command, sizeof(command), R"(rm -r %s)", dirName);
    FILE *GNUplotPipe = popen(command, "w");
    fflush(GNUplotPipe);
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
    //RemoveDir(test_info.name());
  }
};

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::TestEventListeners &listeners = testing::UnitTest::GetInstance()->listeners();
  listeners.Append(new MinimalistPrinter);
  return RUN_ALL_TESTS();
}

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
    input->model.free_surf       = 0;
    input->model.writerSubfolder = mutateInputParams->str1;
    const int matrixPhase        = 0;
    if (mutateInputParams->int2) {
      input->materials.aniso_factor[matrixPhase] = 2.0;
      input->model.aniso                         = 1;
    } else {
      input->materials.aniso_factor[matrixPhase] = 1.0;
      input->model.aniso                         = 0;
    }
    if (mutateInputParams->int4) {
      input->materials.npwl[matrixPhase] = 3.0;
      input->materials.cstv[matrixPhase] = 0;
      input->materials.pwlv[matrixPhase] = 1;
    } else {
      input->materials.npwl[matrixPhase] = 1.0;
      input->materials.cstv[matrixPhase] = 1;
      input->materials.pwlv[matrixPhase] = 0;
    }
  }
};

int getStepsCount(char *hdf5FileName) {
  hid_t file            = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t iterationsGroup = H5Gopen(file, "Iterations", H5P_DEFAULT);
  hid_t numberStepsDataset =
          H5Dopen(iterationsGroup, "NumberSteps", H5P_DEFAULT);
  int numberStepsArray[1] = {-291029};
  H5Dread(numberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          numberStepsArray);
  int stepsCount = numberStepsArray[0];
  printf("Number of iterations: %d\n", stepsCount);
  return stepsCount;
}

TEST_F(ShearTemplate, LinearPureshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 0;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, LinearSimpleshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 1;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, LinearPureshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 0;// shear_style
  mutateInputParams->int2              = 1;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, LinearSimpleshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 1;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, NonLinearPureshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 0;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 1;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

TEST_F(ShearTemplate, NonLinearSimpleshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 1;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 1;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

TEST_F(ShearTemplate, NonLinearPureshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 0;// shear_style
  mutateInputParams->int2              = 1;// matrix aniso
  mutateInputParams->int4              = 1;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 20);
}

TEST_F(ShearTemplate, NonLinearSimpleshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int1              = 1;// shear_style
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 1;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("ShearTemplate.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

class SimpleRifting : public ::testing::Test {
  protected:
  MdoodzSetup setup;

  void        SetUp() {
    setup = {
                   .BuildInitialTopography = new BuildInitialTopography_ff{
                           .SetSurfaceZCoord = SetSurfaceZCoord,
            },
                   .SetParticles = new SetParticles_ff{
                           .SetHorizontalVelocity = SetHorizontalVelocity,
                           .SetVerticalVelocity   = SetVerticalVelocity,
                           .SetPhase              = SetPhase,
                           .SetGrainSize          = SetGrainSize,
            },
                   .SetBCs = new SetBCs_ff{
                           .SetBCVx    = SetPureShearBCVx,
                           .SetBCVz    = SetPureShearBCVz,
                           .SetBCT     = SetBCT,
                           .SetBCPType = SetBCPType,
                           .SetBCTNew  = SetBCTNew,
            },
                   .MutateInput = MutateInput,
    };
  }

  static double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    const double TopoLevel = -0.0e3 / instance->scaling.L;
    const double h_pert    = instance->model.user3 / instance->scaling.L;
    return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
    const bool   isBelowLithosphere   = coordinates.z < -lithosphereThickness;

    if (isBelowLithosphere) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    const int asthenospherePhase = 3;
    return instance->materials.gs_ref[asthenospherePhase];
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

  static void MutateInput(MdoodzInput *input, MutateInputParams *mutateInputParams) {
    input->model.writerSubfolder = mutateInputParams->str1;
    const int crustPhase         = 0;
    if (mutateInputParams->int2) { //  aniso
      input->materials.aniso_factor[crustPhase] = 2.0;
      input->model.aniso                        = 1;
      input->model.fstrain                      = 1;
      input->materials.pwlv[crustPhase]         = 0;
    } else {
      input->materials.aniso_factor[crustPhase] = 1.0;
      input->model.aniso                        = 0;
      input->model.fstrain                      = 0;
      input->materials.pwlv[crustPhase]         = 0;
    }
    if (mutateInputParams->int4) { // non-linear
      input->materials.npwl[crustPhase] = 3.0;
      input->materials.cstv[crustPhase] = 0;
      input->materials.pwlv[crustPhase] = 1;
    } else {
      input->materials.npwl[crustPhase] = 0;
      input->materials.cstv[crustPhase] = 1;
      input->materials.pwlv[crustPhase] = 0;
    }
  }
};

TEST_F(SimpleRifting, FreeSurfaceLinearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("RiftingChenin.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(SimpleRifting, FreeSurfaceLinearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int2              = 1; //  aniso
  mutateInputParams->int4              = 0; // non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("RiftingChenin.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(SimpleRifting, FreeSurfaceNonlinearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int2              = 0;// matrix aniso
  mutateInputParams->int4              = 1;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("RiftingChenin.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

TEST_F(SimpleRifting, FreeSurfaceNonlinearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  MutateInputParams *mutateInputParams = (MutateInputParams *) malloc(sizeof(MutateInputParams));
  setup.mutateInputParams              = mutateInputParams;
  mutateInputParams->int2              = 1;// matrix aniso
  mutateInputParams->int4              = 0;// non-linear
  mutateInputParams->str1              = testName;
  RunMDOODZ("RiftingChenin.txt", &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}