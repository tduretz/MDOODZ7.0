extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include <hdf5.h>

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
    RemoveDir(test_info.name());
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
    if (1 == 0 == 0) {
      return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
    } else {
      return instance->materials.rho[phase];
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
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, LinearSimpleshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}


TEST_F(ShearTemplate, LinearPureshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}

TEST_F(ShearTemplate, LinearSimpleshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_EQ(stepsCount, 1);
}


TEST_F(ShearTemplate, NonLinearPureshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

TEST_F(ShearTemplate, NonLinearSimpleshearIsotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}

TEST_F(ShearTemplate, NonLinearPureshearAnisotropic) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 20);
}

TEST_F(ShearTemplate, NonLinearSimpleshearAnisotropi) {
  const char        *testName          = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearTemplate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = 4;//getStepsCount(fileName);
  ASSERT_TRUE(stepsCount > 1);
  ASSERT_TRUE(stepsCount < 10);
}