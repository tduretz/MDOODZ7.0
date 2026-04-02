extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class ConvergenceRate : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase       = SetPhase,
                    .SetTemperature = SetTemperature,
                    .SetDensity     = SetDensity,
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

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }
};

TEST_F(ConvergenceRate, NewtonPwl) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ConvergenceRate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int newtonSteps = getStepsCount(fileName);
  printf("Newton iterations for power-law: %d\n", newtonSteps);

  // Newton should converge (0 means immediate convergence, which is valid)
  ASSERT_GE(newtonSteps, 0);
  // Newton should converge quickly (quadratic convergence, typically < 15)
  EXPECT_LT(newtonSteps, 15);

  // Check that residuals are below the non-linear tolerance
  double finalRes = getFinalResidual(fileName, "rx_abs");
  printf("Newton final residual: %e\n", finalRes);
  EXPECT_LT(finalRes, 1e-8);

  free(inputName);
  free(fileName);
}

TEST_F(ConvergenceRate, PicardPwl) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ConvergenceRate/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int picardSteps = getStepsCount(fileName);
  printf("Picard iterations for power-law: %d\n", picardSteps);

  // Picard should converge (0 means immediate convergence, valid for well-conditioned problems)
  ASSERT_GE(picardSteps, 0);
  // Picard should converge (linearly) within reasonable iterations
  EXPECT_LT(picardSteps, 40);

  // Check final residual
  double finalRes = getFinalResidual(fileName, "rx_abs");
  printf("Picard final residual: %e\n", finalRes);
  EXPECT_LT(finalRes, 1e-8);

  free(inputName);
  free(fileName);
}

TEST_F(ConvergenceRate, NewtonFasterThanPicard) {
  // This test requires both Newton and Picard to have already run
  // and compares their iteration counts from the HDF5 output
  char *newtonFile;
  asprintf(&newtonFile, "%s/%s", "NewtonPwl", "Output00001.gzip.h5");
  char *picardFile;
  asprintf(&picardFile, "%s/%s", "PicardPwl", "Output00001.gzip.h5");

  int newtonSteps = getStepsCount(newtonFile);
  int picardSteps = getStepsCount(picardFile);
  printf("Newton: %d iterations, Picard: %d iterations\n", newtonSteps, picardSteps);

  // Newton (quadratic convergence) should use fewer iterations than Picard (linear)
  EXPECT_LE(newtonSteps, picardSteps);

  free(newtonFile);
  free(picardFile);
}
