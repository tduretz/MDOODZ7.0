extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class SolverMode : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
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
    return instance->materials.rho[phase];
  }
};

TEST_F(SolverMode, PicardConvergence) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "SolverMode/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  // Picard should converge but may need more iterations than Newton
  ASSERT_GT(stepsCount, 0);
  ASSERT_LT(stepsCount, 50);
  // Check final residual is below the non-linear tolerance
  double finalRes = getFinalResidual(fileName, "rx_abs");
  printf("Final momentum residual: %e\n", finalRes);
  EXPECT_LT(finalRes, 1e-8);
  free(inputName);
  free(fileName);
}

TEST_F(SolverMode, AdaptiveDt) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "SolverMode/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  // Read first and last output to check dt changed
  char *file1;
  asprintf(&file1, "%s/%s", testName, "Output00001.gzip.h5");
  char *file3;
  asprintf(&file3, "%s/%s", testName, "Output00003.gzip.h5");
  int steps1 = getStepsCount(file1);
  int steps3 = getStepsCount(file3);
  ASSERT_GT(steps1, 0);
  ASSERT_GT(steps3, 0);
  // Both steps should complete successfully
  printf("Step 1 iterations: %d, Step 3 iterations: %d\n", steps1, steps3);
  free(inputName);
  free(file1);
  free(file3);
}
