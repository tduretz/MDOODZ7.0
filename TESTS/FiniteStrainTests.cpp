extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class FiniteStrain : public ::testing::Test {
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

TEST_F(FiniteStrain, PureShearStrain) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "FiniteStrain/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  // Read output after 5 steps
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00005.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // Finite strain tensor F = [[Fxx, Fxz],[Fzx, Fzz]]
  // In pure shear ALE mode, grid deforms with material so F tracks
  // cumulative deformation. After 5 steps of pure shear:
  double maxFxx = getMaxFieldValue(fileName, "Centers", "Fxx");
  double minFzz = getMinFieldValue(fileName, "Centers", "Fzz");
  double meanFxx = getMeanFieldValue(fileName, "Centers", "Fxx");
  double meanFzz = getMeanFieldValue(fileName, "Centers", "Fzz");
  printf("Max Fxx: %e, Min Fzz: %e\n", maxFxx, minFzz);
  printf("Mean Fxx: %e, Mean Fzz: %e\n", meanFxx, meanFzz);

  // Fxx and Fzz should be positive (physical requirement)
  EXPECT_GT(maxFxx, 0.0);
  EXPECT_GT(minFzz, 0.0);

  // Check area conservation: det(F) = Fxx*Fzz - Fxz*Fzx ~ 1 for incompressible flow
  double detF = meanFxx * meanFzz;
  printf("det(F) ~ %e\n", detF);
  // det(F) should be close to 1 (volume preserving)
  EXPECT_NEAR(detF, 1.0, 0.2);

  free(inputName);
  free(fileName);
}
