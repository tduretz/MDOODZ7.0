extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class Compressibility : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = SetPhase,
                    .SetDensity = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureShearBCVx,
                    .SetBCVz = SetPureShearBCVz,
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

TEST_F(Compressibility, DilatantFlow) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Compressibility/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // With compressible=1, elastic compressibility (bet > 0) allows
  // small velocity divergence in response to pressure changes
  // Check that the simulation completes and produces valid fields
  double maxP = getMaxFieldValue(fileName, "Centers", "P");
  double minP = getMinFieldValue(fileName, "Centers", "P");
  printf("Pressure range: min=%e, max=%e\n", minP, maxP);
  EXPECT_LT(fabs(maxP), 1e40);

  // Viscosity should be positive and bounded
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("Viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);

  free(inputName);
  free(fileName);
}
