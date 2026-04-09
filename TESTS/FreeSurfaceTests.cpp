extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class FreeSurface : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .BuildInitialTopography = new BuildInitialTopography_ff{
                    .SetSurfaceZCoord = SetSurfaceZCoord,
            },
            .SetParticles = new SetParticles_ff{
                    .SetPhase       = SetPhase,
                    .SetTemperature = SetTemperature,
                    .SetDensity     = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetNoSlipBCVx,
                    .SetBCVz = SetNoSlipBCVz,
            },
    };
  }

  static double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
    return 0.0;
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    // Place a dense block in the upper middle
    const double z0 = 0.5 * (instance->model.zmax - instance->model.zmin) * 0.3;
    if (coordinates.x * coordinates.x + (coordinates.z - z0) * (coordinates.z - z0) < radius * radius) {
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

  static SetBC SetNoSlipBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == N || position == NW || position == NE ||
        position == S || position == SW || position == SE ||
        position == W || position == E) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static SetBC SetNoSlipBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == W || position == E ||
        position == SW || position == SE || position == NW || position == NE ||
        position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }
};

TEST_F(FreeSurface, SinkingBlock) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "FreeSurface/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  // Read output after 3 steps
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00003.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // With free surface, a denser block (rho=3300 vs 3000) should sink under gravity
  // This creates vertical velocities and topographic deflection
  double maxVz = getMaxFieldValue(fileName, "VzNodes", "Vz");
  double minVz = getMinFieldValue(fileName, "VzNodes", "Vz");
  printf("Vz range: [%e, %e]\n", minVz, maxVz);

  // The block should generate non-zero vertical velocity (sinking)
  double maxAbsVz = fmax(fabs(maxVz), fabs(minVz));
  EXPECT_GT(maxAbsVz, 0.0);

  // Pressure should be positive (lithostatic-like)
  double maxP = getMaxFieldValue(fileName, "Centers", "P");
  printf("Max pressure: %e\n", maxP);
  EXPECT_GT(maxP, 0.0);

  free(inputName);
  free(fileName);
}
