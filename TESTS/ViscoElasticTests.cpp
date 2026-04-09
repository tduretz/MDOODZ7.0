extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class ViscoElastic : public ::testing::Test {
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

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }
};

TEST_F(ViscoElastic, StressAccumulation) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ViscoElastic/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  // Read first and last output
  char *file1;
  asprintf(&file1, "%s/%s", testName, "Output00001.gzip.h5");
  char *file5;
  asprintf(&file5, "%s/%s", testName, "Output00005.gzip.h5");

  int steps1 = getStepsCount(file1);
  int steps5 = getStepsCount(file5);
  ASSERT_GE(steps1, 0);
  ASSERT_GE(steps5, 0);

  // Elastic stress accumulation: deviatoric stress should grow over time steps
  // In Maxwell visco-elastic model, stress approaches 2*eta*eps_dot for large t/t_Maxwell
  double maxSxxd_1 = getMaxFieldValue(file1, "Centers", "sxxd");
  double maxSxxd_5 = getMaxFieldValue(file5, "Centers", "sxxd");
  printf("Max sxxd step 1: %e, step 5: %e\n", maxSxxd_1, maxSxxd_5);

  // Stress should increase (or at least not decrease) as elastic stress builds up  
  EXPECT_GE(fabs(maxSxxd_5), fabs(maxSxxd_1) * 0.9);  // allow small numerical tolerance

  // Check elastic strain rate is non-zero
  double maxEiiEl = getMaxFieldValue(file1, "Centers", "eII_el");
  printf("Max elastic strain rate invariant: %e\n", maxEiiEl);
  EXPECT_GT(maxEiiEl, 0.0);

  free(inputName);
  free(file1);
  free(file5);
}
