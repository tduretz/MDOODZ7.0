extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class NeumannBC : public ::testing::Test {
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
    return 0;
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }
};

TEST_F(NeumannBC, StressBCPureShear) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "NeumannBC/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // The built-in SetPureShearBCVx uses type=13 (constant_shear_stress) on N/S
  // This implements free-slip: zero shear stress on horizontal boundaries
  // Verify that the shear stress is near zero on vertices (which includes boundaries)
  double maxSxz = getMaxFieldValue(fileName, "Vertices", "sxz");
  double minSxz = getMinFieldValue(fileName, "Vertices", "sxz");
  printf("Shear stress range: min=%e, max=%e\n", minSxz, maxSxz);

  // In pure shear with free-slip, shear stress should be small everywhere
  // (it's exactly zero on boundaries and small in the interior for cstv)
  EXPECT_LT(fabs(maxSxz), 1.0);
  EXPECT_LT(fabs(minSxz), 1.0);

  // Deviatoric normal stress should be non-zero from pure shear deformation
  double maxSxxd = getMaxFieldValue(fileName, "Centers", "sxxd");
  printf("Max sxxd: %e\n", maxSxxd);
  EXPECT_GT(fabs(maxSxxd), 0.0);

  // Velocity should be non-zero (driven by bkg_strain_rate)
  double maxVx = getMaxFieldValue(fileName, "VxNodes", "Vx");
  printf("Max Vx: %e\n", maxVx);
  EXPECT_GT(fabs(maxVx), 0.0);

  free(inputName);
  free(fileName);
}
