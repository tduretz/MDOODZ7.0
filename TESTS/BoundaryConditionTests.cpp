extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class BoundaryCondition : public ::testing::Test {
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

  static SetBC SetNoSlipBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == N || position == NW || position == NE ||
        position == S || position == SW || position == SE) {
      bc.value = 0.0;
      bc.type  = 0;
    } else if (position == W || position == E) {
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
        position == SW || position == SE || position == NW || position == NE) {
      bc.value = 0.0;
      bc.type  = 0;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }
};

TEST_F(BoundaryCondition, PeriodicSimpleShear) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "BoundaryCondition/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GT(stepsCount, 0);
  // For periodic BCs, pressure should be well-defined
  double maxP = getMaxFieldValue(fileName, "Centers", "P");
  double minP = getMinFieldValue(fileName, "Centers", "P");
  printf("Pressure range: min=%e, max=%e\n", minP, maxP);
  // Pressure should be finite
  EXPECT_LT(maxP, 1e40);
  EXPECT_GT(minP, -1e40);
  free(inputName);
  free(fileName);
}

TEST_F(BoundaryCondition, FreeSlipPureShear) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "BoundaryCondition/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GT(stepsCount, 0);
  // Free slip: shear stress should be near zero at boundaries
  // Check that the solution converged and stresses are bounded
  double maxSxz = getMaxFieldValue(fileName, "Vertices", "sxz");
  double minSxz = getMinFieldValue(fileName, "Vertices", "sxz");
  printf("Shear stress range: min=%e, max=%e\n", minSxz, maxSxz);
  EXPECT_LT(maxSxz, 1e40);
  EXPECT_GT(minSxz, -1e40);
  // Velocity field check: Vx should be anti-symmetric (extensional in x)
  double maxVx = getMaxFieldValue(fileName, "VxNodes", "Vx");
  double minVx = getMinFieldValue(fileName, "VxNodes", "Vx");
  printf("Vx range: [%e, %e]\n", minVx, maxVx);
  EXPECT_GT(maxVx, 0.0);
  EXPECT_LT(minVx, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(BoundaryCondition, NoSlip) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "BoundaryCondition/%s.txt", testName);

  // Override with no-slip BCs
  MdoodzSetup noSlipSetup = {
          .SetParticles = new SetParticles_ff{
                  .SetPhase   = SetPhase,
                  .SetDensity = SetDensity,
          },
          .SetBCs = new SetBCs_ff{
                  .SetBCVx = SetNoSlipBCVx,
                  .SetBCVz = SetNoSlipBCVz,
          },
  };

  RunMDOODZ(inputName, &noSlipSetup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GT(stepsCount, 0);
  // Solution should converge
  double maxP = getMaxFieldValue(fileName, "Centers", "P");
  printf("Max pressure: %e\n", maxP);
  EXPECT_LT(maxP, 1e40);
  // With no-slip BCs and gravity, velocities exist due to buoyancy but are bounded
  double maxAbsVx = getMaxAbsFieldValue(fileName, "VxNodes", "Vx");
  double maxAbsVz = getMaxAbsFieldValue(fileName, "VzNodes", "Vz");
  printf("Max |Vx|: %e, Max |Vz|: %e\n", maxAbsVx, maxAbsVz);
  EXPECT_LT(maxAbsVx, 1e40);
  EXPECT_LT(maxAbsVz, 1e40);
  free(inputName);
  free(fileName);
}
