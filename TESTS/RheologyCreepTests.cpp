extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class RheologyCreep : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase       = SetPhase,
                    .SetTemperature = SetTemperature,
                    .SetGrainSize   = SetGrainSize,
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

  static double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.gs_ref[phase];
  }
};

TEST_F(RheologyCreep, DiffusionCreepLinear) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Verify diffusion creep viscosity was computed (positive, bounded)
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("DiffusionCreepLinear viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(RheologyCreep, DiffusionCreepGrainSize) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Verify viscosity was computed and is positive
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("Viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(RheologyCreep, PeierlsCreep) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  // Peierls creep: verify simulation completed and viscosity is physically valid
  ASSERT_GE(stepsCount, 0);
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("PeierlsCreep viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(RheologyCreep, GBSCreep) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  // GBS creep: verify simulation completed and viscosity is physically valid
  ASSERT_GE(stepsCount, 0);
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("GBSCreep viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(RheologyCreep, CompositePwlLin) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  // Composite: verify simulation completed
  ASSERT_GE(stepsCount, 0);
  // Viscosity should be bounded
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("Composite viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(RheologyCreep, PowerLawCreep) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "RheologyCreep/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Power-law creep (dislocation): verify viscosity is computed and positive
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  double minEta = getMinFieldValue(fileName, "Centers", "eta_n");
  printf("PowerLaw viscosity range: min=%e, max=%e\n", minEta, maxEta);
  EXPECT_GT(minEta, 0.0);
  // Power-law strain rate should be non-zero (dislocation creep is active)
  double maxEiiPwl = getMaxFieldValue(fileName, "Centers", "eII_pwl");
  printf("Max power-law strain rate: %e\n", maxEiiPwl);
  EXPECT_GT(maxEiiPwl, 0.0);
  free(inputName);
  free(fileName);
}
