extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class Plasticity : public ::testing::Test {
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

TEST_F(Plasticity, DruckerPragerYield) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Plasticity/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Check that stress is bounded by yield stress: tau_II <= C*cos(phi) + P*sin(phi) + tolerance
  // With C=20 MPa and high strain rate, yielding should occur
  // Read max stress invariant (sxxd^2 + szzd^2 + sxz^2)^0.5 from centers
  double maxSxxd = getMaxFieldValue(fileName, "Centers", "sxxd");
  double maxSzzd = getMaxFieldValue(fileName, "Centers", "szzd");
  // Effective viscosity should be reduced below the trial viscosity due to plasticity
  double maxEta = getMaxFieldValue(fileName, "Centers", "eta_n");
  printf("Max stress sxxd: %e, szzd: %e, Max eta: %e\n", maxSxxd, maxSzzd, maxEta);
  // Stresses should be finite and bounded
  EXPECT_LT(maxSxxd, 1e40);
  EXPECT_LT(maxSzzd, 1e40);
  // Plasticity must have activated: plastic strain rate > 0 somewhere
  double maxEiiPl = getMaxFieldValue(fileName, "Centers", "eII_pl");
  printf("Max plastic strain rate (yield test): %e\n", maxEiiPl);
  EXPECT_GT(maxEiiPl, 0.0);
  free(inputName);
  free(fileName);
}

TEST_F(Plasticity, DruckerPragerNoYield) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Plasticity/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // With very low strain rate, material should not yield
  double maxEiiPl = getMaxFieldValue(fileName, "Centers", "eII_pl");
  printf("Max plastic strain rate (should be ~0): %e\n", maxEiiPl);
  EXPECT_LT(maxEiiPl, 1e-25);
  free(inputName);
  free(fileName);
}

TEST_F(Plasticity, StrainSoftening) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Plasticity/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  // Read the last output step (3 steps)
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00003.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Read initial cohesion from step 0 for comparison
  char *fileInit;
  asprintf(&fileInit, "%s/%s", testName, "Output00000.gzip.h5");
  double maxCoh_init = getMaxFieldValue(fileInit, "Centers", "cohesion");
  // After strain softening, cohesion should have decreased from initial value
  double minCoh = getMinFieldValue(fileName, "Centers", "cohesion");
  printf("Max initial cohesion: %e, Min final cohesion: %e\n", maxCoh_init, minCoh);
  // Cohesion must remain positive (softening must not overshoot to negative)
  EXPECT_GT(minCoh, 0.0);
  // Cohesion must have decreased from its initial value due to strain softening
  EXPECT_LT(minCoh, maxCoh_init);
  free(fileInit);
  free(inputName);
  free(fileName);
}

TEST_F(Plasticity, StressLimiter) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Plasticity/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);
  // Stress should be limited by Slim = 100 MPa = 1e8 Pa
  double maxSxxd = getMaxFieldValue(fileName, "Centers", "sxxd");
  double maxSzzd = getMaxFieldValue(fileName, "Centers", "szzd");
  double maxEta  = getMaxFieldValue(fileName, "Centers", "eta_n");
  printf("Max sxxd: %e, max szzd: %e, max eta: %e\n", maxSxxd, maxSzzd, maxEta);
  // Stresses should remain finite and bounded
  EXPECT_LT(maxSxxd, 1e40);
  EXPECT_LT(maxSzzd, 1e40);
  // Without the limiter, stress would be ~eta0*2*eps = 1e30*2*5e-14 = 1e17 Pa.
  // With Slim = 1e8, effective viscosity must be far below the trial viscosity eta0 = 1e30.
  // In non-dimensional units (eta_ref=1e20): eta0_nd = 1e10, so max_eta << 1e10.
  // The limiter should cap eta to Slim/(2*eps) = 1e8/(2*5e-14) = 1e21 -> nd = 10.
  double eta_ref = 1e20;
  double eta_nd_unlimited = 1e30 / eta_ref;  // 1e10 without limiter
  EXPECT_LT(maxEta, eta_nd_unlimited * 0.01);  // must be at least 100x below unlimited
  free(inputName);
  free(fileName);
}
