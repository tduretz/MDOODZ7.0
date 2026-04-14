extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class InterpEquiv : public ::testing::Test {
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
                    .SetBCVx  = SetPureShearBCVx,
                    .SetBCVz  = SetPureShearBCVz,
                    .SetBCT   = SetBCT,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    return 0;
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    const double Tbkg = (instance->model.user0 + zeroC) / instance->scaling.T;
    const double dT   = instance->model.user2 / instance->scaling.T;
    const double rad  = instance->model.user3 / instance->scaling.L;
    const double r2   = coordinates.x * coordinates.x + coordinates.z * coordinates.z;
    if (rad > 0.0 && r2 < rad * rad) {
      return Tbkg + dT;
    }
    return Tbkg;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates, double gridTemperature) {
    SetBC bc;
    double Ttop = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (position == N || position == NE || position == NW) {
      bc.value = Ttop;
      bc.type  = 1;
    } else {
      bc.value = 0.0;
      bc.type  = 0;
    }
    return bc;
  }
};

TEST_F(InterpEquiv, Mode0vsMode3) {
  // Run mode 0 (standard interpolation)
  RunMDOODZ("InterpEquiv/InterpEquiv_Mode0.txt", &setup);
  // Run mode 3 (fused P2Mastah interpolation)
  RunMDOODZ("InterpEquiv/InterpEquiv_Mode3.txt", &setup);

  const char *file0 = "InterpEquiv_Mode0/Output00005.gzip.h5";
  const char *file3 = "InterpEquiv_Mode3/Output00005.gzip.h5";

  // Compare temperature fields
  auto T0 = readFieldAsArray(file0, "Centers", "T");
  auto T3 = readFieldAsArray(file3, "Centers", "T");
  ASSERT_EQ(T0.size(), T3.size());
  double maxDiff_T = 0.0;
  for (size_t i = 0; i < T0.size(); i++) {
    double diff = fabs(T0[i] - T3[i]);
    if (diff > maxDiff_T) maxDiff_T = diff;
  }
  printf("T max|diff|: %e\n", maxDiff_T);
  EXPECT_LT(maxDiff_T, 1e-6);

  // Compare pressure fields
  auto P0 = readFieldAsArray(file0, "Centers", "P");
  auto P3 = readFieldAsArray(file3, "Centers", "P");
  ASSERT_EQ(P0.size(), P3.size());
  double maxDiff_P = 0.0;
  for (size_t i = 0; i < P0.size(); i++) {
    double diff = fabs(P0[i] - P3[i]);
    if (diff > maxDiff_P) maxDiff_P = diff;
  }
  printf("P max|diff|: %e\n", maxDiff_P);
  EXPECT_LT(maxDiff_P, 1e-6);

  // Compare viscosity fields (eta_n on centres)
  auto eta0 = readFieldAsArray(file0, "Centers", "eta_n");
  auto eta3 = readFieldAsArray(file3, "Centers", "eta_n");
  ASSERT_EQ(eta0.size(), eta3.size());
  double maxDiff_eta = 0.0;
  for (size_t i = 0; i < eta0.size(); i++) {
    double diff = fabs(eta0[i] - eta3[i]);
    if (diff > maxDiff_eta) maxDiff_eta = diff;
  }
  printf("eta_n max|diff|: %e\n", maxDiff_eta);
  EXPECT_LT(maxDiff_eta, 1e-6);
}
